! Orchestration of the land model. Per grid cell and time step:
!   1. each tile runs its canopy (interception of water and heat)
!   2. tile fluxes are aggregated, weighted by tile fraction, onto the shared soil column
!   3. soil water (Richards) then soil heat (diffusion) advance the column
! All steps are functions of their arguments: no state is held in the module.
module land_mod
  use kinds_mod, only: wp
  use constants_mod, only: water_density, latent_heat_vaporisation
  use types_mod, only: SoilType, SoilColumn, CanopyParameters, CanopyState, TileFluxes, &
                       CanopyStepResult, Tile, CellDiagnostics, GridCell, AtmosphericForcing
  use soil_properties_mod, only: water_content
  use soil_water_mod, only: SoilWaterResult, soil_water_step, column_water_storage
  use soil_heat_mod, only: soil_heat_step
  use canopy_mod, only: canopy_step
  implicit none
  private

  public :: geometric_layer_thicknesses, build_grid_cell, aggregate_tile_fluxes
  public :: cell_step, land_step

contains

  pure function geometric_layer_thicknesses(n_layers, top_layer_thickness, growth_factor) result(thicknesses)
    integer, intent(in) :: n_layers
    real(wp), intent(in) :: top_layer_thickness, growth_factor
    real(wp) :: thicknesses(n_layers)
    integer :: layer
    do layer = 1, n_layers
      thicknesses(layer) = top_layer_thickness * growth_factor**(layer - 1)
    end do
  end function geometric_layer_thicknesses

  ! Create a cell in hydrostatic equilibrium with a water table at the given depth.
  ! One tile is created per entry of tile_fractions, with tile_type_index equal to its position.
  function build_grid_cell(soil_type_index, soil, tile_fractions, layer_thickness, water_table_depth, &
                           initial_temperature) result(cell)
    integer, intent(in) :: soil_type_index
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: tile_fractions(:), layer_thickness(:), water_table_depth, initial_temperature
    type(GridCell) :: cell

    real(wp), parameter :: fraction_tolerance = 1.0e-8_wp
    real(wp) :: layer_centre_depth(size(layer_thickness))
    integer :: layer_index
    integer :: tile_index

    if (abs(sum(tile_fractions) - 1.0_wp) > fraction_tolerance) then
      error stop 'build_grid_cell: tile fractions must sum to 1'
    end if

    cell%soil_type_index = soil_type_index
    layer_centre_depth = cumulative_sum(layer_thickness) - 0.5_wp * layer_thickness
    cell%column%layer_thickness = layer_thickness
    cell%column%matric_head = layer_centre_depth - water_table_depth
    cell%column%theta = water_content(soil, cell%column%matric_head)
    cell%column%temperature = [(initial_temperature, layer_index = 1, size(layer_thickness))]

    allocate(cell%tiles(size(tile_fractions)))
    do tile_index = 1, size(tile_fractions)
      cell%tiles(tile_index)%tile_type_index = tile_index
      cell%tiles(tile_index)%fraction = tile_fractions(tile_index)
      cell%tiles(tile_index)%canopy_state%temperature = initial_temperature
    end do

    cell%diagnostics%soil_water_storage = column_water_storage(soil, layer_thickness, cell%column%matric_head)
  end function build_grid_cell

  pure function cumulative_sum(values) result(running_total)
    real(wp), intent(in) :: values(:)
    real(wp) :: running_total(size(values))
    integer :: index
    running_total(1) = values(1)
    do index = 2, size(values)
      running_total(index) = running_total(index - 1) + values(index)
    end do
  end function cumulative_sum

  ! Fraction-weighted mean of the tile fluxes.
  pure function aggregate_tile_fluxes(tiles) result(mean_fluxes)
    type(Tile), intent(in) :: tiles(:)
    type(TileFluxes) :: mean_fluxes
    integer :: tile_index
    do tile_index = 1, size(tiles)
      mean_fluxes%throughfall = mean_fluxes%throughfall &
          + tiles(tile_index)%fraction * tiles(tile_index)%fluxes%throughfall
      mean_fluxes%canopy_evaporation = mean_fluxes%canopy_evaporation &
          + tiles(tile_index)%fraction * tiles(tile_index)%fluxes%canopy_evaporation
      mean_fluxes%soil_evaporation_demand = mean_fluxes%soil_evaporation_demand &
          + tiles(tile_index)%fraction * tiles(tile_index)%fluxes%soil_evaporation_demand
      mean_fluxes%ground_heat_flux = mean_fluxes%ground_heat_flux &
          + tiles(tile_index)%fraction * tiles(tile_index)%fluxes%ground_heat_flux
    end do
  end function aggregate_tile_fluxes

  ! Soil evaporation (kg m-2 s-1): demand scaled by top-layer wetness and capped so that a step
  ! never removes more than half of the water above residual in the top layer.
  pure function limit_soil_evaporation(soil, column, evaporation_demand, dt) result(evaporation)
    type(SoilType), intent(in) :: soil
    type(SoilColumn), intent(in) :: column
    real(wp), intent(in) :: evaporation_demand, dt
    real(wp) :: evaporation

    real(wp) :: available_water, wetness
    available_water = max(column%theta(1) - soil%theta_res, 0.0_wp)
    wetness = min(available_water / (soil%theta_sat - soil%theta_res), 1.0_wp)
    evaporation = min(evaporation_demand * wetness, &
                      0.5_wp * available_water * column%layer_thickness(1) * water_density / dt)
  end function limit_soil_evaporation

  pure function cell_step(cell, soil, tile_library, forcing, dt) result(new_cell)
    type(GridCell), intent(in) :: cell
    type(SoilType), intent(in) :: soil
    type(CanopyParameters), intent(in) :: tile_library(:)
    type(AtmosphericForcing), intent(in) :: forcing
    real(wp), intent(in) :: dt
    type(GridCell) :: new_cell

    type(CanopyStepResult) :: canopy_result
    type(TileFluxes) :: cell_fluxes
    type(SoilWaterResult) :: water_result
    real(wp) :: soil_evaporation, heat_flux_into_soil
    real(wp) :: canopy_storage_change, soil_storage_old, soil_storage_new
    integer :: tile_index

    new_cell = cell

    ! 1. Canopy, per tile. All tiles see the same soil surface temperature.
    canopy_storage_change = 0.0_wp
    do tile_index = 1, size(cell%tiles)
      canopy_result = canopy_step(tile_library(cell%tiles(tile_index)%tile_type_index), &
                                  cell%tiles(tile_index)%canopy_state, forcing, &
                                  cell%column%temperature(1), dt)
      new_cell%tiles(tile_index)%canopy_state = canopy_result%state
      new_cell%tiles(tile_index)%fluxes = canopy_result%fluxes
      canopy_storage_change = canopy_storage_change + cell%tiles(tile_index)%fraction &
          * (canopy_result%state%water_store - cell%tiles(tile_index)%canopy_state%water_store)
    end do

    ! 2. Aggregate onto the soil column
    cell_fluxes = aggregate_tile_fluxes(new_cell%tiles)
    soil_evaporation = limit_soil_evaporation(soil, cell%column, cell_fluxes%soil_evaporation_demand, dt)

    ! 3. Soil water, then soil heat using the updated water content
    water_result = soil_water_step(soil, cell%column%layer_thickness, cell%column%matric_head, &
                                   cell_fluxes%throughfall, soil_evaporation, dt)
    new_cell%column%matric_head = water_result%matric_head
    new_cell%column%theta = water_result%theta
    heat_flux_into_soil = cell_fluxes%ground_heat_flux - latent_heat_vaporisation * soil_evaporation
    new_cell%column%temperature = soil_heat_step(soil, cell%column%layer_thickness, water_result%theta, &
                                                 cell%column%temperature, heat_flux_into_soil, dt)

    ! Diagnostics, including a closed water budget: P = d(canopy) + E_canopy + E_soil + runoff + drainage + d(soil)
    soil_storage_old = column_water_storage(soil, cell%column%layer_thickness, cell%column%matric_head)
    soil_storage_new = column_water_storage(soil, cell%column%layer_thickness, water_result%matric_head)
    new_cell%diagnostics = CellDiagnostics( &
        precipitation = forcing%precipitation, &
        throughfall = cell_fluxes%throughfall, &
        canopy_evaporation = cell_fluxes%canopy_evaporation, &
        soil_evaporation = soil_evaporation, &
        infiltration = water_result%infiltration_rate, &
        surface_runoff = water_result%surface_runoff_rate, &
        drainage = water_result%drainage_rate, &
        ground_heat_flux = heat_flux_into_soil, &
        canopy_water_storage = sum(new_cell%tiles%fraction * new_cell%tiles%canopy_state%water_store), &
        soil_water_storage = soil_storage_new, &
        water_balance_error = forcing%precipitation * dt - canopy_storage_change &
            - (cell_fluxes%canopy_evaporation + soil_evaporation + water_result%surface_runoff_rate &
               + water_result%drainage_rate) * dt - (soil_storage_new - soil_storage_old), &
        solver_converged = water_result%converged)
  end function cell_step

  pure function land_step(cells, soil_library, tile_library, forcing, dt) result(new_cells)
    type(GridCell), intent(in) :: cells(:)
    type(SoilType), intent(in) :: soil_library(:)
    type(CanopyParameters), intent(in) :: tile_library(:)
    type(AtmosphericForcing), intent(in) :: forcing
    real(wp), intent(in) :: dt
    type(GridCell) :: new_cells(size(cells))

    integer :: cell_index
    do cell_index = 1, size(cells)
      new_cells(cell_index) = cell_step(cells(cell_index), soil_library(cells(cell_index)%soil_type_index), &
                                        tile_library, forcing, dt)
    end do
  end function land_step

end module land_mod
