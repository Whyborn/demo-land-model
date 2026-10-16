program test_land
  use kinds_mod, only: wp
  use types_mod, only: SoilType, CanopyParameters, GridCell, Tile, TileFluxes, AtmosphericForcing
  use land_mod, only: geometric_layer_thicknesses, build_grid_cell, aggregate_tile_fluxes, land_step
  use constants_mod, only: pi, seconds_per_day
  use test_utilities_mod, only: assert_true, assert_close, example_soil_library, example_tile_library
  implicit none

  call test_flux_aggregation()
  call test_multi_day_run()

contains

  subroutine test_flux_aggregation()
    type(Tile) :: tiles(2)
    type(TileFluxes) :: mean_fluxes

    tiles(1) = Tile(tile_type_index=1, fraction=0.25_wp)
    tiles(2) = Tile(tile_type_index=2, fraction=0.75_wp)
    tiles(1)%fluxes = TileFluxes(throughfall=4.0_wp, canopy_evaporation=1.0_wp, &
                                 soil_evaporation_demand=2.0_wp, ground_heat_flux=100.0_wp)
    tiles(2)%fluxes = TileFluxes(throughfall=8.0_wp, canopy_evaporation=3.0_wp, &
                                 soil_evaporation_demand=6.0_wp, ground_heat_flux=-20.0_wp)
    mean_fluxes = aggregate_tile_fluxes(tiles)
    call assert_close(mean_fluxes%throughfall, 7.0_wp, 1.0e-12_wp, 'aggregated throughfall is fraction-weighted')
    call assert_close(mean_fluxes%canopy_evaporation, 2.5_wp, 1.0e-12_wp, 'aggregated canopy evaporation')
    call assert_close(mean_fluxes%soil_evaporation_demand, 5.0_wp, 1.0e-12_wp, 'aggregated soil evaporation demand')
    call assert_close(mean_fluxes%ground_heat_flux, 10.0_wp, 1.0e-12_wp, 'aggregated ground heat flux')
  end subroutine test_flux_aggregation

  ! Three cells (different soils and tile mixes) over three days including two rain events.
  subroutine test_multi_day_run()
    type(SoilType) :: soil_library(3)
    type(CanopyParameters) :: tile_library(3)
    type(GridCell) :: cells(3)
    type(AtmosphericForcing) :: forcing
    real(wp), parameter :: dt = 1800.0_wp
    real(wp) :: layer_thickness(10), worst_balance_error, lowest_temperature, highest_temperature
    real(wp) :: lowest_wetness, highest_wetness
    logical :: all_converged
    integer :: step_index, cell_index

    soil_library = example_soil_library()
    tile_library = example_tile_library()
    layer_thickness = geometric_layer_thicknesses(10, 0.05_wp, 1.3_wp)
    cells(1) = build_grid_cell(1, soil_library(1), [0.6_wp, 0.3_wp, 0.1_wp], layer_thickness, 2.0_wp, 285.0_wp)
    cells(2) = build_grid_cell(2, soil_library(2), [0.2_wp, 0.7_wp, 0.1_wp], layer_thickness, 1.5_wp, 285.0_wp)
    cells(3) = build_grid_cell(3, soil_library(3), [0.0_wp, 0.0_wp, 1.0_wp], layer_thickness, 1.0_wp, 285.0_wp)

    worst_balance_error = 0.0_wp
    lowest_temperature = huge(1.0_wp)
    highest_temperature = -huge(1.0_wp)
    lowest_wetness = huge(1.0_wp)
    highest_wetness = -huge(1.0_wp)
    all_converged = .true.
    do step_index = 1, 144
      forcing = diurnal_forcing((step_index - 1) * dt)
      cells = land_step(cells, soil_library, tile_library, spread(forcing, 1, size(cells)), dt)
      do cell_index = 1, size(cells)
        worst_balance_error = max(worst_balance_error, abs(cells(cell_index)%diagnostics%water_balance_error))
        lowest_temperature = min(lowest_temperature, minval(cells(cell_index)%column%temperature))
        highest_temperature = max(highest_temperature, maxval(cells(cell_index)%column%temperature))
        lowest_wetness = min(lowest_wetness, &
            minval(cells(cell_index)%column%theta) - soil_library(cells(cell_index)%soil_type_index)%theta_res)
        highest_wetness = max(highest_wetness, &
            maxval(cells(cell_index)%column%theta) - soil_library(cells(cell_index)%soil_type_index)%theta_sat)
        all_converged = all_converged .and. cells(cell_index)%diagnostics%solver_converged
      end do
    end do

    call assert_true(worst_balance_error < 1.0e-6_wp, 'cell water budget closes (error < 1e-6 kg m-2 per step)')
    call assert_true(lowest_temperature > 260.0_wp .and. highest_temperature < 320.0_wp, 'soil temperatures stay physical')
    call assert_true(lowest_wetness > 0.0_wp .and. highest_wetness < 1.0e-3_wp, &
                     'water content between residual and saturation')
    call assert_true(all_converged, 'Richards solver converged on every step')
  end subroutine test_multi_day_run

  ! Diurnal cycle plus a 3 h rain event on day 2 and a 2 h event on day 4 (same for every cell).
  pure function diurnal_forcing(time_seconds) result(forcing)
    real(wp), intent(in) :: time_seconds
    type(AtmosphericForcing) :: forcing
    real(wp), parameter :: seconds_per_hour = 3600.0_wp
    real(wp) :: day_phase, daylight_factor

    day_phase = 2.0_wp * pi * (time_seconds - 6.0_wp * seconds_per_hour) / seconds_per_day
    daylight_factor = max(sin(day_phase), 0.0_wp)
    forcing%air_temperature = 288.0_wp + 6.0_wp * sin(day_phase - pi / 4.0_wp)
    forcing%net_radiation = -50.0_wp + 450.0_wp * daylight_factor
    forcing%potential_evaporation = 6.0e-5_wp * daylight_factor
    forcing%precipitation = 0.0_wp
    if (time_seconds >= 30.0_wp * seconds_per_hour .and. time_seconds < 33.0_wp * seconds_per_hour) then
      forcing%precipitation = 10.0_wp / seconds_per_hour
    else if (time_seconds >= 78.0_wp * seconds_per_hour .and. time_seconds < 80.0_wp * seconds_per_hour) then
      forcing%precipitation = 4.0_wp / seconds_per_hour
    end if
  end function diurnal_forcing

end program test_land
