! Prescribed canopy interception of water and heat for one tile.
!   water: a fraction of rain over the canopy is caught in a bounded store; the overflow
!          falls through. The store evaporates in proportion to its wetness.
!   heat:  the canopy absorbs a Beer's-law fraction of net radiation and relaxes towards air
!          temperature; the remainder reaches the ground, together with a canopy-ground exchange.
! A tile with zero LAI and zero storage behaves as bare soil.
module canopy_mod
  use kinds_mod, only: wp
  use types_mod, only: CanopyParameters, CanopyState, TileFluxes, CanopyStepResult, AtmosphericForcing
  implicit none
  private

  public :: canopy_cover_fraction, canopy_step, default_tile_library

contains

  elemental function canopy_cover_fraction(parameters) result(cover_fraction)
    type(CanopyParameters), intent(in) :: parameters
    real(wp) :: cover_fraction
    cover_fraction = 1.0_wp - exp(-parameters%radiation_extinction * parameters%leaf_area_index)
  end function canopy_cover_fraction

  pure function canopy_step(parameters, state, forcing, soil_surface_temperature, dt) result(step_result)
    type(CanopyParameters), intent(in) :: parameters
    type(CanopyState), intent(in) :: state
    type(AtmosphericForcing), intent(in) :: forcing
    real(wp), intent(in) :: soil_surface_temperature, dt
    type(CanopyStepResult) :: step_result

    real(wp) :: cover_fraction, precipitation_mass, intercepted_mass, store_after_rain
    real(wp) :: wet_fraction, evaporation_mass, absorbed_radiation, transmitted_radiation
    real(wp) :: canopy_temperature_new

    cover_fraction = canopy_cover_fraction(parameters)

    ! Water
    precipitation_mass = forcing%precipitation * dt
    intercepted_mass = min(parameters%interception_fraction * cover_fraction * precipitation_mass, &
                           max(parameters%max_water_storage - state%water_store, 0.0_wp))
    store_after_rain = state%water_store + intercepted_mass
    if (parameters%max_water_storage > 0.0_wp) then
      wet_fraction = store_after_rain / parameters%max_water_storage
    else
      wet_fraction = 0.0_wp
    end if
    evaporation_mass = min(store_after_rain, forcing%potential_evaporation * cover_fraction * wet_fraction * dt)

    step_result%state%water_store = store_after_rain - evaporation_mass
    step_result%fluxes%throughfall = (precipitation_mass - intercepted_mass) / dt
    step_result%fluxes%canopy_evaporation = evaporation_mass / dt
    step_result%fluxes%soil_evaporation_demand = forcing%potential_evaporation * (1.0_wp - cover_fraction)

    ! Heat (implicit relaxation of canopy temperature towards air temperature)
    absorbed_radiation = cover_fraction * forcing%net_radiation
    transmitted_radiation = forcing%net_radiation - absorbed_radiation
    canopy_temperature_new = (parameters%heat_capacity * state%temperature &
                              + dt * (absorbed_radiation + parameters%canopy_air_exchange * forcing%air_temperature)) &
                             / (parameters%heat_capacity + dt * parameters%canopy_air_exchange)
    step_result%state%temperature = canopy_temperature_new
    step_result%fluxes%ground_heat_flux = transmitted_radiation &
        + parameters%canopy_ground_exchange * (canopy_temperature_new - soil_surface_temperature)
  end function canopy_step

  ! Broadleaf forest, grass and bare soil.
  pure function default_tile_library() result(library)
    type(CanopyParameters) :: library(3)
    library(1) = CanopyParameters(name='broadleaf', leaf_area_index=5.0_wp, radiation_extinction=0.5_wp, &
                                  max_water_storage=1.0_wp, interception_fraction=0.7_wp, &
                                  heat_capacity=2.0e4_wp, canopy_air_exchange=20.0_wp, &
                                  canopy_ground_exchange=5.0_wp)
    library(2) = CanopyParameters(name='grass', leaf_area_index=2.0_wp, radiation_extinction=0.5_wp, &
                                  max_water_storage=0.4_wp, interception_fraction=0.5_wp, &
                                  heat_capacity=5.0e3_wp, canopy_air_exchange=15.0_wp, &
                                  canopy_ground_exchange=8.0_wp)
    library(3) = CanopyParameters(name='bare_soil', leaf_area_index=0.0_wp, radiation_extinction=0.5_wp, &
                                  max_water_storage=0.0_wp, interception_fraction=0.0_wp, &
                                  heat_capacity=0.0_wp, canopy_air_exchange=10.0_wp, &
                                  canopy_ground_exchange=10.0_wp)
  end function default_tile_library

end module canopy_mod
