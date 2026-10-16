program test_canopy
  use kinds_mod, only: wp
  use types_mod, only: CanopyParameters, CanopyState, CanopyStepResult, AtmosphericForcing
  use canopy_mod, only: canopy_step, default_tile_library
  use test_utilities_mod, only: assert_true, assert_close
  implicit none

  type(CanopyParameters) :: tile_library(3)

  tile_library = default_tile_library()
  call test_water_conservation_and_bounds(tile_library(1))
  call test_bare_soil_passes_everything(tile_library(3))

contains

  subroutine test_water_conservation_and_bounds(parameters)
    type(CanopyParameters), intent(in) :: parameters
    type(CanopyState) :: state
    type(CanopyStepResult) :: step_result
    type(AtmosphericForcing) :: forcing
    real(wp), parameter :: dt = 600.0_wp
    real(wp) :: worst_balance_error, balance_error, lowest_store, highest_store, lowest_throughfall
    integer :: step_index

    state = CanopyState(water_store=0.0_wp, temperature=285.0_wp)
    worst_balance_error = 0.0_wp
    lowest_store = 0.0_wp
    highest_store = 0.0_wp
    lowest_throughfall = 0.0_wp
    do step_index = 1, 400
      forcing%precipitation = merge(1.0e-3_wp, 0.0_wp, mod(step_index, 100) < 20)
      forcing%potential_evaporation = 4.0e-5_wp
      forcing%air_temperature = 288.0_wp
      forcing%net_radiation = 200.0_wp
      step_result = canopy_step(parameters, state, forcing, 286.0_wp, dt)
      balance_error = (forcing%precipitation - step_result%fluxes%throughfall &
                       - step_result%fluxes%canopy_evaporation) * dt &
                      - (step_result%state%water_store - state%water_store)
      worst_balance_error = max(worst_balance_error, abs(balance_error))
      lowest_store = min(lowest_store, step_result%state%water_store)
      highest_store = max(highest_store, step_result%state%water_store)
      lowest_throughfall = min(lowest_throughfall, step_result%fluxes%throughfall)
      state = step_result%state
    end do
    call assert_true(worst_balance_error < 1.0e-12_wp, 'canopy water is conserved each step')
    call assert_true(lowest_store >= 0.0_wp .and. highest_store <= parameters%max_water_storage + 1.0e-12_wp, &
                     'canopy store stays within [0, max storage]')
    call assert_true(lowest_throughfall >= 0.0_wp, 'throughfall is never negative')
    call assert_true(highest_store > 0.5_wp * parameters%max_water_storage, 'canopy store actually fills')
  end subroutine test_water_conservation_and_bounds

  subroutine test_bare_soil_passes_everything(parameters)
    type(CanopyParameters), intent(in) :: parameters
    type(CanopyStepResult) :: step_result
    type(AtmosphericForcing) :: forcing

    forcing = AtmosphericForcing(precipitation=1.0e-3_wp, air_temperature=290.0_wp, net_radiation=300.0_wp, &
                                 potential_evaporation=3.0e-5_wp)
    step_result = canopy_step(parameters, CanopyState(water_store=0.0_wp, temperature=280.0_wp), forcing, &
                              285.0_wp, 600.0_wp)
    call assert_close(step_result%fluxes%throughfall, forcing%precipitation, 1.0e-12_wp, 'bare soil: all rain reaches ground')
    call assert_close(step_result%state%temperature, forcing%air_temperature, 1.0e-12_wp, 'bare soil: surface air temperature')
    call assert_close(step_result%fluxes%soil_evaporation_demand, forcing%potential_evaporation, 1.0e-12_wp, &
                      'bare soil: full evaporative demand on soil')
  end subroutine test_bare_soil_passes_everything

end program test_canopy
