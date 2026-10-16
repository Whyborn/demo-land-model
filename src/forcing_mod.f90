! Synthetic atmospheric forcing: diurnal cycles plus two rain events.
module forcing_mod
  use kinds_mod, only: wp
  use constants_mod, only: pi, seconds_per_day
  use types_mod, only: AtmosphericForcing
  implicit none
  private

  public :: synthetic_forcing

contains

  pure function synthetic_forcing(time_seconds) result(forcing)
    real(wp), intent(in) :: time_seconds
    type(AtmosphericForcing) :: forcing

    real(wp), parameter :: seconds_per_hour = 3600.0_wp
    real(wp) :: day_phase, daylight_factor

    day_phase = 2.0_wp * pi * (time_seconds - 6.0_wp * seconds_per_hour) / seconds_per_day
    daylight_factor = max(sin(day_phase), 0.0_wp)

    forcing%air_temperature = 288.0_wp + 6.0_wp * sin(day_phase - pi / 4.0_wp)    ! peaks at 15:00
    forcing%net_radiation = -50.0_wp + 450.0_wp * daylight_factor
    forcing%potential_evaporation = 6.0e-5_wp * daylight_factor

    forcing%precipitation = 0.0_wp
    if (time_seconds >= 30.0_wp * seconds_per_hour .and. time_seconds < 33.0_wp * seconds_per_hour) then
      forcing%precipitation = 10.0_wp / seconds_per_hour      ! 10 mm h-1 for 3 h
    else if (time_seconds >= 78.0_wp * seconds_per_hour .and. time_seconds < 80.0_wp * seconds_per_hour) then
      forcing%precipitation = 4.0_wp / seconds_per_hour       ! 4 mm h-1 for 2 h
    end if
  end function synthetic_forcing

end module forcing_mod
