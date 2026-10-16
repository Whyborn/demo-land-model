module constants_mod
  use kinds_mod, only: wp
  implicit none
  private

  real(wp), parameter, public :: water_density = 1000.0_wp              ! kg m-3
  real(wp), parameter, public :: water_specific_heat = 4186.0_wp        ! J kg-1 K-1
  real(wp), parameter, public :: latent_heat_vaporisation = 2.5e6_wp    ! J kg-1
  real(wp), parameter, public :: pi = 3.14159265358979323846_wp
  real(wp), parameter, public :: seconds_per_day = 86400.0_wp
end module constants_mod
