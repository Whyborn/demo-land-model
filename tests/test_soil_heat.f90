program test_soil_heat
  use kinds_mod, only: wp
  use constants_mod, only: pi
  use types_mod, only: SoilType
  use soil_properties_mod, only: thermal_conductivity, volumetric_heat_capacity
  use soil_heat_mod, only: soil_heat_step, column_heat_content
  use test_utilities_mod, only: assert_close, example_soil_library
  implicit none

  type(SoilType) :: soil, soil_library(3)

  soil_library = example_soil_library()
  soil = soil_library(2)
  call test_energy_conservation(soil)
  call test_constant_flux_half_space(soil)

contains

  subroutine test_energy_conservation(soil)
    type(SoilType), intent(in) :: soil
    integer, parameter :: n_layers = 10
    real(wp) :: layer_thickness(n_layers), theta(n_layers), temperature(n_layers)
    real(wp) :: energy_before
    real(wp), parameter :: dt = 900.0_wp, heat_flux = 80.0_wp
    integer :: layer, step_index

    layer_thickness = 0.1_wp
    theta = 0.25_wp
    do layer = 1, n_layers
      temperature(layer) = 280.0_wp + real(layer, wp)
    end do

    energy_before = column_heat_content(soil, layer_thickness, theta, temperature)
    do step_index = 1, 50
      temperature = soil_heat_step(soil, layer_thickness, theta, temperature, 0.0_wp, dt)
    end do
    call assert_close(column_heat_content(soil, layer_thickness, theta, temperature), energy_before, &
                      1.0e-12_wp, 'zero-flux heat step conserves energy')

    do step_index = 1, 50
      temperature = soil_heat_step(soil, layer_thickness, theta, temperature, heat_flux, dt)
    end do
    call assert_close(column_heat_content(soil, layer_thickness, theta, temperature), &
                      energy_before + heat_flux * 50.0_wp * dt, 1.0e-12_wp, 'heat content rises by flux * time')
  end subroutine test_energy_conservation

  ! Semi-infinite solid with constant surface flux G:
  !   T(z,t) - T0 = (2 G / lambda) sqrt(kappa t) * ierfc(z / (2 sqrt(kappa t)))
  subroutine test_constant_flux_half_space(soil)
    type(SoilType), intent(in) :: soil
    integer, parameter :: n_layers = 100, probe_layer = 5
    real(wp), parameter :: layer_size = 0.02_wp, dt = 60.0_wp, heat_flux = 100.0_wp, initial_temperature = 285.0_wp
    real(wp) :: layer_thickness(n_layers), theta(n_layers), temperature(n_layers)
    real(wp) :: conductivity, diffusivity, elapsed, probe_depth, similarity, analytic_rise
    integer :: step_index, n_steps

    layer_thickness = layer_size
    theta = 0.2_wp
    temperature = initial_temperature
    n_steps = 1440
    do step_index = 1, n_steps
      temperature = soil_heat_step(soil, layer_thickness, theta, temperature, heat_flux, dt)
    end do

    conductivity = thermal_conductivity(soil, theta(1))
    diffusivity = conductivity / volumetric_heat_capacity(soil, theta(1))
    elapsed = n_steps * dt
    probe_depth = (real(probe_layer, wp) - 0.5_wp) * layer_size
    similarity = probe_depth / (2.0_wp * sqrt(diffusivity * elapsed))
    analytic_rise = 2.0_wp * heat_flux / conductivity * sqrt(diffusivity * elapsed) &
                    * (exp(-similarity**2) / sqrt(pi) - similarity * erfc(similarity))
    call assert_close(temperature(probe_layer) - initial_temperature, analytic_rise, 3.0e-2_wp, &
                      'heat diffusion matches half-space analytic solution')
  end subroutine test_constant_flux_half_space

end program test_soil_heat
