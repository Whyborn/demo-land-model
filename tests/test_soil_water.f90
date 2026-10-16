program test_soil_water
  use kinds_mod, only: wp
  use constants_mod, only: water_density
  use types_mod, only: SoilType
  use soil_properties_mod, only: water_content, matric_head_from_water_content, hydraulic_conductivity
  use soil_water_mod, only: SoilWaterResult, soil_water_step, column_water_storage
  use test_utilities_mod, only: assert_true, assert_close, example_soil_library
  implicit none

  type(SoilType) :: soil_library(3)

  soil_library = example_soil_library()
  call test_retention_inverse(soil_library(2))
  call test_mass_balance(soil_library(2))
  call test_steady_infiltration(soil_library(1))

contains

  subroutine test_retention_inverse(soil)
    type(SoilType), intent(in) :: soil
    real(wp) :: matric_head
    matric_head = -3.7_wp
    call assert_close(matric_head_from_water_content(soil, water_content(soil, matric_head)), matric_head, &
                      1.0e-8_wp, 'matric head inverts water content')
    call assert_close(hydraulic_conductivity(soil, 0.0_wp), soil%k_sat, 1.0e-12_wp, 'K(0) equals K_sat')
  end subroutine test_retention_inverse

  ! Over many steps with rain, evaporation and drainage, the storage change matches the net flux.
  subroutine test_mass_balance(soil)
    type(SoilType), intent(in) :: soil
    integer, parameter :: n_layers = 12
    real(wp) :: layer_thickness(n_layers), matric_head(n_layers), depth
    real(wp) :: storage_before, worst_error, step_error, supply_rate, evaporation_rate
    real(wp), parameter :: dt = 600.0_wp
    type(SoilWaterResult) :: step_result
    integer :: layer, step_index

    layer_thickness = 0.1_wp
    depth = 0.0_wp
    do layer = 1, n_layers
      depth = depth + layer_thickness(layer)
      matric_head(layer) = (depth - 0.5_wp * layer_thickness(layer)) - 1.0_wp
    end do

    worst_error = 0.0_wp
    do step_index = 1, 288
      supply_rate = merge(2.0e-3_wp, 0.0_wp, step_index > 40 .and. step_index < 70)   ! ~7 mm/h pulse
      evaporation_rate = 2.0e-5_wp
      storage_before = column_water_storage(soil, layer_thickness, matric_head)
      step_result = soil_water_step(soil, layer_thickness, matric_head, supply_rate, evaporation_rate, dt)
      matric_head = step_result%matric_head
      step_error = (column_water_storage(soil, layer_thickness, matric_head) - storage_before) &
                   - (step_result%infiltration_rate - evaporation_rate - step_result%drainage_rate) * dt
      worst_error = max(worst_error, abs(step_error))
    end do
    call assert_true(worst_error < 1.0e-6_wp, 'Richards step conserves mass (worst per-step error < 1e-6 kg m-2)')
  end subroutine test_mass_balance

  ! Constant supply below K_sat: drainage converges to the supply once the wetting front passes.
  subroutine test_steady_infiltration(soil)
    type(SoilType), intent(in) :: soil
    integer, parameter :: n_layers = 20
    real(wp) :: layer_thickness(n_layers), matric_head(n_layers), supply_rate
    real(wp), parameter :: dt = 300.0_wp
    type(SoilWaterResult) :: step_result
    integer :: step_index

    layer_thickness = 0.05_wp
    matric_head = -5.0_wp
    supply_rate = 0.3_wp * soil%k_sat * water_density
    do step_index = 1, 288    ! one day
      step_result = soil_water_step(soil, layer_thickness, matric_head, supply_rate, 0.0_wp, dt)
      matric_head = step_result%matric_head
    end do
    call assert_close(step_result%drainage_rate, supply_rate, 1.0e-2_wp, 'steady drainage equals supply')
    call assert_close(step_result%surface_runoff_rate, 0.0_wp, 1.0e-12_wp, 'no runoff below K_sat')
  end subroutine test_steady_infiltration

end program test_soil_water
