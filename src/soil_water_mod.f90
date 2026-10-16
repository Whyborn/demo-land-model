! 1D Richards equation (mixed form), cell-centred finite volumes, implicit Euler,
! modified Picard iteration (Celia et al. 1990). The prognostic variable is matric head.
!
! Depth increases downward; fluxes are positive downward.
!   top boundary:    flux = infiltration - soil evaporation (infiltration limited by capacity)
!   bottom boundary: free drainage (unit hydraulic gradient)
! If Picard iteration fails to converge, the step is recursively split in two halves.
module soil_water_mod
  use kinds_mod, only: wp
  use constants_mod, only: water_density
  use types_mod, only: SoilType
  use soil_properties_mod, only: water_content, hydraulic_conductivity, specific_moisture_capacity
  use tridiagonal_mod, only: solve_tridiagonal
  implicit none
  private

  public :: SoilWaterResult, soil_water_step, column_water_storage

  ! Rates are time-averaged over the step, in kg m-2 s-1 (positive downward / out of the column).
  type :: SoilWaterResult
    real(wp), allocatable :: matric_head(:)
    real(wp), allocatable :: theta(:)
    real(wp) :: infiltration_rate = 0.0_wp
    real(wp) :: surface_runoff_rate = 0.0_wp
    real(wp) :: drainage_rate = 0.0_wp
    logical :: converged = .false.
  end type SoilWaterResult

  integer, parameter :: max_picard_iterations = 40
  integer, parameter :: max_step_halvings = 8
  real(wp), parameter :: water_content_tolerance = 1.0e-10_wp
  ! Each Picard update may change the head by at most max(minimum, fraction * |head|). This keeps the
  ! iteration stable in very dry soil, where the capacity is tiny and the raw update is huge.
  real(wp), parameter :: minimum_head_change_limit = 1.0_wp
  real(wp), parameter :: relative_head_change_limit = 0.5_wp

contains

  ! Total water stored in the column (kg m-2).
  pure function column_water_storage(soil, layer_thickness, matric_head) result(storage)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: layer_thickness(:), matric_head(:)
    real(wp) :: storage
    storage = water_density * sum(layer_thickness * water_content(soil, matric_head))
  end function column_water_storage

  ! Advance the column by dt. Water supply and evaporation are rates in kg m-2 s-1.
  pure function soil_water_step(soil, layer_thickness, matric_head, water_supply_rate, &
                                evaporation_rate, dt) result(result_state)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: layer_thickness(:), matric_head(:)
    real(wp), intent(in) :: water_supply_rate, evaporation_rate, dt
    type(SoilWaterResult) :: result_state

    result_state = advance_with_halving(soil, layer_thickness, matric_head, water_supply_rate, &
                                        evaporation_rate, dt, 0)
  end function soil_water_step

  recursive pure function advance_with_halving(soil, layer_thickness, matric_head, water_supply_rate, &
                                               evaporation_rate, dt, halving_depth) result(result_state)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: layer_thickness(:), matric_head(:)
    real(wp), intent(in) :: water_supply_rate, evaporation_rate, dt
    integer, intent(in) :: halving_depth
    type(SoilWaterResult) :: result_state

    type(SoilWaterResult) :: first_half, second_half

    result_state = picard_attempt(soil, layer_thickness, matric_head, water_supply_rate, &
                                  evaporation_rate, dt)
    if (result_state%converged .or. halving_depth >= max_step_halvings) return

    first_half = advance_with_halving(soil, layer_thickness, matric_head, water_supply_rate, &
                                      evaporation_rate, 0.5_wp * dt, halving_depth + 1)
    second_half = advance_with_halving(soil, layer_thickness, first_half%matric_head, water_supply_rate, &
                                       evaporation_rate, 0.5_wp * dt, halving_depth + 1)
    result_state%matric_head = second_half%matric_head
    result_state%theta = second_half%theta
    result_state%infiltration_rate = 0.5_wp * (first_half%infiltration_rate + second_half%infiltration_rate)
    result_state%surface_runoff_rate = 0.5_wp * (first_half%surface_runoff_rate + second_half%surface_runoff_rate)
    result_state%drainage_rate = 0.5_wp * (first_half%drainage_rate + second_half%drainage_rate)
    result_state%converged = first_half%converged .and. second_half%converged
  end function advance_with_halving

  pure function picard_attempt(soil, layer_thickness, matric_head_old, water_supply_rate, &
                               evaporation_rate, dt) result(result_state)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: layer_thickness(:), matric_head_old(:)
    real(wp), intent(in) :: water_supply_rate, evaporation_rate, dt
    type(SoilWaterResult) :: result_state

    integer :: n_layers, layer, iteration
    real(wp), dimension(size(layer_thickness)) :: theta_old, matric_head_iterate, theta_iterate
    real(wp), dimension(size(layer_thickness)) :: capacity, conductivity, matrix_lower, matrix_diagonal
    real(wp), dimension(size(layer_thickness)) :: matrix_upper, rhs, matric_head_next
    real(wp), dimension(size(layer_thickness)) :: interface_conductivity, interface_conductance, head_change_limit
    real(wp) :: infiltration_capacity, infiltration, top_flux, storage_coefficient, theta_change
    real(wp) :: drainage_flux
    logical :: update_was_limited

    n_layers = size(layer_thickness)
    theta_old = water_content(soil, matric_head_old)
    matric_head_iterate = matric_head_old

    ! Infiltration is limited by the flux the surface soil can take at zero surface head.
    infiltration_capacity = max(soil%k_sat * (1.0_wp - matric_head_old(1) / (0.5_wp * layer_thickness(1))), &
                                0.0_wp)
    infiltration = min(water_supply_rate / water_density, infiltration_capacity)
    top_flux = infiltration - evaporation_rate / water_density
    drainage_flux = 0.0_wp
    result_state%converged = .false.

    do iteration = 1, max_picard_iterations
      theta_iterate = water_content(soil, matric_head_iterate)
      capacity = specific_moisture_capacity(soil, matric_head_iterate)
      conductivity = hydraulic_conductivity(soil, matric_head_iterate)

      interface_conductivity = 0.0_wp
      interface_conductance = 0.0_wp
      do layer = 1, n_layers - 1
        interface_conductivity(layer) = 0.5_wp * (conductivity(layer) + conductivity(layer + 1))
        interface_conductance(layer) = interface_conductivity(layer) &
                                       / (0.5_wp * (layer_thickness(layer) + layer_thickness(layer + 1)))
      end do

      matrix_lower = 0.0_wp
      matrix_upper = 0.0_wp
      do layer = 1, n_layers
        storage_coefficient = layer_thickness(layer) / dt * capacity(layer)
        matrix_diagonal(layer) = storage_coefficient
        rhs(layer) = storage_coefficient * matric_head_iterate(layer) &
                     - layer_thickness(layer) / dt * (theta_iterate(layer) - theta_old(layer))
        if (layer > 1) then
          matrix_diagonal(layer) = matrix_diagonal(layer) + interface_conductance(layer - 1)
          matrix_lower(layer) = -interface_conductance(layer - 1)
          rhs(layer) = rhs(layer) + interface_conductivity(layer - 1)
        else
          rhs(layer) = rhs(layer) + top_flux
        end if
        if (layer < n_layers) then
          matrix_diagonal(layer) = matrix_diagonal(layer) + interface_conductance(layer)
          matrix_upper(layer) = -interface_conductance(layer)
          rhs(layer) = rhs(layer) - interface_conductivity(layer)
        else
          rhs(layer) = rhs(layer) - conductivity(layer)
        end if
      end do
      drainage_flux = conductivity(n_layers)

      matric_head_next = solve_tridiagonal(matrix_lower, matrix_diagonal, matrix_upper, rhs)
      head_change_limit = max(minimum_head_change_limit, relative_head_change_limit * abs(matric_head_iterate))
      update_was_limited = any(abs(matric_head_next - matric_head_iterate) > head_change_limit)
      matric_head_next = matric_head_iterate &
                         + max(min(matric_head_next - matric_head_iterate, head_change_limit), -head_change_limit)
      theta_change = maxval(abs(water_content(soil, matric_head_next) - theta_iterate))
      matric_head_iterate = matric_head_next
      if (theta_change < water_content_tolerance .and. .not. update_was_limited) then
        result_state%converged = .true.
        exit
      end if
    end do

    result_state%matric_head = matric_head_iterate
    result_state%theta = water_content(soil, matric_head_iterate)
    result_state%infiltration_rate = infiltration * water_density
    result_state%surface_runoff_rate = water_supply_rate - infiltration * water_density
    result_state%drainage_rate = drainage_flux * water_density
  end function picard_attempt

end module soil_water_mod
