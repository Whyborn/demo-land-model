! 1D heat diffusion, cell-centred finite volumes, implicit Euler.
!   top boundary:    prescribed heat flux into the soil (W m-2)
!   bottom boundary: zero flux
! Heat capacity and conductivity depend on water content but no phase change is modelled.
module soil_heat_mod
  use kinds_mod, only: wp
  use types_mod, only: SoilType
  use soil_properties_mod, only: thermal_conductivity, volumetric_heat_capacity
  use tridiagonal_mod, only: solve_tridiagonal
  implicit none
  private

  public :: soil_heat_step, column_heat_content

contains

  ! Total heat content relative to 0 K (J m-2).
  pure function column_heat_content(soil, layer_thickness, theta, temperature) result(heat_content)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: layer_thickness(:), theta(:), temperature(:)
    real(wp) :: heat_content
    heat_content = sum(volumetric_heat_capacity(soil, theta) * layer_thickness * temperature)
  end function column_heat_content

  pure function soil_heat_step(soil, layer_thickness, theta, temperature, ground_heat_flux, dt) &
      result(temperature_new)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: layer_thickness(:), theta(:), temperature(:)
    real(wp), intent(in) :: ground_heat_flux, dt
    real(wp) :: temperature_new(size(temperature))

    integer :: n_layers, layer
    real(wp), dimension(size(temperature)) :: conductivity, heat_capacity_per_area
    real(wp), dimension(size(temperature)) :: matrix_lower, matrix_diagonal, matrix_upper, rhs
    real(wp), dimension(size(temperature)) :: interface_conductance

    n_layers = size(temperature)
    conductivity = thermal_conductivity(soil, theta)
    heat_capacity_per_area = volumetric_heat_capacity(soil, theta) * layer_thickness

    interface_conductance = 0.0_wp
    do layer = 1, n_layers - 1
      interface_conductance(layer) = 1.0_wp / (0.5_wp * layer_thickness(layer) / conductivity(layer) &
                                               + 0.5_wp * layer_thickness(layer + 1) / conductivity(layer + 1))
    end do

    matrix_lower = 0.0_wp
    matrix_upper = 0.0_wp
    do layer = 1, n_layers
      matrix_diagonal(layer) = heat_capacity_per_area(layer) / dt
      rhs(layer) = heat_capacity_per_area(layer) / dt * temperature(layer)
      if (layer > 1) then
        matrix_diagonal(layer) = matrix_diagonal(layer) + interface_conductance(layer - 1)
        matrix_lower(layer) = -interface_conductance(layer - 1)
      else
        rhs(layer) = rhs(layer) + ground_heat_flux
      end if
      if (layer < n_layers) then
        matrix_diagonal(layer) = matrix_diagonal(layer) + interface_conductance(layer)
        matrix_upper(layer) = -interface_conductance(layer)
      end if
    end do

    temperature_new = solve_tridiagonal(matrix_lower, matrix_diagonal, matrix_upper, rhs)
  end function soil_heat_step

end module soil_heat_mod
