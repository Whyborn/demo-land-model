! Hydraulic (van Genuchten-Mualem) and thermal soil property functions.
! Above saturation (matric_head >= 0) water content grows linearly via specific storage,
! so the Richards solver conserves mass even for a saturated column.
module soil_properties_mod
  use kinds_mod, only: wp
  use constants_mod, only: water_density, water_specific_heat
  use types_mod, only: SoilType
  implicit none
  private

  public :: water_content, matric_head_from_water_content, hydraulic_conductivity
  public :: specific_moisture_capacity, thermal_conductivity, volumetric_heat_capacity
  public :: default_soil_library

  real(wp), parameter :: relative_conductivity_floor = 1.0e-12_wp

contains

  pure function van_genuchten_m(soil) result(m_exponent)
    type(SoilType), intent(in) :: soil
    real(wp) :: m_exponent
    m_exponent = 1.0_wp - 1.0_wp / soil%van_genuchten_n
  end function van_genuchten_m

  ! Effective saturation for a negative matric head.
  elemental function effective_saturation(soil, matric_head) result(saturation)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: matric_head
    real(wp) :: saturation
    saturation = (1.0_wp + (soil%van_genuchten_alpha * abs(matric_head))**soil%van_genuchten_n) &
                 **(-van_genuchten_m(soil))
  end function effective_saturation

  elemental function water_content(soil, matric_head) result(theta)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: matric_head
    real(wp) :: theta
    if (matric_head >= 0.0_wp) then
      theta = soil%theta_sat + soil%specific_storage * matric_head
    else
      theta = soil%theta_res + (soil%theta_sat - soil%theta_res) * effective_saturation(soil, matric_head)
    end if
  end function water_content

  ! Inverse of water_content.
  elemental function matric_head_from_water_content(soil, theta) result(matric_head)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: theta
    real(wp) :: matric_head
    real(wp) :: saturation, m_exponent
    if (theta >= soil%theta_sat) then
      matric_head = (theta - soil%theta_sat) / soil%specific_storage
    else
      m_exponent = van_genuchten_m(soil)
      saturation = max((theta - soil%theta_res) / (soil%theta_sat - soil%theta_res), 1.0e-12_wp)
      matric_head = -(saturation**(-1.0_wp / m_exponent) - 1.0_wp)**(1.0_wp / soil%van_genuchten_n) &
                    / soil%van_genuchten_alpha
    end if
  end function matric_head_from_water_content

  elemental function hydraulic_conductivity(soil, matric_head) result(conductivity)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: matric_head
    real(wp) :: conductivity
    real(wp) :: saturation, m_exponent
    if (matric_head >= 0.0_wp) then
      conductivity = soil%k_sat
    else
      m_exponent = van_genuchten_m(soil)
      saturation = effective_saturation(soil, matric_head)
      conductivity = soil%k_sat * sqrt(saturation) &
                     * (1.0_wp - (1.0_wp - saturation**(1.0_wp / m_exponent))**m_exponent)**2
      conductivity = max(conductivity, relative_conductivity_floor * soil%k_sat)
    end if
  end function hydraulic_conductivity

  ! d(theta)/d(matric_head); the specific storage applies when saturated.
  elemental function specific_moisture_capacity(soil, matric_head) result(capacity)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: matric_head
    real(wp) :: capacity
    real(wp) :: m_exponent, n_exponent, scaled_head_power
    if (matric_head >= 0.0_wp) then
      capacity = soil%specific_storage
    else
      m_exponent = van_genuchten_m(soil)
      n_exponent = soil%van_genuchten_n
      scaled_head_power = (soil%van_genuchten_alpha * abs(matric_head))**n_exponent
      capacity = (soil%theta_sat - soil%theta_res) * m_exponent * n_exponent * scaled_head_power &
                 / abs(matric_head) * (1.0_wp + scaled_head_power)**(-m_exponent - 1.0_wp)
    end if
  end function specific_moisture_capacity

  elemental function thermal_conductivity(soil, theta) result(conductivity)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: theta
    real(wp) :: conductivity
    real(wp) :: relative_saturation
    relative_saturation = min(max(theta / soil%theta_sat, 0.0_wp), 1.0_wp)
    conductivity = soil%conductivity_dry &
                   + (soil%conductivity_sat - soil%conductivity_dry) * sqrt(relative_saturation)
  end function thermal_conductivity

  elemental function volumetric_heat_capacity(soil, theta) result(heat_capacity)
    type(SoilType), intent(in) :: soil
    real(wp), intent(in) :: theta
    real(wp) :: heat_capacity
    heat_capacity = soil%heat_capacity_dry + theta * water_density * water_specific_heat
  end function volumetric_heat_capacity

  ! Three representative soils (approximate Carsel & Parrish 1988 values).
  pure function default_soil_library() result(library)
    type(SoilType) :: library(3)
    library(1) = SoilType(name='sand', theta_sat=0.43_wp, theta_res=0.045_wp, &
                          van_genuchten_alpha=14.5_wp, van_genuchten_n=2.68_wp, k_sat=8.25e-5_wp, &
                          specific_storage=1.0e-4_wp, heat_capacity_dry=1.3e6_wp, &
                          conductivity_dry=0.3_wp, conductivity_sat=2.2_wp)
    library(2) = SoilType(name='loam', theta_sat=0.43_wp, theta_res=0.078_wp, &
                          van_genuchten_alpha=3.6_wp, van_genuchten_n=1.56_wp, k_sat=2.89e-6_wp, &
                          specific_storage=1.0e-4_wp, heat_capacity_dry=1.2e6_wp, &
                          conductivity_dry=0.25_wp, conductivity_sat=1.6_wp)
    library(3) = SoilType(name='clay_loam', theta_sat=0.41_wp, theta_res=0.095_wp, &
                          van_genuchten_alpha=1.9_wp, van_genuchten_n=1.31_wp, k_sat=6.4e-7_wp, &
                          specific_storage=1.0e-4_wp, heat_capacity_dry=1.1e6_wp, &
                          conductivity_dry=0.2_wp, conductivity_sat=1.4_wp)
  end function default_soil_library

end module soil_properties_mod
