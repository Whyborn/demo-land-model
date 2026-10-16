module test_utilities_mod
  use kinds_mod, only: wp
  use types_mod, only: SoilType, CanopyParameters
  implicit none
  private

  public :: assert_true, assert_close, write_text_file, example_soil_library, example_tile_library
  public :: expect_failure_of_subprocess

contains

  subroutine assert_true(condition, description)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: description
    if (.not. condition) then
      write(*, '(2a)') 'FAIL: ', description
      error stop 1
    end if
    write(*, '(2a)') 'pass: ', description
  end subroutine assert_true

  ! Relative tolerance when expected is non-zero, absolute otherwise.
  subroutine assert_close(actual, expected, tolerance, description)
    real(wp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: description
    real(wp) :: scale
    scale = merge(abs(expected), 1.0_wp, expected /= 0.0_wp)
    if (abs(actual - expected) > tolerance * scale) then
      write(*, '(2a,2(a,es12.5))') 'FAIL: ', description, ' actual=', actual, ' expected=', expected
      error stop 1
    end if
    write(*, '(2a)') 'pass: ', description
  end subroutine assert_close

  ! Write text to a file; lines are separated by achar(10) in the text.
  subroutine write_text_file(file_path, text)
    character(len=*), intent(in) :: file_path, text
    integer :: unit_number
    open(newunit=unit_number, file=file_path, status='replace', access='stream', form='unformatted', action='write')
    write(unit_number) text
    close(unit_number)
  end subroutine write_text_file

  ! Run this same test executable again with one argument and check that it exits with an error.
  ! The test program is expected to dispatch on its first command-line argument.
  subroutine expect_failure_of_subprocess(case_name, description)
    character(len=*), intent(in) :: case_name, description
    character(len=512) :: executable_path
    integer :: exit_status, command_status
    call get_command_argument(0, executable_path)
    call execute_command_line(trim(executable_path) // ' ' // case_name // ' > /dev/null 2>&1', &
                              wait=.true., exitstat=exit_status, cmdstat=command_status)
    call assert_true(command_status == 0 .and. exit_status /= 0, description)
  end subroutine expect_failure_of_subprocess

  ! Small soil and tile sets for tests (IDs 1-3 mirror the example configuration).
  pure function example_soil_library() result(library)
    type(SoilType) :: library(3)
    library(1) = SoilType(id=1, name='sand', theta_sat=0.43_wp, theta_res=0.045_wp, &
                          van_genuchten_alpha=14.5_wp, van_genuchten_n=2.68_wp, k_sat=8.25e-5_wp, &
                          specific_storage=1.0e-4_wp, heat_capacity_dry=1.3e6_wp, &
                          conductivity_dry=0.3_wp, conductivity_sat=2.2_wp)
    library(2) = SoilType(id=2, name='loam', theta_sat=0.43_wp, theta_res=0.078_wp, &
                          van_genuchten_alpha=3.6_wp, van_genuchten_n=1.56_wp, k_sat=2.89e-6_wp, &
                          specific_storage=1.0e-4_wp, heat_capacity_dry=1.2e6_wp, &
                          conductivity_dry=0.25_wp, conductivity_sat=1.6_wp)
    library(3) = SoilType(id=3, name='clay_loam', theta_sat=0.41_wp, theta_res=0.095_wp, &
                          van_genuchten_alpha=1.9_wp, van_genuchten_n=1.31_wp, k_sat=6.4e-7_wp, &
                          specific_storage=1.0e-4_wp, heat_capacity_dry=1.1e6_wp, &
                          conductivity_dry=0.2_wp, conductivity_sat=1.4_wp)
  end function example_soil_library

  pure function example_tile_library() result(library)
    type(CanopyParameters) :: library(3)
    library(1) = CanopyParameters(id=1, name='broadleaf', leaf_area_index=5.0_wp, radiation_extinction=0.5_wp, &
                                  max_water_storage=1.0_wp, interception_fraction=0.7_wp, &
                                  heat_capacity=2.0e4_wp, canopy_air_exchange=20.0_wp, &
                                  canopy_ground_exchange=5.0_wp)
    library(2) = CanopyParameters(id=2, name='grass', leaf_area_index=2.0_wp, radiation_extinction=0.5_wp, &
                                  max_water_storage=0.4_wp, interception_fraction=0.5_wp, &
                                  heat_capacity=5.0e3_wp, canopy_air_exchange=15.0_wp, &
                                  canopy_ground_exchange=8.0_wp)
    library(3) = CanopyParameters(id=3, name='bare_soil', leaf_area_index=0.0_wp, radiation_extinction=0.5_wp, &
                                  max_water_storage=0.0_wp, interception_fraction=0.0_wp, &
                                  heat_capacity=0.0_wp, canopy_air_exchange=10.0_wp, &
                                  canopy_ground_exchange=10.0_wp)
  end function example_tile_library

end module test_utilities_mod
