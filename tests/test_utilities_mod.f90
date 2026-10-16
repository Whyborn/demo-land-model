module test_utilities_mod
  use kinds_mod, only: wp
  implicit none
  private

  public :: assert_true, assert_close

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

end module test_utilities_mod
