module errors_mod
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private

  public :: fail_with_message

contains

  ! Report a fatal input or configuration error and stop with a non-zero exit status.
  subroutine fail_with_message(message)
    character(len=*), intent(in) :: message
    write(error_unit, '(2a)') 'ERROR: ', message
    error stop 1
  end subroutine fail_with_message

end module errors_mod
