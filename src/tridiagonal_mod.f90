module tridiagonal_mod
  use kinds_mod, only: wp
  implicit none
  private

  public :: solve_tridiagonal

contains

  ! Thomas algorithm. lower(1) and upper(n) are never used.
  pure function solve_tridiagonal(lower, diagonal, upper, rhs) result(solution)
    real(wp), intent(in) :: lower(:), diagonal(:), upper(:), rhs(:)
    real(wp) :: solution(size(rhs))

    real(wp) :: scaled_upper(size(rhs)), scaled_rhs(size(rhs)), pivot
    integer :: n, row

    n = size(rhs)
    pivot = diagonal(1)
    scaled_upper(1) = upper(1) / pivot
    scaled_rhs(1) = rhs(1) / pivot
    do row = 2, n
      pivot = diagonal(row) - lower(row) * scaled_upper(row - 1)
      scaled_upper(row) = upper(row) / pivot
      scaled_rhs(row) = (rhs(row) - lower(row) * scaled_rhs(row - 1)) / pivot
    end do

    solution(n) = scaled_rhs(n)
    do row = n - 1, 1, -1
      solution(row) = scaled_rhs(row) - scaled_upper(row) * solution(row + 1)
    end do
  end function solve_tridiagonal

end module tridiagonal_mod
