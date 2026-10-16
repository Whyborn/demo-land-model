program test_tridiagonal
  use kinds_mod, only: wp
  use tridiagonal_mod, only: solve_tridiagonal
  use test_utilities_mod, only: assert_close
  implicit none

  real(wp) :: lower(4), diagonal(4), upper(4), rhs(4), solution(4)

  ! Matrix with 2 on the diagonal and -1 off-diagonal; solution [1,2,3,4] gives rhs [0,0,0,5].
  lower = [0.0_wp, -1.0_wp, -1.0_wp, -1.0_wp]
  diagonal = 2.0_wp
  upper = [-1.0_wp, -1.0_wp, -1.0_wp, 0.0_wp]
  rhs = [0.0_wp, 0.0_wp, 0.0_wp, 5.0_wp]

  solution = solve_tridiagonal(lower, diagonal, upper, rhs)
  call assert_close(solution(1), 1.0_wp, 1.0e-12_wp, 'tridiagonal solution(1)')
  call assert_close(solution(2), 2.0_wp, 1.0e-12_wp, 'tridiagonal solution(2)')
  call assert_close(solution(3), 3.0_wp, 1.0e-12_wp, 'tridiagonal solution(3)')
  call assert_close(solution(4), 4.0_wp, 1.0e-12_wp, 'tridiagonal solution(4)')
end program test_tridiagonal
