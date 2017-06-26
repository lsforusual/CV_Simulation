function C0FORN (u, x, n, G)
! Computes u(0) for points u(0:n-1) at positions x(0:n-1), arbitrarily
! spaced, given the gradient G at x = 0.
  use STUFF;   implicit none
  integer   :: n
  real(dbl) :: C0FORN, u(0:n-1), x(0:n-1), G

  real(dbl),allocatable :: c(:)
  real(dbl) :: sum, b(15), deriv, G0FORN ! b are the weights

  ALLOCATE (C(0:n-1)) 
  C = u
  C(0) = 0 ! To trick G0FORN to compute the sum excluding at X = 0.
  sum = G0FORN (C, X, n)
  call FORNBERG (1, x(0), x, u, n, deriv, b)
  C0FORN = (G - sum) / b(1)
end function C0FORN
