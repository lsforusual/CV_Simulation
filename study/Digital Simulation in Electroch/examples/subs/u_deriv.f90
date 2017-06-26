subroutine U_DERIV (u, x, i, n, which, deriv, coeffs)
! To calculate any u_x(i,n) (which=1) or u_xx(i,n) (which=2)
! and to return at the same time the coefficients for it.
! Because of the limited accuracy in the results for n > 12,
! this is the maximum allowed. This is an outdated routine,
! and should be replaced by FORN.

! U(1:N) is the array of function values;
! X(1:N) is the array of positions at which the U lie;
! I is the point to which the derivative holds;
! N is the length of the arrays U and X
! WHICH determines which derivative is to be computed (1 or 2);
! DERIV is the computed derivative;
! COEFFS(1:N) is the coefficient array that yields DERIV.
  use STUFF;   implicit none
  integer        :: i, n, which
  real(kind=dbl) :: u(1:n), x(1:n), deriv, coeffs(1:n)

  integer        :: j, k, factorial
  real(kind=dbl) :: gpow
  real(kind=dbl),allocatable :: g(:), H(:,:)

  if (n == 2) then
    if (which == 1) then
      deriv = (u(2) - u(1)) / (x(2) - x(1))
      coeffs(1) = -1 / (x(2) - x(1))
      coeffs(2) =  1 / (x(2) - x(1))
    else
      deriv = 0;  print *, " 2-pt 2nd deriv impossible."
    endif
  else if (n > 12) then
    print '(" > 12 pts not allowed in U_DERIV.")'
  else
    ALLOCATE (g(n-1), H(n-1,n-1))
!   Establishing the row of h's:
   j = 0
    do k = 1, n
      if (k /= i) then
        j = j + 1
        g(j) = x(k) - x(i)
      endif
    enddo
!   The factorial factor for later:
    factorial = 1
    do j = 2, which;  factorial = factorial * j;  enddo
!   Setting up the matrix to be inverted:
    do j = 1, n-1
      gpow = 1
      do k = 1, n-1
        gpow = gpow * g(j)
        H(j,k) = gpow
      enddo
    enddo
!   Inversion:
    call MATINV (H(1:n-1,1:n-1), n-1)
!   Row no "which" now generates the coeffs:
    j = 0;  coeffs = 0
    do k = 1, n
      if (k /= i) then
        j = j + 1
        coeffs(k) = H(which,j) * factorial
      endif
    enddo
    coeffs(i) = - SUM(coeffs(1:n)) ! (coeffs(i) was zero!)
!   Putting the derivative together:
    deriv = SUM (coeffs(1:n) * u(1:n))

    DEALLOCATE (g, H)
  endif
end subroutine U_DERIV
