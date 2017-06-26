subroutine JCOBI(ND,N,N0,N1, al,be, dif1,dif2,dif3, root)
! F-90 adaptation of the subroutine published in Villadsen and
! Michelsen (1978) "Solution of Differential Equation Models
! by Polynomial Approximation". This computes the roots of the
! shifted Jacobian polynomials, as used in orthogonal
! collocation.

! ND is the dimension of the arrays, and should here be N+1.
! N is the number of internal points.
! N0 is passed as zero to exclude the point at X = 0, else 1.
! N1 as for X0, but for X = 1
! al, be are polynomial coefficients, in the present context
!          set to zero.
! difx(..): xth derivative values at the collocation points.
! root(..): the roots, i.e. the collocation positions.

  use STUFF;  implicit none

  integer        :: ND, N, N0, N1
  real(kind=dbl) :: al,be, dif1(ND),dif2(ND),dif3(ND), root(ND)

  integer        :: i, j, nt
  real(kind=dbl) :: ab, ad, ap, x, xd, xd1, xp, xp1, xn, xn1, &
                    z, z1, zc, y

  ab = al + be;   ad = be - al;   ap = be * al
  dif1(1) = (ad/(ab+2) + 1) / 2;  dif2(1) = 0
  if (N >= 2) then
    do i = 2, N
      z1 = i - 1;  z = ab + 2*z1
      dif1(i) = (ab * ad / z / (z+2) + 1) / 2
      if (i == 2) then
        dif2(i) = (ab + ap + z1) / z**2 / (z+1)
        cycle
      endif
      z = z**2;  y = z1 * (ab + z1);  y = y * (ap + y)
      dif2(i) = y / z / (z-1)
    enddo
  endif
! Root determination by Newton's method with supression
!    of previously determined roots
  x = 0
  do i = 1, N
    do ! old label 25
      xd = 0; xn = 1; xd1 = 0; xn1 = 0
      do j = 1, N
        xp = (dif1(j)-x)*xn - dif2(j)*xd
        xp1 = (dif1(j)-x)*xn1 - dif2(j)*xd1 - xn
        xd = xn;  xd1 = xn1;  xn = xp;  xn1 = xp1
      enddo
      zc = 1;  z = xn / xn1
      if (i /= 1) then
        do j = 2, i
          zc = zc - z / (x-root(j-1))
        enddo
      endif
      z = z / zc;  x = x - z
      if (ABS(z) .lt. 1.0E-09_dbl) exit
    enddo
    root(i) = x
    x = x + 0.0001
  enddo
! Add possible interpolation points at X = 0 or X = 1
  NT = N + N0 + N1
  if (n /= 0) then
    do i = 1, N
      j = N + 1 - i
      root(j+1) = root(j)
    enddo
    root(1) = 0
  endif
  if (N1 == 1) root(NT)= 1
! Now evaluate derivatives of polynomial
  do i = 1, NT
    x = root(i)
    dif1(i) = 1;  dif2(i) = 0;  dif3(i) = 0
    do j = 1, NT
      if (j /= i) then
        y = x - root(j)
        dif3(i) = y * dif3(i) + 3*dif2(i)
        dif2(i) = y * dif2(i) + 2*dif1(i)
        dif1(i) = y * dif1(i)
      endif
    enddo
  enddo
end subroutine JCOBI
