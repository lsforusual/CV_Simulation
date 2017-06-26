subroutine FORNBERG (degree, x0, x, u, N, derivative, w)
! Returns a derivative at any point X0 given a number N of positions X
! and values U at these X-points, using the Fornberg algorithm,
! B. Fornberg, Math. Comp. 51 (1988) 699-706. The degree can be
! 0, 1, 2, .., where zero means an interpolation. X0 need not be among
! the X-positions, but should preferably be. X0 may coincide with any
! x-position.
! This routine calls FORNBERG_SUB, and returns the weights which
! produce the derivatives, as well as the actual derivative.
! The arrays are received here numbered 1:N but in FORNBERG_SUB, they
! are assumed indexed as 0:N-1. This is the main reason for this front-end
! to the actual algorithm.
  use STUFF;   implicit none
  integer        :: degree, N
  real(kind=dbl) :: x0, x(N), u(N), derivative, w(N)

  call FORNBERG_SUB (degree, N-1, x0, X(1:N), w(1:N))
  derivative = SUM (w(1:N) * u(1:N))
end subroutine FORNBERG

subroutine FORNBERG_SUB (M, N, x0, x, w)
! The Fortran version of Fornberg's algorithm.
  use STUFF;   implicit none
  integer        :: M, N
  real(kind=dbl) :: x0, x(0:N), w(0:N)

  integer        :: mm, nn, nu
  real(kind=dbl) :: c1, c2, c3
  real(kind=dbl),allocatable :: delta(:,:,:)

  ALLOCATE (delta(0:M,0:N,0:N))
  delta = 0;  delta (0,0,0) = 1
  c1 = 1
  do nn = 1, N
    c2 = 1
    do nu = 0, nn-1
      c3 = x(nn) - x(nu)
      c2 = c2 * c3
      if (nn <= M) delta(nn, nn-1,nu) = 0
      delta (0,nn,nu) = (x(nn)-x0) * delta(0,nn-1,nu) / c3 ! mm = 0
      do mm = 1, MIN(nn,M)
        delta(mm,nn,nu) = ((x(nn)-x0) * delta(mm,nn-1,nu)     &
                            - mm * delta(mm-1,nn-1,nu))   / c3
      enddo
    enddo
    delta(0,nn,nn) =  - c1/c2 * (x(nn-1)-x0) * delta(0,nn-1,nn-1) ! mm = 0
    do mm = 1, MIN(nn,M)
      delta(mm,nn,nn) = c1/c2 * (mm * delta(mm-1,nn-1,nn-1)           &
                                 - (x(nn-1)-x0) * delta(mm,nn-1,nn-1))
    enddo
    c1 = c2
  enddo
  w(0:N) = delta (M, N, 0:N)
  DEALLOCATE (delta)
end subroutine FORNBERG_SUB
