subroutine FORN (M, N, x0, x, w)
! The Fortran version of Fornberg's algorithm. It returns only
! the coefficients w for a given approximation to a derivative
! of degree M (or for an interpolation if M = 0) at X0 on N points
! whose location is in X(0:N-1). Note that N is the actual number
! of points, but they are indexed here as 0:N-1.
  use STUFF;   implicit none
  integer        :: M, N
  real(kind=dbl) :: x0, x(0:N-1), w(0:N-1)

  integer        :: mm, nn, nu
  real(kind=dbl) :: c1, c2, c3
  real(kind=dbl),allocatable :: delta(:,:,:)

  ALLOCATE (delta(0:M,0:N-1,0:N-1))
  delta = 0;  delta (0,0,0) = 1
  c1 = 1
  do nn = 1, N-1
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
  w(0:N-1) = delta (M, N-1, 0:N-1)
  DEALLOCATE (delta)
end subroutine FORN
