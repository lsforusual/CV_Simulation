subroutine DAMPED_EXPANSION (H_1, N, L, gamma_1, x, ierr)
! To find that gamma_1 value which satisfies both H_1 and X(N) = L
! and to generate the resulting X values. The damping function
! has been chosen so as to be optimal but might be fiddled with.
! The damping is wrt X.
  use STUFF;   implicit none
  integer   :: N, ierr
  real(dbl) :: X(0:N), L, H_1, gamma_1

  integer   :: i
  real(dbl) :: dX, low, high

  ierr = 0
! 1. Find the right gamma_1 - it must lie between 1 and 2.
  low = 1;  high = 2
  call XN_FIND (N, L, H_1, high, X)
  if (X(N) > L)  then ! OK, gamma_1 lies between 1 and 2.
    do
      gamma_1 = (low + high) / 2
      call XN_FIND (N, L, H_1, gamma_1, X)
      if (X(N) < L) then
        low = gamma_1
      else
        high = gamma_1
      endif
      if (ABS(high - low) < small) exit
    enddo
  else ! Gamma_1 > 1, no good.
    print '(" Gamma_1 > 1, no good, XN =", e12.3)', X(N)
    ierr = 1
  endif

CONTAINS
subroutine XN_FIND (N, L, H_1, gamma_1, X)
! Finds X(N) for the given gamma_1 and at the same time generates
! the X sequence.
  integer   :: N
  real(dbl) :: L, H_1, gamma_1, X(0:N)

  integer   :: i
  integer,parameter :: damp_fac=6 ! Adjust if necessary
  real(dbl)         :: gamma, Hi, z, DAMPFUNC
  DAMPFUNC(z) = EXP(-damp_fac*z) ! One-line damping function

  x(0) = 0;  X(1) = H_1;    Hi = H_1
  do i = 2, N
      gamma = 1 + (gamma_1 - 1) * DAMPFUNC(X(i-1)/L)
      Hi = Hi * gamma
      X(i) = X(i-1)  +  Hi
  enddo
end subroutine XN_FIND
end subroutine DAMPED_EXPANSION
