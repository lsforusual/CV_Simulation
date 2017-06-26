program CHRONO_CN
! To simulate chronopotentiometry, outputting the C(0)
! at doubling intervals and a data file for plotting.
! Crank-Nicolson is used here, with equal intervals, and
! the Pearson start with M subintervals within the first step.
! The error is taken to be that in C(0,T).
! Variable npt for G is used.
  use STUFF;  implicit none
  character(len=40) :: filename
  integer           :: NT, N, npt, M, iT, next_out, Xlim=6
  real(kind=dbl)    :: dT, H, lambda, err, logerr, GH, C0FUNC
  real(kind=dbl), allocatable :: C(:), ad(:), bd(:)

  call FILSPC (4)
  read*, NT, lambda, npt, M
  dT = 1.0_dbl / NT;  H = SQRT(dT/lambda);  N = CEILING(Xlim/H)
  ALLOCATE (C(0:N+1), ad(N), bd(N))
  GH = SQRT(pi) / 2 * H                   ! const G * H
  C = 1;   C(0) = C0FUNC (C, npt, GH)      ! C initialised, n-pt G
  print '(" CN chronopot. simulation.")'  ! Data input summary:
  print '(" NT     =", i6)', NT
  print '(" Lambda =", f9.2)', lambda
  print '(" npt    =", i6)', npt
  print '(" N      =", i6, " pts along X.")', N
  print '(i3, " Pearson substeps within first step.")', M
  print '(/4x, "iT", 7x, "T", 7x, "C(0)", 4x, "log10(err)")'
  write (4,'("# Chronopot. simulation by CN.")') ! File Headers
  write (4,'("# NT, lambda, M =", i6, f8.3, i6)') NT,lambda,M
  write (4,'("#", 4x, "T", 8x, " C(0)", 4x, "logerr")')
  do iT = 1, M                             ! The Pearson steps
    call CNCHRONO (C, N, npt, lambda/M, ad, bd, GH) ! A small CN step
  enddo
  next_out = 2
  do iT = 2, NT                              ! Grand T-loop
    call  CNCHRONO (C, N, npt, lambda, ad, bd, GH)      ! A CN step
    err = C(0)-1+SQRT(iT*dT); logerr = LOG10(ABS(err))
    write (4,'(2f10.5, f8.3)') it*dT, C(0), logerr !Write out
    if (iT == next_out .OR. iT == NT) then
      print '(i6, 2f10.3, f9.2)', iT, iT*dT, C(0), logerr
      next_out = 2 * next_out
    endif
  enddo
  DEALLOCATE (C, ad, bd)
end program CHRONO_CN


subroutine CNCHRONO (C, N, npt, lambda, ad, bd, GH)
! To solve the Crank-Nicolson system, by the backwards/forwards
! scheme. Boundary conditions for chronopotentiometry assumed.
  use STUFF; implicit none
  integer        ::  N, npt
  real(kind=dbl) :: C(0:N), lambda, ad(N), bd(N), GH

  integer        :: i
  real(kind=dbl) :: a, a1, bi, c1, c2, c3, G, u(0:5), v(0:5), G0FUNC

  a  = - 2/lambda * (1 + lambda)
  a1 = - 2/lambda * (1 - lambda)

! Backwards from C(N), to generate all a' and b' values recursively:
! a'(N) = a(N),  and b'(N) = b(N) - C*  (C* is the bulk conc.).

  ad(N) = a
  bd(N) = -C(N-1) + a1*C(N) - 2*C(N+1)

  C3 = C(N)
  C2 = C(N-1)
  do  i = N-1, 1, -1
    C1 = C(i-1)
    bi = -C1 + a1*C2 - C3
    ad(i) = a - 1/ad(i+1)
    bd(i) = bi - bd(i+1)/ad(i+1)
    C3 = C2
    C2 = C1
  enddo

! C'(0) using the u-v device, and the 6-point G-approximation:
  u(0) = 0;  v(0) = 1
  do i = 1, 5     ! The recursive expressions for u and v:
    u(i) = (bd(i) - u(i-1)) / ad(i);  v(i) = -v(i-1) / ad(i)
  enddo
  C(0) = (GH - G0FUNC(u,npt,1.0_dbl)) / G0FUNC(v,npt,1.0_dbl)

! Forward again, replacing all C with C' values:

  do i = 1, N
    C(i) = (bd(i) - c(i-1)) / ad(i)
  enddo
end subroutine CNCHRONO
