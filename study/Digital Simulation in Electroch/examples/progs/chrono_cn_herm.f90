program CHRONO_CN_HERM
! To simulate chronopotentiometry, outputting the C(0) and
! their errors at doubling intervals and a data file for plotting.
! Crank-Nicolson is used here, with equal intervals, and
! the Pearson start with M subintervals within the first step.
! 2(2) or 2(3) Hermitian G is applied (resp. herm = 2 or 3)
! but an npt GH is used to initialise C_0.
  use STUFF;  implicit none
  character(len=40) :: filename
  integer           :: NT, N, npt, herm, M, iT, next_out, Xlim=6
  real(kind=dbl)    :: dT, H, lambda, err, logerr, GH, C0FUNC
  real(kind=dbl), allocatable :: C(:), ad(:), bd(:)

  call FILSPC (4)
  read*, NT, lambda, npt, herm, M
  dT = 1.0_dbl / NT;  H = SQRT(dT/lambda);  N = CEILING(Xlim/H)
  ALLOCATE (C(0:N+1), ad(0:N), bd(0:N))
  GH = SQRT(pi) / 2 * H                   ! const G * H
  C = 1;   C(0) = C0FUNC (C, npt, GH)      ! C initialised, n-pt G
  print '(" CN chronopot. simulation with Hermitian 2(2) or 2(3).")'
  print '(" NT     =", i6)', NT
  print '(" Lambda =", f9.2)', lambda
  print '(" npt    =", i6)', npt
  print '(" N      =", i6, " pts along X.")', N
  print '(i3, " Pearson substeps within first step.")', M
  print '(" Hermitian G0, the 2(", i1, ") form.")', herm
  print '(/4x, "iT", 7x, "T", 7x, "C(0)", 4x, "log10(err)")'
  write (4,'("# Chronopot. simulation by CN.")') ! File Headers
  write (4,'("# NT, lambda, M =", i6, f8.3, i6)') NT,lambda,M
  write (4,'("#", 4x, "T", 8x, " C(0)", 4x, "logerr")')
  do iT = 1, M                             ! The Pearson steps
    call CNCHRONO (C, N, npt, herm, lambda/M, ad, bd, GH) ! A small CN step
  enddo
  next_out = 2
  do iT = 2, NT                              ! Grand T-loop
    call  CNCHRONO (C, N, npt, herm, lambda, ad, bd, GH)      ! A CN step
    err = C(0)-1+SQRT(iT*dT); logerr = LOG10(ABS(err))
    write (4,'(2f10.5, f8.3)') it*dT, C(0), logerr !Write out
    if (iT == next_out .OR. iT == NT) then
      print '(i6, 2f10.3, f9.2)', iT, iT*dT, C(0), logerr
      next_out = 2 * next_out
    endif
  enddo
  DEALLOCATE (C, ad, bd)
end program CHRONO_CN_HERM

subroutine CNCHRONO (C, N, npt, herm, lambda, ad, bd, GH)
! To solve the Crank-Nicolson system, by the backwards/forwards
! scheme. Boundary conditions for chronopotentiometry assumed.
  use STUFF; implicit none
  integer        ::  N, npt, herm
  real(kind=dbl) :: C(0:N), lambda, ad(0:N), bd(0:N), GH

  integer        :: i
  real(kind=dbl) :: a, b, a1, bi, c1, c2, c3, G, u(0:5), v(0:5), &
                    G0FUNC

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

! C'(0) using the 2(2) or 2(3) Hermitian approximation for G:
  if (herm == 2) then
    a = - (2*lambda + 1) / (2*lambda)
    ad(0) = 1 / a
    bd(0) = (GH - C(0)/(2*lambda)) / a
    C(0) = (bd(0) - ad(0)*bd(1)/ad(1)) / (1 - ad(0)/ad(1))
  else if (herm == 3) then
    a = - (3*lambda + 1) / (3*lambda)
    b = (6*lambda - 1) / (6*lambda)
    ad(0) = b / a
    bd(0) = (GH - C(0)/(3*lambda) - C(1)/(6*lambda)) / a
    C(0) = (bd(0) - ad(0)*bd(1)/ad(1)) / (1 - ad(0)/ad(1))
  else
    STOP 'Illegal Hermitian order.'
  endif

! Forward again, replacing all C with C' values:
  do i = 1, N
    C(i) = (bd(i) - c(i-1)) / ad(i)
  enddo
end subroutine CNCHRONO
