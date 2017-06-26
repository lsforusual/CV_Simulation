program COTT_CN
! To simulate the Cottrell experiment, outputting the current
! at doubling intervals and a data file for plotting.
! Crank-Nicolson is used here, with equal intervals, and
! the Pearson start with M subintervals within the first step.
  use STUFF;  implicit none
  character(len=40)          :: filename
  integer                    :: NT, N, M, iT, next_out, Xlim=6
  real(kind=dbl)             :: dT, H, lambda, G, Ganal, &
                                err, logerr, G0FUNC
  real(kind=dbl),allocatable :: C(:), ad(:), bd(:)

  print '(" Output data file name?")'; read '(a40)', filename
  OPEN (unit=1, file=TRIM(filename), status='new')  ! Data file
  print '(" NT, lambda, M?")';   read*, NT, lambda, M
  dT = 1.0_dbl / NT;  H = SQRT(dT/lambda);  N = CEILING(Xlim/H)
  ALLOCATE (C(0:N+1), ad(N), bd(N))
  C = 1;  C(0) = 0    ! C initialised
  print '(" CN Cottrell simulation.")'  ! Data input summary:
  print '(" NT     =", i6)', NT
  print '(" Lambda =", f9.2)', lambda
  print '(" N      =", i6, " pts along X.")', N
  print '(i3, " Pearson substeps within first step.")', M
  print '(/4x, "iT", 7x, "T", 7x, "Gsim", 4x, "log10(err)")'
  write (1,'("# Cottrell simulation by CN.")')   ! File Headers
  write (1,'("# NT, lambda, M =", i6, f8.3, i6)') NT,lambda,M
  write (1,'("#",4x,"T",8x," Gsim",4x,"logerr")')
  do iT = 1, M                             ! The Pearson steps
    call CNCOTT (C, N, lambda/M, ad, bd)   ! A small CN step
  enddo
  next_out = 2
  do iT = 2, NT                              ! Grand T-loop
    call  CNCOTT (C, N, lambda, ad, bd)         ! A CN step
    G = G0FUNC(C,6,H);  Ganal = 1/SQRT(pi*iT*dT)! 6-pt approx
    err = (G-Ganal)/Ganal; logerr = LOG10(ABS(err))
    write (1,'(2f10.5, f8.3)') it*dT, G, logerr ! Write to file
    if (iT == next_out .OR. iT == NT) then
      print '(i6, 2f10.3, f10.2)', iT, iT*dT, G, logerr
      next_out = 2 * next_out
    endif
  enddo
  DEALLOCATE (C, ad, bd)
end program COTT_CN


subroutine CNCOTT (C, N, lambda, ad, bd)
! To solve the Crank-Nicolson system, by the backwards/forwards
! scheme. The Cottrell boundary condition is assumed, C(0) = 0.
  use STUFF; implicit none
  integer        ::  N
  real(kind=dbl) :: C(0:N), lambda, ad(N), bd(N)

  integer        :: i
  real(kind=dbl) :: a, a1, bi, c1, c2, c3

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

! Forward again, replacing all C with C' values:

  C(0) = 0    ! Cottrell boundary condition
  do i = 1, N
    C(i) = (bd(i) - c(i-1)) / ad(i)
  enddo
end subroutine CNCOTT
