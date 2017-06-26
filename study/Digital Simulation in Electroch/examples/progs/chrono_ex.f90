program CHRONOEX
! To simulate chronopotentiometry, outputting the boundary
! concentration at doubling intervals and plot data.
  use STUFF;  implicit none
  character(len=40)          :: filename
  integer                    :: NT, N, iT, iX, next_out, Xlim=6
  real(kind=dbl)             :: dT, H, lambda, G, C1, C2, C3, &
                                err, logerr, C0FUNC
  real(kind=dbl),allocatable :: C(:)

  print '(" Output data file name?")'; read '(a40)', filename
  OPEN (unit=1, file=TRIM(filename), status='new')  ! Data file
  print '(" NT, lambda?")';   read*, NT, lambda
  dT = 1.0_dbl / NT;  H = SQRT(dT/lambda);  N = CEILING(Xlim/H)
  G = SQRT(pi) / 2                          ! constant G value
  ALLOCATE (C(0:N+1));  C = 1               ! C initialised
  print '(4x, "iT", 5x, "T", 5x, "C(0)", 2x, "C(0)(analyt)", 2x, "log10(err)")'
  write (1,'("# Cottrell simulation.")')     ! Header lines
  write (1,'("# NT, lambda =", i6, f8.3)') NT, lambda
  next_out = 1
  C(0) = C0FUNC(C, 7, G*H)                  ! 7-pt G approx.
  do iT = 1, NT                             ! Grand T-loop
    C1 = C(0);  C2 = C(1)                    ! Running scalars
    do iX = 1, N                             ! X-loop
      C3 = C(iX+1)                            ! Next scalar
      C(iX) = C2 + lambda * (C1 - 2*C2 + C3)  ! new C in array
      C1 = C2;  C2 = C3                       ! Scalars shifted
    enddo
    C(0) = C0FUNC(C, 7, G*H)                  ! 7-pt G approx.
    err = C(0)-1+SQRT(iT*dT); logerr = LOG10(ABS(err))
    write (1,'(2f10.5)') it*dT, C(0)         ! Write to file
    if (iT == next_out .OR. iT == NT) then
      print '(i6, 3f8.3, f14.2)', iT, iT*dT, C(0), 1-SQRT(iT*dT), logerr
      next_out = 2 * next_out
    endif
  enddo
  DEALLOCATE (C)
end program CHRONOEX


