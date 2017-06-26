program COTTEX
! To simulate the Cottrell experiment, outputting the current
! at doubling intervals and a data file for plotting.
  use STUFF;  implicit none
  character(len=40)          :: filename
  integer                    :: NT, N, iT, iX, next_out, Xlim=6
  real(kind=dbl)             :: dT, H, lambda, G, C1, C2, C3, Ganal, logerr
  real(kind=dbl),allocatable :: C(:)

  print '(" Output data file name?")'; read '(a40)', filename
  OPEN (unit=1, file=TRIM(filename), status='new')  ! Data file
  print '(" NT, lambda?")';   read*, NT, lambda
  dT = 1.0_dbl / NT;  H = SQRT(dT/lambda);  N = CEILING(Xlim/H)
  ALLOCATE (C(0:N+1));  C = 1;  C(0) = 0    ! C initialised
  print '(4x, "iT", 5x, "T", 7x, "G", 3x, "log10(error)")' !Header
  write (1,'("# Cottrell simulation.")')     ! Header lines
  write (1,'("# NT, lambda =", i6, f8.3)') NT, lambda
  next_out = 1
  do iT = 1, NT                               ! Grand T-loop
    C1 = C(0);  C2 = C(1)                     ! Running scalars
    do iX = 1, N                              ! X-loop
      C3 = C(iX+1)                            ! Next scalar
      C(iX) = C2 + lambda * (C1 - 2*C2 + C3)  ! new C in array
      C1 = C2;  C2 = C3                       ! Scalars shifted
    enddo
    G = (-3*C(0) + 4*C(1) - C(2)) / (2*H)     ! 3-point G-approx
    Ganal = 1 / SQRT(pi*iT*dT)
    logerr = LOG10(ABS((G-Ganal)/Ganal))
    write (1,'(3f10.5)') it*dT, G, Ganal      ! Write to file
    if (iT == next_out .OR. iT == NT) then
      print '(i6, 2f8.3, f9.2)', iT, iT*dT, G, logerr
      next_out = 2 * next_out
    endif
  enddo
  DEALLOCATE (C)
end program COTTEX
