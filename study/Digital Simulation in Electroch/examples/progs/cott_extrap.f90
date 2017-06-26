program COTT_EXTRAP
! To simulate the Cottrell experiment, outputting the current
! at doubling intervals and a data file for plotting.
! 2nd-order BI-based extrapolation is used here, with unequal
! intervals. Three-point G-approximation and the direct four-
! point formula for the second spatial derivative are used.
! D is the temp. concentration array for the two half steps.
  use STUFF;  implicit none
  integer                    :: NT, N, iT, next_out, i
  real(kind=dbl)             :: dT, H1, Xlim, G, Ganal,T,    &
                                gamma, err, logerr,          &
                                EE_FAC, G0FORN
  real(kind=dbl),allocatable :: X(:), C(:), D(:),            &
                                a1(:), a2(:), a3(:), a4(:),  &
                                b1(:), b2(:), b3(:), b4(:),  &
                                ad(:), bd(:)

  call FILSPC (4)
  print '(" NT, X(1), N?")';   read*, NT, H1, N
  dT = 1.0_dbl / NT
  ALLOCATE (C(0:N+1), D(0:N+1), X(0:N+2), &
            a1(N), a2(N), a3(N), a4(N),   &
            b1(N), b2(N), b3(N), b4(N),   &
            ad(N), bd(N))
  C = 1; C(0) = 0                            ! C's initialised 
  x(0) = 0;   X(1) = H1;  Xlim = 6                  ! x-values:
  gamma = EE_FAC (X(1), N, Xlim)
  do i = 2, N+2;  X(i) = X(1) * (gamma**i-1) / (gamma-1); enddo
  print '(" BI/Extrap (4-pt 2nd der) Cottrell simulation.")'
  print '(" NT     =", i6)', NT
  print '(" X(1)   =", f9.2)', X(1)
  print '(" N      =", i6, " pts along X.")', N
  print '(" gamma  =", f12.5, " (found by iteration)")', gamma
  print '(/4x, "iT", 7x, "T", 7x, "Gsim", 5x, "log(err)")'
  write (4,'("# Cottrell sim. by extrap-2 (4-pt 2nd deriv).")')
  write (4,'("# NT, X(1), N =", i6, f8.4, i6)') NT, H1, N
  write (4,'("# gamma =",f8.5)') gamma  ! Last 2 header lines
  write (4,'("#",4x,"T",8x," Gsim",4x,"logerr")')
  call COEFFS (X, N, dT,   a1, a2, a3, a4)
  call COEFFS (X, N, dT/2, b1, b2, b3, b4)
  next_out = 1
  do iT = 1, NT                              ! Grand T-loop
    D = C    ! Copy the present concs into D.
    call STEP (C, N, a1, a2, a3, a4, ad, bd)   ! A whole step
    call STEP (D, N, b1, b2, b3, b4, ad, bd)   ! A half step
    call STEP (D, N, b1, b2, b3, b4, ad, bd)   ! ..and another.
    C = 2*D - C                     ! The extrapolation
    T = iT * dT
    G = G0FORN (C, X, 3);     Ganal = 1/SQRT(pi*T)
    err = (G-Ganal)/Ganal; logerr = LOG10(ABS(err))
    write (4,'(2f10.5, f8.3, f12.6)') T, G, logerr, G-Ganal
    if (iT == next_out .OR. iT == NT) then
      print '(i6, 3f10.3)', iT, T, G, logerr
      next_out = 2 * next_out
    endif
  enddo
  DEALLOCATE (C, D, a1, a3, b1, b3, ad, bd)
end program COTT_EXTRAP


subroutine COEFFS (X, N, dT, a1, a2, a3, a4)
! To precalculate the constants.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: X(0:N+1), K, a1(N), a2(N),  &
                    a3(N), a4(N), dT
  integer        :: i
  real(kind=dbl) :: alpha1, alpha2, alpha3, alpha4, &
                    h1, h2, h3, det

  do i = 1, N ! The four-point formulas:
    h1 = x(i)-x(i-1);  h2 = x(i+1)-x(i);  h3 = x(i+2)-x(i)
    det = (- h1*h2**2*h3**3 + h1*h2**3*h3**2 - h1**2*h2*h3**3 &
           + h1**2*h2**3*h3 - h1**3*h2*h3**2 + h1**3*h2**2*h3)&
          / 12
    alpha1 = (-h2*h3**3 + h2**3*h3) / (6*det)
    alpha2 = (  h2*h3**3 - h2**3*h3 + h1*h3**3 - h1**3*h3 &
              - h1*h2**3 + h1**3*h2) / (6*det)
    alpha3 = (-h1*h3**3 + h1**3*h3) / (6*det)
    alpha4 = (h1*h2**3 - h1**3 * h2) / (6*det)
    a1(i) =   (alpha2 - 1/dT) / alpha1
    a2(i) = alpha3 / alpha1
    a3(i) = alpha4 / alpha1
    a4(i) = -1 / dT / alpha1
  enddo
end subroutine COEFFS

subroutine STEP (C, N, a1, a2, a3, a4, ad, bd)
! To solve the BI system, by the backwards/forwards
! scheme. The Cottrell boundary condition is assumed,
! all C(0) = 0.
  use STUFF; implicit none
  integer        ::  N
  real(kind=dbl) :: C(0:N+1), a1(N), a2(N), a3(N), a4(N), &
                    ad(N), bd(N)
  integer        :: i
  real(kind=dbl) :: bi

! Backwards from C(N) to generate a' and b' values recursively:
  ad(N) = a1(N)
  bi = a4(N) * C(N); bd(N) = bi - (a2(N)+a3(N))*C(N+1)
  ad(N-1) = a1(N-1) - a2(N-1)/ad(N)
  bi = a4(N-1) * C(N-1)
  bd(N-1) = bi - a3(N-1)*C(N+1) - a2(N-1)*bd(N)/ad(N)
  do i = N-2, 1, -1
    ad(i) = a1(i) - (a2(i) - a3(i)/ad(i+2)) / ad(i+1)
    bi = a4(i)*C(i)
    bd(i) = bi - a2(i) * bd(i+1)/ad(i+1) &
               - a3(i) * (bd(i+2) - bd(i+1)/ad(i+1)) / ad(i+2)
  enddo

! Forward again, replacing all C with C' values:
  do i = 1, N
    C(i) = (bd(i) - C(i-1)) / ad(i)
  enddo
end subroutine STEP
