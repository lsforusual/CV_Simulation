program CV_EX
! EX CV simulation. A quasireversible reaction A + e --> B is
! assumed, taken as reversible if K0 > 1000.
! A data file is output of the current (G) vs potential p,
! outputting only when a change > 0.001 has occurred to minimise
! output. Peak and trough current and -potential are printed.
! Data input:
! 1. p_start, p_rev, nT/punit, lambda, K0
! 2. Current output file name.
  use STUFF;  implicit none
  integer                    :: N, nTper, nT, it, nout
  real(kind=dbl),parameter   :: xlim=6.0_dbl, punit=0.025691_dbl
  real(kind=dbl),allocatable :: CA(:), CB(:)
  real(kind=dbl)             :: lambda, K0, p_start, p_rev, p, &
                                Xlimit, H, dT, dp, g1, g2, g3, &
                                gmax, pmax, G_peak, p_peak,    &
                                G_trough, p_trough, gsave, G0FUNC
  print '(" p(start), p(reverse), nT/p, lambda, K0?")'
  read *, p_start, p_rev, nTper, lambda, K0
  call FILSPC (4)

  dT = 1.0_dbl / nTper;   H = SQRT(dT/lambda)
  Xlimit = Xlim * SQRT(ABS(p_start-p_rev)*2)
  N = Xlimit / H;   ALLOCATE (CA(0:N+1), CB(0:N+1))
  write (4,'("# EX CV simulation.")')   ! Plotting file header:
  write (4,'("# Start & reverse pot:", 2f8.2)') p_start, p_rev
  write (4,'("# nT/p_unit, lambda, K0:", i6, 2f8.2)') &
                nTper,     lambda, K0
  print '(" CV_EX")'
  print '(" Lambda         =", f10.3)', lambda
  print '(" H              =", f10.3)', H
  print '(" nT per p unit  =", i6)', nTper
  print '(" N              =", i6)', N
  print '(" P_start, P_rev =", 2f10.3)', p_start, p_rev 
  print '(" K0             =", f10.3)', K0
! Assume a scan from p_start, initially going negative to p_rev,
! and back again. Note that |dp| = dT.
  nT = 2 * nTper * NINT(ABS(p_start-p_rev))
  p = p_start;  dp = -dT            ! Start going negative
  CA = 1;  CB = 0                   ! Initial CA & CB at p_start
  g1 = 0;  g2 = 0; gsave = 0;  nout = 0
  write (4,'(2f10.5)') p_start,0.0 ! Starting point for plot
  do iT = 1, nT
    call ONESTEP (CA, CB, N, p+dp, H, lambda, K0)
    g3 = G0FUNC (CA, 5, H)              ! Current, 5-pt approx.
    if (-2.0 < p  .and. p < 2.0) then   ! Eliminate noise peaks
      if ((g2-g1)*(g3-g2) <= 0) then          ! A peak/trough ?
        call MINMAX (g1, g2, g3, gmax, pmax) ! Find it
        if (dp < 0) then
          g_peak = gmax; p_peak = p + dp*pmax
        else
          g_trough = gmax; p_trough = p + dp*pmax
        endif
      endif
    endif
    if (ABS(g3-gsave) > 0.001) then
      write (4,'(2f10.5)') p+dp, g3; gsave = g3; nout = nout + 1
    endif
    if (iT == nT/2) dp = -dp         ! Reverse scan starts now
    p = p + dp;  g1 = g2;  g2 = g3
  enddo
  print '(i6, " points written into plot file.")', nout
  print '(" Top current and -p =", 2f8.4)', g_peak, p_peak
  print '(" Bot current and -p =", 2f8.4)', g_trough, p_trough
end program CV_EX


subroutine ONESTEP (CA, CB, N, p, H, lambda, K0)
! Does a single step forward to p, calculating new concs.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: CA(0:N+1), CB(0:N+1), p, H, lambda, K0

  integer        :: iX
  real(kind=dbl) :: CA1, CA2, CA3, CB1, CB2, CB3, &
                    a11, a12, a21, a22, b1, b2, G0FUNC, G0BETA

! It is assumed that the boundary concs are correct.
  CA1 = CA(0);  CA2 = CA(1);  CB1 = CB(0); CB2 = CB(1)
  do ix = 1, N
    CA3 = CA(ix+1);  CB3 = CB(ix+1)
    CA(ix) = CA2 + lambda * (CA1 - 2*CA2 + CA3)
    CB(ix) = CB2 + lambda * (CB1 - 2*CB2 + CB3)
    CA1 = CA2;  CA2 = CA3;  CB1 = CB2;  CB2 = CB3
  enddo
  
! Now adjust the C(0)'s to conform with the new profiles:
  CA(0) = 0;  CB(0) = 0    ! Device to enable G0FUNC(1..n-1)
  if (K0 > 1000) then      ! Reversible case:
    CA(0) = -(G0FUNC(CA,5,1.0_dbl) + G0FUNC(CB,5,1.0_dbl)) &
            / (G0BETA(0,5)*(1 + EXP(-p)))
    CB(0) = CA(0) * EXP(-p)
  else
    a11 = K0 * H * EXP(-0.5*p) - G0BETA(0,5) ! B/V, alpha = 0.5
    a12 = - K0 * H * EXP(0.5*p)              ! Kb *H
    a21 = G0BETA (0,5)
    a22 = a21
    b1 = G0FUNC (CA,5,1.0_dbl)
    b2 = -b1 - G0FUNC(CB,5,1.0_dbl)
    CA(0) = (b1*a22-b2*a12) / (a11*a22-a12*a21)
    CB(0) = (b2 - a21*CA(0)) / a22
  endif
end subroutine ONESTEP
