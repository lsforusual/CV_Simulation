program UMDE_VB
! UMDE simulation using Verbrugge/Baker, with m-point approximations
! across theta rows and up Gamma columns, same for both directions.
! Nth sets the number of intervals along theta, Nga those along
! gamma. Index convention: i for gamma, j for theta, both
! starting at 0. Theta runs from 0 to theta_max, while gamma
! runs from 0 to Gamma_max, to be calculated from Tmax, the maximum T to
! which the simulation is to run.
! The array of concentrations is mapped onto a large matrix, using KMAP.
! Three-point BDF is used, started with a single BI step.
! C_1 is the previous concentration array at T-dT, C the
! current one. As soon as there is a solution, C replaces C_1,
! the new one replaces C. The whole grid, including the points
! along the electrode and the bulk points, is included as
! "unknowns", but they are not all included in the derivative coeffs,
! except the first derivative wrt Gamma at Gamma = 0, for the current.
! The Harwell suite MA28 solves the sparse system at each step.
! For the first BI, and the first BDF step, MA28 provides an LU
! decomposition; thereafter only back-substitution is needed.
! All the various coefficients are wrapped up in the two w's. These
! are two-dimensional, the first index indicating the index i or j at which
! the set applies, and the second the number of values in that set.
! The windows' start and end indices are i1, i2 along gamma, and j1, j2
! along theta. F is the function value at (i,j).
! dT is left out of the coefficients, and must be inserted in A_SET.
  use STUFF;   implicit none
  integer           :: Nth, Nga, N, m, nT, fac
  real(dbl)         :: dT, Xlim

  call FILSPC (4)
  print '(" Nth, Nga, m?")';           read *, Nth, Nga, m
  print '(" dT, nT?")';                read *, dT, nT
  print '(" Xlim?")';                  read *, Xlim
  print '(" Factor for LICN, LRCN?")'; read *, fac
! The factor is for MA28. It might need increasing beyond the MA28
! recommendations. Values as high as 100 were needed, increasing for
! higher N. We hope for a modest N*N however, and thus a smaller factor.
  print '(" m-pt sparse BI/BDF3 UMDE simul., VB grid.")'
  print '(" Nth, Nga                =", 2i10)', Nth, Nga
  print '(" m, derivative windows   =", i10)', m
  print '(" dT, nT, Tmax            =", f10.3, i10, f10.2)', dT, nT, nT*dT
  print '(" Xlim                    =", f12.1)', Xlim
  print '(" MA28 Workspace factor   =", i10)', fac

  call RUN (Nth, Nga, m, dT, nT, Xlim, fac)
end program UMDE_VB

subroutine RUN (Nth, Nga, m, dT, nT, Xlim, MAfac)
! The weights are declared as (0:Nth,nTh+1) and (0:NGa,nGa+1) which
! is needed for OC; if not OC, only (:,1:m) are actually used.
  use STUFF;    implicit none
  character(len=1) :: distr
  logical          :: output
  integer          :: Nth, Nga, N, m, nT, MAfac
  real(dbl)        :: dT, Xlim

  integer          :: iT, i, j, inext
  real(dbl)        :: dth, dga, theta_max, Gamma_max, T, Tmax, &
                      curr, L, curr_exact, err
  real(dbl),allocatable :: C_1(:,:), C(:,:), theta(:), Gamma(:), &
                           x(:), wth(:,:), wga(:,:), F(:,:)
! MA28-specific declarations:
  integer,parameter     :: solve_direct = 1
  integer               :: NZ28, fac, licn, lirn, iflag
  integer, allocatable  :: irn(:), icn(:), ikeep(:), iw(:)
  real(dbl)             :: u=1.0_dbl
  real(dbl),allocatable :: A(:), W(:), rhs(:)

  fac = MAfac
  ALLOCATE (C_1(0:Nga,0:Nth), C(0:Nga,0:Nth), theta(0:Nth), Gamma(0:NGa), &
            F(Nga-1,Nth-1), wth(0:Nth,m), wga(0:Nga,m))
  Tmax = nT * dT;   L = Xlim * SQRT(Tmax)
  Gamma_max = ACOSH(1+L)
  Gamma_max = Gamma_max / (1 + Gamma_max) ! VB-type Gamma
  theta_max = pi / 2                      ! VB theta limit
  dth = theta_max / Nth;  dga = Gamma_max / Nga
  C = 1;  C(0,:) = 0;  C_1 = C            ! Initialisation
! The grid:
  do j = 0, Nth
    theta(j) = j * dth
  enddo
  do i = 0, NGa
    Gamma(i) = i * dGa
  enddo
  call COEFFS (theta, Gamma, Nth, Nga, m, dT, wth, wga, F)
  N = (Nth+1) * (Nga+1) ! No. of unknowns, grid points.
  print '(" Large NxN matrix, N     =", i10)', N
  print '(" Gamma_max               =", f13.6)', Gamma_max
  call A_PRELIM (NZ28, N, m, Nth, Nga)         ! Count the nonzeroes.
  print '(" NZ28 nonzero elements:   ", i10)', NZ28
  lirn = fac*NZ28;  licn = fac*NZ28            ! Working space for MA28.
  ALLOCATE (A(lirn), icn(licn), irn(lirn), RHS(N), &
                          ikeep(5*N), iw(8*N), W(N))
  write (4,'(/"# MINI_UMDE, Nth, Nga, m, m =", 4i6)') Nth,Nga,m,m
  print '(5x, "iT", 8x, "T", 9x, "Curr", 4x, "%rel.err")'
  write(4,'("#", 6x, "T", 8x,"%rel.err")')

  inext = nT / 10
  do iT = 1, NT
    if (iT < 3) then              ! iT is either 1 (BI) or > 1 (BDF)
      call A_SET (A, NZ28, icn, irn, N, m, Nth, Nga, dth, dga, &
                  F, wth, wga, iT)
      call MA28AD (N, NZ28, A, licn, irn, lirn, icn, u, ikeep,  &
                                 iw, w, iflag) ! LUD for BI and BDF
      if (iflag /= 0) STOP "Something wrong in MA28AD. Abort 1."
    endif
    call B_SET (RHS, N, Nth, Nga, C_1, C, F, iT) ! RHS, for BI or BDF

    call MA28CD (N, A, licn, icn, ikeep, RHS, w, solve_direct)
    call SHIFTEM (C_1, C, RHS, Nth, Nga, N)      ! Shift concs down.
    call CURRENT (C, theta, Nth, Nga, m, wga(0,1:m), curr)
    T = iT*dT
    call UMDE_ERROR (T, curr, curr_exact, err)
    write(4,'(3f12.6)') T, 100*err/curr_exact
    if (iT == inext .OR. iT == NT) then
      print '(i8, f12.6, f10.5, f9.4)', iT, T, curr, 100*err/curr_exact
      inext = inext + nT / 10
    endif
  enddo
  DEALLOCATE (C_1, C, theta, Gamma, F, wth, wga, irn, icn, ikeep, iw, W, RHS)
end subroutine RUN

subroutine COEFFS (theta, Gamma, Nth, Nga, m, dT, wth, wga, F)
! Precomputes the coeffs for the discretisations.
  use STUFF;   implicit none
  integer          :: Nth, Nga, m
  real(dbl)        :: theta(0:Nth), Gamma(0:NGa), dT, wth(0:Nth,Nth+1), &
                      wga(0:Nga,NGa+1), F(Nga-1,Nth-1)

  integer   :: i, j, k, i1, i2, j1, j2
  real(dbl),allocatable :: al(:), be(:)
  real(dbl) :: b_theta, a_Gamma, b_Gamma

  ALLOCATE (al(m), be(m))
! The F's, a 2D array, incorporating dT:
  do i = 1, Nga-1;  do j = 1, Nth-1
    F(i,j) = SIN(theta(j))**2 + SINH(Gamma(i)/(1-Gamma(i)))**2
  enddo;  enddo
  F = F / dT
! Theta weights:
! LH edge,  dC/dth = 0, only the betas
  call FORN (1, m, theta(0), theta(0:m-1), wth(0,1:m))
! Bulk points
  do j = 1, Nth-1
    call I1I2 (j, Nth, m, j1, j2)
    call FORN (2, m, theta(j), theta(j1:j2), al(1:m))
    call FORN (1, m, theta(j), theta(j1:j2), be(1:m))
    b_theta = - TAN(theta(j))
    wth(j,1:m) = al(1:m)  +  b_theta * be(1:m)
  enddo
!  RH edge, only betas for dC/dth = 0
  call I1I2 (Nth, Nth, m, j1, j2)
  call FORN (1, m, theta(Nth), theta(j1:j2), wth(Nth,1:m))

! Gamma weights
! Bottom, i = 0, beta's only in wGa for current dC(0)/dGa:
  call FORN (1, m, Gamma(0), Gamma(0:m-1), wga(0,1:m))
! The bulk points, full wGa's:
  do i = 1, NGa-1
    call I1I2 (i, Nga, m, i1, i2)
    call FORN (2, m, Gamma(i), Gamma(i1:i2), al(1:m))
    call FORN (1, m, Gamma(i), Gamma(i1:i2), be(1:m))
    a_gamma = (1 - Gamma(i))**4
    b_Gamma = (1-Gamma(i))**2 * TANH(Gamma(i)/(1-Gamma(i))) &
              - 2*(1-Gamma(i))**3
    wga(i,1:m) = a_Gamma * al(1:m)  +  b_Gamma * be(1:m)
  enddo
  DEALLOCATE (al, be)
end subroutine COEFFS

subroutine A_PRELIM (NZ28, N, m, Nth, Nga)
! Preliminary run through A, to find no. of unknowns.
  use STUFF;     implicit none
  integer :: N, m, NZ28, Nth, Nga

  integer :: i, j, k, l, i1, i2, j1, j2, NZ

  NZ = 0 ! Count of nonzero elements to be.
! 1. The base line, electrode, here all C(0,j) = 0:
  do j = 0, Nth
    NZ = NZ + 1
  enddo
! 2. The diffusion region, a row at a time:
  do i = 1, Nga-1
    call I1I2 (i, Nga, m, i1, i2)
!   First the LHS edge, j=0: here we have dC/dth = 0, forward diff., wth only
    do k = 0, m-1
      NZ = NZ + 1
    enddo
!   Now the bulk across theta
    do j = 1, Nth-1
!     The diagonal point (i,j) first:
      NZ = NZ + 1
!     Now the others in the stretch of m points, leaving out (i,j)
      call I1I2 (j, Nth, m, j1, j2)
      do k = j1, j2
        if (k == j) cycle
        NZ = NZ + 1
      enddo
!     The points up the column stretch, leaving out (i,j):
      do k = i1, i2
        if (k == i) cycle
        NZ = NZ + 1
      enddo
    enddo
!   Last, the RHS edge, dC/dth = 0, backwards diff.
    call I1I2 (Nth, Nth, m, j1, j2)
    do k = j1, j2
      NZ = NZ + 1
    enddo
  enddo
!  3. Top edge, all C = 1
  do j = 0, Nth
    NZ = NZ + 1
  enddo
  NZ28 = NZ
end subroutine A_PRELIM

subroutine A_SET (A, NZ28, icn, irn, N, m, Nth, Nga, dth, dga, &
                                              F, wth, wga, iT)
! Sets the large matrix A, using the passed coeffs.
! If iT = 1, A is set for BI, otherwise for 3-pt BDF.
  use STUFF;    implicit none
  integer   :: N, NZ28, m, irn(NZ28), icn(NZ28), Nth, Nga, iT
  real(dbl) :: A(NZ28), dth, dGa, F(Nga-1,Nth-1), &
                      wth(0:Nth,m), wga(0:NGa,m)

  integer   :: i, j, k, i1, i2, j1, j2, iw, jw, krow, kcol, NZ, KMAP

! Positions krow and kcol are the mapped indices.
! k is a loop counter.
! i1, i2 and j1, j2 are the indices for the ends of the discretisation
! array stretches.
! iw and jw are indices along a given disretisation stretch.
! NZ is the count of non-zeroes.

  A = 0; NZ = 0 ! Count of nonzero elements to be.
! 1. The base line, electrode, here all C(0,j) = 0:
  do j = 0, Nth
    krow = KMAP (0, j, Nth)
    NZ = NZ + 1;  irn(NZ) = krow;  icn(NZ) = krow;  A(NZ) = 1
  enddo
 
! 2. The diffusion region, a row at a time:
  do i = 1, Nga-1
    call I1I2 (i, Nga, m, i1, i2)
!   First the LHS edge, j=0: where dC/dth = 0, forward difference
    krow = KMAP (i, 0, Nth)
    do k = 0, m-1
      NZ = NZ + 1;  kcol = KMAP (i, k, Nth);  jw = k + 1
      irn(NZ) = krow;  icn(NZ) = kcol;  A(NZ) = wth(0,jw)
    enddo
!   Now the bulk across theta
    do j = 1, Nth-1
      call I1I2 (j, Nth, m, j1, j2)
      krow = KMAP (i, j, Nth)
!     The diagonal point (i,j) first:
      NZ = NZ + 1;   irn(NZ) = krow;  icn(NZ) = krow
      iw = i - i1 + 1;  jw = j - j1 + 1
      A(NZ) = wth(j,jw) + wga(i,iw) - F(i,j) ! BI case
      if (iT > 1) A(NZ) = A(NZ) - F(i,j)/2 ! Must be BDF
!     Now the others in the stretch of m points, leaving out (i,j)
      do k = j1, j2
        if (k == j) cycle
        kcol = KMAP (i, k, Nth);  jw = k - j1 + 1
        NZ = NZ + 1;   irn(NZ) = krow;  icn(NZ) = kcol;  A(NZ) = wth(j,jw)
      enddo
!     The points up the column stretch, leaving out (i,j):
      do k = i1, i2
        if (k == i) cycle
        kcol = KMAP (k, j, Nth);  iw = k - i1 + 1
        NZ = NZ + 1;   irn(NZ) = krow;  icn(NZ) = kcol;  A(NZ) = wGa(i,iw)
      enddo
    enddo
!   Last, the RHS edge, where again dC/dth = 0, backwards difference
    krow = KMAP (i, Nth, Nth)
    call I1I2 (Nth, Nth, m, j1, j2)
    do k = j1, j2
      NZ = NZ + 1;  kcol = KMAP (i, k, Nth);  jw = k - j1 + 1
      irn(NZ) = krow;  icn(NZ) = kcol;  A(NZ) = wth(Nth,jw)
    enddo
  enddo

!  3. Top edge, all C = 1
  do j = 0, Nth
    krow = KMAP (NGa, j, Nth);  kcol = KMAP (NGa, j, Nth)
    NZ = NZ + 1;  irn(NZ) = krow;  icn(NZ) = kcol;  A(NZ) = 1
  enddo
! Check the number of unknowns against those found in A_PRELIM:
  if (NZ /= NZ28) STOP ! Problem in A_SET with no. nonzeroes"
end subroutine A_SET


subroutine B_SET (B, N, Nth, Nga, C_1, C, F, iT)
! Sets the RHS vector B, according to whether this is BI (iT = 1)
! or BDF (iT>1). C_{bulk} = 1 is assumed. The sequence through the
! grid is the same as for SET_A.
  use STUFF;    implicit none
  integer   :: Nth, Nga, N, iT
  real(dbl) :: C_1(0:Nga,0:Nth), C(0:Nga,0:Nth), B(N), F(Nga-1,Nth-1)

  integer        :: i, j, k, KMAP

  B = 0
! 1. The disk, bottom edge, all = 0
  do j = 0, Nth
    k = KMAP(0,j,Nth);  B(k) = 0
  enddo
! 2. The bulk region, diffusion equation:
  do i = 1, Nga-1
!   First the LHS edge: here we have dC/dth = 0, forward diff.
    k = KMAP (i, 0, Nth)
    B(k) = 0
!   Now the bulk across theta
    do j = 1, Nth-1
      k = KMAP (i, j, Nth)
      B(k) = -F(i,j) * C(i,j)
      if (iT > 1) B(k) = F(i,j)/2 * (C_1(i,j)-4*C(i,j))
    enddo
!   Last, the RHS edge, dC/dth = 0, backwards diff.
    k = KMAP (i, Nth, Nth)
    B(k) = 0
  enddo
!  3. Top edge, all C = 1
  do j = 0, Nth
    k = KMAP (Nga, j, Nth)
    B(k) = 1
  enddo
end subroutine B_SET


subroutine I1I2 (i, N, m, i1, i2)
! Computes i1 and i2 for the stretch of m points along N, so that
! they center i where possible, but make an asymmetric pair where not.
  implicit none
  integer :: i, N, m, i1, i2

  integer :: mm

  mm = m / 2
  i1 = MAX(0, i-mm);  i2 = i1 + m - 1 ! Bulk case
  if (i2 > N) then;  i2 = N;  i1 = N - m + 1;   endif
end subroutine I1I2


subroutine SHIFTEM (C_1, C, B, Nth, Nga, N)
! To shift the new concs C down to past C_1, and fold B into C.
  use STUFF;   implicit none

  integer        :: Nth, Nga, N
  real(kind=dbl) :: C_1(0:Nga,0:Nth), C(0:Nga,0:Nth), B(N)

  integer        :: i, j, k
  C_1 = C
  do k = 1, N
    call KUNMAP (k, i, j, nth)
    C(i,j) = B(k)
  enddo
end subroutine SHIFTEM


integer function KMAP (i, j, Nth)
! Maps (i,j) in space to k'th row of the big matrix.
  implicit none
  integer :: i, j, Nth

  if (j > Nth) STOP "Bad j value in KMAP. Abort."
  KMAP = j + 1 + i*(Nth+1)
end function KMAP

subroutine KUNMAP (k, i, j, Nth)
! Unmaps k back to (i,j) in space.
  implicit none
  integer :: k, i, j, Nth

  i = (k-1) / (Nth+1)
  j = MOD(k,Nth+1) - 1;   if(j < 0) j = Nth

end subroutine KUNMAP


subroutine CURRENT (C, theta, Nth, Nga, m, wga, curr)
! Calculates the current, by Simpson integration. The derivatives
! are approximated using the precomputed m-point coeffs.
! The Simpson-like function U_SIMP is used, even though the theta
! intervals are equal. This allows the use of even or odd NTh.
  use STUFF;     implicit none
  integer   :: Nth, Nga, m
  real(dbl) :: C(0:Nga,0:Nth), theta(0:Nth), wga(m), curr

  integer               :: i, j
  real(dbl),allocatable :: u(:)
  real(dbl)             :: deriv, U_SIMP

  ALLOCATE (u(0:Nth))

  do j = 0, Nth
    deriv = SUM (wga(1:m) * C(0:m-1,j))
    u(j) = COS(theta(j)) * deriv * pi / 2
  enddo
  curr = U_SIMP (theta(0:Nth), u(0:Nth), Nth)
  DEALLOCATE (u)
end subroutine CURRENT
