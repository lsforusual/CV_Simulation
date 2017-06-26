program UMDE_DIRECT
! UME simulation using direct m-point discretisation on an unequal grid,
! with exponentially expanding intervals. BI/BDF(3-pt).
! nA is the number of points R(0:nA) where R(nA) = 1, the disk edge.
! nR is the total number of points in R, R(nR) being in the bulk.
! Thus there are nR-nA+1 points R(nA+1:nR) between the edge and bulk.
! nZ is the number of intervals along Z, Z(0:nZ) (pt NZ in the bulk).
! Index convention: i for Z, j for R.
! The array of concs is mapped onto a large matrix, using index k.
! Along the insulating plane, an m-pt approximation to dC/dZ = 0 is used
! as boundary condition.
! C_1 is the previous concs at T-dT, C the present one. 
! Three-point BDF starting with BI is used.
! The whole grid, including the outer points to nR+1, nZ+1, are taken
! as "unknowns". The package routines out of MA28 are used to solve for
! each step.
! The current and its error are written out into a file for plotting.
  use STUFF;   implicit none
  integer           :: nZ, nA, nR, m, nT, N, fac
  real(kind=dbl)    :: dR1, dZ1, dT, Xlim

  call FILSPC (4)
  read *, nZ, nA, nR
  read *, m
  read *, dZ1, dR1
  read *, dT, nT
  read *, Xlim
  read *, fac   ! Factor for MA28

  print '(" m-pt sparse BI/BDF UMDE sim., unequal intervals,", &
          & " direct on (R,Z).")'
  print '(" nZ, nA, nR           =", 3i10)', nZ, nA, nR
  print '(" m                    =", i10)', m
  print '(" Large NxN matrix, N  =", i10)', (nZ+1)*(nR+1)
  print '(" dZ1, dR1             =", 2es12.1)',  dZ1, dR1
  print '(" Xlim                 =", f10.3)', Xlim
  Xlim = Xlim * SQRT(dT*NT)
  print '(" Zmax, Rmax           =", 2f10.3)', Xlim, 1+Xlim
  print '(" dT, nT               =", f10.3, i10)', dT, nT
  print '(" Tmax                 =", f10.3)', dT*NT
  print '(" MA28 factor          =", i10)', fac
  call RUN (nZ, nA, nR, m, dZ1, dR1, Xlim, dT, NT, fac)
end program UMDE_DIRECT

subroutine RUN (nZ, nA, nR, m, dZ1, dR1, Xlim, dT, NT, fac)
  use STUFF;    implicit none
  logical   :: output
  integer   :: nZ, nA, nR, m, nT, fac
  real      :: t1, t2, t3
  real(dbl) :: dR1, dZ1, dT, T, Xlim

  integer               :: N, iT, out_interval
  real(dbl),allocatable :: Z(:), R(:), C_1(:,:), C(:,:),  &
                           alZ(:,:), alR(:,:), beR(:), beZ(:)
  real(dbl)             :: gamR1, gamR2, gamZ, curr, curr_exact, err
! MA28-specific declarations:
  integer               :: Nelem, NZ28, ncol, nrow, ier
  integer, allocatable  :: irow(:), jcol(:), ikeep(:), iw(:)
  real(dbl)             :: u=1.0_dbl
  real(dbl),allocatable :: A(:), W(:), RHS(:)

  N = (nR+1) * (nZ+1)
  ALLOCATE (C_1(0:nZ,0:nR), C(0:nZ,0:nR), Z(0:nZ), &
            R(0:nR), alZ(m,nZ), alR(m,nR), beR(0:m-1), beZ(0:m-1))
  C = 1;  C(0,0:nA) = 0;  C_1 = C
  call GRID (Z, R, dZ1, dR1, Xlim, nZ, nA, nR, gamZ, gamR1, gamR2)
  call COEFFS (Z, R, nZ, nR, Xlim, m, alZ, alR, beR, beZ)
  NZ28 = (2*m-1) * N ! Conservative initial guess, A_SET will give exact value.
  ncol = fac * NZ28;  nrow = fac * NZ28
  print '(" Expansion gamZ       =",  f10.3)',  gamZ
  print '(" Expansion gamR1, -2  =", 2f10.3)',  gamR1, gamR2
  print '(" Max lambda in Z-dir  =", e12.3)', dT/dZ1**2
  print '(" Max lambda in R-dir  =", e12.3)', dT/dR1**2
  print '(/" NZ28 starts as  :", i10)', NZ28
  ALLOCATE (A(ncol), jcol(ncol), irow(nrow), rhs(N), ikeep(5*N), iw(8*N), W(N))
  call A_SET_PRELIM (Nelem, NZ, nA, nR, m)
  print '(" NZ28 is actually:", i10)', Nelem
  write (4,'(/"# MINI_DIRECT, NZ, NA, NR, m =", 4i6)') NZ, NA, NR, m
  print '(6x, "iT", 8x, "T", 7x, "curr", 5x, "%err")'
  write(4,'("#", 6x, "T", 6x, "curr", 4x, "%rel.err")')

  out_interval = NT / 10
  do iT = 1, NT
    T = iT*dT
    if (iT <= 2) then             ! First and second steps only
      call A_SET (A, NZ28, jcol,irow, N, nZ,nA,nR, m, alZ,alR,beR,beZ, dT, iT)
      if (NZ28 /= Nelem) STOP " Problem with NZ28"
      call MA28AD (N, NZ28, A, ncol,irow,nrow,jcol, u, ikeep, iw, w, ier)!LUD
      if (ier /= 0) STOP "Problem in MA28AD. Abort."
    endif
    call B_SET (RHS, N, nZ, nA, nR, m, C_1, C, dT, iT)  ! Setting RHS vector
    call MA28CD (N, A, ncol, jcol, ikeep, rhs, w, 1)    ! Solve
    call SHIFTEM (C_1, C, rhs, nZ, nR, N) ! Shift concs down, unfold B.
    call CURRENT (C, nZ, nA, nR, m, Z, R, curr)
    call UMDE_ERROR (T, curr, curr_exact, err)
    write(4,'(3e14.6)') T, curr, 100*err/curr_exact
    if (MOD(iT,out_interval) == 0 .or. iT==nT) then
      print '(i8, f12.6, 2(f10.4, 2f8.3))', &
              iT, T, curr, 100*err/curr_exact
    endif
  enddo
end subroutine RUN

subroutine GRID (Z, R, dZ1, dR1, Xlim, nZ, nA, nR, gamZ, gamR1, gamR2)
! To compute the grid, i.e the R and Z values. The Z values expand
! exponentially from Z = 0 upwards, whereas the R values expand exp.
! away from the R = 1 point at R(nA) in both directions. This gives
! two different gamma values for the expansion along R.
  use STUFF;  implicit none
  integer        :: nZ, nA, nR
  real(kind=dbl) :: Z(0:nZ), R(0:nR), dZ1, dR1, Xlim, gamZ, gamR1, gamR2

  integer        :: i, j
  real(kind=dbl) :: EE_FAC

  gamZ = EE_FAC (dZ1, nZ, Xlim)   ! The Z values, simple
  Z(0) = 0
  do i = 1, nZ
    Z(i) = dZ1 * (gamZ**i - 1) / (gamZ - 1)
  enddo

  gamR1 = EE_FAC (dR1, nA, 1.0_dbl)  ! The R's, over the disk, index 0:nA
  R(0) = 0; R(nA) = 1                ! axis and edge
  do j = 1, nA -1                    ! Backwards, in from the edge:
    R(nA-j) = 1 - dR1*(gamR1**j - 1) / (gamR1 - 1)
  enddo
  gamR2 = EE_FAC (dr1, nR-nA, Xlim)  ! Outside the disk on plane
  do j = 1, nR-nA                    ! Over the insulating plane:
    R(nA+j) = 1 + dR1*(gamR2**j-1)/(gamR2-1)
  enddo
end subroutine GRID

subroutine COEFFS (Z, R, nZ, nR, Xlim, m, alZ, alR, beR, beZ)
! Precomputes the coeffs for the discretisations.
! The alR are composites for the 2nd and 1st derivatives
! combined (the first divided by R).
  use STUFF;    implicit none
  integer        :: nZ, nR, m
  real(kind=dbl) :: Z(0:nZ), R(0:nR), Xlim, alZ(m,nZ), alR(m,nR), &
                    beR(0:m-1), beZ(0:m-1)

  integer        :: i, j, k, i1, i2, j1, j2
  real(kind=dbl),allocatable :: al(:), be(:)

  ALLOCATE (al(m), be(m))
  alZ = 0;  alR = 0
! The coeffs along Z; very simple, FORN does it directly:
  call FORN (1, m, Z(0), Z(0:m-1), beZ(0:m-1)) ! 1st deriv for bc
  do i = 1, nZ-1
    call I1I2 (i, NZ, m, i1, i2)
    call FORN (2, m, Z(i), Z(i1:i2), alZ(1:m,i))
  enddo
! The coeffs along R:
  call FORN (1, m, R(0), R(0:m-1), beR(0:m-1)) ! 1st deriv for axis "bc"
  do j = 1, nR-1
    call I1I2 (j, NR, m, j1, j2)
    call FORN (1, m, R(j), R(j1:j2), be(1:m))
    call FORN (2, m, R(j), R(j1:j2), al(1:m))
    alR(1:m,j) = al(1:m) + be(1:m)/R(j)
  enddo
  DEALLOCATE (al, be)
end subroutine COEFFS

subroutine A_SET_PRELIM (Nelem, NZ, nA, nR, m)
! A prelim. run-through setting A, only to count its elements, Nelem.
  use STUFF;     implicit none
  integer :: Nelem, NZ, nA, nR, m

  integer :: k, i, j, i1, i2, j1, j2, ii, jj

  k = 0
! The disk itself, all C = 0, but must be included:
  do j = 0, nA;     k = k + 1;   enddo
! The insulating plane:
  do j = nA+1, nR-1
    do i = 0, m-1;       k = k + 1;     enddo
  enddo
! The internal points that vary
  do i = 1, nZ-1
    call I1I2(i,NZ,m,i1,i2)
    do j = 0, m-1          ! The axis approximation
      k = k + 1
    enddo
    do j = 1, nR-1
      call I1I2(j,NR,m,j1,j2)
      k = k + 1 ! Diagonal
!     Now the others in the stretch of m points along R, leaving out (i,j)
      do jj = j1, j2
        if (jj == j) cycle
        k = k + 1
      enddo
!     ...and the others in the stretch of m points along Z, leaving out (i,j)
      do ii = i1, i2
        if (ii == i) cycle
        k = k + 1
      enddo
    enddo
  enddo
! The top row, all = 1, entry 1 in matrix:
  do j = 0, nR;    k = k + 1;  enddo
! The right-most column, same thing:
  do i = 0, nZ-1;    k = k + 1;  enddo
  Nelem = k
end subroutine A_SET_PRELIM

subroutine A_SET (A, NZ28, jcol, irow, N, nZ, nA, nR, m, alZ, alR, &
                  beR, beZ, dT, iT)
! Sets the large matrix A, using the passed coeffs. This depends on iT;
! if it equals 1 it is for BI, if > 1, for BDF.
  use STUFF;     implicit none
  integer        :: N, NZ28, irow(NZ28), jcol(NZ28), nZ, nA, nR, m, iT
  real(kind=dbl) :: A(NZ28), alZ(1:m,nZ),alR(1:m,nR), beR(0:m-1), &
                    beZ(0:m-1), dT

  integer        :: i,ii, j,jj, k, i1,i2, j1,j2, iw,jw, krow,kcol, KMAP

  A = 0; k = 0
! The disk itself, all C = 0, but must be included:
  do j = 0, nA
    krow = KMAP(0,j,nR)
    k = k + 1;  irow(k) = krow;  jcol(k) = krow;  A(k) = 1
  enddo
! The insulating plane:
  do j = nA+1, nR-1
    krow = KMAP(0,j,nR)
    do i = 0, m-1
      k = k + 1;  kcol = KMAP(i,j,NR)
      irow(k) = krow;  jcol(k) = kcol;  A(k) = beZ(i)
    enddo
  enddo
! The internal points that vary
  do i = 1, nZ-1
    call I1I2 (i, NZ, m, i1, i2)
    krow = KMAP (i, 0, NR) ! Axis point itself
    do j = 0, m-1          ! The axis approximation
      k = k + 1;  kcol = krow + j
      irow(k) = krow;  jcol(k) = kcol; A(k) = beR(j)
    enddo
    do j = 1, nR-1
      krow = KMAP (i, j, NR)
      call I1I2 (j, NR, m, j1, j2)
!     The diagonal point (i,j) first:
      k = k + 1;   irow(k) = krow;  jcol(k) = krow
      iw = i - i1 + 1;  jw = j - j1 + 1
      A(k) = alR(jw,j) + alZ(iw,i) - 1/dT       ! BI case
      if (iT > 1) A(k) = A(k) - 0.5_dbl/dT      ! Must be BDF
!     Now the others in the stretch of m points along R, leaving out (i,j)
      do jj = j1, j2
        if (jj == j) cycle
        kcol = KMAP (i, jj, NR);  jw = jj - j1 + 1
        k = k + 1;   irow(k) = krow;  jcol(k) = kcol;  A(k) = alR(jw,j)
      enddo
!     ...and the others in the stretch of m points along Z, leaving out (i,j)
      do ii = i1, i2
        if (ii == i) cycle
        kcol = KMAP (ii, j, NR);  iw = ii - i1 + 1
        k = k + 1;   irow(k) = krow;  jcol(k) = kcol;  A(k) = alZ(iw,i)
      enddo
    enddo
  enddo
! Top and RHS edges, all = 1
  do j = 0, nR
    krow = KMAP(NZ,j,nR)
    k = k + 1;  irow(k) = krow;  jcol(k) = krow;  A(k) = 1
  enddo
  do i = 0, nZ-1       ! (top already done)
    krow = KMAP(i,NR,NR)
    k = k + 1;  irow(k) = krow;  jcol(k) = krow;  A(k) = 1
  enddo
  NZ28 = k
end subroutine A_SET

subroutine B_SET (B, N, nZ, nA, nR, m, C_1, C, dT, iT)
! Sets the constants vector B, according to whether this is
! BI or BDF. C_{bulk} = 1 is assumed
  use STUFF;     implicit none
  integer        :: nZ, nA, nR, N, m, iT
  real(kind=dbl) :: C_1(0:nZ,0:nR), C(0:nZ,0:nR), B(N), dT

  integer        :: i, j, k, KMAP

  B = 0
! The disk and the insulating plane, all = 0
  do j = 0, nR-1
    k = KMAP(0,j,nR);  B(k) = 0
  enddo
! Variable points in diffusion space
  do i = 1, nZ-1
    k = KMAP(i,0,NR);  B(k) = 0  ! Axis point, zero dC/dR
    do j = 1, nR-1               ! Internal points
      k = KMAP(i,j,nR)
      if (iT == 1) then          ! BI
        B(k) = - C(i,j) / dT
      else                       ! BDF
        B(k) = C_1(i,j)/(2*dT) - 2*C(i,j)/dT
      endif
    enddo
  enddo
! The top row, all = 1:
  do j = 0, nR
    k = KMAP(NZ,j,nR);  B(k) = 1
  enddo
! The right-most column, same thing:
  do i = 0, nZ-1         ! (top already done)
    k = KMAP(i,NR,NR);   B(k) = 1
  enddo
end subroutine B_SET

subroutine SHIFTEM (C_1, C, B, nZ, nR, N)
! To shift the new concs C down to past C_1, and fold B into C.
  use STUFF;   implicit none
  integer        :: nZ, nR, N
  real(kind=dbl) :: C_1(0:nZ,0:nR), C(0:nZ,0:nR), B(N)

  integer        :: i, j, k
  C_1 = C
  do k = 1, N
    call KUNMAP (k, i, j, nR)
    C(i,j) = B(k)
  enddo
end subroutine SHIFTEM

integer function KMAP (i, j, nR)
! Maps (i,j) in space to k'th row of the big matrix.
  implicit none
  integer :: i, j, nR

  if (j > nR) STOP "Bad j value in KMAP. Abort."
  KMAP = j + 1 + i*(nR+1)
end function KMAP

subroutine KUNMAP (k, i, j, nR)
! Unmaps k back to (i,j) in space.
  implicit none
  integer :: k, i, j, nR

  i = (k-1) / (nR+1)
  j = MOD(k,nR+1) - 1;   if(j < 0) j = nR
end subroutine KUNMAP

subroutine I1I2 (i, N, m, i1, i2)
! Computes i1 and i2 for the stretch of m points along N, so that
! they center i where possible, but make an asymmetric pair where not.
  implicit none
  integer :: i, N, m, i1, i2

  integer :: mm

  mm = (m+1) / 2
  i1 = MAX(0, i-mm+1);  i2 = i1 + m - 1  ! Bulk
  if (i2 > N) then;  i2 = N;  i1 = N - m + 1;  endif
end subroutine I1I2

subroutine CURRENT (C, nZ, nA, nR, m, Z, R, curr)
! Calculates the currents, by the Simpson-like method.
! The derivatives are approximated using m points.
  use STUFF;    implicit none
  integer   :: nZ, nA, nR, m
  real(dbl) :: C(0:nZ,0:nR), Z(0:nZ), R(0:nR), curr

  integer   :: j
  real(dbl),allocatable :: u(:), w(:)
  real(dbl) :: deriv, G0FORN, U_SIMP

  ALLOCATE (u(0:NZ), w(NZ))

  do j = 0, nA          ! m-point local current
     u(j) = R(j) * G0FORN (C(0:m-1,j), Z, m)
  enddo
  curr = U_SIMP (R(0:nA), u(0:nA), nA) * pi / 2
  DEALLOCATE (u, w)
end subroutine CURRENT
