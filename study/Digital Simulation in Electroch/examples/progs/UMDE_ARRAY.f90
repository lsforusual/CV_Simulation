program UMDE_ARRAY
! UMDE array simulation of a linear potential scan (LSV) for
! a reversible redox reaction using the diffusion zone approach.  
! Direct discretisation using m-point approximations on an unequal grid, 
! with exponentially expanding intervals.
! BDF (3-pt) is used, started with a single BI step to avoid start-up
! complications.
!
! Input parameter:
! filename the potential-current data will be written into this file. 
! Pstart   starting potential, dimensionless.
! Pstop    final potential, dimensionless.
! nT       number of time intervals per p-unit.
! pscan    scan rate parameter, dimensionless.
! nZ       sets the number of points along Z, Z(0:nZ) including the
!          outer bulk.   
! Rd       Radius of the diffusion zone, dimensionless. 
! nR       sets the total number of points R(0:nR), the point R(nR) 
!          being the point of the diffusion zone edge.
! nA       sets the number of disk points R(0:nA) where R(nA) = 1,
!          the disk edge.
! dZ1      the first space interval in Z-direction.
! dR1      the first space interval at both sides of the the electrode 
!          edge, R=1.
! m        number of points of the finite difference formulas. m=3
!          (5 point stencil in 2-D) or m=5 (9 points stencil) are recommended.
! fac      A factor to estimate the required working space (memory) 
!          of the MA28 sparse matrix solver. A value of fac=20 is 
!          sufficient in most cases.
! output: file "filename" with potential-current data 
!          file "cprofile" with concentration data at the final potential  
!  
! Index convention: i for Z, j for R.
! The array of concs is mapped onto a large matrix, using KMAP.
! C_1 is the previous concentration array at T-dT, C the
! current one. As soon as there is a solution, C replaces C_1. 
! The Harwell suite MA28 solves the sparse system at each step.
! For the first BI, and the first BDF step, this provides an LU
! decomposition; thereafter only back-substitution is needed.

  use STUFF;  implicit none
  integer        :: nZ, nA, nR, N, nT, nTot, ier, j, L, iT, m, i
  real(kind=dbl) :: dR1, dZ1, Rd, Zmax, dT, T, &
                    gamR1, gamR2, gamZ, curr, P, dP, &
                    Pstart, Pstop, pscan
  real(kind=dbl),allocatable :: Z(:), R(:), C_1(:,:), C(:,:), &
                                alZ(:,:), alR(:,:)
! MA28-specific declarations:
  integer,parameter          :: solve_direct = 1
  integer                    :: NZ28, fac, licn, lirn, iflag
  integer, allocatable       :: irn(:), icn(:), ikeep(:), iw(:)
  real(kind=dbl)             :: u=1.0_dbl
  real(kind=dbl),allocatable :: A(:), W(:), rhs(:)

  call FILSPC (4)   ! Output file for current
  call FILSPC (14)  ! .. for the final concs.
  read*, Pstart, Pstop, nT 
  read*, pscan, Rd 
  read*, nZ, nA, nR 
  read*, dZ1, dR1
  read*, m 
  read*, fac

  dT = 1.0_dbl / nT ! time step per p unit
  nTot =  NINT(nT * ABS(pstop-pstart)) ! Total no. of time steps 
  N = (nR+1) * (nZ+1) ! Total no. of grid points
  Zmax = 6 * SQRT(ABS(pstop-pstart))/pscan ! Maximum value of Z
  ALLOCATE (C_1(0:nZ,0:nR), C(0:nZ,0:nR), Z(0:nZ),  &
            R(0:nR), alZ(1:m,0:nZ), alR(1:m,0:nR))
  C = 1;  do j = 0, nA;  C(0,j) = 0;  enddo;  ! Initialisation

  call GRID (Z, R, dZ1, dR1, nZ, nA, nR, Zmax, Rd, &
             gamZ, gamR1, gamR2)
  call COEFFS (Z, R, nZ, nR, alZ, alR, m, pscan)

  print '(" m-pt BDF3 UMDE array LSV sim., unequal intervals, &
            direct on (R,Z).")'
  print '(" m                      =", i10)', m   
  print '(" nZ, nA, nR             =", 3i10)', nZ, nA, nR
  print '(" Large NxN matrix, N    =", i10)', (nZ+1)*(nR+1)
  print '(" dZ1, dR1               =", 2es10.2)',  dZ1, dR1
  print '(" Zmax, Rd               =", 2f10.3)', Zmax, Rd
  print '(" Expansion gamZ         =",  f10.3)',  gamZ
  print '(" Expansion gamR1, gamR2 =", 2f10.3)',  gamR1, gamR2
  print '(" Pstart, Pstop, pscan   =", 3f10.3)', Pstart, Pstop, pscan
  print '(" dT, nT, nTot           =", es10.2, 2i10)', dT, nT, nTot

  call  A_PRELIM (NZ28, N, nZ, nA, nR, m)
  print '(/" NZ28                 :", i10)', NZ28
  print '(" MA28 Workspace factor=", i10)', fac
  licn = fac*NZ28;  lirn = fac*NZ28 ! Working space for MA28.
  ALLOCATE (A(licn), icn(licn), irn(lirn), rhs(N), ikeep(5*N),&
            iw(8*N), W(N))
  write (4,'("# m-pt BDF3 UMDE array LSV sim., unequal intervals, &
            direct on (R,Z).")')
  write (4,'("# UMDE_ARRAY NZ, NA, NR, m =", 4i6)') NZ, NA, NR, m
  write (4,'("# Rd, pscan =",2f10.4)') Rd, pscan
  print '(/6x,"P", 12x, "curr")'
  write(4,'("#", 6x, "P", 12x, "curr")')
   
! if (iflag /= 0) STOP "Something wrong in MA28AD. Abort 2."
  P = Pstart; dP=dT
  do iT = 1, nTot  
    P = P - dP
    if (iT <= 2) then    ! First and second steps only
      call A_SET (A, NZ28, icn, irn, N, nZ, nA, nR, alZ, alR, dT, iT, m)
      call MA28AD (N, NZ28, A, licn, irn, lirn, icn, u, ikeep, iw, w, &
                  iflag) !LUD
      if (iflag /= 0) then
        print*, "Something wrong in MA28AD, iflag, it:", iflag, it 
        print*, "Abort." 
        STOP       
      endif
    endif
    call B_SET (rhs, N, nZ, nA, nR, C_1, C, dT, iT, p)  ! RHS
    call MA28CD (N, A, licn, icn, ikeep, rhs, w, solve_direct) !Solve
    call C_UPDATE(C_1, C, rhs, nZ, nR, N) 
    call CURRENT (C, nZ, nA, nR, m, Z, R, curr)
    write(4,'(3e14.6)') P, curr
    if (mod(iT,8) == 0) then !every 8th cycle print the output on the screen 
      print '(2f12.6)', P, curr
    endif
  enddo
  call CPROFILE (C, Z, R, nZ, nR)
  DEALLOCATE (C_1, C, Z, R, alZ, alR)
  DEALLOCATE (A, icn, irn, rhs, ikeep, iw, W)
end program UMDE_ARRAY


subroutine GRID (Z, R, dZ1, dR1, nZ, nA, nR, Zmax, Rmax, &
                                        gamZ, gamR1, gamR2)
! To construct the grid, i.e. the R and Z values. The Z values
! expand exponentially from Z = 0 upwards, whereas the R values
! expand exponentially away from the R = 1 point at R(nA) in
! both directions. This gives 3 different gamma values.
  use STUFF;   implicit none
  integer        :: nZ, nA, nR
  real(kind=dbl) :: Z(0:nZ), R(0:nR), dZ1, dR1, Zmax, &
                    Rmax, gamZ, gamR1, gamR2

  integer        :: i, j
  real(kind=dbl) :: EE_FAC

! Z values, first the variable points, expanding:
  gamZ = EE_FAC (dZ1, nZ, Zmax)   ! The Z values, simple
  Z(0) = 0
  do i = 1, nZ
    Z(i) = dZ1 * (gamZ**i - 1) / (gamZ - 1)
  enddo
! R-values, first over the disk, expanding towards the axis:
  gamR1 = EE_FAC (dR1, nA, 1.0_dbl)
  R(0) = 0; R(nA) = 1
  do j = 1, nA -1
    R(nA-j) = 1 - dR1*(gamR1**j - 1) / (gamR1 - 1)
  enddo
! The insulating plane, expanding outwards:
  gamR2 = EE_FAC (dr1, nR-nA, Rmax-1)
  do j = 1, nR-nA
    R(nA+j) = 1 + dR1*(gamR2**j-1)/(gamR2-1)
  enddo 
end subroutine GRID


subroutine COEFFS (Z, R, nZ, nR, alZ, alR, m, pscan)
! Precomputes the coefficients for the discretisations.
! They are alpha_R(1:m) and alpha_Z(1:m) for all rows and cols.
! The alpha_R are composites for the 2nd and 1st derivatives
! combined (the first divided by R).
  use STUFF;   implicit none
  integer        :: nZ, nR, m
  real(kind=dbl) :: Z(0:nZ), R(0:nR), alZ(1:m,0:nZ), &
                    alR(1:m,0:nR), pscan

  integer        :: i, j, k, i1, i2, j1, j2
  real(kind=dbl) :: alpha(1:m), beta(1:m)

  alZ = 0;  alR = 0
! The coeffs along Z; inner space
  do i = 1, nZ-1
    call I1I2(i, nZ, m, i1, i2)
    call FORN(2, m, Z(i), Z(i1:i2), alpha(1:m))
    alZ(1:m,i) = alpha(1:m) / pscan**2 
  enddo
! The coeffs along R; inner space
  do j=1, nR-1
    call I1I2(j, nR, m, j1, j2)
    call FORN(2, m, R(j), R(j1:j2), alpha(1:m)) 
    call FORN(1, m, R(j), R(j1:j2),  beta(1:m))
    alR(1:m,j) = (alpha(1:m) + beta(1:m)/R(j)) / pscan**2 
  enddo 
! no flux at R(0)=0 and R(nR)=Rd
  j=0
  call FORN(1, m, R(j), R(j:j+m-1), beta(1:m))
  alR(1:m,j) = beta(1:m)
  j=nR
  call FORN(1, m, R(j), R(j-m+1:j), beta(1:m))
  alR(1:m,j) = beta(1:m)
! no-flux at insulator, Z=0, 1<R<Rd
  i=0
  call FORN(1, m, Z(i), Z(i:i+m-1), beta(1:m))
  alZ(1:m,i) = beta(1:m)    
end subroutine COEFFS


subroutine A_PRELIM (NZ28, N, nZ, nA, nR, m)
! Goes through the same motions as A_PRELIM, but only to count
! the non-zero elements in A. It was constructed by removing
! stuff from A_SET, keeping only the k increases.
  use STUFF;   implicit none
  integer        :: N, NZ28, nZ, nA, nR, m

  integer        :: i, j, k, mp

  k = 0
! The disk itself, all C = 0, entry 1 in matrix:
  do j = 0, nA; k = k + 1;  enddo

! The top side, entry 1 in matrix:
  do j = 0, nR; k = k + 1; enddo  

! The insulator, m pts
  do j=nA+1, nR
    do mp=1, m
      k = k + 1
    enddo
  enddo 

! The left and ride side, m pts
  do j=1, 2  
    do i=1, nZ-1; do mp=1, m; k = k+1; enddo; enddo
  enddo

! Inner points, m points in each direction, (m+m-1) star
  do i=1, nZ-1; do j=1, nR-1
    do mp=1, m+m-1; k=k+1; enddo
  enddo; enddo   

  NZ28 = k
end subroutine A_PRELIM


subroutine A_SET (A, NZ28, icn, irn, N, nZ, nA, nR, alZ, alR, &
                                         dT, iT, m)
! Sets the large matrix A, using the passed coeffs.
! If TIME_INT = 1 (BI), A is set for BI, otherwise for BDF(3pt).
  use STUFF;    implicit none
  integer        :: N, NZ28, irn(NZ28), icn(NZ28), nZ, nA, nR,&
                    iT, m
  real(kind=dbl) :: A(NZ28), alZ(1:m,0:nZ), alR(1:m,0:nR), dT

  integer        :: i, j, k, krow, kcol, im, i1, i2, j1, j2, ia, &
                    ja, KMAP             

  A = 0; k = 0 ! k is the count of nonzero elements to be.

! the top  
  do j=0, nR; krow = KMAP(nZ, j, nR)
      k = k + 1;  irn(k) = krow;  icn(k) = krow;  A(k) = 1
 enddo

! The disk itself, all C = 0, entry 1 in matrix:
  do j = 0, nA
    krow = KMAP(0,j,nR)
    k = k + 1;  irn(k) = krow;  icn(k) = krow;  A(k) = 1
  enddo

! left and right boundaries  
  do i=1, nZ-1
     call I1I2(i, nZ, m, i1, i2)
! left boundary, j=0, dC/dR=0, m pts. forward differences      
     krow = KMAP(i,0,nR)
     do im=0, m-1
        k = k+1; kcol=KMAP(i,im,nR); ja = im+1
        irn(k) = krow; icn(k) = kcol; A(k) = alR(ja,0)
     enddo
! right boundary, j=nR, dC/dR=0, m pts. backward differences
     krow = KMAP(i,nR,nR)
     call I1I2(nR, nR, m, j1, j2)    
     do im=j1, j2
        k=k+1; kcol=KMAP(i,im,nR); ja=im-j1+1
        irn(k)=krow; icn(k)=kcol; A(k)=alR(ja,nR)
     enddo
  enddo

!the insulating plane, i=0, dC/dZ=0, m pts forward differences
  do j=na+1, nR
    krow = KMAP(0,j,nR)
    call I1I2(0, nZ, m, i1, i2)
    do im=0, m-1
       k=k+1; ja=im+1; kcol=KMAP(im,j,nR)
       irn(k) = krow; icn(k) = kcol; A(k) = alZ(ja,0)
    enddo
  enddo   

!the inner points
  do i=1, nZ-1
    call I1I2(i, nZ, m, i1, i2)
    do j=1, nR-1
      call I1I2(j, nR, m, j1, j2)
      krow = KMAP(i, j, nR)
! The central point (i,j) first
        k=k+1; irn(k)=krow; icn(k)=krow
        ia=i-i1+1; ja=j-j1+1
        A(k) = alR(ja,j) + alZ(ia,i) - 1/dT    ! BI
        if (iT > 1) A(k) = A(k) - 0.5/dT       ! BDF correction 
!  Now the others in the stretch of m points, leaving out (i,j)
        do im=j1, j2
           if (im==j) cycle
           kcol=KMAP(i,im,nR); ja=im-j1+1
           k=k+1; irn(k)=krow; icn(k)=kcol; A(k) = alR(ja,j)
        enddo
!  The points up the column stretch, leaving out (i,j)         
        do im=i1, i2
           if(im==i) cycle
           kcol=KMAP(im,j,nR); ia=im-i1+1
           k=k+1; irn(k)=krow; icn(k)=kcol; A(k)=alZ(ia,i)
        enddo
     enddo        
   enddo
   
   if(k /= NZ28) then 
      print*, "k    = ", k  
      print*, "NZ28 = ", NZ28
      stop 'Something wrong with number of elements in A_Set'  
   endif  
end subroutine A_SET


subroutine B_SET (B, N, nZ, nA, nR, C_1, C, dT, iT, p)
! Sets the RHS vector B, according to whether this is BI or
! BDF-3. C_{bulk} = 1 is assumed. The sequence is the same as
! for SET_A.
  use STUFF;   implicit none
  integer        :: nZ, nA, nR, N, iT
  real(kind=dbl) :: C_1(0:nZ,0:nR), C(0:nZ,0:nR), B(N), dT, p
  integer        :: i, j, k, kk, KMAP

  B = 0; k = 0

! The disk and the insulating plane
  do j = 0, nR
    k = k+1
    kk = KMAP(0,j,nR)
    if(j <= nA) then
       B(kk) = exp(p) / (1 + exp(p)) !Nernst eq.
    else   
       B(kk) = 0  
    endif
  enddo

!left and right side
 j=0
 do i=1, nZ-1
     k=k+1; kk = KMAP(i,0,nR); B(kk) = 0
     k=k+1; kk = KMAP(i,nR,nR); B(kk) = 0
 enddo 

! The internal variable points, -C_{i,j}/dT or BDF equivalent:
  do i = 1, nZ-1;    do j = 1, nR-1
    kk = KMAP (i,j,nR)
    k=k+1
    if (iT == 1) then
      B(kk) = - C(i,j) / dT
    else
      B(kk) = C_1(i,j)/(2*dT) - 2*C(i,j)/dT
    endif
  enddo;  enddo

! The top row, all = 1:
  do j = 0, nR
    k=k+1; kk = KMAP(nZ,j,nR);  B(kk) = 1
  enddo

   if(k /= N) then 
      print*, "k    = ", k  
      print*, "N    = ", N
      stop 'Something wrong with number of elements in B_Set'  
   endif 
end subroutine B_SET

subroutine C_UPDATE(C_1, C, B, nZ, nR, N)
! To shift the new concs C down to past C_1, and fold B into C
 use stuff; implicit none
 integer :: nZ, nR, N
 real(kind=dbl), intent(in) :: B(N)
 real(kind=dbl), intent(out) :: C(0:nZ,0:nR), C_1(0:nZ,0:nR)
 integer :: i, j, k, KMAP
 C_1 = C
 do i=0, nZ
   do j=0, nR
      k=KMAP(i, j, nR)
      C(i,j) = B(k)
   enddo
 enddo
end subroutine C_UPDATE


integer function KMAP (i, j, nR)
! Maps (i,j) in space to k'th row of the big matrix.
  implicit none

  integer :: i, j, nR

  if (j > nR) STOP "Bad j value in KMAP. Abort."
  KMAP = j + 1 + i*(nR+1)
end function KMAP


subroutine CURRENT (C, nZ, nA, nR, m, Z, R, curr)
! Calculates the current, by a Simpson-like integration on the unequal grid. 
! Derivatives are approxinmated by m points.
  use STUFF;   implicit none

  integer        :: nZ, nA, nR, m
  real(kind=dbl) :: C(0:nZ,0:nR), Z(0:nZ), R(0:nR), curr

  integer        :: j
  real(kind=dbl),allocatable :: u(:), w(:)
  real(kind=dbl) :: deriv, U_SIMP

  ALLOCATE (u(0:na), w(na))
  u(0) = 0                 ! Setting up the current slices:
  do j = 1, nA
    call FORNBERG (1, Z(0), Z, C(0:m-1,j), m, deriv, w)
    deriv = (C(1,j) - C(0,j)) / Z(1)
    u(j) = R(j) * deriv
  enddo
  curr = U_SIMP (R(0:nA), u(0:nA), nA) * pi / 2 
  DEALLOCATE (u, w)
end subroutine CURRENT

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

subroutine CPROFILE(C, Z, R, nZ, nR)
 use stuff; implicit none
 integer, intent(in) :: nZ, nR
 real(kind=dbl), intent(in) :: C(0:nZ,0:nR), Z(0:nZ), R(0:nR)

 integer i, j

 do i = 0, nZ
    write(14,*) ' '
    do j=0, nR
       write(14,'(2i5, 3f16.8)') i, j, Z(i), R(j), C(i,j)
    enddo
 enddo
end subroutine cprofile 

