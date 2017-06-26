function ERFC(X)
! Modified, Feb-2001, to Fortran90. Original IBM name: DERFC

!     THIS DOUBLE-PRECISION FUNCTION IS DESIGNED TO CALCULATE
!     ERFC(X) FOR ALL VALUES OF X.

!**********************************************************************

  use STUFF;  implicit none

  real(kind=dbl) :: ERFC, X
  real(kind=dbl) :: A0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, &
                    b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, &
                    b12, b13, b14, b15, b16, b17, b18, c0, c1, c2, c3,&
                    c4, c5, c6, c7, d1, d2, d3, d4, d5, d6, d7,       &
                    xx, pa, terf, terfc, z, test, f

  if (ABS(x) < TINY(1.0_dbl)) then
      ERFC = 1
  else
!     DEFINE CONSTANTS USED IN THE EXPANSIONS.
    A0=1.128379167095513_dbl;   A1=-.3761263890318352_dbl
    A2=0.1128379167094419_dbl;  A3=-.0268661706431117_dbl
    A4=0.0052239776061200_dbl;  A5=-.0008548325929373_dbl
    A6=0.0001205529357828_dbl;  A7=-.0000149247123228_dbl
    A8=0.0000016447131768_dbl;  A9=-.0000001620631488_dbl
    A10=0.0000000137109838_dbl; A11=-.0000000007779473_dbl
    B0=0.0156249828805124_dbl;  B1=-.0607180339828094_dbl
    B2=0.1037958188348712_dbl;  B3=-.0980513980932520_dbl
    B4=0.0492094856954216_dbl;  B5=-.0042334886622643_dbl
    B6=-.0107101851479655_dbl;  B7=0.0062390502665986_dbl
    B8=-.0003713317915058_dbl;  B9=-.0010720864304375_dbl
    B10=0.0004325546823446_dbl; B11=0.0000409917610466_dbl
    B12=-.0000772102070427_dbl; B13=0.0000145149904597_dbl
    B14=0.0000065839110963_dbl; B15=-.0000033359210544_dbl
    B16=0.0000000465742265_dbl; B17=0.0000005954826731_dbl
    B18=0.0000000838046807_dbl
    C0=0.5641895835477550_dbl;  C1=1.499999999360903_dbl
    C2=3.499994562238586_dbl;   C3=5.497021993223130_dbl
    C4=7.290752372050481_dbl;   C5=7.047138950515183_dbl
    C6=3.897999536575896_dbl;   C7=1.222637858242353_dbl
    D1=-.2820947917731498_dbl;  D2=-1.499999909587046_dbl
    D3=-4.999659851121717_dbl;  D4=-10.40367962659367_dbl
    D5=-14.51757254291838_dbl;  D6=-7.698257129012999_dbl
    D7=-1.311612589053097_dbl

!     METHOD OF COMPUTATION FOR ERF AND ERFC DEPENDS ON THE
!     RANGE IN WHICH THE ARGUMENT LIES.
!***********************************************************************
!                         X < 0
!     FOR NEGATIVE VALUES USE ERF(-X)=-ERF(X) AND ERFC(-X)=2-ERFC(X)

    XX = ABS(X)

    if (xx < 1.0) then
      PA=A0+A1*XX**2+A2*XX**4+A3*XX**6+A4*XX**8+A5*XX**10+A6*XX**12+   &
      A7*XX**14+A8*XX**16+A9*XX**18+A10*XX**20+A11*XX**22
      TERF=XX*PA
      TERFC = 1 - TERF
    else if (xx < 2.040452) then   ! In this range and THE NEXT, USE THE
!                                    RELATIONSHIP ERF(XX)=1-ERFC(XX)
      Z=XX-1.70947265625_dbl
      TERFC=B0+B1*Z+B2*Z**2+B3*Z**3+B4*Z**4+B5*Z**5+B6*Z**6+B7*Z**7+   &
            B8*Z**8+B9*Z**9+B10*Z**10+B11*Z**11+B12*Z**12+B13*Z**13+   &
            B14*Z**14+B15*Z**15+B16*Z**16+B17*Z**17+B18*Z**18

    else IF(XX < 13.306) then      ! 2.040452 <= XX < 13.306:
      Z=XX*XX
      F=C0+D1/(Z+C1+D2/(Z+C2+D3/(Z+C3+D4/(Z+C4+D5/(Z+C5+               &
        D6/(Z+C6+D7/(Z+C7)))))))
      TERFC= EXP(-Z)*F/XX
    else                           ! XX >= 13.306:
      TERFC = 0
    endif

    IF (x < 0) then
       TERFC = 2 - TERFC
    endif

    ERFC = TERFC
  endif
end function ERFC
