 program CHANNEL_BAND
! To simulate the steady state at a band in a channel with
! laminar flow going through it. The equation then reduces to
! a 1D parabolic form, and a march along X can be used.
! There are two independent variables apart from simulation
! variables: Pe, the Peclet number and H, the channel height in
! units of the electrode length. Pe determines
! Y_max = 6 / SQRT(Pe), and H must be > Y_max, so that the
! bulk concentration boundary condition at Y_m can be applied.
! Simulation parameters are NX, the number of intervals in X,
! and NY, those between 0 and Y_max.
! Current is integrated by the trapezium rule.
  use STUFF;   implicit none
  integer        :: NX, NY, i, j
  real(kind=dbl) :: H, Pe, Y_max, Isum, Ilocal, bi, dX, dY
  real(kind=dbl),allocatable :: C(:), ai(:),a2i(:), ad(:),bd(:)

  print '(" Pe, H?")';  read *, Pe, H
  Y_max = 6 / SQRT(Pe)
  if (Y_max > H) STOP "H too small. Abort."
  print '(" NX, NY?")';  read *, NX, NY
  ALLOCATE (C(0:NY+1), ai(NY), a2i(NY), ad(NY), bd(NY))
  dX = 1.0_dbl / NX;   dY = Y_max / NY
  call PRECALC (ai, a2i, NX, NY, H, Pe, Y_max)
  print '(" CHANNEL_BAND")'
  print '(" Pe    =", f10.2)', Pe
  print '(" H     =", f10.2)', H
  print '(" NX    =", i7)', NX
  print '(" NY    =", i7)', NY
  print '(" Y_max =", f11.3)', Y_max

  Isum = 0 ! Current integration start with I = 0.
  C = 1;  C(0) = 0
  do j = 1, NX
    ad(NY) = ai(NY)
    bi = a2i(NY)*C(NY);  bd(NY) = bi - C(NY+1)
    do i = NY-1, 1, -1
      ad(i) = ai(i) - 1/ad(i+1)
      bi = a2i(i) * C(i)
      bd(i) = bi - bd(i+1)/ad(i+1)
    enddo
    C(0) = 0
    do i = 1, NY
      C(i) = (bd(i) - C(i-1)) / ad(i)
    enddo
    Ilocal = (-3*C(0) + 4*C(1) - C(2)) / (2*dY)
    Isum = Isum + Ilocal ! trapezium integration process
  enddo
  Isum = Isum - Ilocal/2 ! Correction to last integration point
  Isum = Isum * dX ! Final current
  print '(/" Current       =", f8.2)', Isum
  print '(" I / Pe**(1/3) =", f8.2)', Isum / Pe**0.3333333
end program CHANNEL_BAND

subroutine PRECALC (ai, a2i, NX, NY, H, Pe, Y_max)
! Precomputes the two constant arrays ai and a2i.
  use STUFF;   implicit none
  integer        :: NX, NY
  real(kind=dbl) :: ai(NY), a2i(NY), Pe, H, Y_max, lami

  integer        :: i
  real(kind=dbl) :: Y, VX, dY, dX

  dY = Y_max / NY ;  dX = 1.0_dbl / NX
  do i = 1, NY
    Y = i * dY
    VX = 1.5_dbl * (1 - ((Y-H)/H)**2)
    lami = dX / (VX * Pe * dY**2)
    ai(i) = - 2 -  1/lami
    a2i(i) = - 1 / lami
  enddo
end subroutine PRECALC

