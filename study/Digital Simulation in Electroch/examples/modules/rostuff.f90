module ROSTUFF
! Module for the Rosenbrock coefficients.
  use STUFF;   implicit none
  real(kind=dbl) :: gamma, gamma1, gamma2, gamma3, &
                    alpha1,  alpha2,  alpha3,      &
                    a21, a31, a32,                 &
                    c21, c31, c32, m1, m2, m3
CONTAINS
  subroutine ROCOEFFS (order)
  ! Sets the Rosenbrock coeffs for orders 2 or 3,
  ! accessing the module ROSTUFF where they are.
  ! Order 2 is ROS2, order 3 is ROWDA3.
  use STUFF;  implicit none
  integer :: order

  gamma1  = 0;   gamma2 = 0;   gamma3 = 0 ! Zero defaults
  alpha1  = 0;   alpha2 = 0;   alpha3 = 0
  a21 = 0;   a31 = 0;   a32 = 0
  c21 = 0;   c31 = 0;   c32 = 0
  m1 = 0
  select case (order)
    case (2)
      gamma = 1.707106781186547_dbl
      gamma1 = 0;   gamma2 = - gamma
      alpha2 = 1
      a21 = 0.5857864376269050_dbl
      c21 = - 1.171572875253810_dbl
      m1  = 0.8786796564403575_dbl
      m2  = 0.2928932188134525_dbl
    case (3)
      gamma = 0.435866521508459_dbl
      gamma1 = gamma;  gamma2 = 0.6044552840655588_dbl
                       gamma3 = 6.3797887993448800_dbl
      alpha2 = 0.7_dbl;  alpha3 = 0.7_dbl
      a21 = 1.605996252195329_dbl
      a31 = a21;  a32 = 0
      c21 = 0.8874044410657823_dbl
      c31 = 23.98747971635035_dbl
      c32 = 5.263722371562130_dbl
      m1  = 2.236727045296589_dbl
      m2  = 2.250067730969645_dbl
      m3  = -0.209251404439032_dbl
    end select
  end subroutine ROCOEFFS
end module ROSTUFF
