function G0FUNC (C, n, H)
! Computes the n-point gradient of C wrt X at X=0, given the
! coefficients from G0BETA. N-values of 2..7 only are allowed,
! else zero is returned.
  use STUFF;  implicit none
  integer        :: n, i
  real(kind=dbl) :: C(0:*), H, G0FUNC

  integer        :: i
  real(kind=dbl) :: sum, G0BETA

  if (n>1 .and. n<8) then
    sum = 0
    do i = 0, n-1
      sum = sum + G0BETA(i,n)*C(i)
    enddo
    G0FUNC = sum / H
  else
    G0FUNC = 0
  endif
end function G0FUNC
