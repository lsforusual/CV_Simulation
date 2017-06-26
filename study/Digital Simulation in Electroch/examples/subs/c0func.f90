function C0FUNC (C, n, GH)
! Computes C0, given G*H and the n-point gradient of C wrt X
! at X=0, given the coeffícients from G0BETA. N-values of 2..7
! only are allowed, else zero is returned.
  use STUFF;  implicit none
  integer        :: n, i
  real(kind=dbl) :: C(0:*), GH, C0FUNC

  real(kind=dbl) :: sum, G0BETA

  if (n>1 .and. n<8) then
    sum = 0
    do i = 1, n-1
      sum = sum + G0BETA(i,n)*C(i)
    enddo
    C0FUNC = (GH - sum) / G0BETA(0,n)
  else
    C0FUNC = 0
  endif
end function C0FUNC
