module STUFF

! General-purpose module.

  implicit none

  integer, parameter         :: sgl=selected_real_kind(6),  &
                                dbl=selected_real_kind(14), &
                                qud=selected_real_kind(30)

! CPU-time measurement: CPU_T0 is start time, -1 is present.
! These are set with a call to CPUNUL and CPUOUT resp.
  real                       :: cpu_t0, cpu_t1

  real(kind=dbl), parameter  :: small = 1.0E-08

  real(kind=dbl), parameter  :: pi  = 3.14159265358979_dbl
  real(kind=qud), parameter  :: qpi = 3.141592653589793238462643383279503_qud

end module STUFF
