subroutine FILSPC (lun)
! Reads in file name and opens the file, for LUNs ending with 3 or 4.
! If the LUN ends in 3, it is an input file, if 4, an output file.
! We allow a maximum of 10 files (FILNUM's) for each LUN-ending (FILTYPE).
! Thus, you can have input LUNs 3, 13, 23, .. and output LUNs 4, 14, 24...
  implicit none
  character(LEN=40) :: filname
  integer           :: filnum, filtype, lun

  filnum = lun / 10
  filtype = MOD (lun, 10)

! STOP if LUN is out of range:
  if (lun<1  .or. lun>99)  stop 'Bad LUN'
! Range OK, go on:

  read '(a)', filname

  select case (filtype)
  case (3)             ! Input file:
      open (unit=lun, file=filname, status='old',action="read")
  case (4)             ! Output file:
      open (unit=lun, file=filname, status='new',action="write")
  case default
      STOP " Bad LUN"
  end select
end subroutine FILSPC
