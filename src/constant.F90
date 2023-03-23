!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Tue 30 Mar 2021 11:14:31 AM CST
!-------------------------------------------------------------------------------

module cons
  implicit none

  integer, parameter :: MK = 8

  integer, parameter :: LSS = 32, LMS = 128, LLS = 512

  integer, parameter :: NOFATALERR = 0
  integer, parameter :: FAIL2CHECK = 1
  integer, parameter :: FAIL4FILE = 10, FAIL2OPEN = 11, FAIL2CLOSE = 12, &
    & FAIL2READ = 13, FAIL2WRITE = 14

end module cons

! vim:ft=fortran:ts=4:sw=2:et:ai:
