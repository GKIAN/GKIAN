!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Tue 22 Dec 2020 10:12:11 AM CST
!-------------------------------------------------------------------------------

program main

  !$ use omp_lib
  use cons
  use math
  use paraMod
  use grtcMod
  implicit none

  character(len = LLS) :: middleName, kernelFile, dKernelFile
  character(len = LSS) :: fmtStr
  integer :: fileID, i, j
  !$ real(kind = MK) :: tbeg, tend

  call paraGetArguments()

  call mathInitialize()
  call paraInitialize()
  call grtcInitialize()

  middleName = inputFile(index(inputFile, '/', .true.) + 1 &
    & :index(inputFile, '.nml') - 1)
  if(len_trim(middleName) == 0) then
    call paraErrorExcept(NOFATALERR, &
      & "Failed to automatically format outputfile name.")
    kernelFile = 'gPS02.ker'
  else
    kernelFile = 'gPS02.' // trim(adjustl(middleName)) // '.ker'
  end if
  dKernelFile = 'd' // trim(adjustl(kernelFile))

  !$ tbeg = omp_get_wtime()
  call grtcEvaluate()
  !$ tend = omp_get_wtime()
  !$ write(*, '(A, F10.6, A)') 'OpenMP elapsed time is ', tend - tbeg, ' s.'

  write(fmtStr, '(A, G0, A, A)') '(G0, ', cNum, '(2X, G0)', ')'
  write(*, '(A)') '> Output kernel spectrum file: [' &
    & // trim(adjustl(kernelFile)) // '] ...'
  open(newunit = fileID, file = kernelFile)
    write(fileID, fmtStr) 0, (c(j), j = 1, cNum)
    do i = 1, fNum
      write(fileID, fmtStr) f(i), (kerMat(i, j), j = 1, cNum)
    end do
  close(fileID)
  if(hasDeriv) then
    write(*, '(A)') '> Output kernel derivative spectrum file: [' &
      & // trim(adjustl(dKernelFile)) // '] ...'
    open(newunit = fileID, file = dKernelFile, access = 'stream')
      write(fileID) dKerMat
    close(fileID)
  end if

  call grtcFinalize()
  call paraFinalize()

end program

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
