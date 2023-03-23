!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Tue 22 Dec 2020 10:37:01 AM CST
!-------------------------------------------------------------------------------

module paraMod

  use cons
  use math
  implicit none
  public

  character(len = LLS) :: inputFile = 'input.nml', modelFile = 'model.dat'

  integer :: fNum, cNum
  real(kind = MK), private :: fMin = 0.1_MK, fMax =  10.0_MK, fStep = 0.1_MK
  real(kind = MK), private :: cMin = 1.0_MK, cMax = 100.0_MK, cStep = 1.0_MK
  real(kind = MK), private :: wAdjust = 3.0_MK
  real(kind = MK) :: eps

  real(kind = MK), allocatable :: f(:), c(:)

  integer :: nLayer = 0
  real(kind = MK), allocatable :: h(:), z(:), rho(:), beta(:), alpha(:), &
    & Qs(:), Qp(:)

  logical :: hasDeriv = .false.

  interface paraInput
    module procedure paraInputReal1D
    module procedure paraInputReal2D
  end interface paraInput

  namelist /para/ fMin, fMax, fNum, cMin, cMax, cNum, wAdjust
  private para

  contains

    subroutine paraGetArguments(input, model, deriv)
      character(len = *), intent(in), optional :: input, model
      logical, intent(in), optional :: deriv
      !f2py character(len = *) optional, intent(in) :: input = 'input.nml'
      !f2py character(len = *) optional, intent(in) :: model = 'model.dat'
      !f2py logical optional, intent(in) :: deriv = 0
      integer :: nArgs, pivot = 0
      character(len = LSS) :: argStr = 'NULL'
      nArgs = command_argument_count()
      if(nArgs >= 1) then
        call get_command_argument(1, argStr)
        if(argStr == '-h' .or. argStr == '--help') then
          call paraPrintHelp()
          stop
        else if(argStr == '-d' .or. argStr == '--deriv') then
          hasDeriv = .true.
          pivot = 1
        end if
        if(nArgs >= 1 + pivot) call get_command_argument(1 + pivot, inputFile)
        if(nArgs >= 2 + pivot) call get_command_argument(2 + pivot, modelFile)
        if(nArgs >= 3 + pivot) &
          & call paraErrorExcept(NOFATALERR, 'Some invalid arguments found.')
      else
        if(present(input)) inputFile = input
        if(present(model)) modelFile = model
        if(present(deriv)) hasDeriv = deriv
      end if

#ifdef DEBUG
      write(*, '(A)') 'paraGetArguments: inputFile = ' // trim(inputFile)
      write(*, '(A)') 'paraGetArguments: modelFile = ' // trim(modelFile)
      write(*, '(A, G0)') 'paraGetArguments: hasDeriv = ', hasDeriv
#endif
    end subroutine paraGetArguments

    subroutine paraPrintHelp()
      character(len = LSS) :: exeName
      call get_command_argument(0, exeName)
      write(*, *)
      write(*, '(A)') 'Usage:'

      write(*, '(A)') '  ' // trim(adjustl(exeName)) &
        & // ' [-d] [ inputFile [modelFile] ]'
      write(*, '(A)') "               calculate Green's function integration " &
        & // 'kernel and/or its derivative for [gPS0]_2.'
      write(*, *)
      write(*, '(A)') '  -d           to calculate derivative'
      write(*, '(A)') '  inputFile    the input file, [default: input.nml]'
      write(*, '(A)') '  modelFile    the model file, [default: model.dat]'
      write(*, *)
    end subroutine paraPrintHelp

    subroutine paraInitialize(wadj)
      real(kind = MK), intent(in), optional :: wadj
      !f2py real(kind = MK) optional, intent(in) :: wadj = 0.0
      integer :: fileID, ioStatus
      logical :: readNML = .true.
      integer :: i, iTemp

      if(present(wadj)) then
        if(wadj /= 0.0_MK) readNML = .false.
      end if

      if(readNML) then
        ! read input parameters from file
        open(newunit = fileID, file = inputFile, status = 'old')
          read(fileID, nml = para)
        close(fileID)
      else
        wAdjust = wadj
      end if

      ! read model parameters from file
      !=> get number of lines of file
      open(newunit = fileID, file = modelFile, status = 'old')
        nLayer = - 1
        read(fileID, *, iostat = ioStatus)
        do while(ioStatus == 0)
          nLayer = nLayer + 1
          read(fileID, *, iostat = ioStatus) iTemp
        end do
      close(fileID)
      !=> allocate model parameters
      allocate(h(nLayer))
      allocate(z(0:nLayer))
      allocate(rho(nLayer))
      allocate(beta(nLayer), alpha(nLayer))
      allocate(Qs(nLayer), Qp(nLayer))
      !=> read model parameters
      call paraInputModel(modelFile)

      ! other parameters
      eps = pi / wAdjust
      fStep = (fMax - fMin)/(fNum - 1)
      cStep = (cMax - cMin)/(cNum - 1)
      allocate(f(fNum), c(cNum))
      do i = 1, fNum
        f(i) = fMin + (i - 1) * fStep
      end do
      do i = 1, cNum
        c(i) = cMin + (i - 1) * cStep
      end do

#ifdef DEBUG
      write(*, '(A, 3(1X, G0))') 'paraInitialize: f =', fMin, fMax, fNum
      write(*, '(A, 3(1X, G0))') 'paraInitialize: c =', cMin, cMax, cNum
#endif
    end subroutine paraInitialize

    subroutine paraReconfigure(nLayer_, fNum_, cNum_, wAdjust_)
      integer, intent(in), optional :: nLayer_, fNum_, cNum_
      real(kind = MK), intent(in), optional :: wAdjust_
      !f2py integer optional, intent(in) :: nLayer_ = 0, fNum_ = 0, cNum_ = 0
      !f2py real(kind = MK) optional, intent(in) :: wAdjust_ = 3.0
      if(present(nLayer_)) then
        nLayer = nLayer_
        if(allocated(h)) deallocate(h)
        if(allocated(z)) deallocate(z)
        if(allocated(rho  )) deallocate(rho  )
        if(allocated(beta )) deallocate(beta )
        if(allocated(alpha)) deallocate(alpha)
        if(allocated(Qs)) deallocate(Qs)
        if(allocated(Qp)) deallocate(Qp)
        allocate(h(nLayer)); h = 0.0_MK
        allocate(z(0:nLayer)); z = 0.0_MK
        allocate(rho  (nLayer)); rho   = 0.0_MK
        allocate(beta (nLayer)); beta  = 0.0_MK
        allocate(alpha(nLayer)); alpha = 0.0_MK
        allocate(Qs(nLayer)); Qs = 0.0_MK
        allocate(Qp(nLayer)); Qp = 0.0_MK
      end if
      if(present(fNum_)) then
        fNum = fNum_
        if(allocated(f)) deallocate(f)
        allocate(f(fNum)); f = 0.0_MK
      end if
      if(present(cNum_)) then
        cNum = cNum_
        if(allocated(c)) deallocate(c)
        allocate(c(cNum)); c = 0.0_MK
      end if
      if(present(wAdjust_)) wAdjust = wAdjust_
    end subroutine paraReconfigure

    subroutine paraInputModel(fileName)
      character(len = *), intent(in) :: fileName
      integer :: fileID, i, iTemp
      open(newunit = fileID, file = fileName, status = 'old')
        read(fileID, *)
        do i = 1, nLayer
#ifndef WITHQ
          read(fileID, *) iTemp, z(i - 1), rho(i), beta(i), alpha(i)
#else
          read(fileID, *) iTemp, z(i - 1), rho(i), beta(i), alpha(i), &
            & Qs(i), Qp(i)
#endif
        end do
      close(fileID)
      z(nLayer) = inf
      h = z(1:nLayer) - z(0:nLayer - 1)
    end subroutine paraInputModel

    subroutine paraInputReal1D(fileName, realVar)
      character(len = *), intent(in) :: fileName
      real(kind = MK), intent(out) :: realVar(:)
      integer :: fileID
      open(newunit = fileID, file = fileName, status = 'old')
        read(fileID, *) realVar
      close(fileID)
    end subroutine paraInputReal1D
    subroutine paraInputReal2D(fileName, realVar)
      character(len = *), intent(in) :: fileName
      real(kind = MK), intent(out) :: realVar(:, :)
      integer :: fileID
      open(newunit = fileID, file = fileName, status = 'old')
        read(fileID, *) realVar
      close(fileID)
    end subroutine paraInputReal2D

    subroutine paraFinalize()
      if(allocated(f)) deallocate(f)
      if(allocated(c)) deallocate(c)
      if(allocated(h)) deallocate(h)
      if(allocated(z)) deallocate(z)
      if(allocated(rho)) deallocate(rho)
      if(allocated(beta)) deallocate(beta)
      if(allocated(alpha)) deallocate(alpha)
      if(allocated(Qs)) deallocate(Qs)
      if(allocated(Qp)) deallocate(Qp)
    end subroutine paraFinalize

    subroutine paraErrorExcept(errorCode, suppInfo)
      integer, intent(in) :: errorCode
      character(len = *), intent(in) :: suppInfo
      if(errorCode == NOFATALERR) then
#ifndef COLORPRINT
        write(*, '(A)') repeat('-', len(suppInfo) + 13)
        write(*, '(A)') '>>> Warning: ' // suppInfo
        write(*, '(A)') repeat('-', len(suppInfo) + 13)
#else
        write(*, '(A)') char(27) // '[00;93;100mWarning:' &
          & // char(27) // '[0m ' // suppInfo
#endif
      else
#ifndef COLORPRINT
        write(*, '(A)') repeat('=', len(suppInfo) + 14)
        write(*, '(A)') '>>>>>> ERROR: ' // suppInfo
        write(*, '(A)') repeat('=', len(suppInfo) + 14)
#else
        write(*, '(A)') char(27) // '[01;91;100mERROR:' &
          & // char(27) // '[04;39;49m ' // suppInfo // char(27) // '[0m'
#endif
        stop errorCode
      end if
    end subroutine paraErrorExcept

end module paraMod

! vim:ft=fortran tw=80 ts=2 sw=2 et ai 
