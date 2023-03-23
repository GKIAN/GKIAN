!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Wed 23 Dec 2020 09:27:16 AM CST
!-------------------------------------------------------------------------------

#define _BETA    11
#define _ALPHA   12
#define _RHO     10

#define CALCK(fi, cj) kwn = 2.0_MK * pi * fi / cj

module grtcMod

  use cons
  use math
  use paraMod, only: nLayer, h, rho, beta, alpha, Qp, Qs, fNum, cNum, eps, &
    & f, c, hasDeriv
  implicit none
  private

  integer :: N

  complex(kind = MK) :: omg, kwn
  !$OMP THREADPRIVATE(omg, kwn)

  complex(kind = MK), allocatable :: gRho(:), gBeta(:), gAlpha(:)
  complex(kind = MK), allocatable :: nu(:), gam(:), mu(:), chi(:)
  !$OMP THREADPRIVATE(gRho, gBeta, gAlpha, nu, gam, mu, chi)

  complex(kind = MK), allocatable :: E(:, :, :), iE(:, :, :)
  complex(kind = MK), allocatable :: gTd(:, :, :), gRud(:, :, :), gRdu(:, :, :)
  complex(kind = MK) :: Eus(2, 2)
  complex(kind = MK) :: gPS02
  !$OMP THREADPRIVATE(E, iE, gTd, gRud, gRdu, Eus, gPS02)

  complex(kind = MK), allocatable :: dnu(:), dgam(:), dmu(:)
  !$OMP THREADPRIVATE(dnu, dgam, dmu)

  complex(kind = MK), allocatable :: dE(:, :, :)
  complex(kind = MK), allocatable :: dgTd(:, :, :, :), dgRud(:, :, :, :), &
    & dgRdu(:, :, :, :)
  complex(kind = MK), allocatable :: dgPS02(:)
  !$OMP THREADPRIVATE(dE, dgTd, dgRud, dgRdu, dgPS02)

  real(kind = MK), allocatable, public :: kerMat(:, :), dKerMat(:, :, :)

  public grtcInitialize, grtcSpectrum, grtcDerivative, grtcEvaluate, &
    & grtcStream, grtcFinalize

  public grtcSetMedia, grtcSetEigen, grtcCoefficient, grtcKernel

  contains

    subroutine grtcInitialize()

      N = nLayer - 1

      allocate(kerMat(fNum, cNum)); kerMat = 0.0_MK

      !$OMP PARALLEL
      Eus = (0.0_MK, 0.0_MK)

      allocate(gRho  (nLayer)); gRho   = (0.0_MK, 0.0_MK)
      allocate(gBeta (nLayer)); gBeta  = (0.0_MK, 0.0_MK)
      allocate(gAlpha(nLayer)); gAlpha = (0.0_MK, 0.0_MK)
      allocate(nu (nLayer)); nu  = (0.0_MK, 0.0_MK)
      allocate(gam(nLayer)); gam = (0.0_MK, 0.0_MK)
      allocate(mu (nLayer)); mu  = (0.0_MK, 0.0_MK)
      allocate(chi(nLayer)); chi = (0.0_MK, 0.0_MK)

      allocate( E(4, 4, nLayer));  E = (0.0_MK, 0.0_MK)
      allocate(iE(4, 4, nLayer)); iE = (0.0_MK, 0.0_MK)

      allocate(gRud(2, 2, 0:0)); gRud = (0.0_MK, 0.0_MK)
      allocate(gTd (2, 2, 1:N));     gTd  = (0.0_MK, 0.0_MK)
      allocate(gRdu(2, 2, 1:N + 1)); gRdu = (0.0_MK, 0.0_MK)
      !$OMP END PARALLEL

      if(.not. hasDeriv) return

      allocate(dKerMat(fNum, cNum, nLayer)); dKerMat = 0.0_MK

      !$OMP PARALLEL
      allocate(dnu (nLayer)); dnu  = (0.0_MK, 0.0_MK)
      allocate(dgam(nLayer)); dgam = (0.0_MK, 0.0_MK)
      allocate(dmu (nLayer)); dmu  = (0.0_MK, 0.0_MK)

      allocate(dE(4, 4, nLayer)); dE = (0.0_MK, 0.0_MK)

      allocate(dgRud(2, 2, 0:0, 1:1));     dgRud = (0.0_MK, 0.0_MK)
      allocate(dgTd (2, 2, 1:N, 1:N + 1));     dgTd  = (0.0_MK, 0.0_MK)
      allocate(dgRdu(2, 2, 1:N + 1, 1:N + 1)); dgRdu = (0.0_MK, 0.0_MK)

      allocate(dgPS02(nLayer)); dgPS02 = (0.0_MK, 0.0_MK)
      !$OMP END PARALLEL
    end subroutine grtcInitialize

    subroutine grtcSetMedia(freq)
      real(kind = MK), intent(in) :: freq
      omg = 2.0_MK * pi * freq - eps * (0.0_MK, 1.0_MK)

      gRho = rho
#ifndef WITHQ
      gBeta = beta
      gAlpha = alpha
#else
      gBeta  = beta  / ( 1.0_MK + log(omg /(2.0_MK * pi)) / (pi * Qs) &
        & + (0.0_MK, 1.0_MK) / (2.0_MK * Qs) )
      gAlpha = alpha / ( 1.0_MK + log(omg /(2.0_MK * pi)) / (pi * Qp) &
        & + (0.0_MK, 1.0_MK) / (2.0_MK * Qp) )
#endif
      mu = gRho * gBeta * gBeta
    end subroutine grtcSetMedia

    subroutine grtcSetEigen(toDeriv)
      logical, intent(in), optional :: toDeriv
      complex(kind = MK) :: tau(nLayer)
      logical :: notDeriv

      if(present(toDeriv)) then
        notDeriv = (.not. (hasDeriv .and. toDeriv))
      else
        notDeriv = (.not. hasDeriv)
      end if

      nu = sqrt(kwn * kwn - (omg / gBeta) ** 2)

      chi = kwn * kwn + nu * nu
      gam = sqrt(kwn * kwn - (omg / gAlpha) ** 2)
      E(1, 1, :) = gAlpha * kwn
      E(2, 1, :) = gAlpha * gam
      E(3, 1, :) = - 2.0_MK * gAlpha * mu * kwn * gam
      E(4, 1, :) = - gAlpha * mu * chi
      E(1, 2, :) = gBeta * nu
      E(2, 2, :) = gBeta * kwn
      E(3, 2, :) = - gBeta * mu * chi
      E(4, 2, :) = - 2.0_MK * gBeta * mu * kwn * nu
      E(1, 3:4, :) =   E(1, 1:2, :)
      E(2, 3:4, :) = - E(2, 1:2, :)
      E(3, 3:4, :) = - E(3, 1:2, :)
      E(4, 3:4, :) =   E(4, 1:2, :)
      E = E / omg

      tau =  gBeta / (2.0_MK * gAlpha * mu * gam * nu * omg)
      iE(1, 1, :) =   tau * 2.0_MK * gBeta * mu * kwn * gam * nu
      iE(1, 2, :) = - tau * gBeta * mu * nu * chi
      iE(1, 3, :) = - tau * gBeta * kwn * nu
      iE(1, 4, :) =   tau * gBeta * gam * nu
      iE(2, 1, :) = - tau * gAlpha * mu * gam * chi
      iE(2, 2, :) =   tau * 2.0_MK * gAlpha * mu * kwn * gam * nu
      iE(2, 3, :) =   tau * gAlpha * gam * nu
      iE(2, 4, :) = - tau * gAlpha * kwn * gam
      iE(3:4, 1, :) =   iE(1:2, 1, :)
      iE(3:4, 2, :) = - iE(1:2, 2, :)
      iE(3:4, 3, :) = - iE(1:2, 3, :)
      iE(3:4, 4, :) =   iE(1:2, 4, :)

      if(notDeriv) return

#if DERIV == _BETA
      dnu  = omg * omg / (nu * gBeta ** 3)
      dmu  = 2.0_MK * gRho * gBeta
      dE(3, 1, :) = - 4.0_MK * gAlpha * kwn * gam * gRho * gBeta
      dE(4, 1, :) = - 4.0_MK * kwn * kwn * gAlpha * gRho * gBeta
      dE(1, 2, :) = nu + gBeta * dnu
      dE(2, 2, :) = kwn
      dE(3, 2, :) = - mu * ( chi + 4.0_MK * kwn * kwn )
      dE(4, 2, :) = - 2.0_MK * kwn * mu * ( nu + chi / nu )
      dE(1,   4, :) =   dE(1,   2, :)
      dE(2,   4, :) = - dE(2,   2, :)
      dE(3, 3:4, :) = - dE(3, 1:2, :)
      dE(4, 3:4, :) =   dE(4, 1:2, :)
      dE = dE / omg
#elif DERIV == _ALPHA
      dgam = omg * omg / (gam * gAlpha ** 3)
      dE(1, 1, :) = kwn
      dE(2, 1, :) = gam + gAlpha * dgam
      dE(3, 1, :) = - 2.0_MK * mu * kwn * ( gam + gAlpha * dgam )
      dE(4, 1, :) = - mu * chi
      dE(1, 3, :) =   dE(1, 1, :)
      dE(2, 3, :) = - dE(2, 1, :)
      dE(3, 3, :) = - dE(3, 1, :)
      dE(4, 3, :) =   dE(4, 1, :)
      dE = dE / omg
#elif DERIV == _RHO
      dmu = gBeta * gBeta
      dE(3, 1, :) = - dmu * 2.0_MK * gAlpha * kwn * gam
      dE(4, 1, :) = - dmu * gAlpha * chi
      dE(3, 2, :) = - dmu * gBeta * chi
      dE(4, 2, :) = - dmu * 2.0_MK * gBeta * kwn * nu
      dE(3, 3:4, :) = - dE(3, 1:2, :)
      dE(4, 3:4, :) =   dE(4, 1:2, :)
      dE = dE / omg
#endif
    end subroutine grtcSetEigen

    subroutine grtcCoefficient(toDeriv)
      logical, intent(in), optional :: toDeriv
      complex(kind = MK) :: L(2, 2, 0:1) = (0.0_MK, 0.0_MK)
      complex(kind = MK) :: iE21(2, 2), MT(2, 2), MR(2, 2), W(4, 4)
      complex(kind = MK) :: dL(2, 2, 0:1) = (0.0_MK, 0.0_MK)
      complex(kind = MK) :: dM(2, 2), dW(4, 4, 0:1)
      !$OMP THREADPRIVATE(L, dL)
      logical :: notDeriv
      integer i, ii

      if(present(toDeriv)) then
        notDeriv = (.not. (hasDeriv .and. toDeriv))
      else
        notDeriv = (.not. hasDeriv)
      end if

      !=> for the free surface:
      iE21 = MatInv22(E(3:4, 1:2, 1))
      gRud(:, :, 0) = - matmul( iE21, E(3:4, 3:4, 1) )

      if(hasDeriv) then
        dM = matmul( dE(3:4, 1:2, 1), gRud(:, :, 0) ) + dE(3:4, 3:4, 1)
        dgRud(:, :, 0, 1) = - matmul(iE21, dM)
      end if

      !=> for j = s, s + 1, ..., N:
      L(1, 1, 0) = exp( - gam(N + 1) * h(N + 1) )
      L(2, 2, 0) = exp( -  nu(N + 1) * h(N + 1) )
      dL(1, 1, 0) = (0.0_MK, 0.0_MK)
      dL(2, 2, 0) = (0.0_MK, 0.0_MK)

      do i = N, 1, - 1
        L(:, :, 1) = L(:, :, 0)
        L(1, 1, 0) = exp( - gam(i) * h(i) )
        L(2, 2, 0) = exp( -  nu(i) * h(i) )
        W = matmul(iE(:, :, i), E(:, :, i + 1))
        MT = MatInv22( W(1:2, 1:2) + matmul( matmul(W(1:2, 3:4), &
          & L(:, :, 1)), gRdu(:, :, i + 1) ) )
        gTd(:, :, i) = matmul(MT, L(:, :, 0))
        MR = W(3:4, 1:2) + matmul( matmul(W(3:4, 3:4), L(:, :, 1)), &
          & gRdu(:, :, i + 1) )
        gRdu(:, :, i) = matmul(MR, gTd(:, :, i))

        if(notDeriv) cycle

        dL(:, :, 1) = dL(:, :, 0)
        dL(1, 1, 0) = - h(i) * dgam(i) * L(1, 1, 0)
        dL(2, 2, 0) = - h(i) * dnu (i) * L(2, 2, 0)
        dW(:, :, 0) = - matmul( matmul(iE(:, :, i), dE(:, :, i)), W )
        dW(:, :, 1) = matmul(iE(:, :, i), dE(:, :, i + 1))

        do ii = N + 1, i + 1, - 1
          dM = matmul( matmul(W(1:2, 3:4), L(:, :, 1)), dgRdu(:, :, i + 1, ii) )
          dgTd(:, :, i, ii) = - matmul( matmul(MT, dM), gTd(:, :, i) )
        end do
        dgTd(:, :, i, i) = matmul(MT, dL(:, :, 0))
        do ii = 0, 1
          dM = dW(1:2, 1:2, ii) + matmul( matmul(dW(1:2, 3:4, ii), &
            & L(:, :, 1)), gRdu(:, :, i + 1) )
          dgTd(:, :, i, i + ii) = dgTd(:, :, i, i + ii) &
            & - matmul( matmul(MT, dM), gTd(:, :, i) )
        end do
        dM = matmul( matmul(W(1:2, 3:4), dL(:, :, 1)), gRdu(:, :, i + 1) )
        dgTd(:, :, i, i + 1) = dgTd(:, :, i, i + 1) &
          & - matmul( matmul(MT, dM), gTd(:, :, i) )

        do ii = N + 1, i + 1, - 1
          dM = matmul( matmul(W(3:4, 3:4), L(:, :, 1)), dgRdu(:, :, i + 1, ii) )
          dgRdu(:, :, i, ii) = matmul(dM, gTd(:, :, i)) &
            & + matmul(MR, dgTd(:, :, i, ii))
        end do
        dgRdu(:, :, i, i) = matmul(MR, dgTd(:, :, i, i))
        do ii = 0, 1
          dM = dW(3:4, 1:2, ii) + matmul( matmul(dW(3:4, 3:4, ii), &
            & L(:, :, 1)), gRdu(:, :, i + 1) )
          dgRdu(:, :, i, i + ii) = dgRdu(:, :, i, i + ii) &
            & + matmul(dM, gTd(:, :, i))
        end do
        dM = matmul( matmul(W(3:4, 3:4), dL(:, :, 1)), gRdu(:, :, i + 1) )
        dgRdu(:, :, i, i + 1) = dgRdu(:, :, i, i + 1) + matmul(dM, gTd(:, :, i))
      end do
    end subroutine grtcCoefficient

    subroutine grtcKernel(toDeriv)
      logical, intent(in), optional :: toDeriv
      complex(kind = MK) :: Mj(2), Ms(2, 2)
      complex(kind = MK) :: YPS(2), f0(2), P0(2, 2)
      complex(kind = MK) :: dM(2, 2), dYPS(2), dP(2, 2)
      complex(kind = MK) :: dEus(2, 2) = (0.0_MK, 0.0_MK)
      complex(kind = MK) :: df0(2) = (0.0_MK, 0.0_MK)
      !$OMP THREADPRIVATE(dEus, df0)
      logical :: notDeriv
      integer :: i

      if(present(toDeriv)) then
        notDeriv = (.not. (hasDeriv .and. toDeriv))
      else
        notDeriv = (.not. hasDeriv)
      end if

      Eus(1, 1) = exp( - gam(1) * h(1) )
      Eus(2, 2) = exp( -  nu(1) * h(1) )

      Mj = matmul( E(2, 1:2, 1), gRud(:, :, 0) ) + E(2, 3:4, 1)
      Ms = MatInv22( I22 - matmul( matmul( Eus, gRdu(:, :, 1) ), &
        & gRud(:, :, 0) ) )
      YPS = matmul(Mj, Ms)

      P0 = matmul(Eus, gRdu(:, :, 1)) - I22
      f0(1) = 0.5_MK / (gRho(1) * gAlpha(1) * omg)
      f0(2) = - kwn / (2.0_MK * gRho(1) * gBeta(1) * omg * nu(1))
      gPS02 = sum( matmul(YPS, P0) * f0 )

      if(notDeriv) return

      dEus(1, 1) = - h(1) * dgam(1) * Eus(1, 1)
      dEus(2, 2) = - h(1) * dnu (1) * Eus(2, 2)

      !=> for dYPS:
      do i = N + 1, 1, - 1
        dM = - matmul( matmul(Eus, dgRdu(:, :, 1, i)), gRud(:, :, 0) )
        dYPS = - matmul( matmul(YPS, dM), Ms )
        dgPS02(i) = sum( matmul(dYPS, P0) * f0 )
      end do
      dM = - matmul( matmul(Eus, gRdu(:, :, 1)), dgRud(:, :, 0, 1) ) &
        & - matmul( matmul(dEus, gRdu(:, :, 1)), gRud(:, :, 0) )
      dYPS = - matmul( matmul(YPS, dM), Ms )
      dgPS02(1) = dgPS02(1) + sum( matmul(dYPS, P0) * f0 )

      dM = matmul( E(1:2, 1:2, 1), dgRud(:, :, 0, 1) ) &
        & + matmul( dE(1:2, 1:2, 1), gRud(:, :, 0) ) + dE(1:2, 3:4, 1)
      dYPS = matmul(dM(2, :), Ms)
      dgPS02(1) = dgPS02(1) + sum( matmul(dYPS, P0) * f0 )

      !=> for dP:
      do i = 1, N + 1
        dP = matmul(Eus, dgRdu(:, :, 1, i))
        dgPS02(i) = dgPS02(i) + sum( matmul(YPS, dP) * f0 )
      end do
      dP = matmul(dEus, gRdu(:, :, 1))
      dgPS02(1) = dgPS02(1) + sum( matmul(YPS, dP) * f0 )

      !=> for df0:
#if DERIV == _BETA
      df0(2) = - (kwn / nu(1)) ** 2 / gBeta(1) * f0(2)
#elif DERIV == _ALPHA
      df0(1) = - f0(1) / gAlpha(1)
#elif DERIV == _RHO
      df0 = - f0 / gRho(1)
#endif
      dgPS02(1) = dgPS02(1) + sum( matmul(YPS, P0) * df0 )
    end subroutine grtcKernel

    subroutine grtcSpectrum()
      integer i, j
      !$OMP PARALLEL DO PRIVATE(i, j)
      do i = 1, fNum
        call grtcSetMedia(f(i))
        do j = 1, cNum
          CALCK(f(i), c(j))
          call grtcSetEigen(.false.)
          call grtcCoefficient(.false.)
          call grtcKernel(.false.)
          kerMat(i, j) = imag(gPS02)
        end do
      end do
      !$OMP END PARALLEL DO
    end subroutine grtcSpectrum

    subroutine grtcDerivative()
      integer i, j
      if(.not. hasDeriv) return
      !$OMP PARALLEL DO PRIVATE(i, j)
      do i = 1, fNum
        call grtcSetMedia(f(i))
        do j = 1, cNum
          CALCK(f(i), c(j))
          call grtcSetEigen()
          call grtcCoefficient()
          call grtcKernel()
          dKerMat(i, j, :) = imag(dgPS02(:))
        end do
      end do
      !$OMP END PARALLEL DO
    end subroutine grtcDerivative

    subroutine grtcEvaluate()
      integer i, j
      !$OMP PARALLEL DO PRIVATE(i, j)
      do i = 1, fNum
        call grtcSetMedia(f(i))
        do j = 1, cNum
          CALCK(f(i), c(j))
          call grtcSetEigen()
          call grtcCoefficient()
          call grtcKernel()
          kerMat(i, j) = imag(gPS02)
          if(hasDeriv) dKerMat(i, j, :) = imag(dgPS02(:))
        end do
      end do
      !$OMP END PARALLEL DO
    end subroutine grtcEvaluate

    subroutine grtcStream()
      integer i
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1, cNum
        call grtcSetMedia(f(i))
        CALCK(f(i), c(i))
        call grtcSetEigen()
        call grtcCoefficient()
        call grtcKernel()
        kerMat(1, i) = imag(gPS02)
        if(hasDeriv) dKerMat(1, i, :) = imag(dgPS02(:))
      end do
      !$OMP END PARALLEL DO
    end subroutine grtcStream

    subroutine grtcFinalize()
      if(allocated(kerMat)) deallocate(kerMat)

      !$OMP PARALLEL
      if(allocated(gRho  )) deallocate(gRho  )
      if(allocated(gBeta )) deallocate(gBeta )
      if(allocated(gAlpha)) deallocate(gAlpha)
      if(allocated(nu )) deallocate(nu )
      if(allocated(gam)) deallocate(gam)
      if(allocated(mu )) deallocate(mu )
      if(allocated(chi)) deallocate(chi)

      if(allocated( E)) deallocate( E)
      if(allocated(iE)) deallocate(iE)

      if(allocated(gTd )) deallocate(gTd )
      if(allocated(gRud)) deallocate(gRud)
      if(allocated(gRdu)) deallocate(gRdu)
      !$OMP END PARALLEL

      if(.not. hasDeriv) return

      if(allocated(dKerMat)) deallocate(dKerMat)

      !$OMP PARALLEL
      if(allocated(dnu )) deallocate(dnu )
      if(allocated(dgam)) deallocate(dgam)
      if(allocated(dmu )) deallocate(dmu )

      if(allocated(dE)) deallocate(dE)

      if(allocated(dgTd )) deallocate(dgTd )
      if(allocated(dgRud)) deallocate(dgRud)
      if(allocated(dgRdu)) deallocate(dgRdu)

      if(allocated(dgPS02)) deallocate(dgPS02)
      !$OMP END PARALLEL
    end subroutine grtcFinalize

end module grtcMod

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
