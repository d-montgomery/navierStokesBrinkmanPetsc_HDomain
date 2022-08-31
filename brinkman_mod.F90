      MODULE brinkman_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys
      USE utilities_mod
      USE params_mod

      IMPLICIT NONE

CONTAINS ! Subroutines for adding in the Brinkman term

! -----------------------------------------------------------------------
! SUBROUTINE: Set alpha using Carman-Kozeny relation
! -----------------------------------------------------------------------
      SUBROUTINE carmanKozeny(t,a)
      PetscScalar, DIMENSION(:,:), INTENT(IN) :: t
      PetscScalar, DIMENSION(:,:), INTENT(OUT) :: a
      PetscScalar :: C
      
      C = 1D12 ! 1 / m^2
      a = C * (0.6 * t)**2 / (1 - 0.6*t)**3
      END SUBROUTINE carmanKozeny   

! -----------------------------------------------------------------------
! SUBROUTINE: Set platelet fraction theta for current time step
! -----------------------------------------------------------------------
      FUNCTION plateletFrac(t,Tf) RESULT (theta)
      PetscScalar :: t, Tf, t_stp, theta
    
      t_stp = 0.667*Tf
      theta = t / t_stp + 0.5D0*(SIGN(1.D0,t-t_stp) + 1)*(1 - t / t_stp)

      END FUNCTION plateletFrac
      
! -----------------------------------------------------------------------
! SUBROUTINE: A test Brinkman term that simulates a clot in the injury
!             channel with radius r, and steepness parameter S
! -----------------------------------------------------------------------
      SUBROUTINE test_B(theta,S,r,XP,B)
      PetscScalar, DIMENSION(:,:), INTENT(IN) :: theta, XP
      PetscScalar :: S, r, mid
      PetscScalar, DIMENSION(:,:), INTENT(OUT) :: B
      
      CALL carmanKozeny(theta,B)
      
      mid = Xchar * (WL + Wmid + WR) / 2D0
      B = B * 0.25D0 * (1. + ERF(S*(XP - (mid - r)))) & 
          * (1 - ERF(S*(XP- (mid + r)))) / Bchar
      END SUBROUTINE test_B

     
! -----------------------------------------------------------------------
! SUBROUTINE: Interpolate the Brinkman term to U-momentum cell center
! -----------------------------------------------------------------------
      SUBROUTINE Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                           XU,XV,YV,B,Bp)
      PetscInt, INTENT(IN) :: i, j, rowStartU, rowStartV, rowStartP
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV, YV
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B
      PetscScalar, INTENT(OUT) :: Bp
      PetscScalar :: dx, XP, xe, xw, Be, Bw

      dx = XV(i,j+1) - XV(i,j)
      XP = XU(i,j)
      xe = XV(i,j+1) 
      xw = XV(i,j)
      Be = B(i-1,j); Bw = B(i-1,j-1)
      Bp = ( Be*(XP - xw)-Bw*(XP - xe) ) / dx
      END SUBROUTINE Bp_Ucoeff

! -----------------------------------------------------------------------
! SUBROUTINE: Interpolate the Brinkman term to V-momentum cell center
! -----------------------------------------------------------------------
      SUBROUTINE Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                           XU,YU,YV,B,Bp)
      PetscInt, INTENT(IN) :: i, j, rowStartU, rowStartV, rowStartP
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU, YU
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: YV
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B
      PetscScalar, INTENT(OUT) :: Bp
      PetscScalar :: dy, YP, yn, ys, Bn, Bs

      dy = YU(i+1,j)-YU(i,j)
      YP = YV(i,j)
      yn = YU(i+1,j)
      ys = YU(i,j)
      Bn = B(i,j-1)
      Bs = B(i-1,j-1)
      Bp = (Bn*(YP-ys) - Bs*(YP-yn))/dy

      END SUBROUTINE Bp_Vcoeff
      
      END MODULE brinkman_mod
