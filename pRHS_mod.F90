      MODULE pRHS_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
#include <petsc/finclude/petscvec.h>  
      USE petscsys
      USE petscvec

! Necessary F90 Modules
      USE params_mod
      IMPLICIT NONE

CONTAINS

! -----------------------------------------------------------------------
! Subroutine to build the RHS vector in pressure-Poisson equation 
! -----------------------------------------------------------------------
      SUBROUTINE pRHS(N,NNXY,rowStartU,rowStartV,i_locStart,i_locEnd,&
                      U,V,P0,XU,YV,yp,dt,f,ierr)
      IMPLICIT NONE
      ! Inputs
      PetscInt, DIMENSION(:), INTENT(IN) :: N !N=[NxL,NxM,Nx,NyL,NyM,Ny]
      PetscInt, INTENT(IN) :: NNXY ! Au is NNXY x NNXY matrix
      PetscInt, INTENT(IN) :: rowStartU, rowStartV, i_locStart, i_locEnd
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: U
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: V
      PetscScalar, DIMENSION(i_locStart-2:,:), INTENT(IN) :: P0  
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: YV  
      PetscScalar, DIMENSION(:), INTENT(IN) :: yp
      PetscScalar, INTENT(IN) :: dt 
      Vec, INTENT(INOUT) :: f ! vector
      PetscErrorCode, INTENT(INOUT) :: ierr
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny, Nsml
      PetscInt :: i,j,j_ctr,row
      PetscScalar :: dx,dy,val,Pb,BC,as,an

      
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
      
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL) + (Nx - NxM)

! --- Fill Inner Cells: Left Lower Legs of H ----------------------------  
      DO i = MAX(1,i_locStart), MIN(NyL, i_locEnd)
      DO j = 1,NxL
            row = Nsml*(i-1)+j-1

            dx = XU(i+1,j+1) - XU(i+1,j)  
            dy = YV(i+1,j+1) - YV(i,j+1)
            
            BC = 0
!            IF (i == 1) THEN ! Apply fixed pressure at boundary
!                  as = 1/(yp(i+1)-yp(i))/dy
!                  Pb = (3*P0(i,j) - P0(i+1,j))/2
!                  BC = -2*as*( Pw_out - Pb)
!            END IF
            
            val = 1.D0/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx &
                              +  ( V(i+1,j+1) - V(i,j+1) )/dy ) + BC

            CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
      END DO
      END DO

! --- Fill Inner Cells: Right Lower Legs of H ----------------------------
      DO i = MAX(1,i_locStart),MIN(NyL,i_locEnd)
            j_ctr = NxL
      DO j = NxM+1,Nx

            j_ctr = j_ctr + 1
            row = Nsml*(i-1)+j_ctr-1
            dx = XU(i+1,j+1) - XU(i+1,j)  
            dy = YV(i+1,j+1) - YV(i,j+1)

            BC = 0
!            IF (i == 1) THEN ! Apply fixed pressure at boundary
!                  as = 1/(yp(i+1)-yp(i))/dy
!                  Pb = (3*P0(i,j) - P0(i+1,j))/2
!                  BC = -2*as*( Pb_out - Pb)
!            END IF
            
            val = 1.D0/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx &
                              +  ( V(i+1,j+1) - V(i,j+1) )/dy ) + BC

            CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
      END DO
      END DO

! --- Fill Inner Cells: Injury Channel -------------------------------------
      DO i = MAX(NyL+1,i_locStart), MIN(NyM, i_locEnd)
      DO j = 1, Nx
            row = Nsml*(NyL) + (i - NyL -1) * Nx + j - 1

            dx = XU(i+1,j+1) - XU(i+1,j)  
            dy = YV(i+1,j+1) - YV(i,j+1)

            val = 1.D0/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx &
                              +  ( V(i+1,j+1) - V(i,j+1) )/dy )

            CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
      END DO
      END DO

! --- Fill Inner Cells: Left Upper Legs of H ----------------------------
      DO i = MAX(NyM+1, i_locStart), MIN(Ny, i_locEnd)
      DO j = 1,NxL
            row = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j - 1

            dx = XU(i+1,j+1) - XU(i+1,j)  
            dy = YV(i+1,j+1) - YV(i,j+1)

            BC = 0
!            IF (i == Ny) THEN ! Apply fixed pressure at boundary
!                  an = 1/(yp(i+2)-yp(i+1))/dy
!                  Pb = (3*P0(i-1,j) - P0(i,j))/2
!                  BC = -2*an*( Pw_in - Pb)
!            END IF
            
            val = 1.D0/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx &
                              +  ( V(i+1,j+1) - V(i,j+1) )/dy ) + BC

            CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
      END DO
      END DO

! --- Fill Inner Cells: Right Upper Legs of H ----------------------------
      DO i = MAX(NyM+1, i_locStart), MIN(Ny, i_locEnd)
            j_ctr = NxL
      DO j = NxM+1,Nx
            j_ctr = j_ctr + 1
            row = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j_ctr-1

            dx = XU(i+1,j+1) - XU(i+1,j)  
            dy = YV(i+1,j+1) - YV(i,j+1)
            
            BC = 0
!            IF (i == Ny) THEN ! Apply fixed pressure at boundary
!                  an = 1/(yp(i+2)-yp(i+1))/dy
!                  Pb = (3*P0(i-1,j) - P0(i,j))/2
!                  BC = -2*an*( Pb_in - Pb)
!            END IF
            
            val = 1.D0/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx &
                              +  ( V(i+1,j+1) - V(i,j+1) )/dy ) + BC

            CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
      END DO
      END DO

!! ---  Augmentation ------------------------------------------------------ 
!      IF (i_locEnd == Ny) THEN 
!            CALL VecSetValues(f,1,NNXY-1,0.D0,INSERT_VALUES,ierr)
!      END IF

      CALL VecAssemblyBegin(f,ierr)   
      CALL VecAssemblyEnd(f,ierr)
      END SUBROUTINE pRHS
      
      END MODULE pRHS_mod
