      MODULE Bv_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
#include <petsc/finclude/petscmat.h>  
      USE petscsys
      USE petscmat

! Necessary F90 Modules
      USE params_mod
      USE brinkman_mod
      IMPLICIT NONE

CONTAINS

! -----------------------------------------------------------------------
! Subroutine to build the Bv matrix in v-momentum equation 
! -----------------------------------------------------------------------
      SUBROUTINE build_Bv(rank,N,NNXY,locDim,rowStartU,rowStartV,&
                          rowStartP,XU,YU,YV,B,A,ierr)
      IMPLICIT NONE
      ! Inputs
      PetscInt, INTENT(IN) :: rank !rank in communicator
      PetscInt, DIMENSION(:), INTENT(IN) :: N !N=[NxL,NxM,Nx,NyL,NyM,Ny]
      PetscInt, INTENT(IN) :: rowStartU,rowStartV,rowStartP
      PetscInt, INTENT(IN) :: NNXY ! Au is NNXY x NNXY matrix
      PetscInt, INTENT(IN) :: locDim !number of rows ea. core allocates
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU,YU 
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: YV 
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B
      Mat, INTENT(INOUT) :: A ! matrix
      PetscErrorCode, INTENT(INOUT) :: ierr
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny, Nsml
      PetscInt :: i,j,j_ctr,iStart,iEnd,jEnd,i_locStart,i_locEnd,row
      PetscScalar :: dx, dy, Bp

      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
      
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL + 2) + (Nx + 2 - NxM)
      
!      CALL MatSetSizes(A, locdim, locdim, NNXY, NNXY, ierr)
!      CALL MatSetFromOptions(A, ierr)
!      CALL MatSeqAIJSetPreallocation(A, 5, PETSC_NULL_INTEGER, ierr)
!      CALL MatMPIAIJSetPreallocation(A, 5, PETSC_NULL_INTEGER, 4, &
!                PETSC_NULL_INTEGER, ierr)

      CALL MatGetOwnershipRange(A, i_locStart, i_locEnd, ierr)

!! --- Boundary Conditions -----------------------------------------------
!      ! Top Boundaries (\partial \Omega_i, i = 1,3)
!      iStart = NNXY - Nsml
!      IF (i_locEnd .GE. iStart) THEN
!            DO i = MAX(i_locStart,iStart), i_locEnd-1
!                  CALL MatSetValue(A,i,i,0.D0,INSERT_VALUES,ierr)
!            END DO
!      END IF
!
!      ! Bottom Boundaries (\partial \Omega_i, i = 2,4)
!      iEnd = Nsml -1
!      IF (i_locStart .LE. iEnd) THEN
!            DO i = i_locStart, MIN(i_locEnd,iEnd)
!                  CALL MatSetValue(A,i,i,0.D0,INSERT_VALUES,ierr) 
!            END DO
!      END IF
!
!      ! Left/Right Boundary - Lower Legs (\partial \Omega_i, i=5,6)
!      iStart = Nsml; 
!      iEnd = NyL*Nsml - 1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = 2,NyL
!            ! Left Boundary
!            row = (i-1)*Nsml
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            row = i*Nsml - 1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!
!      ! Left/Right Boundary - Middle of H (\partial \Omega_i, i = 5,6)
!      iStart = NyL*Nsml; 
!      iEnd = NyL*Nsml + (NyM+1)*(Nx+2)-1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = NyL+1,NyM+1
!            ! Left Boundary
!            row = NyL*Nsml + (i - NyL - 1) * (Nx+2)
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            row = NyL*Nsml + (i - NyL) * (Nx+2) - 1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!
!      ! Left/Right Boundary - Upper Legs (\partial \Omega_i, i = 5,6)
!      iStart = NyL*Nsml + (NyM+1 - NyL)*(Nx+2)
!      iEnd = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (Ny-NyM - 1)*Nsml-1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = NyM+2, Ny
!            ! Left Boundary
!            row = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (i-NyM - 2)*Nsml
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            row = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (i-NyM - 1)*Nsml-1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!      
!      ! Inner Left/Right Walls - Upper Legs (\Gamma_i, i=1,3)
!      iStart = NyL*Nsml + (NyM+1-NyL)*(Nx+2) + NxL+1
!      iEnd   = NyL*Nsml + (NyM+1-NyL)*(Nx+2) + (Ny-NyM-2)*Nsml + NxL+2
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = NyM+2, Ny
!            ! Left Boundary
!            row = NyL*Nsml + (NyM+1-NyL)*(Nx+2) + (i-NyM-2)*Nsml + NxL+1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            row = NyL*Nsml + (NyM+1-NyL)*(Nx+2) + (i-NyM-2)*Nsml + NxL+2
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!      
!      ! Inner Left/Right Walls - Lower Legs (\Gamma_i, i=2,4)
!      iStart = Nsml + NxL + 1 
!      iEnd = (NyL-1)*Nsml + NxL+2
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = 2,NyL
!            ! Left Boundary
!            row = (i-1)*Nsml + NxL+1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            
!            ! Right Boundary
!            row = row+1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!
!
!      ! Upper Mid-H Wall (\Gamma_5)
!      iStart = NyL*Nsml + (NyM-NyL) * (Nx+2) + NxL + 1
!      iEnd   = NyL*Nsml + (NyM-NyL) * (Nx+2) + NxM
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN
!            DO j = NxL + 2, NxM + 1
!              row = NyL*Nsml + (NyM-NyL) * (Nx+2) + j - 1
!              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)
!              END IF
!            END DO
!      END IF
!
!      ! Lower Mid-H Wall (\Gamma_6)
!      iStart = NyL*Nsml + NxL+1
!      iEnd = NyL*Nsml + NxM-1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN
!            DO j = NxL+2, NxM+1
!              row = NyL*Nsml + j - 1
!              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0.D0,INSERT_VALUES,ierr)
!              END IF
!            END DO
!      END IF
!
!! --- End of Boundary Conditions ----------------------------------------

! --- Fill Inner Cells: Left Lower Legs of H ----------------------------
      iStart = Nsml+1
      iEnd = Nsml*NyL + NxL
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
      
            DO i = 2,NyL+1
            DO j = 2,NxL+1
            row = Nsml*(i-1)+j-1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell Height and Width
                dx = XU(i,j) - XU(i,j-1)
                dy = YU(i+1,j)-YU(i,j)

                ! Get Brinkman term on v-momentum cell centers
                CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                               XU,YU,YV,B,Bp)
                Bp = Bp * dx * dy / 2
                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
            END IF 
            END DO
            END DO
      END IF

! --- Fill Inner Cells: Right Lower Legs of H ----------------------------
      iStart = Nsml+NxL+3
      iEnd = Nsml*(NyL-1) + NxL+4 + (Nx-NxM) - 1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = 2,NyL
            j_ctr = NxL+3 ! Starting value for index on right of H
              DO j = NxM+2,Nx+1
                j_ctr = j_ctr + 1
                row = Nsml*(i-1) + j_ctr - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell Height and Width
                dx = XU(i,j) - XU(i,j-1)
                dy = YU(i+1,j)-YU(i,j)

                ! Get Brinkman term on v-momentum cell centers
                CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                               XU,YU,YV,B,Bp)
                Bp = Bp * dx * dy / 2
                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
                END IF 
              END DO
            END DO
      END IF

! --- Fill Inner Cells: Right at Lower Interface with Injury Channel -----
      iStart = Nsml*NyL + NxM+1
      iEnd = Nsml*NyL + (Nx+1)
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            i = NyL+1
              DO j = NxM+2,Nx+1
                row = Nsml*(i-1) + j - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell Height and Width
                dx = XU(i,j) - XU(i,j-1)
                dy = YU(i+1,j)-YU(i,j)

                ! Get Brinkman term on v-momentum cell centers
                CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                               XU,YU,YV,B,Bp)
                Bp = Bp * dx * dy / 2
                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
                END IF 
              END DO
      END IF

! --- Fill Inner Cells: Injury Channel ----------------------------
      iStart = Nsml*NyL + (Nx+2) + 1
      iEnd = Nsml*NyL + (NyM - NyL)*(Nx+2) + NxL+1 - 1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = NyL+2, NyM+1
              IF (i == NyM+1) THEN
                jEnd = NxL+1 
              ELSE 
                jEnd = Nx+1
              END IF

              DO j = 2,jEnd
                row = Nsml*NyL + (i-NyL-1)*(Nx+2) + j - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell Height and Width
                dx = XU(i,j) - XU(i,j-1)
                dy = YU(i+1,j)-YU(i,j)

                ! Get Brinkman term on v-momentum cell centers
                CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                               XU,YU,YV,B,Bp)
                Bp = Bp * dx * dy / 2
                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
                END IF 
              END DO
            END DO
      END IF

! --- Fill Inner Cells: Right at Upper Interface with Injury Channel -----
      iStart = Nsml*NyL + (NyM-NyL)*(Nx+2) + NxM+2 - 1
      iEnd   = Nsml*NyL + (NyM-NyL)*(Nx+2) + Nx
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            i = NyM+1
              DO j = NxM+2,Nx+1
                row = Nsml*NyL + (i-1-NyL)*(Nx+2) + j - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell Height and Width
                dx = XU(i,j) - XU(i,j-1)
                dy = YU(i+1,j)-YU(i,j)

                ! Get Brinkman term on v-momentum cell centers
                CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                               XU,YU,YV,B,Bp)
                Bp = Bp * dx * dy / 2
                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
                END IF 
              END DO
      END IF

! --- Fill Inner Cells: Left Upper Legs of H ----------------------------
      iStart = Nsml*NyL + (NyM+1-NyL)*(Nx+2) + 1 
      iEnd   = Nsml*NyL + (NyM+1-NyL)*(Nx+2) + (Ny-NyM-2)*Nsml + NxL
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = NyM+2, Ny
            DO j = 2,NxL+1
            row = Nsml*NyL + (NyM+1-NyL)*(Nx+2) + (i-NyM-2)*Nsml + j - 1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell Height and Width
                dx = XU(i,j) - XU(i,j-1)
                dy = YU(i+1,j)-YU(i,j)

                ! Get Brinkman term on v-momentum cell centers
                CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                               XU,YU,YV,B,Bp)
                Bp = Bp * dx * dy / 2
                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
            END IF 
            END DO
            END DO
      END IF

! --- Fill Inner Cells: Right Upper Legs of H ----------------------------
      iStart = Nsml*NyL + (NyM+1-NyL)*(Nx+2) + NxL + 2 
      iEnd = NNXY - Nsml - 2
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = NyM+2, Ny
              j_ctr = NxL+3
              DO j = NxM+2,Nx+1
              j_ctr = j_ctr + 1
              row = Nsml*NyL+(NyM+1-NyL)*(Nx+2) + (i-NyM-2)*Nsml+j_ctr-1
      
              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell Height and Width
                dx = XU(i,j) - XU(i,j-1)
                dy = YU(i+1,j)-YU(i,j)

                ! Get Brinkman term on v-momentum cell centers
                CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,&
                               XU,YU,YV,B,Bp)
                Bp = Bp * dx * dy / 2
                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
              END IF 
              END DO
            END DO
      END IF

      CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)   
      CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      END SUBROUTINE build_Bv

      END MODULE Bv_mod
