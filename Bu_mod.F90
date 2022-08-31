      MODULE Bu_mod
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
! Subroutine to build the Bu matrix in u-momentum equation 
! -----------------------------------------------------------------------
      SUBROUTINE build_Bu(rank,N,NNXY,locDim,rowStartU,rowStartV,&
                          rowStartP,XU,XV,YV,B,A,ierr)
      IMPLICIT NONE
      ! Inputs
      PetscInt, INTENT(IN) :: rank !rank in communicator
      PetscInt, DIMENSION(:), INTENT(IN) :: N !N=[NxL,NxM,Nx,NyL,NyM,Ny]
      PetscInt, INTENT(IN) :: rowStartU,rowStartV,rowStartP
      PetscInt, INTENT(IN) :: NNXY ! Bu is NNXY x NNXY matrix
      PetscInt, INTENT(IN) :: locDim !number of rows ea. core allocates
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU 
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV,YV 
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B
      Mat, INTENT(INOUT) :: A ! matrix
      PetscErrorCode, INTENT(INOUT) :: ierr
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny, Nsml
      PetscInt :: i,j,j_ctr,iStart,iEnd,jEnd,i_locStart,i_locEnd,row
      PetscScalar :: dx,dy,Bp

      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
      
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL + 1) + (Nx + 1 - NxM)
      
!      CALL MatSetSizes(A, locdim, locdim, NNXY, NNXY, ierr)
!      CALL MatSetFromOptions(A, ierr)
!      CALL MatSeqAIJSetPreallocation(A, 5, PETSC_NULL_INTEGER, ierr)
!      CALL MatMPIAIJSetPreallocation(A, 5, PETSC_NULL_INTEGER, 2, &
!                PETSC_NULL_INTEGER, ierr)

      CALL MatGetOwnershipRange(A, i_locStart, i_locEnd, ierr)

!! --- Boundary Conditions -----------------------------------------------
!      ! Top Boundaries (\partial \Omega_i, i = 1,3)
!      iStart = NNXY - Nsml
!      IF (i_locEnd .GE. iStart) THEN
!            DO i = MAX(i_locStart,iStart), i_locEnd-1
!                  CALL MatSetValue(A,i,i,0D0,INSERT_VALUES,ierr)
!            END DO
!      END IF
!
!      ! Bottom Boundaries (\partial \Omega_i, i = 2,4)
!      iEnd = Nsml -1
!      IF (i_locStart .LE. iEnd) THEN
!            DO i = i_locStart, MIN(i_locEnd,iEnd)
!                  CALL MatSetValue(A,i,i,0D0,INSERT_VALUES,ierr) 
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
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            row = i*Nsml - 1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!
!      ! Left/Right Boundary - Middle of H (\partial \Omega_i, i = 5,6)
!      iStart = NyL*Nsml; 
!      iEnd = NyL*Nsml + (NyM+2)*(Nx+1)-1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = NyL+1,NyM+2
!            ! Left Boundary
!            row = NyL*Nsml + (i - NyL - 1) * (Nx+1)
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            row = NyL*Nsml + (i - NyL) * (Nx+1) - 1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!
!      ! Left/Right Boundary - Upper Legs (\partial \Omega_i, i = 5,6)
!      iStart = NyL*Nsml + (NyM+2 - NyL)*(Nx+1)
!      iEnd = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (Ny+1 - NyM - 2)*Nsml-1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = NyM+3, Ny+1
!            ! Left Boundary
!            row = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (i-NyM - 3)*Nsml
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            row = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (i-NyM - 2)*Nsml-1
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!      
!      ! Inner Left/Right Walls - Lower Legs (\Gamma_i, i=2,4)
!      iStart = Nsml + NxL 
!      iEnd = NyL*Nsml + NxL+1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = 2,NyL+1
!            ! Left Boundary
!            row = (i-1)*Nsml + NxL
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            IF (i < NyL + 1) THEN
!                  row = (i-1)*Nsml + NxL + 1
!            ELSE
!                  row = (i-1)*Nsml + NxM
!            END IF
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!
!      ! Inner Left/Right Walls - Upper Legs (\Gamma_i, i=2,4)
!      iStart = NyL*Nsml + (NyM+1 - NyL)*(Nx+1) + NxL
!      iEnd = NyL*Nsml + (NyM+2-NyL)*(Nx+1) + (Ny-NyM-2)*Nsml + NxL + 2
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
!            DO i = NyM+2, Ny+1
!            ! Left Boundary
!            IF (i == NyM + 2) THEN
!              row = NyL*Nsml + (NyM+1 - NyL)*(Nx+1) + NxL 
!            ELSE
!              row = NyL*Nsml + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + NxL
!            END IF
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            ! Right Boundary
!            IF (i == NyM+2) THEN
!                  row = NyL*Nsml + (NyM+1 - NyL)*(Nx+1) + NxM
!            ELSE
!                  row = NyL*Nsml + (NyM+2-NyL)*(Nx+1) & 
!                        + (i-NyM-3)*Nsml + NxL + 1
!            END IF
!            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)  
!            END IF
!            END DO
!      END IF
!
!      ! Upper Mid-H Wall (\Gamma_5)
!      iStart = NyL*Nsml + (NyM + 1 - NyL) * (Nx+1) + NxL + 1
!      iEnd = NyL*Nsml + (NyM + 1 - NyL) * (Nx+1) + NxM - 1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN
!            DO j = NxL + 2, NxM
!              row = NyL*Nsml + (NyM + 1 - NyL) * (Nx+1) + j - 1
!              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)
!              END IF
!            END DO
!      END IF
!
!      ! Lower Mid-H Wall (\Gamma_6)
!      iStart = NyL*Nsml + NxL+1
!      iEnd = NyL*Nsml + NxM-1
!      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN
!            DO j = NxL + 2, NxM
!              row = NyL*Nsml + j - 1
!              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
!                  CALL MatSetValue(A,row,row,0D0,INSERT_VALUES,ierr)
!              END IF
!            END DO
!      END IF
!
!! --- End of Boundary Conditions ----------------------------------------

! --- Fill Inner Cells: Left Lower Legs of H ----------------------------
      iStart = Nsml+1
      iEnd = Nsml*NyL + NxL-1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
      
            DO i = 2,NyL+1
            DO j = 2,NxL
            row = Nsml*(i-1)+j-1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell height and width
                dy = YV(i,j)-YV(i-1,j)
                dx = XV(i,j+1) - XV(i,j)

                ! Get coeffiicient Bp
                CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                                           XU,XV,YV,B,Bp)
                Bp = dx*dy*Bp/2D0

                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)
            END IF 
            END DO
            END DO
      END IF

! --- Fill Inner Cells: Right Lower Legs of H ----------------------------
      iStart = Nsml+NxL+2
      iEnd = Nsml*NyL-2
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = 2,NyL
            j_ctr = NxL+2 ! Starting value for index on right of H
              DO j = NxM+2,Nx
                j_ctr = j_ctr + 1
                row = Nsml*(i-1) + j_ctr - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell height and width
                dy = YV(i,j)-YV(i-1,j)
                dx = XV(i,j+1) - XV(i,j)

                ! Get coeffiicient Bp
                CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                                           XU,XV,YV,B,Bp)
                Bp = dx*dy*Bp/2D0

                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)
                END IF 
              END DO
            END DO
      END IF

! --- Fill Inner Cells: Right at Lower Interface with Injury Channel -----
      iStart = Nsml*NyL + NxM+1
      iEnd = Nsml*NyL + (Nx)
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            i = NyL+1
              DO j = NxM+2,Nx
                row = Nsml*(i-1) + j - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell height and width
                dy = YV(i,j)-YV(i-1,j)
                dx = XV(i,j+1) - XV(i,j)

                ! Get coeffiicient Bp
                CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                                           XU,XV,YV,B,Bp)
                Bp = dx*dy*Bp/2D0

                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
                END IF 
              END DO
      END IF

! --- Fill Inner Cells: Injury Channel ----------------------------
      iStart = Nsml*NyL + (Nx+1) + 1
      iEnd = Nsml*NyL + (NyM+1 - NyL)*(Nx+1) + Nx - 1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = NyL+2, NyM+2
              IF (i == NyM+2) THEN
                jEnd = NxL  
              ELSE 
                jEnd = Nx
              END IF

              DO j = 2,jEnd
                row = Nsml*NyL + (i-NyL-1)*(Nx+1) + j - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell height and width
                dy = YV(i,j)-YV(i-1,j)
                dx = XV(i,j+1) - XV(i,j)

                ! Get coeffiicient Bp
                CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                                           XU,XV,YV,B,Bp)
                Bp = dx*dy*Bp/2D0

                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
                END IF 
              END DO
            END DO
      END IF

! --- Fill Inner Cells: Right at Upper Interface with Injury Channel -----
      iStart = Nsml*NyL + (NyM+2-NyL-1)*(Nx+1) + NxM+2 - 1
      iEnd = Nsml*NyL + (NyM+2-NyL-1)*(Nx+1) + Nx - 1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            i = NyM+2
              DO j = NxM+2,Nx
                row = Nsml*NyL + (i-NyL-1)*(Nx+1) + j - 1
      
                IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell height and width
                dy = YV(i,j)-YV(i-1,j)
                dx = XV(i,j+1) - XV(i,j)

                ! Get coeffiicient Bp
                CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                                           XU,XV,YV,B,Bp)
                Bp = dx*dy*Bp/2D0

                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
                END IF 
              END DO
      END IF

! --- Fill Inner Cells: Left Upper Legs of H ----------------------------
      iStart = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + 1 
      iEnd = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (Ny-NyM-2)*Nsml + NxL - 1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = NyM+3, Ny+1
            DO j = 2,NxL
            row = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j - 1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell height and width
                dy = YV(i,j)-YV(i-1,j)
                dx = XV(i,j+1) - XV(i,j)

                ! Get coeffiicient Bp
                CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                                           XU,XV,YV,B,Bp)
                Bp = dx*dy*Bp/2D0

                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
            END IF 
            END DO
            END DO
      END IF

! --- Fill Inner Cells: Right Upper Legs of H ----------------------------
      iStart = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + NxL + 1 
      iEnd = NNXY - Nsml - 2
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            DO i = NyM+3, Ny+1
              j_ctr = NxL+2
              DO j = NxM+2,Nx
              j_ctr = j_ctr + 1
              row = Nsml*NyL+(NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml+j_ctr-1
      
              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
                ! Cell height and width
                dy = YV(i,j)-YV(i-1,j)
                dx = XV(i,j+1) - XV(i,j)

                ! Get coeffiicient Bp
                CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,&
                                           XU,XV,YV,B,Bp)
                Bp = dx*dy*Bp/2D0

                ! Principal
                CALL MatSetValue(A,row,row,Bp,INSERT_VALUES,ierr)  
              END IF 
              END DO
            END DO
      END IF

      CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)   
      CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      END SUBROUTINE build_Bu

      END MODULE Bu_mod
