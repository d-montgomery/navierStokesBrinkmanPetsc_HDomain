      MODULE indexing_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

      IMPLICIT NONE

CONTAINS
! -----------------------------------------------------------------------
! Subroutine: Round-robin distribution of n values to P processors
! -----------------------------------------------------------------------
      SUBROUTINE dim_dist(n,P,locDim)
      PetscInt, INTENT(IN) :: n, P
      PetscInt, DIMENSION(0:), INTENT(OUT) :: locDim
      PetscInt :: rem, i
      
      locDim = n/P
      rem = MOD(n,P)
      DO i = 0, rem - 1
            locDim(i) = locDim(i) + 1
      END DO

      END SUBROUTINE dim_dist

! -----------------------------------------------------------------------
! Subroutine: Determine start/stop indexing given a vector locDim of 
!             a distribution of n points to P processors
! -----------------------------------------------------------------------
      SUBROUTINE ind_dist(n,P,rank,locDim,rowStart,rowStop)
      PetscInt, INTENT(IN) :: n, P, rank
      PetscInt, DIMENSION(0:), INTENT(IN) :: locDim
      PetscInt, INTENT(OUT) :: rowStart, rowStop
      PetscInt :: i 

      rowStart = 1; rowStop = locDim(0)
      DO i = 1, rank
            rowStart = rowStart + locDim(i-1)
            rowStop = rowStart + locDim(i)-1
      END DO
      END SUBROUTINE ind_dist

! -----------------------------------------------------------------------
! Subroutine: Distribution of Cells for U
! -----------------------------------------------------------------------
      SUBROUTINE u_indx(P,k,NyL,NyM,Ny,dimU,rowStartU,rowEndU)
      PetscInt, INTENT(IN) :: P, k, NyL, NyM, Ny
      PetscInt, DIMENSION(:), INTENT(OUT) :: dimU, rowStartU, rowEndU
      ! Variables
      PetscInt :: i
      PetscInt, ALLOCATABLE, DIMENSION(:) :: dimLowU, dimLowMidU, &
                                             dimUpMidU, dimUpperU
      IF (k > 0) THEN
      
      ALLOCATE(dimLowU(k), dimLowMidU(k), dimUpMidU(k), dimUpperU(k) )
            CALL dim_dist(NyL,k,dimLowU)
            CALL dim_dist(FLOOR(REAL(NyM+2-NyL)/2.0),k,dimLowMidU)
            CALL dim_dist(CEILING(REAL(NyM+2-NyL)/2.0),k,dimUpMidU)
            CALL dim_dist(Ny-NyM, k, dimUpperU)
            dimU = (/ dimLowU, dimLowMidU, dimUpMidU, dimUpperU /)
      DEALLOCATE(dimLowU, dimLowMidU, dimUpMidU, dimUpperU )
      
      ELSE
            dimU = Ny+2
      END IF

      rowStartU(1) = 1
      DO i = 2,P
          rowStartU(i) = rowStartU(i-1) + dimU(i-1)
      END DO

      DO i = 1,P-1
         rowEndU(i) = rowStartU(i+1) - 1
      END DO
      rowEndU(P) = Ny+2
      
      END SUBROUTINE u_indx

! -----------------------------------------------------------------------
! Subroutine: Distribution of Cells for V
! -----------------------------------------------------------------------
      SUBROUTINE v_indx(P,k,NyL,NyM,Ny,dimV,rowStartV,rowEndV)
      PetscInt, INTENT(IN) :: P, k, NyL, NyM, Ny
      PetscInt, DIMENSION(:), INTENT(OUT) :: dimV, rowStartV, rowEndV
      ! Variables
      PetscInt :: i
      PetscInt, ALLOCATABLE, DIMENSION(:) :: dimLowV, dimLowMidV, &
                                             dimUpMidV, dimUpperV
      IF (k > 0) THEN
      
      ALLOCATE(dimLowV(k), dimLowMidV(k), dimUpMidV(k), dimUpperV(k) )
            CALL dim_dist(NyL,k,dimLowV)
            CALL dim_dist(FLOOR(REAL(NyM+1-NyL)/2.0),k,dimLowMidV)
            CALL dim_dist(CEILING(REAL(NyM+1-NyL)/2.0),k,dimUpMidV)
            CALL dim_dist(Ny-NyM,k,dimUpperV)
            dimV = (/ dimLowV, dimLowMidV, dimUpMidV, dimUpperV /)
      DEALLOCATE(dimLowV, dimLowMidV, dimUpMidV, dimUpperV )
      
      ELSE
            dimV = Ny+1
      END IF

      rowStartV(1) = 1
      DO i = 2,P
          rowStartV(i) = rowStartV(i-1) + dimV(i-1)
      END DO

      DO i = 1,P-1
         rowEndV(i) = rowStartV(i+1) - 1
      END DO
      rowEndV(P) = Ny+1
      
      END SUBROUTINE v_indx

! -----------------------------------------------------------------------
! Subroutine: Distribution of Cells for P
! -----------------------------------------------------------------------
      SUBROUTINE p_indx(P,k,NyL,NyM,Ny,dimP,rowStartP,rowEndP)
      PetscInt, INTENT(IN) :: P, k, NyL, NyM, Ny
      PetscInt, DIMENSION(:), INTENT(OUT) :: dimP, rowStartP, rowEndP
      ! Variables
      PetscInt :: i
      PetscInt, ALLOCATABLE, DIMENSION(:) :: dimLowP, dimLowMidP, &
                                             dimUpMidP, dimUpperP
      IF (k > 0) THEN
      
      ALLOCATE(dimLowP(k), dimLowMidP(k), dimUpMidP(k), dimUpperP(k) ) 
            CALL dim_dist(NyL-1,k,dimLowP)
            CALL dim_dist(FLOOR(REAL(NyM+2-NyL)/2.0),k,dimLowMidP)
            CALL dim_dist(CEILING(REAL(NyM+2-NyL)/2.0),k,dimUpMidP)
            CALL dim_dist(Ny-NyM-1,k,dimUpperP)
            dimP = (/ dimLowP, dimLowMidP, dimUpMidP, dimUpperP /)
      DEALLOCATE(dimLowP, dimLowMidP, dimUpMidP, dimUpperP ) 
      
      ELSE
            dimP = Ny
      END IF

      rowStartP(1) = 1
      DO i = 2,P
          rowStartP(i) = rowStartP(i-1) + dimP(i-1)
      END DO

      DO i = 1,P-1
         rowEndP(i) = rowStartP(i+1) - 1
      END DO
      rowEndP(P) = Ny
      
      END SUBROUTINE p_indx

! -----------------------------------------------------------------------
! Subroutine: Distribution of rows for U & V Linear Systems
! -----------------------------------------------------------------------
      SUBROUTINE dim_LinSysUV(Nsml,Nbig,Nsys,k,locDim,locDimSys)
      PetscInt, INTENT(IN) :: Nsml,Nbig,Nsys,k
      PetscInt, DIMENSION(0:), INTENT(IN) :: locDim
      PetscInt, DIMENSION(0:), INTENT(OUT) :: locDimSys
      PetscInt :: j
      
      IF (k == 0) THEN 
            locDimSys = Nsys
      ELSE
            DO j = 0, k-1
                  locDimSys(j) = locDim(j)*Nsml         !Lower Legs
                  locDimSys(k+j) = locDim(k+j)*Nbig     !Lower Injury
                  locDimSys(2*k+j) = locDim(2*k+j)*Nbig !Upper Injury
                  locDimSys(3*k+j) = locDim(3*k+j)*Nsml !Upper Legs
            END DO
      END IF
      END SUBROUTINE dim_LinSysUV

! -----------------------------------------------------------------------
! Subroutine: Distribution of rows for P Linear System
! -----------------------------------------------------------------------
      SUBROUTINE dim_LinSysP(Nsml,Nbig,Nsys,k,locDim,locDimSys)
      PetscInt, INTENT(IN) :: Nsml,Nbig,Nsys,k
      PetscInt, DIMENSION(0:), INTENT(IN) :: locDim
      PetscInt, DIMENSION(0:), INTENT(OUT) :: locDimSys
      PetscInt :: j, val1, val2
      
      IF (k == 0) THEN 
            locDimSys = Nsys
      ELSE
            DO j = 0, k-1
                  !Lower Legs
                  locDimSys(j) = locDim(j)*Nsml       
                  
                  ! Lower Injury Channel
                  IF (j == 0) THEN
                    val1 = Nsml
                    val2 = Nbig*(locDim(k+j)-1)
                    locDimSys(k+j) = val1 + val2     !Lower Injury interface
                  ELSE
                    locDimSys(k+j) = locDim(k+j)*Nbig !Lower Injury
                  END IF
                  
                  ! Upper Injury Channel
                  IF (j == k-1) THEN
                    val1 = Nsml
                    val2 = Nbig*(locDim(2*k+j)-1)
                    locDimSys(2*k+j) = val1 + val2  !Upper Injury Interface
                  ELSE
                    locDimSys(2*k+j) = locDim(2*k+j)*Nbig !Upper Injury
                  END IF
                  
                  !Upper Legs
                  locDimSys(3*k+j) = locDim(3*k+j)*Nsml 
            END DO
      END IF
      END SUBROUTINE dim_LinSysP

      END MODULE indexing_mod
