      MODULE reshape_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
#include <petsc/finclude/petscvec.h>  
      USE petscsys
      USE petscvec
! Necessary F90 Modules
      USE params_mod
      
      IMPLICIT NONE

CONTAINS
! -------------------------------------------------------------------------
! Subroutine that reshapes the U solution vector into a 2D array
! -------------------------------------------------------------------------
      SUBROUTINE reshapeU(N,locStart,locEnd,sysStart,u_vecF90,U)
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscInt, INTENT(IN) :: locStart, locEnd, sysStart
      PetscScalar, DIMENSION(sysStart+1:) :: u_vecF90
      PetscScalar, DIMENSION(locStart-1:,:), INTENT(OUT) :: U
      PetscInt :: Nsml,NxL,NxM,Nx,NyL,NyM,Ny,i, j, row 
      
      U = 0.D0
      
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
    
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL + 1) + (Nx + 1 - NxM)

      ! Lower Legs of H
      IF (locStart .LE. NyL) THEN
            ! Lower Legs - Left
            DO j = 1, NxL+1
            DO i = MAX(1,locStart), MIN(NyL,locEnd)
                row = j + (i-1)*Nsml
                U(i,j) = u_vecF90( row )
            END DO
            END DO

            ! Lower Legs - Right
            DO i = MAX(1,locStart), MIN(NyL, locEnd)
                row = NxL+1 + (i-1)*Nsml
            DO j = NxM+1, Nx+1
                row = row + 1 
                U(i,j) = u_vecF90( row )
            END DO
            END DO
      END IF

      ! Middle of H
      IF ((locStart .LE. NyM+2) .AND. (locEnd .GE. NyL+1)) THEN
            DO j = 1, Nx+1
            DO i = MAX(NyL+1,locStart), MIN(NyM+2,locEnd)
                row = NyL*Nsml + j + (i-NyL-1)*(Nx+1)
                U(i,j) = u_vecF90( row )
            END DO
            END DO
      END IF

      ! Upper Legs of H
      IF (locEnd .GE. NyM+3) THEN
            ! Upper Legs - Left
            DO j = 1, NxL+1
            DO i = MAX(NyM+3,locStart), MIN(Ny+2,locEnd)
                row = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) &
                      + j + (i-NyM-3)*Nsml
                U(i,j) = u_vecF90( row )
            END DO
            END DO

            ! Upper Legs - Right
            DO i = MAX(NyM+3,locStart), MIN(Ny+2,locEnd)
                row = NxL+1 + NyL*Nsml + (NyM+2 - NyL)*(Nx+1) & 
                      + (i-NyM-3)*Nsml
            DO j = NxM+1, Nx+1
                row = row + 1
                U(i,j) = u_vecF90( row )
            END DO
            END DO
      END IF

      END SUBROUTINE reshapeU

! -------------------------------------------------------------------------
! Subroutine that reshapes the V solution vector into a 2D array
! -------------------------------------------------------------------------
      SUBROUTINE reshapeV(N,locStart,locEnd,sysStart,v_vecF90,V)
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscInt, INTENT(IN) :: locStart, locEnd, sysStart
      PetscScalar, DIMENSION(sysStart+1:) :: v_vecF90
      PetscScalar, DIMENSION(locStart-1:,:), INTENT(OUT) :: V
      PetscInt :: Nsml,NxL,NxM,Nx,NyL,NyM,Ny,i, j, row 
      
      V = 0.D0
      
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
    
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL + 2) + (Nx + 2 - NxM)

      ! Lower Legs of H
      IF (locStart .LE. NyL) THEN
            ! Lower Legs - Left
            DO j = 1, NxL+2
            DO i = MAX(1,locStart), MIN(NyL,locEnd)
                row = j + (i-1)*Nsml
                V(i,j) = v_vecF90( row )
            END DO
            END DO

            ! Lower Legs - Right
            DO i = MAX(1,locStart), MIN(NyL, locEnd)
                row = NxL+2 + (i-1)*Nsml
            DO j = NxM+1, Nx+2
                row = row + 1 
                V(i,j) = v_vecF90( row )
            END DO
            END DO
      END IF

      ! Middle of H
      IF ((locStart .LE. NyM+1) .AND. (locEnd .GE. NyL+1)) THEN
            DO j = 1, Nx+2
            DO i = MAX(NyL+1,locStart), MIN(NyM+1,locEnd)
                row = NyL*Nsml + j + (i-NyL-1)*(Nx+2)
                V(i,j) = v_vecF90( row )
            END DO
            END DO
      END IF

      ! Upper Legs of H
      IF (locEnd .GE. NyM+2) THEN
            ! Upper Legs - Left
            DO j = 1, NxL+2
            DO i = MAX(NyM+2,locStart), MIN(Ny+1,locEnd)
                row = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) &
                      + j + (i-NyM-2)*Nsml
                V(i,j) = v_vecF90( row )
            END DO
            END DO

            ! Upper Legs - Right
            DO i = MAX(NyM+2,locStart), MIN(Ny+1,locEnd)
                row = NxL+2 + NyL*Nsml + (NyM+1 - NyL)*(Nx+2) & 
                      + (i-NyM-2)*Nsml
            DO j = NxM+1, Nx+2
                row = row + 1
                V(i,j) = v_vecF90( row )
            END DO
            END DO
      END IF

      END SUBROUTINE reshapeV

! -------------------------------------------------------------------------
! Subroutine that reshapes the PHI solution vector into a 2D array
! -------------------------------------------------------------------------
      SUBROUTINE reshapePHI(N,locStart,locEnd,sysStart,p_vecF90,P)
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscInt, INTENT(IN) :: locStart, locEnd, sysStart
      PetscScalar, DIMENSION(sysStart+1:) :: p_vecF90
      PetscScalar, DIMENSION(locStart-2:,:), INTENT(OUT) :: P
      PetscInt :: Nsml,NxL,NxM,Nx,NyL,NyM,Ny,i, j, row 
      
      P = 0.D0
      
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
    
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL) + (Nx - NxM)

      ! Lower Legs of H
      IF (locStart .LE. NyL) THEN
            ! Lower Legs - Left
            DO j = 1, NxL
            DO i = MAX(1,locStart), MIN(NyL,locEnd)
                row = j + (i-1)*Nsml
                P(i,j) = p_vecF90( row )
            END DO
            END DO

            ! Lower Legs - Right
            DO i = MAX(1,locStart), MIN(NyL, locEnd)
                row = NxL + (i-1)*Nsml
            DO j = NxM+1, Nx
                row = row + 1 
                P(i,j) = p_vecF90( row )
            END DO
            END DO
      END IF

      ! Middle of H
      IF ((locStart .LE. NyM) .AND. (locEnd .GE. NyL+1)) THEN
            DO j = 1, Nx
            DO i = MAX(NyL+1,locStart), MIN(NyM,locEnd)
                row = NyL*Nsml + j + (i-NyL-1)*(Nx)
                P(i,j) = p_vecF90( row )
            END DO
            END DO
      END IF

      ! Upper Legs of H
      IF (locEnd .GE. NyM+1) THEN
            ! Upper Legs - Left
            DO j = 1, NxL
            DO i = MAX(NyM+1,locStart), MIN(Ny,locEnd)
                row = NyL*Nsml + (NyM - NyL)*(Nx) &
                      + j + (i-NyM-1)*Nsml
                P(i,j) = p_vecF90( row )
            END DO
            END DO

            ! Upper Legs - Right
            DO i = MAX(NyM+1,locStart), MIN(Ny,locEnd)
                row = NxL + NyL*Nsml + (NyM - NyL)*(Nx) & 
                      + (i-NyM-1)*Nsml
            DO j = NxM+1, Nx
                row = row + 1
                P(i,j) = p_vecF90( row )
            END DO
            END DO
      END IF

      END SUBROUTINE reshapePHI

      END MODULE reshape_mod
