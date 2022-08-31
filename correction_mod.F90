      MODULE correction_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

      IMPLICIT NONE

CONTAINS
! -----------------------------------------------------------------------
! Subroutine: 
! -----------------------------------------------------------------------
      SUBROUTINE correction(N,dt,xp,yp,rowStartU,rowEndU,rowStartV,&
                            rowEndV,rowStartP,rowEndP,P0,U,V,P,BC_in)
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscScalar, INTENT(IN) :: dt
      PetscScalar, DIMENSION(:), INTENT(IN) :: xp, yp
      PetscInt, INTENT(IN) :: rowStartU,rowEndU,rowStartV,rowEndV,&
                             rowStartP,rowEndP
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: P0
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(INOUT) :: U
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(INOUT) :: V
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(INOUT) :: P
      PetscScalar, DIMENSION(:), INTENT(IN) :: BC_in
      ! Variables
      PetscInt :: NxL,NxM,Nx,NyL,NyM,Ny,i,j

      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);

      ! update U - Wash Channel
      DO j = 2, NxL
      DO i = MAX(2,rowStartU-1), MIN(Ny+1,rowEndU+1)
      U(i,j) = U(i,j) - dt * (P(i-1,j) - P(i-1,j-1)) / (xp(j+1)-xp(j))
      END DO
      END DO
      IF (rowStartU < 2) THEN
            U(1,1:NxL+1) = U(2,1:NxL+1); ! du/dy = 0 on bottom boundary
      END IF
      IF (rowEndU > Ny) THEN
            U(Ny+2,2:NxL) = 0 ! u = 0 on top boundary
      END IF

      ! update U - Injury Channel
      DO i = MAX(NyL+2,rowStartU-1), MIN(NyM+1,rowEndU+1)
      DO j = NxL+1, NxM+1
      U(i,j) = U(i,j) - dt * ( P(i-1,j) - P(i-1,j-1)) / (xp(j+1)-xp(j))
      END DO
      END DO

      ! update U - Blood Channel
      DO j = NxM+2, Nx
      DO i = MAX(2,rowStartU-1), MIN(Ny+1,rowEndU+1)
      U(i,j) = U(i,j) - dt * (P(i-1,j) - P(i-1,j-1)) / (xp(j+1)-xp(j) )
      END DO
      END DO
      IF (rowStartU < 2) THEN
            U(1,NxM+1:Nx+1) = U(2,NxM+1:Nx+1) ! du/dy = 0 on bottom boundary
      END IF
      IF (rowEndU > Ny) THEN
            U(Ny+2, NxM+2:Nx) = 0 ! u = 0 on top boundary
      END IF

      ! update V - Wash Channel
      DO j = 2, NxL+1
      DO i = MAX(2,rowStartV-1), MIN(Ny,rowEndV+1)
      V(i,j) = V(i,j) - dt * (P(i,j-1) - P(i-1,j-1)) / (yp(i+1)-yp(i))
      END DO
      END DO
      IF (rowStartV < 2) THEN
            V(1,1:NxL+2) = V(2,1:NxL+2) !dv/dy = 0 on bottom boundary
      END IF
      IF (rowEndV > Ny) THEN
            V(Ny+1,2:NxL+1) = BC_in(2:NxL+1) !v = inlet condition on top
      END IF

      ! update V - Injury Channel
      DO i = MAX(NyL+2,rowStartV-1), MIN(NyM,rowEndV+1)
      DO j = NxL+2, NxM+1
      V(i,j) = V(i,j) - dt * (P(i,j-1) - P(i-1,j-1)) / (yp(i+1)-yp(i))
      END DO
      END DO

      ! update V - Blood Channel
      DO j = NxM+2, Nx+1
      DO i = MAX(2,rowStartV-1), MIN(Ny,rowEndV+1)
      V(i,j) = V(i,j) - dt * (P(i,j-1) - P(i-1,j-1)) / (yp(i+1)-yp(i) )
      END DO
      END DO
      IF (rowStartV < 2) THEN
            V(1,NxM+1:Nx+2) = V(2,NxM+1:Nx+2) !dv/dy = 0 on bottom boundary
      END IF
      IF (rowEndV > Ny) THEN
            V(Ny+1,NxM+2:Nx+1) = BC_in(NxM+2:Nx+1) !v = inlet condition
      END IF

      ! Update Pressure
      P = P0 + P

      ! Zero out interior of H
      IF (rowStartP .LE. NyL) THEN
            P(rowStartP-1:MIN(NyL,rowEndP+1),NxL+1:NxM) = 0.D0
      ELSE IF (rowEndP .GE. NyM+1) THEN
            P(MAX(NyM+1,rowStartP-1):rowEndP+1,NxL+1:NxM) = 0.D0
      END IF
      END SUBROUTINE correction

      END MODULE correction_mod
