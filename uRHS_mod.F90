      MODULE uRHS_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
#include <petsc/finclude/petscvec.h>  
      USE petscsys
      USE petscvec

! Necessary F90 Modules
      USE params_mod
      USE brinkman_mod
      IMPLICIT NONE

CONTAINS
! -----------------------------------------------------------------------
! Subroutine for inner cells in RHS of u-momentum equation 
! -----------------------------------------------------------------------
      SUBROUTINE uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
      IMPLICIT NONE
      ! Inputs
      PetscInt, INTENT(IN) :: i,j
      PetscInt, INTENT(IN) :: rowStartU, rowStartV, rowStartP
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU, YU
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV, YV  
      PetscScalar, INTENT(IN) :: t,dt,Re ! current t, dt, & Re number
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: U0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: P0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B0
      PetscScalar, DIMENSION(rowStartU:,:), INTENT(IN) :: Q,Qo,NLU,NLU0
      ! Output
      PetscScalar, INTENT(OUT) :: val
      ! Variables
      PetscScalar :: dx, dy, XE, XP, XW, YN, YP, YS
      PetscScalar :: UP, UE, UW, UN, US
      PetscScalar :: dw, de, ds, dn, D
      PetscScalar :: Pe, Pw, PR, Bp

      ! Cell height/width
      dx = XV(i,j+1) - XV(i,j);
      dy = YV(i,j) - YV(i-1,j);

      ! Location of Nodes
      XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
      YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);

      ! Velocities           
      UP = U0(i,j);
      UE = U0(i,j+1); UN = U0(i+1,j); 
      UW = U0(i,j-1); US = U0(i-1,j);


      ! Diffusion Terms
      dw = - 1.D0/Re * (UP - UW)*dy/(XP-XW)/2.D0;
      de = 1.D0/Re * (UE - UP)*dy/(XE-XP)/2.D0;
      ds = -1.D0/Re * (UP - US)*dx/(YP-YS)/2.D0;
      dn = 1.D0/Re * (UN - UP)*dx/(YN-YP)/2.D0;
      D = de + dn + dw + ds;

      ! Pressure Terms
      Pe = P0(i-1,j); Pw = P0(i-1,j-1);
      PR =  (Pe-Pw)* dy;

      ! Brinkman Term
      CALL Bp_Ucoeff(i,j,rowStartU,rowStartV,rowStartP,XU,XV,YV,B0,Bp)

      val = (1.D0/dt - Bp/2)*U0(i,j)*dx*dy - DBLE(1.5)*NLU(i,j) & 
            + DBLE(0.5)*NLU0(i,j) + D - PR &
            + (Q(i,j) + Qo(i,j))*dx*dy/2.D0
      END SUBROUTINE uRHS_inner
! -----------------------------------------------------------------------
! Subroutine to build the RHS vector in u-momentum equation 
! -----------------------------------------------------------------------
      SUBROUTINE uRHS(N,NNXY,rowStartU,rowStartV,rowStartP,XU,YU,XV,YV,&
                      t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,f,ierr)
      IMPLICIT NONE
      ! Inputs
      PetscInt, DIMENSION(:), INTENT(IN) :: N !N=[NxL,NxM,Nx,NyL,NyM,Ny]
      PetscInt, INTENT(IN) :: NNXY ! Au is NNXY x NNXY matrix
      PetscInt, INTENT(IN) :: rowStartU, rowStartV, rowStartP
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU, YU
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV, YV  
      PetscScalar, INTENT(IN) :: t,dt,Re ! current t, dt, & Re number
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: U0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: P0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B0
      PetscScalar, DIMENSION(rowStartU:,:), INTENT(IN) :: Q,Qo,NLU,NLU0
      Vec, INTENT(INOUT) :: f ! vector
      PetscErrorCode, INTENT(INOUT) :: ierr
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny, Nsml
      PetscInt :: i,j,j_ctr,iStart,iEnd,jEnd,i_locStart,i_locEnd,row
      PetscScalar :: val

      
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
      
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL + 1) + (Nx + 1 - NxM)

      ! Get local start/stop index
      CALL VecGetOwnershipRange(f, i_locStart, i_locEnd, ierr)
      
! --- Boundary Conditions are all set to zero ---------------------------

! --- Fill Inner Cells: Left Lower Legs of H ----------------------------
      iStart = Nsml+1
      iEnd = Nsml*NyL + NxL-1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
      
            DO i = 2,NyL+1
            DO j = 2,NxL
            row = Nsml*(i-1)+j-1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
              CALL uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
              
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
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
              CALL uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
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
              CALL uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
            END IF 
            END DO
      END IF

! --- Fill Inner Cells: Injury Channel ----------------------------
      iStart = Nsml*NyL + (Nx+1) + 1
      iEnd = Nsml*NyL + (NyM - NyL)*(Nx+1) + Nx - 1
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
              CALL uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
            END IF 
            END DO
            END DO
      END IF

! --- Fill Inner Cells: Right at Upper Interface with Injury Channel ----
      iStart = Nsml*NyL + (NyM+2-NyL-1)*(Nx+1) + NxM+2 - 1
      iEnd = Nsml*NyL + (NyM+2-NyL-1)*(Nx+1) + Nx - 1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            i = NyM+2
            DO j = NxM+2,Nx
            row = Nsml*NyL + (i-NyL-1)*(Nx+1) + j - 1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
              CALL uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
                
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
              CALL uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr) 
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
              CALL uRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,U0,P0,B0,Q,Qo,NLU,NLU0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr) 
              
            END IF 
            END DO
            END DO
      END IF

      CALL VecAssemblyBegin(f,ierr)   
      CALL VecAssemblyEnd(f,ierr)
      END SUBROUTINE uRHS
      
      END MODULE uRHS_mod
