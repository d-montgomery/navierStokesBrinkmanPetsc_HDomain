      MODULE vRHS_mod
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
! Subroutine for the inlet BC of the v-momentum equation
! h = h1 + h3, where h1 and h3 are parabolic flows from inlet 1 and 3
! -----------------------------------------------------------------------
      SUBROUTINE inletBC_v(XV,h)
      PetscScalar, DIMENSION(:), INTENT(IN) :: XV
      PetscScalar, DIMENSION(:), INTENT(OUT) :: h 
      PetscScalar :: alpha, beta, W
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: h1, h3

      ALLOCATE(h1(SIZE(XV)), h3(SIZE(XV)))
      alpha = 6.D0*Q1/(Xchar*WL)**3
      beta = 6.D0*Q3/(Xchar*WR)**3
      W = WL + Wmid + WR

      ! Top Left Inlet \partial \omega_1
      h1 = MIN( alpha*(Xchar*XV - Xchar*WL/2.D0)**2 &
                 - alpha*(Xchar*WL/2)**2 , 0.D0) / Uchar
      ! Top Right Inlet \partial \omega_3
      h3 = MIN( beta*(Xchar*XV - Xchar*W + Xchar*WR/2.D0)**2 &
                     - beta*(Xchar*WR/2)**2, 0.D0) / Uchar 
      h = h1 + h3
      
      DEALLOCATE(h1,h3)
      END SUBROUTINE inletBC_v

! -----------------------------------------------------------------------
! Subroutine for inner cells in RHS of v-momentum equation 
! -----------------------------------------------------------------------
      SUBROUTINE vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
      IMPLICIT NONE
      ! Inputs
      PetscInt, INTENT(IN) :: i,j
      PetscInt, INTENT(IN) :: rowStartU, rowStartV, rowStartP
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU, YU
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV, YV  
      PetscScalar, INTENT(IN) :: t,dt,Re ! current t, dt, & Re number
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: V0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: P0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B0
      PetscScalar, DIMENSION(rowStartV:,:), INTENT(IN) :: Q,Qo,NLV,NLV0
      ! Output
      PetscScalar, INTENT(OUT) :: val
      ! Variables
      PetscScalar :: dx, dy, XE, XP, XW, YN, YP, YS
      PetscScalar :: VP, VE, VW, VN, VS
      PetscScalar :: dw, de, ds, dn, D
      PetscScalar :: Pn, Ps, PR, Bp

      ! Cell height/width
      dx = XU(i,j) - XU(i,j-1)
      dy = YU(i+1,j)-YU(i,j)

      ! Location of Nodes
      XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
      YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);

      ! Velocities
      VP = V0(i,j)
      VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);

      ! Diffusion Terms
      dw = -1.D0/Re * ( VP - VW )*dy/(XP-XW)/2.D0
      de = 1.D0/Re * ( VE - VP )*dy/(XE-XP)/2.D0
      ds = -1.D0/Re * ( VP - VS )*dx/(YP-YS)/2.D0
      dn = 1.D0/Re * ( VN - VP )*dx/(YN-YP)/2.D0
      D = dw + de + ds + dn

      ! Pressure
      Pn = P0(i,j-1); Ps = P0(i-1,j-1)
      PR = (Pn-Ps) * dx

      ! Brinkman Term
      CALL Bp_Vcoeff(i,j,rowStartU,rowStartV,rowStartP,XU,YU,YV,B0,Bp)

      val = (1.D0/dt - Bp/2.D0)*V0(i,j)*dx*dy - DBLE(1.5)*NLV(i,j) &
                              + DBLE(0.5)*NLV0(i,j) + D - PR &
                              + (Q(i,j) + Qo(i,j))*dx*dy/2.D0
      END SUBROUTINE vRHS_inner
! -----------------------------------------------------------------------
! Subroutine to build the RHS vector in v-momentum equation 
! -----------------------------------------------------------------------
      SUBROUTINE vRHS(N,NNXY,rowStartU,rowStartV,rowStartP,XU,YU,XV,YV,&
                      t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,BC,f,ierr)
      IMPLICIT NONE
      ! Inputs
      PetscInt, DIMENSION(:), INTENT(IN) :: N !N=[NxL,NxM,Nx,NyL,NyM,Ny]
      PetscInt, INTENT(IN) :: NNXY ! Au is NNXY x NNXY matrix
      PetscInt, INTENT(IN) :: rowStartU, rowStartV, rowStartP
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU, YU
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV, YV  
      PetscScalar, INTENT(IN) :: t,dt,Re ! current t, dt, & Re number
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: V0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: P0
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(IN) :: B0
      PetscScalar, DIMENSION(rowStartV:,:), INTENT(IN) :: Q,Qo,NLV,NLV0
      PetscScalar, DIMENSION(:), INTENT(IN) :: BC
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
      Nsml = (NxL + 2) + (Nx + 2 - NxM)

      ! Get local start/stop index
      CALL VecGetOwnershipRange(f, i_locStart, i_locEnd, ierr)
      i_locEnd = i_locEnd-1

! --- Boundary Conditions -----------------------------------------------
      ! Top Left BC (\partial \Omega_1)
      iStart = NNXY - Nsml
      IF (i_locEnd .GE. iStart) THEN
            DO j = 1, NxL+2
              row = NNXY - Nsml + j - 1
              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
                CALL VecSetValues(f,1,row,BC(j),INSERT_VALUES,ierr)
              END IF
            END DO
      END IF
      
!      ! Bottom Left BC (\partial \Omega_2)
!      iEnd = NxL+1
!      IF (i_locStart .LE. iEND) THEN
!            DO row = i_locStart, MIN(i_locEnd,iEnd)
!              j = row+1
!              CALL VecSetValues(f,1,row,0.D0,INSERT_VALUES,ierr)
!            END DO
!      END IF

      ! Top Right BC (\partial \Omega_3)
      iStart = NNXY - Nsml + NxL+1
      IF (i_locEnd .GE. iStart) THEN
            j_ctr = NxL + 2
            DO j = NxM+1, Nx+2
              j_ctr = j_ctr + 1
              row = NNXY - Nsml + j_ctr - 1
              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEnd)) THEN
                CALL VecSetValues(f,1,row,BC(j),INSERT_VALUES,ierr)
              END IF
            END DO
      END IF

!      ! Bottom Right BC (\partial \Omega_4)
!      iStart = NxL+2
!      iEnd = iStart + Nx+1 - NxM
!      IF ((i_locStart .LE. iEND) .AND. (i_locEnd .GE. iStart)) THEN
!            row = NxL+1
!            DO j = NxM+1, Nx+2
!              row = row + 1
!              IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
!                CALL VecSetValues(f,1,row,0.D0,INSERT_VALUES,ierr)
!              END IF
!            END DO
!      END IF

! --- End of Boundary Conditions ----------------------------------------

! --- Fill Inner Cells: Left Lower Legs of H ----------------------------
      iStart = Nsml+1
      iEnd = Nsml*NyL + NxL
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
      
            DO i = 2,NyL+1
            DO j = 2,NxL+1
            row = Nsml*(i-1)+j-1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
              CALL vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
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
              CALL vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
                END IF 
              END DO
            END DO
      END IF

! --- Fill Inner Cells: Right at Lower Interface with Injury Channel -----
      iStart = Nsml*NyL + NxM+1
      iEnd = Nsml*NyL + (Nx+1)-1
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            i = NyL+1
            DO j = NxM+2,Nx+1
            row = Nsml*(i-1) + j - 1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
              CALL vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
            END IF 
            END DO
      END IF

! --- Fill Inner Cells: Injury Channel ----------------------------
      iStart = Nsml*NyL + (Nx+2) + 1
      iEnd = Nsml*NyL + (NyM - NyL)*(Nx+2) + NxL + 1
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
              CALL vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
            END IF 
            END DO
            END DO
      END IF

! --- Fill Inner Cells: Right at Upper Interface with Injury Channel ----
      iStart = Nsml*NyL + (NyM-NyL)*(Nx+2) + NxM+2 - 1
      iEnd   = Nsml*NyL + (NyM-NyL)*(Nx+2) + Nx 
      IF ((i_locEnd .GE. iStart) .AND. (i_locStart .LE. iEnd)) THEN 
            i = NyM+1
            DO j = NxM+2,Nx+1
            row = Nsml*NyL + (i-1-NyL)*(Nx+2) + j - 1
      
            IF ((row .GE. i_locStart) .AND. (row .LE. i_locEND)) THEN
              CALL vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr)
                
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
              CALL vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr) 
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
              CALL vRHS_inner(i,j,rowStartU,rowStartV,rowStartP,&
                        XU,YU,XV,YV,t,dt,Re,V0,P0,B0,Q,Qo,NLV,NLV0,val)
              CALL VecSetValues(f,1,row,val,INSERT_VALUES,ierr) 
              
            END IF 
            END DO
            END DO
      END IF

      CALL VecAssemblyBegin(f,ierr)   
      CALL VecAssemblyEnd(f,ierr)
      END SUBROUTINE vRHS
      
      END MODULE vRHS_mod
