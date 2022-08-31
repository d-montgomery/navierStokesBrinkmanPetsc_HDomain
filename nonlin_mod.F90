      MODULE nonlin_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys
! Necessary .F90 Modules
      USE params_mod
      IMPLICIT NONE
CONTAINS

! ------------------------------------------------------------------------
! Compute Nonlinear Terms for U-momentum Equation
! ------------------------------------------------------------------------
      SUBROUTINE nonlin_U(N,rowStartU,rowEndU,rowStartV,XU,YU,XV,YV,&
                          U,V,NLU)
      ! Inputs
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscInt, INTENT(IN) :: rowStartU, rowEndU, rowStartV
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU,YU,U
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV,YV,V
      PetscScalar, DIMENSION(rowStartU:,:), INTENT(OUT) :: NLU
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny, i, j
      PetscScalar :: dy,ysv,ynv,YN,YP,YS
      PetscScalar :: XE,XP,XW
      PetscScalar :: ue,un,uw,us,me,mn,mw,ms
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);

      ! --- Compute NLU for x-momentum equation --------------------------
      NLU = 0.D0 

      ! Wash Channel
      DO j = 2, NxL 
      DO i = MAX(2,rowStartU), MIN(Ny+1,rowEndU)
          dy = YV(i,j) - YV(i-1,j);

          ! Top and Bottom Edges of cell
          ysv = YV(i-1,j); ynv = YV(i,j);

          ! Location of Nodes
          XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
          YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);


          ! East
          ue = 0.5*(U(i,j+1) + U(i,j));
          me = 0.5 * dy * (U(i,j) + U(i,j+1));

          ! North
          IF (i == Ny+1) THEN 
              un = U(i+1,j); ! Don't interpolate at boundary
          ELSE
              un = U(i,j) + (ynv - YP) * (U(i+1,j) - U(i,j)) / (YN - YP)
          END IF
          mn = 0.5 * ( V(i,j) *(XP - XW) + V(i,j+1) *(XE - XP) );


          ! West
          uw = 0.5*(U(i,j) + U(i,j-1));
          mw = 0.5 * dy * (U(i,j) + U(i,j-1));

          ! South
          IF (i == 2) THEN
              us = U(i-1,j); ! Don't interpolate at boundary
          ELSE
              us = U(i,j) + (ysv-YP) * ( U(i-1,j) - U(i,j) ) / (YS-YP);
          END IF
          ms = 0.5 * ( V(i-1,j)*(XP-XW) + V(i-1,j+1)*(XE-XP) );

          NLU(i,j) = me*ue +mn*un -mw*uw - ms*us;
      END DO
      END DO

      ! Injury Channel
      IF ((rowEndU .GE. NyL+2) .AND. (rowStartU .LE. NyM+1)) THEN
      DO j = NxL+1, NxM+1
      DO i = MAX(NyL+2,rowStartU), MIN(NyM+1,rowEndU)
          dy = YV(i,j) - YV(i-1,j);

          ! Top and Bottom Edges of cell
          ysv = YV(i-1,j); ynv = YV(i,j);

          ! Location of Nodes
          XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
          YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);


          ! East
          ue = 0.5*(U(i,j+1) + U(i,j));
          me = 0.5 * dy * (U(i,j) + U(i,j+1));

          ! North
          IF ((i == NyM+1) .AND. (j == NxL+1)) THEN
              un = U(i,j) + (ynv - YP) * (U(i+1,j) - U(i,j)) / (YN - YP)
          ELSE IF ((i == NyM+1) .AND. (j == NxM+1)) THEN
              un = U(i,j) + (ynv - YP) * (U(i+1,j) - U(i,j)) / (YN - YP)
          ELSE IF ( i == NyM+1 ) THEN
                  un = U(i+1,j); ! Don't interpolate at boundary
          ELSE
              un = U(i,j) + (ynv - YP) * (U(i+1,j) - U(i,j)) / (YN - YP)
          END IF
          mn = 0.5 * ( V(i,j) *(XP - XW) + V(i,j+1) *(XE - XP) );


          ! West
          uw = 0.5*(U(i,j) + U(i,j-1));
          mw = 0.5 * dy * (U(i,j) + U(i,j-1));

          ! South
          IF (i == NyL+2) THEN
              us = U(i-1,j); ! Don't interpolate at boundary
          ELSE
              us = U(i,j) + (ysv-YP) * (U(i-1,j) - U(i,j)) / (YS-YP);
          END IF
          ms = 0.5 * ( V(i-1,j)*(XP-XW) + V(i-1,j+1)*(XE-XP) );

          NLU(i,j) = me*ue +mn*un -mw*uw - ms*us;
      END DO
      END DO
      END IF

      ! Blood Channel
      DO j = NxM+2, Nx
      DO i = MAX(2,rowStartU), MIN(Ny+1,rowEndU)
          dy = YV(i,j) - YV(i-1,j);

          ! Top and Bottom Edges of cell
          ysv = YV(i-1,j); ynv = YV(i,j);

          ! Location of Nodes
          XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
          YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);


          ! East
          ue = 0.5*(U(i,j+1) + U(i,j));
          me = 0.5 * dy * (U(i,j) + U(i,j+1));

          ! North
          IF (i == Ny+1) THEN 
              un = U(i+1,j); ! Don't interpolate at boundary
          ELSE
              un = U(i,j) + (ynv - YP) * (U(i+1,j) - U(i,j)) / (YN - YP)
          END IF
          mn = 0.5 * ( V(i,j) *(XP - XW) + V(i,j+1) *(XE - XP) );


          ! West
          uw = 0.5*(U(i,j) + U(i,j-1));
          mw = 0.5 * dy * (U(i,j) + U(i,j-1));

          ! South
          IF (i == 2) THEN
              us = U(i-1,j); ! Don't interpolate at boundary
          ELSE
              us = U(i,j) + (ysv-YP) * ( U(i-1,j) - U(i,j) ) / (YS-YP);
          END IF
          ms = 0.5 * ( V(i-1,j)*(XP-XW) + V(i-1,j+1)*(XE-XP) );

          NLU(i,j) = me*ue +mn*un -mw*uw - ms*us;
      END DO
      END DO
      
      END SUBROUTINE nonlin_U

! -------------------------------------------------------------------------
! Compute Nonlinear Terms for V-momentum Equation
! -------------------------------------------------------------------------
      SUBROUTINE nonlin_V(N,rowStartU,rowStartV,rowEndV,XU,YU,XV,YV,&
                          U,V,NLV)
      ! Inputs
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscInt, INTENT(IN) :: rowStartU, rowStartV, rowEndV
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU,YU,U
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV,YV,V
      PetscScalar, DIMENSION(rowStartV:,:), INTENT(OUT) :: NLV
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny, i, j
      PetscScalar :: YN,YP,YS
      PetscScalar :: dx,XE,XP,XW,xwu,xeu
      PetscScalar :: ve,vn,vw,vs,me,mn,mw,ms
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);

      ! --- Compute NLV for y-momentum equation ---------------------------
      NLV = 0.D0 

      ! Wash Channel
      DO j = 2,NxL+1
      DO i = MAX(2,rowStartV), MIN(Ny,rowEndV)

          ! Cell height/width
          dx = XU(i,j) - XU(i,j-1);

          ! Left and Right Edges of cell
          xwu = XU(i,j-1); xeu = XU(i,j);

          ! Location of Nodes
          XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
          YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);

          ! East
          IF ((j == NxL+1) .AND. (i < NyL+1)) THEN
              ve = V(i,j+1); ! don't interpolate on boundary
          ELSE IF ((j == NxL+1) .AND. (i > NyM+1)) THEN
              ve = V(i,j+1); ! don't interpolate on boundary
          ELSE
              ve = V(i,j) + (xeu-XP)* (V(i,j+1) - V(i,j)) / (XE-XP);
          END IF
          me = 0.5*( U(i+1,j) * (YN-YP) + U(i,j) * (YP-YS) );

          ! North
          vn = 0.5*( V(i+1,j) + V(i,j) );
          mn = 0.5* dx * ( V(i+1,j) + V(i,j) );


          ! West
          IF (j == 2) THEN
              vw = V(i,j-1); ! don't interpolate on boundary
          ELSE
              vw = V(i,j) + (xwu-XP) * ( V(i,j-1)-V(i,j) ) / (XW-XP);
          END IF
          mw = 0.5*( U(i+1,j-1) * (YN-YP) + U(i,j-1) * (YP-YS) );

          ! South
          vs = 0.5*( V(i,j) + V(i-1,j) );
          ms = 0.5*dx * ( V(i,j) + V(i-1,j) );

          NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;
      END DO
      END DO

      ! Injury Channel
      IF ((rowEndV .GE. NyL+2) .AND. (rowStartV .LE. NyM)) THEN
      DO j = NxL+2,NxM+1
      DO i = MAX(NyL+2,rowStartV), MIN(NyM,rowEndV)

          ! Cell height/width
          dx = XU(i,j) - XU(i,j-1);

          ! Left and Right Edges of cell
          xwu = XU(i,j-1); xeu = XU(i,j);

          ! Location of Nodes
          XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
          YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);

          ! East
          ve = V(i,j) + (xeu-XP)* (V(i,j+1) - V(i,j)) / (XE-XP);
          me = 0.5*( U(i+1,j) * (YN-YP) + U(i,j) * (YP-YS) );

          ! North
          vn = 0.5*( V(i+1,j) + V(i,j) );
          mn = 0.5* dx * ( V(i+1,j) + V(i,j) );


          ! West
          vw = V(i,j) + (xwu-XP) * ( V(i,j-1)-V(i,j) ) / (XW-XP);
          mw = 0.5*( U(i+1,j-1) * (YN-YP) + U(i,j-1) * (YP-YS) );

          ! South
          vs = 0.5*( V(i,j) + V(i-1,j) );
          ms = 0.5*dx * ( V(i,j) + V(i-1,j) );

          NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;
      END DO
      END DO
      END IF

      ! Blood Channel
      DO j = NxM+2,Nx+1
      DO i = MAX(2,rowStartV),MIN(Ny,rowEndV)

          ! Cell height/width
          dx = XU(i,j) - XU(i,j-1);

          ! Left and Right Edges of cell
          xwu = XU(i,j-1); xeu = XU(i,j);

          ! Location of Nodes
          XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
          YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);

          ! East
          IF (j == Nx+1) THEN
              ve = V(i,j+1); ! don't interpolate on boundary
          ELSE
              ve = V(i,j) + (xeu-XP)* (V(i,j+1) - V(i,j)) / (XE-XP);
          END IF
          me = 0.5*( U(i+1,j) * (YN-YP) + U(i,j) * (YP-YS) );

          ! North
          vn = 0.5*( V(i+1,j) + V(i,j) );
          mn = 0.5* dx * ( V(i+1,j) + V(i,j) );


          ! West
          IF ((j == NxM+2) .AND. (i < NyL+1)) THEN
              vw = V(i,j-1); ! don't interpolate on boundary
          ELSE IF ((j == NxM+2) .AND. (i > NyM+1)) THEN
              vw = V(i,j-1); ! don't interpolate on boundary
          ELSE
              vw = V(i,j) + (xwu-XP) * ( V(i,j-1)-V(i,j) ) / (XW-XP);
          END IF
          mw = 0.5*( U(i+1,j-1) * (YN-YP) + U(i,j-1) * (YP-YS) );

          ! South
          vs = 0.5*( V(i,j) + V(i-1,j) );
          ms = 0.5*dx * ( V(i,j) + V(i-1,j) );

          NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;
      END DO
      END DO

      END SUBROUTINE nonlin_V

      END MODULE nonlin_mod
