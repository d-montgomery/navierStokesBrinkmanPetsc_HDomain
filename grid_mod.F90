      MODULE grid_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

! Necessary F90 Modules
      USE params_mod
      USE utilities_mod 
    
      IMPLICIT NONE

CONTAINS

! -----------------------------------------------------------------------
! Subdivide H given Nx & Ny: Return N =[NxL, NxM, NxR, NyL, NyM, NyU] 
! -----------------------------------------------------------------------
      SUBROUTINE num_spatialPoints(Nx,Ny,GRID,N)
      PetscInt, INTENT(INOUT) :: Nx, Ny
      PetscInt, INTENT(IN) :: GRID
      PetscInt, DIMENSION(:), INTENT(OUT) :: N
      PetscScalar :: W, H
      PetscInt :: NxL, NxMid, NxR, NyLow, NyMid, NyUp
      PetscScalar :: Sx, Sy, WashPercent, BloodPercent

      ! Determine overall width and height of Domain
      W = WL + Wmid + WR
      H = HL + Hmid + HU

      IF (GRID == 0) THEN ! Distribute Points for Uniform Grid
         ! --- Determine number of grid points in each region of H -----
         ! Number of x-grid points for each region
         NxL = NINT(Nx*WL/W) ! x pts in wash channel
         NxR = NINT(Nx*WR/W) ! x pts in blood channels
         NxMid = NINT(Nx*Wmid/W) ! Nx in Injury Channel
         Nx = NxMid + NxL + NxR ! Total # of points in x-direction

         ! Number of y-grid points for each region
         NyLow = NINT(HL/H*Ny) ! Lower Legs of H
         NyMid = NINT(Hmid/H*Ny) ! Injury Channel
         NyUp = NINT(HU/H*Ny) ! Upper Legs of H
         Ny = NyLow + NyMid + NyUp ! Total # of points in y-direction

         N = (/ NxL, NxL + NxMid, Nx, NyLow, NyLow + NyMid, Ny /)
      ELSE IF (GRID == 1) THEN
         ! percentage of x-pts and y-pts in injury channel
         Sx = DBLE(0.55); 
         Sy = DBLE(0.3);
         
         ! Number of x-grid points for each region
         NxMid = NINT(Sx*Nx);
         WashPercent = WL/(W - Wmid)
         NxL = NINT((1.D0-Sx)*DBLE(Nx)*WashPercent)
         BloodPercent = WR/(W - Wmid)
         NxR = NINT((1.D0-Sx)*DBLE(Nx)*BloodPercent)
         
         ! percentage of y-pts in injury channel
         NyMid = NINT(Sy*Ny);
         NyLow = NINT(HL/(H-Hmid)*(1.D0-Sy)*DBLE(Ny)) ! Injury Channel
         NyUp = NINT(HU/(H-Hmid)*(1.D0-Sy)*DBLE(Ny)) ! Upper Legs of H
         
         Nx = NxL + NxMid + NxR
         Ny = NyLow + NyMid + NyUp
     
         N = (/ NxL, NxL + NxMid, Nx, NyLow, NyLow + NyMid, Ny /)
      END IF 
      END SUBROUTINE num_spatialPoints

! -----------------------------------------------------------------------
! Force zeros in middle to form H 
! -----------------------------------------------------------------------
      SUBROUTINE zeroH(NxL,NyL,NxM,Nym,Ny,rowStartU,rowEndU,rowStartV,&
                       rowEndV,rowStartP,rowEndP,U,V,P)
      !Inputs
      PetscInt, INTENT(IN) :: NxL, NyL, NxM, NyM, Ny,rowStartU,rowEndU,&
                              rowStartV,rowEndV,rowStartP,rowEndP
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(INOUT) :: U
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(INOUT) :: V
      PetscScalar, DIMENSION(rowStartP-2:,:), INTENT(INOUT) :: P
      
      ! Zero out U between legs of H
      IF (rowStartU-1 .LE. NyL) THEN
            U(MAX(1,rowStartU-1):MIN(NyL,rowEndU+1),NxL+2:Nxm) = 0
      END IF
      IF (rowEndU+1 .GE. NyM+3) THEN 
            U(MAX(Nym+3,rowStartU-1):MIN(Ny+2,rowEndU+1),NxL+2:Nxm) = 0
      END IF

      ! Zero out V between legs of H
      IF (rowStartV-1 .LE. NyL) THEN 
            V(MAX(1,rowStartV-1):MIN(NyL,rowEndV+1),NxL+3:Nxm) = 0
      END IF
      IF (rowEndV+1 .GE. NyM+2) THEN
            V(MAX(Nym+2,rowStartV-1):MIN(Ny+1,rowEndV+1),NxL+3:Nxm) = 0
      END IF

      ! Zero out P between legs of H
      IF (rowStartP-2 .LE. NyL) THEN
            P(MAX(1,rowStartP-2):MIN(NyL,rowEndP+2),NxL+1:Nxm) = 0
      END IF
      IF (rowEndP+2 .GE. NyM+1) THEN
            P(MAX(Nym+1,rowStartP-2):MIN(Ny,rowEndP+2),NxL+1:Nxm) = 0
      END IF

      END SUBROUTINE zeroH

! -----------------------------------------------------------------------
! Create Staggered Grid for H domain ------------------------------------
! -----------------------------------------------------------------------
      SUBROUTINE H_grid(N,GRID,XXUU,YYUU,XXVV,YYVV,xp,yp,min_h)
      IMPLICIT NONE
      ! Inputs
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscInt, INTENT(IN) :: GRID
      ! Outputs
      PetscScalar, DIMENSION(:,:), INTENT(OUT) :: XXUU, YYUU, XXVV, YYVV
      PetscScalar, DIMENSION(:), INTENT(OUT) :: xp, yp
      PetscScalar, INTENT(OUT) :: min_h
      ! Intermediate Variables
      PetscScalar :: W, H
      PetscInt :: i, NxL, NxMid, NxR, Nx, NyLow, NyMid, NyUp, Ny 
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: dx,x,xW,xI,xB,xu,xv
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: dy,y,yLow,yI,yUp,yu,yv
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: dxw, dyl, yLtmp
      PetscScalar :: p

      ! Determine overall width and height of Domain
      W = WL + Wmid + WR
      H = HL + Hmid + HU

      ! --- Determine number of grid points in each region of H ---------
      ! Number of x-grid points for each region
      NxL = N(1) ! x pts in wash channel
      NxMid = N(2)-N(1) ! Nx in Injury Channel
      NxR = N(3)-N(2) ! x pts in blood channels
      Nx = N(3) ! Total # of points in x-direction

      ! Number of y-grid points for each region
      NyLow = N(4) ! Lower Legs of H
      NyMid = N(5)-N(4) ! Injury Channel
      NyUp = N(6)-N(5) ! Upper Legs of H
      Ny = N(6) ! Total # of points in y-direction
      
      IF (GRID == 0) THEN ! Create Uniform Grid
         ! Wash Channel - x
         ALLOCATE(dx(1), dy(1))
         ALLOCATE(xW(NxL+1), xI(NxMid+1), xB(NxR+1))
         dx = WL / NxL
         xW = (/ (i*dx, i = 0, NxL) /)

         ! Injury Channel - x
         dx = Wmid / NxMid
         xI = (/ (WL + i*dx, i = 0, NxMid) /)

         ! Blood Channel - x
         dx = WR / NxR
         xB = (/ (WL + Wmid + i*dx, i = 0, NxR) /)

         ! Lower - y
         ALLOCATE(yLow(NyLow+1), yI(NyMid+1), yUp(NyUp+1))
         dy = HL / NyLow
         yLow = (/ (i*dy, i = 0, NyLow) /)

         ! Injury - y
         dy = Hmid / NyMid;
         yI = (/ (HL + i*dy, i = 0, NyMid ) /)

         ! Upper - y
         dy = HU / NyUp;
         yUp = (/ (HL + Hmid + i*dy, i = 0, NyUp) /)
         

         ! Build x and y vectors for grid
         ALLOCATE(x(Nx+1), y(Ny+1))
         x = (/ xW, xI(2:NxMid+1), xB(2:NxR+1) /)
         y = (/ yLow, yI(2:NyMid+1), yUp(2:NyUp+1) /)
         DEALLOCATE(dx, dy, xW, xI, xB, yLow, yI, yUp)
      
      ELSE IF (GRID == 1) THEN !Create Nonuniform Grid
         ALLOCATE(dx(1), dy(1))
         ALLOCATE(xW(NxL+1), xI(NxMid+1), xB(NxR+1))
         
         ! Injury Channel - x
         dx = Wmid / NxMid
         xI = (/ (WL + i*dx, i = 0, NxMid) /)
         
         ! Blood Channel - x
         p = LOG(WL/(xI(NxMid+1)-xI(NxMid)))/LOG(DBLE(NxL)) - 1
         xB(1) = WL+Wmid;
         DO i = 2, NxL+1
            dx = DBLE(i-1)**p * WL / DBLE(NxL)**(p+1)
            xB(i) = WL + Wmid + DBLE(i-1)*dx(1)
         END DO

         ! Wash Channel - x
         ALLOCATE(dxw(NxL))
         CALL fliplr( xB(2:NxL+1) - xB(1:NxL), dxw)
         xW(1) = 0.D0
         DO i = 2,NxL+1
             xW(i) = xW(i-1) + dxw(i-1)
         END DO
         DEALLOCATE(dxw)

         ! Injury - y
         ALLOCATE(yLow(NyLow+1), yI(NyMid+1), yUp(NyUp+1))
         dy = Hmid / NyMid
         yI = (/ (HL + i*dy, i = 0,NyMid) /)
         ! Upper - y
         p = LOG(HU/(yI(NyMid+1)-yI(NyMid)))/LOG(DBLE(NyUp)) - 1
         yUp(1) = HL+Hmid
         DO i = 2, NyUp+1
            dy(1) = DBLE(i-1)**p * HU / NyUp**(p+1)
            yUp(i) = HL + Hmid + (i-1)*dy(1)
         END DO

         ! Lower - y
         ALLOCATE(yLtmp(NyLow+1),dyl(NyLow))
         p = LOG(HL/(yI(NyMid+1)-yI(NyMid)))/LOG(DBLE(NyLow)) - 1
         yLtmp(1) = 0.D0
         DO i = 2, NyLow+1
            dy = DBLE(i-1)**p * HL / DBLE(NyLow)**(p+1)
            yLtmp(i) = 0.D0 + (i-1)*dy(1)
         END DO
         CALL fliplr(yLtmp(2:NyLow+1) - yLtmp(1:NyLow), dyl)
         yLow(NyLow+1) = HL
         DO i = NyLow,1,-1
             yLow(i) = yLow(i+1) - dyl(i);
         END DO
         DEALLOCATE(yLtmp, dyl)
         
         x = (/ xW, xI(2:NxMid), xB /)
         y = (/ yLow, yI(2:NyMid), yUp /)

         DEALLOCATE(dx, dy, xW, xI, xB, yLow, yI, yUp)
      
      END IF
      

      ! Length of each cell (general for non-uniform grid)
      ALLOCATE(dx(Nx), dy(Ny))
      dx = x(2:Nx+1) - x(1:Nx)
      dy = y(2:Ny+1) - y(1:Ny)
      min_h = MINVAL(dx) + MINVAL(dy) ! Smallest cell size for CFL 

      ! --- Location of All Points in a Rectangle ---
      ! Location of pressure (mid points of cells)
      xp = (/ x(1)-dx(1)/2, x(1:Nx) + dx/2, x(Nx+1)+dx(Nx)/2 /)
      yp = (/ y(1)-dy(1)/2, y(1:Ny) + dy/2, y(Ny+1)+dy(Ny)/2 /)

      ! Location of u for rectangle
      ALLOCATE(xu(Nx+1), yu(Ny+2))
      xu = x
      yu = (/ y(1) , y(1:Ny) + dy/2, y(Ny+1) /)

      ! Location of v for rectangle
      ALLOCATE(xv(Nx+2), yv(Ny+1))
      xv = (/ x(1), x(1:Nx) + dx/2, x(Nx+1) /)
      yv = y
      
      ! Make Meshgrid of points for u and v
      CALL meshgrid(xu,yu,XXUU,YYUU)
      CALL meshgrid(xv,yv,XXVV,YYVV)

      ! Adjust for boundaries of H
      YYUU(N(4)+1, N(1)+2:N(2)) = YYVV(N(4)+1,N(1)+2:N(2))
      YYUU(N(5)+2, N(1)+2:N(2)) = YYVV(N(5)+1,N(1)+2:N(2))
      
      XXVV(1:N(4), N(1)+2) = XXUU(1:N(4), N(1)+1)
      XXVV(N(5)+2:Ny+1, N(1)+2) = XXUU(N(5)+2:Ny+1, N(1)+1)
      XXVV(1:N(4), N(2)+1) = XXUU(1:N(4), N(2)+1)
      XXVV(N(5)+2:Ny+1, N(2)+1) = XXUU(N(5)+2:Ny+1, N(2)+1)
      
      DEALLOCATE(x,y,dx,dy,xu,yu,xv,yv)

      END SUBROUTINE H_grid
 
      END MODULE grid_mod
