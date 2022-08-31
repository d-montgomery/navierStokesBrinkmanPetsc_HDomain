      MODULE Ap_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
#include <petsc/finclude/petscmat.h>  
      USE petscsys
      USE petscmat

! Necessary F90 Modules
      USE params_mod
      IMPLICIT NONE

CONTAINS

! -----------------------------------------------------------------------
! Subroutine: Coefficients for Inner Cells of Ap 
! -----------------------------------------------------------------------
      SUBROUTINE innerAp(i,j,rowStartU,rowStartV,XU,YU,XV,YV,&
                         xp,yp,dt,Re,ae,an,aw,as,ap)
      PetscInt, INTENT(IN) :: i,j,rowStartU,rowStartV
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU,YU 
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV,YV
      PetscScalar, DIMENSION(:), INTENT(IN) :: xp, yp
      PetscScalar, INTENT(IN) :: dt, Re ! time step & Reynolds number
      ! Variables 
      PetscScalar :: dx,dy,xe,xj,xw,yn,yi,ys
      PetscScalar :: ae,an,aw,as,ap 

      ! Location of grid points (Note xp & yp have locations of ghost nodes)
      xj = xp(j+1)
      xe = xp(j+2)
      xw = xp(j)
      yi = yp(i+1)
      yn = yp(i+2)
      ys = yp(i)
      dx = XU(i+1,j+1)-XU(i+1,j)
      dy = YV(i+1,j+1)-YV(i,j+1)

      ! Coefficients
      ae = 1.D0/(xe-xj)/dx
      aw = 1.D0/(xj-xw)/dx
      an = 1.D0/(yn-yi)/dy
      as = 1.D0/(yi-ys)/dy 
      ap = -1.D0*(ae+aw+an+as)
      END SUBROUTINE innerAp


! -----------------------------------------------------------------------
! Subroutine to build the Ap matrix for Poisson equation in Phi 
! -----------------------------------------------------------------------
      SUBROUTINE build_Ap(rank,N,NNXY,locDim,rowStartU,rowStartV,&
                          rowstartP,rowEndP,XU,YU,XV,YV,xp,yp,dt,Re,&
                          A,ierr)
      IMPLICIT NONE
      ! Inputs
      PetscInt, INTENT(IN) :: rank !rank in communicator
      PetscInt, DIMENSION(:), INTENT(IN) :: N !N=[NxL,NxM,Nx,NyL,NyM,Ny]
      PetscInt, INTENT(IN) :: rowStartU,rowStartV,rowStartP,rowEndP
      PetscInt, INTENT(IN) :: NNXY ! Au is NNXY x NNXY matrix
      PetscInt, INTENT(IN) :: locDim !number of rows ea. core allocates
      PetscScalar, DIMENSION(rowStartU-1:,:), INTENT(IN) :: XU,YU 
      PetscScalar, DIMENSION(rowStartV-1:,:), INTENT(IN) :: XV,YV 
      PetscScalar, DIMENSION(:), INTENT(IN) :: xp, yp
      PetscScalar, INTENT(IN) :: dt, Re ! time step & Reynolds number
      Mat, INTENT(INOUT) :: A ! matrix
      PetscErrorCode, INTENT(INOUT) :: ierr
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny, Nsml
      PetscInt :: i,j,j_ctr,row
      PetscScalar :: ae,an,aw,as,ap,s_BC,n_BC,Lp,Tp,Rp,Bp

      ! Set if Neumann or Dirichlet Condition on South/North Boundaries
      s_BC = -1.D0 ! negative => Dirichlet-PHI, positive => Neumann-PHI
      n_BC = -1.D0 ! negative => Dirichlet-PHI, positive => Neumann-PHI
      
      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
      
      ! Size of sml sub matrices = Nsub x Nsub
      Nsml = (NxL) + (Nx - NxM)
      
      CALL MatSetSizes(A, locdim, locdim, NNXY, NNXY, ierr)
      CALL MatSetFromOptions(A, ierr)
      CALL MatSeqAIJSetPreallocation(A, 6, PETSC_NULL_INTEGER, ierr)
      CALL MatMPIAIJSetPreallocation(A, 6, PETSC_NULL_INTEGER, 6, &
                PETSC_NULL_INTEGER, ierr)

! --- Fill Inner Cells: Left Lower Legs of H ----------------------------
      DO i = MAX(1,rowStartP), MIN(NyL,rowEndP)
      DO j = 1,NxL
            row = Nsml*(i-1)+j-1

            ! Get coefficients ae,an,aw,as,ap
            CALL innerAp(i,j,rowStartU,rowStartV,&
                         XU,YU,XV,YV,xp,yp,dt,Re,ae,an,aw,as,ap)

            ! Extra Coefficients for Ghost Nodes
            Bp = 0; Tp = 0; Lp = 0; Rp = 0;
            IF (i == 1) THEN ! \partial\Omega_ 2 or 4
            Bp = s_BC*as
            END IF
            IF (j == 1) THEN ! \partial\Omega_5
            Lp = aw
            END IF
            IF (j == NxL) THEN ! \Gamma_2
            Rp = ae
            END IF

            ! East
            IF (j < NxL) THEN
            CALL MatSetValue(A,row,row+1,ae,INSERT_VALUES,ierr)  
            END IF

            ! North
            CALL MatSetValue(A,row,row+Nsml,an,INSERT_VALUES,ierr)

            ! West
            IF (j > 1) THEN
            CALL MatSetValue(A,row,row-1,aw,INSERT_VALUES,ierr)  
            END IF

            ! South
            IF (i > 1) THEN
            CALL MatSetValue(A,row,row-Nsml,as,INSERT_VALUES,ierr)
            END IF

            ! Principal
            ap =  ap + Bp + Lp + Rp
            CALL MatSetValue(A,row,row,ap,INSERT_VALUES,ierr)  
      END DO
      END DO

! --- Fill Inner Cells: Right Lower Legs of H ----------------------------
      DO i = MAX(1,rowStartP), MIN(NyL,rowEndP)
            j_ctr = NxL
      DO j = NxM+1,Nx
            j_ctr = j_ctr + 1
            row = Nsml*(i-1)+j_ctr-1

            ! Get coefficients ae,an,aw,as,ap
            CALL innerAp(i,j,rowStartU,rowStartV,&
                         XU,YU,XV,YV,xp,yp,dt,Re,ae,an,aw,as,ap)

            ! Extra Coefficients for Ghost Nodes
            Bp = 0; Tp = 0; Lp = 0; Rp = 0;
            IF (i == 1) THEN ! \partial\Omega_ 2 or 4
            Bp = s_BC*as 
            END IF
            IF (j == NxM+1) THEN ! \partial\Omega_6
            Lp = aw
            END IF
            IF (j == Nx) THEN ! \Gamma_2
            Rp = ae
            END IF

            ! East
            IF (j < Nx) THEN
            CALL MatSetValue(A,row,row+1,ae,INSERT_VALUES,ierr)  
            END IF

            ! North
            IF (i == NyL) THEN ! At Interface
            CALL MatSetValue(A,row,row+Nx,an,INSERT_VALUES,ierr)
            ELSE
            CALL MatSetValue(A,row,row+Nsml,an,INSERT_VALUES,ierr)
            END IF

            ! West
            IF (j > NxM+1) THEN
            CALL MatSetValue(A,row,row-1,aw,INSERT_VALUES,ierr)  
            END IF

            ! South
            IF (i > 1) THEN
            CALL MatSetValue(A,row,row-Nsml,as,INSERT_VALUES,ierr)
            END IF

            ! Principal
            ap =  ap + Bp + Lp + Rp
            CALL MatSetValue(A,row,row,ap,INSERT_VALUES,ierr)  
      END DO
      END DO

! --- Fill Inner Cells: Injury Channel -------------------------------------
      DO i = MAX(NyL+1,rowStartP), MIN(NyM,rowEndP)
      DO j = 1, Nx
            row = Nsml*(NyL) + (i - NyL -1) * Nx + j - 1

            ! Get coefficients ae,an,aw,as,ap
            CALL innerAp(i,j,rowStartU,rowStartV,&
                     XU,YU,XV,YV,xp,yp,dt,Re,ae,an,aw,as,ap)

            ! Extra Coefficients for Ghost Nodes
            Bp = 0; Tp = 0; Lp = 0; Rp = 0;
            IF (i == NyL+1) THEN ! \Gamma_6
            IF ((j > NxL) .AND. (j < NxM+1)) THEN
            Bp = as ! negative => Dirichlet
            END IF
            END IF

            IF (j == 1) THEN ! \partial\Omega_5
            Lp = aw
            END IF

            IF (i == NyM) THEN
            IF ((j > NxL) .AND. (j < NxM+1)) THEN !\gamma_5
            Tp = an
            END IF
            END IF

            ! East
            IF (j < Nx) THEN
            CALL MatSetValue(A,row,row+1,ae,INSERT_VALUES,ierr)  
            END IF

            ! North
            IF (i < NyM) THEN ! In Injury Channel
            CALL MatSetValue(A,row,row+Nx,an,INSERT_VALUES,ierr)
            ELSE IF ((i == NyM) .AND. (j < NxL+1)) THEN ! Left Interface
            CALL MatSetValue(A,row,row+Nx,an,INSERT_VALUES,ierr)
            ELSE IF ((i == NyM) .AND. (j > NxM)) THEN ! Right Interface
            CALL MatSetValue(A,row,row+Nsml,an,INSERT_VALUES,ierr)
            END IF

            ! West
            IF (j > 1) THEN
            CALL MatSetValue(A,row,row-1,aw,INSERT_VALUES,ierr)  
            END IF

            ! South
            IF (i > NyL+1) THEN !in Injury Channel
            CALL MatSetValue(A,row,row-Nx,as,INSERT_VALUES,ierr)
            ELSE IF ((i == NyL+1) .AND. (j < NxL+1)) THEN
            CALL MatSetValue(A,row,row-Nsml,as,INSERT_VALUES,ierr)
            ELSE IF ((i == NyL+1) .AND. (j > NxM)) THEN
            CALL MatSetValue(A,row,row-Nx,as,INSERT_VALUES,ierr)
            END IF

            ! Principal
            ap =  ap + Tp + Bp + Lp + Rp 
            CALL MatSetValue(A,row,row,ap,INSERT_VALUES,ierr)  
      END DO
      END DO

! --- Fill Inner Cells: Left Upper Legs of H ----------------------------
      DO i = MAX(NyM+1,rowStartP), MIN(Ny,rowEndP)
      DO j = 1,NxL
            row = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j - 1

            ! Get coefficients ae,an,aw,as,ap
            CALL innerAp(i,j,rowStartU,rowStartV,&
                       XU,YU,XV,YV,xp,yp,dt,Re,ae,an,aw,as,ap)

            ! Extra Coefficients for Ghost Nodes
            Bp = 0; Tp = 0; Lp = 0; Rp = 0;
            IF (i == Ny) THEN ! \partial\Omega_ 1
            Tp = n_BC*an ! negative => Dirichlet
            END IF
            IF (j == 1) THEN ! \partial\Omega_5
            Lp = aw
            END IF
            IF (j == NxL) THEN ! \Gamma_1
            Rp = ae
            END IF

            ! East
            IF (j < NxL) THEN
            CALL MatSetValue(A,row,row+1,ae,INSERT_VALUES,ierr)  
            END IF

            ! North
            IF (i < Ny) THEN
            CALL MatSetValue(A,row,row+Nsml,an,INSERT_VALUES,ierr)
            END IF

            ! West
            IF (j > 1) THEN
            CALL MatSetValue(A,row,row-1,aw,INSERT_VALUES,ierr)  
            END IF

            ! South
            IF (i == NyM+1) THEN !At Interface
            CALL MatSetValue(A,row,row-Nx,as,INSERT_VALUES,ierr)
            ELSE
            CALL MatSetValue(A,row,row-Nsml,as,INSERT_VALUES,ierr)
            END IF

            ! Principal
            ap =  ap + Tp + Bp + Lp + Rp
            CALL MatSetValue(A,row,row,ap,INSERT_VALUES,ierr)  
      END DO
      END DO

! --- Fill Inner Cells: Right Upper Legs of H ----------------------------
      DO i = MAX(NyM+1,rowStartP), MIN(Ny,rowEndP)
            j_ctr = NxL
      DO j = NxM+1,Nx
            j_ctr = j_ctr + 1
            row = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j_ctr-1

            ! Get coefficients ae,an,aw,as,ap
            CALL innerAp(i,j,rowStartU,rowStartV,&
                       XU,YU,XV,YV,xp,yp,dt,Re,ae,an,aw,as,ap)

            ! Extra Coefficients for Ghost Nodes
            Bp = 0; Tp = 0; Lp = 0; Rp = 0;
            IF (i == Ny) THEN ! \partial\Omega_3
            Tp = n_BC*an 
            END IF
            IF (j == NxM+1) THEN ! \gamma_3
            Lp = aw
            END IF
            IF (j == Nx) THEN ! \partial\Omega_6
            Rp = ae
            END IF

            ! East
            IF (j < Nx) THEN
            CALL MatSetValue(A,row,row+1,ae,INSERT_VALUES,ierr)  
            END IF

            ! North
            IF (i < Ny) THEN 
            CALL MatSetValue(A,row,row+Nsml,an,INSERT_VALUES,ierr)
            END IF

            ! West
            IF (j > NxM+1) THEN
            CALL MatSetValue(A,row,row-1,aw,INSERT_VALUES,ierr)  
            END IF

            ! South
            CALL MatSetValue(A,row,row-Nsml,as,INSERT_VALUES,ierr)

            ! Principal
            ap =  ap + Tp + Bp + Lp + Rp
            CALL MatSetValue(A,row,row,ap,INSERT_VALUES,ierr)  
      END DO
      END DO

!      ! Add Lagrange Multiplier
!      CALL MatGetOwnershipRange(A, i_start, i_end, ierr)
!      i_end = i_end - 1
!      DO i = i_start, i_end
!            CALL MatSetValue(A,i,NNXY-1,1.D0,INSERT_VALUES,ierr)  
!      END DO
!
!      IF (i_end == NNXY-1) THEN
!            CALL MatSetValue(A,NNXY-1,0,1.D0,INSERT_VALUES,ierr)  
!            CALL MatSetValue(A,NNXY-1,NNXY-1,0.D0,INSERT_VALUES,ierr)  
!      END IF
      
      CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)   
      CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      END SUBROUTINE build_Ap

      END MODULE Ap_mod
