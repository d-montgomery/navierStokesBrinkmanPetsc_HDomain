      MODULE fluids_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>          
#include <petsc/finclude/petscksp.h>
      USE petscsys
      USE petscvec           
      USE petscmat
      USE petscksp    
! Necessary .F90 Modules
      USE params_mod
      USE utilities_mod
      USE indexing_mod
      USE grid_mod
      USE soln_mod
      USE Au_mod
      USE Av_mod
      USE Ap_mod
      USE brinkman_mod
      USE Bu_mod
      USE Bv_mod
      USE initial_P_mod
      USE nonlin_mod
      USE uRHS_mod
      USE vRHS_mod
      USE pRHS_mod
      USE reshape_mod
      USE correction_mod
      IMPLICIT NONE
CONTAINS

! -----------------------------------------------------------------------
! FLUIDS SOLVER W/BRINKMAN SUBROUTINE
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
! Second Order Semi-implicit FVM for solving 2D (Nondimensional) 
! Incompressible Navier-Stokes-Brinkman Equations w/Neumann BC on Southern Bndry:
! u_t + u*u_x + v*u_x = - p_x + 1/Re * (u_xx + u_yy) + fu - Bu,
! v_t + u*v_x + v*v_y = - p_y + 1/Re * (v_xx + v_yy) + fv - Bv,
! u_x + v_y = 0.
! -----------------------------------------------------------------------
      SUBROUTINE fluids_solver(Nx,Ny,GRID,Tf,tol,save_soln)

!--- Inputs -------------------------------------------------------------
      PetscInt :: Nx, Ny ! Number of spatial pts in x/y directions
      PetscInt :: GRID !0 -> uniform grid, 1 -> non-uniform grid
      PetscScalar :: Tf ! Temporal parameters
      PetscScalar :: tol
      PetscInt, INTENT(IN) :: save_soln ! if = 1, write soln to file

! --- Variables ---------------------------------------------------------
      ! MPI Related Variables
      PetscErrorCode :: ierr, fail_code
      PetscInt :: my_rank, num_cores, master = 0
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpi_status
      PetscReal :: t_w0, t_wf, t_tot0, t_totF
      
      ! Indexing for U-System - PETSc matrix and vectors
      PetscInt :: sysStartU, sysEndU
      PetscInt, ALLOCATABLE, DIMENSION(:) :: sysDimU, sysIndxU
      PetscInt :: sysStartV, sysEndV
      PetscInt, ALLOCATABLE, DIMENSION(:) :: sysDimV, sysIndxV
      PetscInt :: sysStartP, sysEndP
      PetscInt, ALLOCATABLE, DIMENSION(:) :: sysDimP, sysIndxP
      
      ! Indexing for U,V,P-Arrays
      PetscInt, ALLOCATABLE, DIMENSION(:) :: dimU, rowStartU, rowEndU
      PetscInt, ALLOCATABLE, DIMENSION(:) :: dimV, rowStartV, rowEndV
      PetscInt, ALLOCATABLE, DIMENSION(:) :: dimP, rowStartP, rowEndP

      ! Temporal and Spatial Variables
      PetscScalar :: t, dt, min_h
      PetscInt :: NxL, NxM, NyL, NyM, Nt
      PetscInt :: Nu, Nv, Np
      PetscInt, DIMENSION(6) :: N
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: XU0,YU0,XV0,YV0,XP0,YP0
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: XU,YU,XU_full,YU_full
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: XV,YV,XV_full, YV_full
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: xp, yp
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: XXPP, YYPP
      
      ! Velocity, Pressure, Nonlinear Terms, and Body Forcing Terms
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: U,U0,U_full
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: V,V0,V_full
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: P,P0,P_full
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: u_vecF90,v_vecF90,p_vecF90
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: NLU0, NLU
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: NLV0, NLV
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: Quo, Qu
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: Qvo, Qv
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: v_BC
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: B, B0, theta, B_full
      PetscScalar :: maxU, maxV, maxUV, steady

      ! Matrices
      Mat :: Au, Bu, Av, Bv, Ap
      
      ! Vectors
      Vec :: fu, u_vec
      Vec :: fv, v_vec
      Vec :: fp, p_vec
      
      ! KSP Objects
      KSP :: ksp_U, ksp_V, ksp_P
      
      ! Indexing
      PetscInt :: i,j,k,ii
      CHARACTER (LEN = 40) :: fileName
! --- End of Variable Declaration ---------------------------------------


! -----------------------------------------------------------------------
!  Initilaize MPI/Petsc Communicator
! -----------------------------------------------------------------------
      ! Initilaize MPI/Petsc & get my_rank & num_cores from MPI
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      t_tot0 = MPI_Wtime()

      IF (ierr .NE. 0) THEN
        PRINT *,'Unable to initialize PETSc'
        STOP
      END IF
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, my_rank, ierr)
      CALL MPI_COMM_SIZE(PETSC_COMM_WORLD, num_cores, ierr)
      
      ! Ensure that num_cores is divisible by 4
      IF ( MOD(num_cores,4) == 0 ) THEN
            k = num_cores/4
      ELSE IF (num_cores == 1) THEN
            k = 0
      ELSE IF ((MOD(num_cores,4) .NE. 0) .AND. (my_rank == master)) THEN
            PRINT *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
            PRINT *, 'ERROR: num_cores must be divisible by 4'
            PRINT *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
            CALL MPI_Abort(PETSC_COMM_WORLD,fail_code, ierr)
      END IF

! -----------------------------------------------------------------------
!  Set all physical constants and characteristic scales
! -----------------------------------------------------------------------
      CALL setparams

! -----------------------------------------------------------------------
!  Build the spatial grid
! -----------------------------------------------------------------------
      ! Get number of spatial points in each region, update Nx & Ny 
      IF (my_rank == master) THEN
            CALL num_spatialPoints(Nx,Ny,GRID,N)
            IF ( num_cores > NINT( DBLE(N(6))/4.0 ) )  THEN 
               PRINT *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               PRINT *, 'ERROR: Too many processors for grid'
               PRINT *, '       decomposition. Try decreasing the number'
               PRINT *, '       of processors or use a finer grid.'
               PRINT *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               CALL MPI_Abort(PETSC_COMM_WORLD,fail_code,ierr)
            END IF
      END IF
      CALL MPI_Bcast(N, 6, MPI_INTEGER, master, PETSC_COMM_WORLD, ierr)
      NxL = N(1); NxM = N(2); Nx = N(3); 
      NyL = N(4); NyM = N(5); Ny = N(6);
      
      ! Get Indexing for Domain Decomposition
      ALLOCATE(dimU(0:num_cores-1),rowStartU(0:num_cores-1),&
               rowEndU(0:num_cores-1), dimV(0:num_cores-1), &
               rowStartV(0:num_cores-1),rowEndV(0:num_cores-1),&
               dimP(0:num_cores-1),rowStartP(0:num_cores-1),&
               rowEndP(0:num_cores-1))
      CALL u_indx(num_cores,k,NyL,NyM,Ny,dimU,rowStartU,rowEndU)
      CALL v_indx(num_cores,k,NyL,NyM,Ny,dimV,rowStartV,rowEndV)
      CALL p_indx(num_cores,k,NyL,NyM,Ny,dimP,rowStartP,rowEndP)
      
      ! Allocate local arrays for grid
      ALLOCATE(XU(rowStartU(my_rank)-1:rowEndU(my_rank)+1,Nx+1),&
               YU(rowStartU(my_rank)-1:rowEndU(my_rank)+1,Nx+1),&
               XV(rowStartV(my_rank)-1:rowEndV(my_rank)+1,Nx+2),&
               YV(rowStartV(my_rank)-1:rowEndV(my_rank)+1,Nx+2),&
               XXPP(rowStartP(my_rank)-2:rowEndP(my_rank)+2,Nx),&
               YYPP(rowStartP(my_rank)-2:rowEndP(my_rank)+2,Nx))
      ALLOCATE(xp(Nx+2), yp(Ny+2))
      XU = 0.D0; YU = 0.D0; XV = 0.D0; YV = 0.D0;
      xp = 0.D0; yp = 0.D0; XXPP = 0.D0; YYPP = 0.D0;
      
      !Master creates full mesh and sends specified rows to each processor
      IF (my_rank == master) THEN
            ALLOCATE(XU0(0:Ny+3,Nx+1), YU0(0:Ny+3,Nx+1))
            ALLOCATE(XV0(0:Ny+2,Nx+2), YV0(0:Ny+2,Nx+2))
            ALLOCATE(XP0(-1:Ny+2,Nx),YP0(-1:Ny+2,Nx)) 
            
            XU0 = 0; YU0 = 0; XV0 = 0; YV0 = 0; XP0 = 0; YP0 = 0;
            CALL H_grid(N,GRID,XU0(1:Ny+2,:),YU0(1:Ny+2,:),&
                        XV0(1:Ny+1,:),YV0(1:Ny+1,:),xp,yp,min_h)
            
            ! Create Mesh Grid for Pressure 
            ! ONLY NEEDED WITH EXACT SOLN and Initial Condition
            CALL meshgrid(xp(2:Nx+1), yp, XP0(0:Ny+1,:), YP0(0:Ny+1,:))

            ! Store allocated local values on master
            XU = XU0(rowStartU(my_rank)-1:rowEndU(my_rank)+1,:)
            YU = YU0(rowStartU(my_rank)-1:rowEndU(my_rank)+1,:)
            XV = XV0(rowStartV(my_rank)-1:rowEndV(my_rank)+1,:)
            YV = YV0(rowStartV(my_rank)-1:rowEndV(my_rank)+1,:)
            ! ONLY NEEDED WITH EXACT SOLN and Initial Condition
            XXPP = XP0(rowStartP(my_rank)-2:rowEndP(my_rank)+2,:)
            YYPP = YP0(rowStartP(my_rank)-2:rowEndP(my_rank)+2,:)
            
            ! Send specified rows to other processors
            IF (num_cores > 1) THEN
               DO i = 1,num_cores-1
               CALL MPI_Send(XU0(rowStartU(i)-1:rowEndU(i)+1,:),&
                             (2+dimU(i))*(Nx+1),MPI_DOUBLE_PRECISION, i, &
                             num_cores + i, PETSC_COMM_WORLD, ierr)
               CALL MPI_Send(YU0(rowStartU(i)-1:rowEndU(i)+1,:),&
                             (2+dimU(i))*(Nx+1),MPI_DOUBLE_PRECISION, i, &
                             2*num_cores + i, PETSC_COMM_WORLD, ierr)
               CALL MPI_Send(XV0(rowStartV(i)-1:rowEndV(i)+1,:),&
                             (2+dimV(i))*(Nx+2),MPI_DOUBLE_PRECISION, i, &
                             3*num_cores + i, PETSC_COMM_WORLD, ierr)
               CALL MPI_Send(YV0(rowStartV(i)-1:rowEndV(i)+1,:),&
                             (2+dimV(i))*(Nx+2),MPI_DOUBLE_PRECISION, i, &
                             4*num_cores + i, PETSC_COMM_WORLD, ierr)
               CALL MPI_Send(XP0(rowStartP(i)-2:rowEndP(i)+2,:),&
                             (4+dimP(i))*(Nx),MPI_DOUBLE_PRECISION, i, &
                             5*num_cores + i, PETSC_COMM_WORLD, ierr)
               CALL MPI_Send(YP0(rowStartP(i)-2:rowEndP(i)+2,:),&
                             (4+dimP(i))*(Nx),MPI_DOUBLE_PRECISION, i, &
                             6*num_cores + i, PETSC_COMM_WORLD, ierr)
               END DO
            END IF

            ! Deallocate the temporary arrays used to build grid       
            DEALLOCATE(XU0,YU0,XV0,YV0,XP0,YP0)
      
      ! Receive allocated rows of grid
      ELSE IF (my_rank > master) THEN
         CALL MPI_Recv(XU,(2+dimU(my_rank))*(Nx+1),MPI_DOUBLE_PRECISION,&
                         master,num_cores+my_rank,PETSC_COMM_WORLD,&
                         mpi_status,ierr)
         CALL MPI_Recv(YU,(2+dimU(my_rank))*(Nx+1),MPI_DOUBLE_PRECISION,& 
                       master,2*num_cores+my_rank,PETSC_COMM_WORLD,&
                       mpi_status,ierr)
         CALL MPI_Recv(XV,(2+dimV(my_rank))*(Nx+2),MPI_DOUBLE_PRECISION,& 
                       master,3*num_cores+my_rank,PETSC_COMM_WORLD,&
                       mpi_status,ierr)
         CALL MPI_Recv(YV,(2+dimV(my_rank))*(Nx+2),MPI_DOUBLE_PRECISION,&
                       master,4*num_cores+my_rank,PETSC_COMM_WORLD,&
                       mpi_status,ierr)
         CALL MPI_Recv(XXPP,(4+dimP(my_rank))*(Nx),MPI_DOUBLE_PRECISION,&
                       master,5*num_cores+my_rank,PETSC_COMM_WORLD,&
                       mpi_status,ierr)
         CALL MPI_Recv(YYPP,(4+dimP(my_rank))*(Nx),MPI_DOUBLE_PRECISION,&
                       master,6*num_cores+my_rank,PETSC_COMM_WORLD,&
                       mpi_status,ierr)
      END IF
      
      ! Broadcast xp, yp, and min_h to all cores
      CALL MPI_Bcast(xp, Nx+2, MPI_DOUBLE_PRECISION, master, &
                     PETSC_COMM_WORLD, ierr)
      CALL MPI_Bcast(yp, Ny+2, MPI_DOUBLE_PRECISION, master, &
                     PETSC_COMM_WORLD, ierr)
      CALL MPI_Bcast(min_h, 1, MPI_DOUBLE_PRECISION, master, &
                     PETSC_COMM_WORLD, ierr)

! -----------------------------------------------------------------------
! Determine temporal step size based on CFL condition 
! -----------------------------------------------------------------------
      ! Get Inlet Boundary Condition
      IF (my_rank == num_cores - 1) THEN 
            ALLOCATE(V_BC(Nx+2))
            CALL inletBC_v(XV(Ny+1,:),V_BC)
            maxUV = -1D0*MINVAL(V_BC)
      END IF
      CALL MPI_Bcast(maxUV,1,MPI_DOUBLE_PRECISION, num_cores-1, &
                     PETSC_COMM_WORLD, ierr)
      dt = min_h / maxUV
      t = 0
      Tf = Tf / Tchar
      dt = 1D-4 / Tchar
! -----------------------------------------------------------------------
! Building the Matrices Au,Av,Ap, and some additional info for U,V,P vecs
! -----------------------------------------------------------------------
      ! Number of elements ea. process allocates/stores for matrices 
      ! Au, Av, Ap solution vectors U,V,P, and RHS vectors fu, fv, fp. 
      ! Note Au is NU x NU, U is NU x 1, Av is NV x NV, etc...
      NU = (Nx +1)*(Ny +2) - (NxM - NxL - 1)*(Ny - NyM + NyL)
      NV = (Nx +2)*(Ny + 1) - (Nxm - Nxl - 2)*(Ny - Nym + NyL)
      NP = Nx*Ny - (NxM - NxL)*(Ny - Nym + NyL)
      
      ! Create vectors with the number of rows ea. processor 
      ! stores for each of the the linear systems
      ALLOCATE(sysDimU(0:num_cores-1), &
               sysDimV(0:num_cores-1), &
               sysDimP(0:num_cores-1))
      CALL dim_LinSysUV((NxL+1)+(Nx+1-NxM),Nx+1,NU,k,dimU,sysDimU) 
      CALL dim_LinSysUV((NxL+2)+(Nx+2-NxM),Nx+2,NV,k,dimV,sysDimV) 
      CALL dim_LinSysP((NxL)+(Nx - NxM),Nx,NP,k,dimP,sysDimP) 
      
      ! Build Matrix Au
      CALL MatCreate(PETSC_COMM_WORLD, Au, ierr)
      CALL build_Au(my_rank,N,NU,sysDimU(my_rank),rowStartU(my_rank),&
                    rowStartV(my_rank),XU,YU,XV,YV,dt,Re,Au,ierr)
!CALL MatView(Au,PETSC_VIEWER_STDOUT_WORLD,ierr)
      
      ! Build Matrix Av
      CALL MatCreate(PETSC_COMM_WORLD, Av, ierr)
      CALL build_Av(my_rank,N,NV,sysDimV(my_rank),rowStartU(my_rank),&
                    rowStartV(my_rank),XU,YU,XV,YV,dt,Re,Av,ierr)

      ! Build Matrix Ap
      CALL MatCreate(PETSC_COMM_WORLD, Ap, ierr)
      CALL build_Ap(my_rank,N,NP,sysDimP(my_rank),rowStartU(my_rank),&
                  rowStartV(my_rank),rowStartP(my_rank),rowEndP(my_rank),&
                  XU,YU,XV,YV,xp,yp,dt,Re,Ap,ierr)

! -----------------------------------------------------------------------
! Brinkman Setup 
! -----------------------------------------------------------------------
      ALLOCATE(theta(rowStartP(my_rank)-2:rowEndP(my_rank)+2,Nx),&
               B(rowStartP(my_rank)-2:rowEndP(my_rank)+2,Nx))
      
      ! Initialize Matrix Bu 
      CALL MatCreate(PETSC_COMM_WORLD,Bu,ierr)
      CALL MatSetSizes(Bu,sysDimU(my_rank),sysDimU(my_rank),NU,NU,ierr)
      CALL MatSetFromOptions(Bu,ierr)
      CALL MatSeqAIJSetPreallocation(Bu,5,PETSC_NULL_INTEGER,ierr)
      CALL MatMPIAIJSetPreallocation(Bu,5,PETSC_NULL_INTEGER,2,&
                      PETSC_NULL_INTEGER, ierr)
      CALL MatZeroEntries(Bu,ierr)

      ! Initialize Matrix Bv
      CALL MatCreate(PETSC_COMM_WORLD, Bv, ierr)
      CALL MatSetSizes(Bv,sysDimV(my_rank),sysDimV(my_rank),NV,NV,ierr)
      CALL MatSetFromOptions(Bv,ierr)
      CALL MatSeqAIJSetPreallocation(Bv,5,PETSC_NULL_INTEGER,ierr)
      CALL MatMPIAIJSetPreallocation(Bv,5,PETSC_NULL_INTEGER,2,&
                      PETSC_NULL_INTEGER, ierr)
      CALL MatZeroEntries(Bv,ierr)
! -----------------------------------------------------------------------
! Initial Conditions 
! -----------------------------------------------------------------------
      ALLOCATE(U0(rowStartU(my_rank)-1:rowEndU(my_rank)+1,Nx+1),&
               V0(rowStartV(my_rank)-1:rowEndV(my_rank)+1,Nx+2),&
               P0(rowStartP(my_rank)-2:rowEndP(my_rank)+2,Nx),&
               B0(rowStartP(my_rank)-2:rowEndP(my_rank)+2,Nx))
      U0 = 0.D0        
      V0 = 0.D0  
      P0 = 0.D0 
      B0 = 0.D0
      B  = 0.D0
      
      ! Get Initial Pressure
      CALL initial_P(N,XXPP,YYPP,rowStartP(my_rank),rowEndP(my_rank),P0)

      ! Zero Out innner H
      CALL zeroH(NxL,NyL,NxM,NyM,Ny,rowStartU(my_rank),rowEndU(my_rank),&
                 rowStartV(my_rank),rowEndV(my_rank),rowStartP(my_rank),&
                 rowEndP(my_rank),U0,V0,P0)

! -----------------------------------------------------------------------
! Calculate Nonlinear Terms for One Step of Forward Euler
! -----------------------------------------------------------------------
      ! --- Nonlinear terms for U: NLU and NLU0
      ALLOCATE(NLU0(rowStartU(my_rank):rowEndU(my_rank),Nx),&
                NLU(rowStartU(my_rank):rowEndU(my_rank),Nx))
      NLU = 0.D0
      NLU0 = 0.D0
      CALL nonlin_U(N,rowStartU(my_rank),rowEndU(my_rank),&
                    rowStartV(my_rank),XU,YU,XV,YV,U0,V0,NLU)
      NLU = 2.D0/3.D0 * NLU ! For one step of Forward Euler

      ! --- Nonlinear terms for V: NLV and NLV0
      ALLOCATE(NLV0(rowStartV(my_rank):rowEndV(my_rank),Nx+1),&
                NLV(rowStartV(my_rank):rowEndV(my_rank),Nx+1))
      NLV = 0.D0
      NLV0 = 0.D0
      CALL nonlin_V(N,rowStartU(my_rank),rowStartV(my_rank),&
                    rowEndV(my_rank),XU,YU,XV,YV,U0,V0,NLV)
      NLV = 2.D0/3.D0 * NLV ! For one step of Forward Euler

! -----------------------------------------------------------------------
! Forcing Terms at t = 0 
! -----------------------------------------------------------------------
      ! Forcing Term for U eqt.
      ALLOCATE(Quo(rowStartU(my_rank):rowEndU(my_rank),Nx+1))
      Quo = 0.D0      
      ! Forcing Term for V eqt
      ALLOCATE(Qvo(rowStartV(my_rank):rowEndV(my_rank),Nx+2))
      Qvo = 0.D0

! -----------------------------------------------------------------------
! Prep U solution for time stepping loop
! -----------------------------------------------------------------------
      ! Create KSP Context for U-momentum
      CALL KSPCreate(PETSC_COMM_WORLD, ksp_u, ierr)
      CALL KSPSetOptionsPrefix(ksp_u,"u_sys_",ierr)
      CALL KSPSetFromOptions(ksp_u, ierr)
      CALL KSPSetOperators(ksp_u, Au, Au, ierr)

      ! Initialize u-RHS vector
      CALL VecCreate(PETSC_COMM_WORLD,fu,ierr) 
      CALL VecSetSizes(fu,sysDimU(my_rank),NU,ierr)
      CALL VecSetFromOptions(fu,ierr)
      
      ! Initialize u-soln vector
      CALL VecCreate(PETSC_COMM_WORLD,u_vec,ierr) 
      CALL VecSetSizes(u_vec,sysDimU(my_rank),NU,ierr)
      CALL VecSetFromOptions(u_vec,ierr)
      
      ! Allocate for body forcing Qu, solution U, and temporary u_vecF90
      ALLOCATE(Qu(rowStartU(my_rank):rowEndU(my_rank),Nx+1))
      ALLOCATE(U(rowStartU(my_rank)-1:rowEndU(my_rank)+1,Nx+1))
      ALLOCATE(u_vecF90(sysDimU(my_rank)))
      u_vecF90 = 0.D0; U = 0.D0;

      ! Get indexing used in reshaping u_vec -> u_vecF90 -> U
      CALL VecGetOwnershipRange(u_vec, sysStartU, sysEndU, ierr)
      ALLOCATE(sysIndxU(sysStartU:sysEndU-1))
      sysIndxU = (/ (i, i= sysStartU, sysEndU-1) /)

! -----------------------------------------------------------------------
! Prep V solution for time stepping loop
! -----------------------------------------------------------------------
      ! Create KSP Context for V-momentum
      CALL KSPCreate(PETSC_COMM_WORLD, ksp_v, ierr)
      CALL KSPSetOptionsPrefix(ksp_v,"v_sys_",ierr)
      CALL KSPSetFromOptions(ksp_v, ierr)
      CALL KSPSetOperators(ksp_v, Av, Av, ierr)

      ! Initialize v-RHS vector
      CALL VecCreate(PETSC_COMM_WORLD,fv,ierr) 
      CALL VecSetSizes(fv,sysDimV(my_rank),NV,ierr)
      CALL VecSetFromOptions(fv,ierr)
       
      ! Initialize v-soln vector
      CALL VecCreate(PETSC_COMM_WORLD,v_vec,ierr) 
      CALL VecSetSizes(v_vec,sysDimV(my_rank),NV,ierr)
      CALL VecSetFromOptions(v_vec,ierr)
      
      ! Allocate for body forcing Qv, solution V, and temporary v_vecF90
      ALLOCATE(Qv(rowStartV(my_rank):rowEndV(my_rank),Nx+2))
      ALLOCATE(V(rowStartV(my_rank)-1:rowEndV(my_rank)+1,Nx+2))
      ALLOCATE(v_vecF90(sysDimV(my_rank)))
      v_vecF90 = 0.D0; V = 0.D0;


      ! Get indexing used in reshaping v_vec -> v_vecF90 -> V
      CALL VecGetOwnershipRange(v_vec, sysStartV, sysEndV, ierr)
      ALLOCATE(sysIndxV(sysStartV:sysEndV-1))
      sysIndxV = (/ (i, i= sysStartV, sysEndV-1) /)
      
! -----------------------------------------------------------------------
! Prep P solution for time stepping loop
! -----------------------------------------------------------------------
      ! Create KSP Context for Pressure-Poisson
      CALL KSPCreate(PETSC_COMM_WORLD, ksp_p, ierr)
      CALL KSPSetOperators(ksp_p, Ap, Ap, ierr)
      CALL KSPSetOptionsPrefix(ksp_p,"p_sys_",ierr)
      CALL KSPSetFromOptions(ksp_p, ierr)

      ! Initialize p-RHS vector
      CALL VecCreate(PETSC_COMM_WORLD,fp,ierr) 
      CALL VecSetSizes(fp,sysDimP(my_rank),NP,ierr)
      CALL VecSetFromOptions(fp,ierr)
      
      ! Initialize p-soln vector
      CALL VecCreate(PETSC_COMM_WORLD,p_vec,ierr) 
      CALL VecSetSizes(p_vec,sysDimP(my_rank),NP,ierr)
      CALL VecSetFromOptions(p_vec,ierr)
       
      ! Allocate for solution P, and temporary p_vecF90
      ALLOCATE(P(rowStartP(my_rank)-2:rowEndP(my_rank)+2,Nx))
      ALLOCATE(p_vecF90(sysDimP(my_rank)))
      p_vecF90 = 0.D0; P = 0.D0;

      ! Get indexing used in reshaping p_vec -> p_vecF90 -> P
      CALL VecGetOwnershipRange(p_vec, sysStartP, sysEndP, ierr)
      ALLOCATE(sysIndxP(sysStartP:sysEndP-1))
      sysIndxP = (/ (i, i= sysStartP, sysEndP-1) /)
      
! -----------------------------------------------------------------------
! Print details before entering time stepping loop
! -----------------------------------------------------------------------
      IF (my_rank == master) THEN
      PRINT *,'------------------------------------------------'&
              '------------------------'
            PRINT *, 'Details prior to time stepping:'
            PRINT *, 'Nx = ', Nx, 'Ny = ', Ny
            PRINT *, 'Max Tf = ',Tchar*Tf, 'dt_dim = ', Tchar*dt
      PRINT *,'------------------------------------------------'&
              '------------------------'
      END IF

! ***********************************************************************
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! Begin Time Stepping
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! ***********************************************************************
      
      ! Start Wall Time and set steadystate = 1
      t_w0 = MPI_Wtime()
      steady = 1
      Nt = 0

DO WHILE (t < Tf) 
      ! Update time and count iterations
      t = t + dt
      Nt = Nt + 1

      ! Get B = alpha(theta)
      theta = plateletFrac(t,Tf)
      !theta = 0.D0
      CALL test_B(0*theta, 2D5, 0.6D-4, Xchar*XXPP, B)
!PRINT *, 'i = ', i, 'Nt+1 = ',Nt+1
!PRINT *, 'theta = '
!DO j = 1, Ny
!      PRINT '(20f8.3)', theta(j,:)
!END DO
!
!PRINT *, 'B = '
!DO j = 1, Ny
!      PRINT '(40F10.3)', B(j,:)
!END DO

! -----------------------------------------------------------------------
! --- Solve u-momentum equation Au * u_vec = fu
! -----------------------------------------------------------------------

      ! Forcing terms at current time  
      Qu = 0.D0

      CALL build_Bu(my_rank,N,NU,sysDimU(my_rank),rowStartU(my_rank),&
                 rowStartV(my_rank),rowStartP(my_rank),XU,XV,YV,B,Bu,ierr)
      
      ! Set Bu = Bu + Au
      CALL MatAXPY(Bu,1D0,Au,DIFFERENT_NONZERO_PATTERN,ierr)

      ! Set KSP Operators for u-momentum equation keeeping the identical
      ! preconditioner matrix for all linear solves. This approach is often
      ! effective when the linear systems do not change very much between
      ! succesive steps.
      CALL KSPSetReusePreconditioner(ksp_u,PETSC_TRUE,ierr)
      CALL KSPSetOperators(ksp_u, Bu, Au, ierr)
      
      ! Get RHS vector fu
      CALL uRHS(N,NU,rowStartU(my_rank),rowStartV(my_rank),&
                rowStartP(my_rank),XU,YU,XV,YV,t,dt,Re,U0,P0,B0,&
                Qu,Quo,NLU,NLU0,fu,ierr)
      
      ! Solve linear system for u_vec
      CALL KSPSolve(ksp_u, fu, u_vec, ierr)
      
      ! Reshape soln vector into global array
      CALL VecGetValues(u_vec, sysDimU(my_rank), sysIndxU, &
                        u_vecF90, ierr)
      CALL reshapeU(N,rowStartU(my_rank),rowEndU(my_rank),&
                    sysStartU,u_vecF90,U)
      
      ! Share rows U-soln: Send/Recv rows shared w/neighboring processor
      IF (num_cores > 1) THEN
      IF (my_rank == master) THEN 
      CALL MPI_Send(U(rowEndU(my_rank),:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 8*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      CALL MPI_Recv(U(rowEndU(my_rank)+1,:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 7*num_cores+my_rank+1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)

      ELSE IF ((my_rank > 0) .AND. (my_rank < num_cores-1)) THEN
      CALL MPI_Recv(U(rowStartU(my_rank)-1,:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 8*num_cores+my_rank-1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)
      CALL MPI_Send(U(rowStartU(my_rank),:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 7*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      CALL MPI_Send(U(rowEndU(my_rank),:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 8*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      CALL MPI_Recv(U(rowEndU(my_rank)+1,:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 7*num_cores+my_rank+1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)
      ELSE IF (my_rank == num_cores-1) THEN 
      CALL MPI_Recv(U(rowStartU(my_rank)-1,:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 8*num_cores+my_rank-1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)
      CALL MPI_Send(U(rowStartU(my_rank),:),Nx+1,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 7*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      END IF
      END IF

! -----------------------------------------------------------------------
! --- Solve v-momentum equation Av * v_vec = fv
! -----------------------------------------------------------------------
      ! Forcing terms at current time  
      Qv = 0.D0
 
      ! Brinkman Matrix Bv
      CALL build_Bv(my_rank,N,NV,sysDimV(my_rank),rowStartU(my_rank),&
                 rowStartV(my_rank),rowStartP(my_rank),XU,YU,YV,B,Bv,ierr)
      
      ! Set Bv = Bv + 1*Av
      CALL MatAXPY(Bv,1D0,Av,DIFFERENT_NONZERO_PATTERN,ierr)
      !CALL MatAXPY(Bv,1D0,Av,SAME_NONZERO_PATTERN,ierr)

      ! Set KSP Operators for v-momentum equation
      CALL KSPSetReusePreconditioner(ksp_v,PETSC_TRUE,ierr)
      CALL KSPSetOperators(ksp_v, Bv, Av, ierr)

      ! Get RHS vector fv
      CALL vRHS(N,NV,rowStartU(my_rank),rowStartV(my_rank),&
                rowStartP(my_rank),XU,YU,XV,YV,&
                t,dt,Re,V0,P0,B0,Qv,Qvo,NLV,NLV0,V_BC,fv,ierr)

      ! Solve linear system for v_vec
      CALL KSPSolve(ksp_v, fv, v_vec, ierr)

      ! Reshape soln vector into global array
      CALL VecGetValues(v_vec, sysDimV(my_rank), sysIndxV, &
                        v_vecF90, ierr)
      CALL reshapeV(N,rowStartV(my_rank),rowEndV(my_rank),&
                    sysStartV,v_vecF90,V)
      
      ! Share rows V-soln: Send/Recv rows shared w/neighboring processor
      IF (num_cores > 1) THEN
      IF (my_rank == master) THEN 
      ! Send row to my_rank+1
      CALL MPI_Send(V(rowEndV(my_rank),:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 10*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      ! Recv from my_rank+1
      CALL MPI_Recv(V(rowEndV(my_rank)+1,:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 9*num_cores+my_rank+1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)

      ELSE IF ((my_rank > 0) .AND. (my_rank < num_cores-1)) THEN
      ! Recv row from my_rank-1
      CALL MPI_Recv(V(rowStartV(my_rank)-1,:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 10*num_cores+my_rank-1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)
      ! Send row to my_rank-1
      CALL MPI_Send(V(rowStartV(my_rank),:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 9*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      ! Send row to my_rank+1
      CALL MPI_Send(V(rowEndV(my_rank),:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 10*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      ! Recv row from my_rank+1
      CALL MPI_Recv(V(rowEndV(my_rank)+1,:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank+1, 9*num_cores+my_rank+1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)
      ELSE IF (my_rank == num_cores-1) THEN 
      ! Recv row from my_rank-1
      CALL MPI_Recv(V(rowStartV(my_rank)-1,:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 10*num_cores+my_rank-1,PETSC_COMM_WORLD,&
                    mpi_status,ierr)
      ! Send row to my_rank-1
      CALL MPI_Send(V(rowStartV(my_rank),:),Nx+2,MPI_DOUBLE_PRECISION,&
                    my_rank-1, 9*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      END IF
      END IF

! -----------------------------------------------------------------------
! --- Solve pressure-Poisson equation Ap * p_vec = fp
! -----------------------------------------------------------------------
      ! Get RHS vector fp
      CALL pRHS(N,NP,rowStartU(my_rank),rowStartV(my_rank),&
                rowStartP(my_rank),rowEndP(my_rank),&
                U,V,P0,XU,YV,yp,dt,fp,ierr)

      ! Solve linear system for p_vec
      CALL KSPSolve(ksp_p, fp, p_vec, ierr)

      ! Reshape soln vector into global array
      CALL VecGetValues(p_vec, sysDimP(my_rank), sysIndxP, &
                        p_vecF90, ierr)
      CALL reshapePHI(N,rowStartP(my_rank),rowEndP(my_rank),&
                    sysStartP,p_vecF90,P)
      
      ! Share rows PHI-soln: Send/Recv rows shared w/neighboring processor
      IF (num_cores > 1) THEN
      IF (my_rank == master) THEN 
      ! Send row to my_rank+1
      CALL MPI_Send(P(rowEndP(my_rank)-1:rowEndP(my_rank),:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank+1,12*num_cores+my_rank,&
                    PETSC_COMM_WORLD,ierr)
      ! Recv from my_rank+1
      CALL MPI_Recv(P(rowEndP(my_rank)+1:rowEndP(my_rank)+2,:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank+1, 11*num_cores+my_rank+1,&
                    PETSC_COMM_WORLD,mpi_status,ierr)

      ELSE IF ((my_rank > 0) .AND. (my_rank < num_cores-1)) THEN
      ! Recv row from my_rank-1
      CALL MPI_Recv(P(rowStartP(my_rank)-2:rowStartP(my_rank)-1,:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank-1, 12*num_cores+my_rank-1,&
                    PETSC_COMM_WORLD,mpi_status,ierr)
      ! Send row to my_rank-1
      CALL MPI_Send(P(rowStartP(my_rank):rowStartP(my_rank)+1,:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank-1, 11*num_cores+my_rank,&
                    PETSC_COMM_WORLD,ierr)
      ! Send row to my_rank+1
      CALL MPI_Send(P(rowEndP(my_rank)-1:rowEndP(my_rank),:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank+1, 12*num_cores+my_rank,&
                    PETSC_COMM_WORLD,ierr)
      ! Recv row from my_rank+1
      CALL MPI_Recv(P(rowEndP(my_rank)+1:rowEndP(my_rank)+2,:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank+1, 11*num_cores+my_rank+1,&
                    PETSC_COMM_WORLD,&
                    mpi_status,ierr)
      ELSE IF (my_rank == num_cores-1) THEN 
      ! Recv row from my_rank-1
      CALL MPI_Recv(P(rowStartP(my_rank)-2:rowStartP(my_rank)-1,:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank-1, 12*num_cores+my_rank-1,&
                    PETSC_COMM_WORLD,mpi_status,ierr)
      ! Send row to my_rank-1
      CALL MPI_Send(P(rowStartP(my_rank):rowStartP(my_rank)+1,:),2*Nx,&
                    MPI_DOUBLE_PRECISION,my_rank-1, 11*num_cores+my_rank,&
                    PETSC_COMM_WORLD,ierr)
      END IF
      END IF

! -----------------------------------------------------------------------
!--- Correction for Fluid and Pressure -- Check Steady State 
! -----------------------------------------------------------------------
      CALL correction(N,dt,xp,yp,rowStartU(my_rank),rowEndU(my_rank),&
                      rowStartV(my_rank),rowEndV(my_rank),&
                      rowStartP(my_rank),rowEndP(my_rank),P0,U,V,P,v_BC)
      
      ! Check steady state
      steady = MAX(MAXVAL(U-U0),MAXVAL(V-V0))
      CALL MPI_ALLreduce(steady,steady,1,MPI_DOUBLE_PRECISION,&
                      MPI_MAX,PETSC_COMM_WORLD,ierr)
      IF ((MOD(Nt,1000) == 0) .AND. my_rank == master) THEN
            PRINT *, 'Steadystate:'
            PRINT *,'t_dim = ', Tchar*t,' steady = ', steady
      END IF
      IF (steady <= tol) THEN
            EXIT
      END IF

! -----------------------------------------------------------------------
!--- Updates for Fluid, Pressure, Nonlinear, and Brinkman Terms 
! -----------------------------------------------------------------------
      !Update Prev. Velocities and Pressure
      U0 = U; V0 = V; P0 = P;
      
      ! Update Prev. Forcing Terms
      Quo = Qu; Qvo = Qv; 
      
      ! Update Prev. Nonlinear Terms
      NLU0 = NLU; NLV0 = NLV; ! Update Prev. Nonlinear Terms
      
      ! Calculate new Nonlinear terms
      CALL nonlin_U(N,rowStartU(my_rank),rowEndU(my_rank),&
              rowStartV(my_rank),XU,YU,XV,YV,U0,V0,NLU)
      CALL nonlin_V(N,rowStartU(my_rank),rowStartV(my_rank),&
                    rowEndV(my_rank),XU,YU,XV,YV,U0,V0,NLV)

      ! Update Brinkman term
      B0 = B; 
      CALL MatZeroEntries(Bu,ierr)
      CALL MatZeroEntries(Bv,ierr)

END DO

! ***********************************************************************
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! END Time Stepping
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! ***********************************************************************
      
      ! Final CPU Time
      t_wf = MPI_Wtime()
      t_totF = MPI_Wtime()
      
      ! Print Results
      CALL MPI_Barrier(PETSC_COMM_WORLD,ierr)
      IF (my_rank == master) THEN
      PRINT *,'------------------------------------------------'&
              '------------------------'
      PRINT *, 'Num_cores = ', num_cores
      PRINT *, 'Physical Params: Re = ', Re, ', rho = ', rho,&
               ', mu = ', mu 
      PRINT *, 'Spatial Params: Nx = ',Nx, ', Ny = ',Ny,', GRID = ',GRID
      PRINT *, 'Temporal Params: Tf = ', Tchar*Tf, ', dt = ',&
               Tchar*dt
      PRINT *, 'Total time: CPU time = ',t_totf-t_tot0

      PRINT *, 'Avg CPU/dt = ',(t_wf-t_w0)/Nt
      PRINT *,'------------------------------------------------'&
              '------------------------'
      END IF

! -----------------------------------------------------------------------
! Write Solutions to .txt file
! -----------------------------------------------------------------------
      IF ((save_soln > 0) .AND. (my_rank == master)) THEN
            ! Allocate full vars for saving
            ALLOCATE(U_full(Ny+2,Nx+1),XU_full(Ny+2,Nx+1),&
                    YU_full(Ny+2,Nx+1),&
                    V_full(Ny+1,Nx+2), XV_full(Ny+1,Nx+2),&
                    YV_full(Ny+1,Nx+2),&
                    P_full(Ny,Nx), XP0(Ny,Nx), YP0(Ny,Nx), B_full(Ny,Nx) )
            U_full = 0; XU_full = 0; YU_full = 0;
            V_full = 0; XV_full = 0; YV_full = 0;
            P_full = 0; XP0 = 0;     YP0 = 0;
            
            ! Create full Grid
            CALL H_grid(N,GRID,XU_full,YU_full,&
                        XV_full,YV_full,xp,yp,min_h)
            CALL meshgrid(xp(2:Nx+1), yp(2:Ny+1), XP0, YP0)

            ! Get solution from master core
            U_full(rowStartU(my_rank):rowEndU(my_rank),:) = &
                    U(rowStartU(my_rank):rowEndU(my_rank),:)
            V_full(rowStartV(my_rank):rowEndV(my_rank),:) = &
                    V(rowStartV(my_rank):rowEndV(my_rank),:)
            P_full(rowStartP(my_rank):rowEndP(my_rank),:) = &
                    P(rowStartP(my_rank):rowEndP(my_rank),:)
            B_full(rowStartP(my_rank):rowEndP(my_rank),:) = &
                    B(rowStartP(my_rank):rowEndP(my_rank),:)
            
            ! Receive Solutions from other cores
            IF (num_cores > 1) THEN
            DO i = 1, num_cores-1
            CALL MPI_Recv(U_full(rowStartU(i):rowEndU(i),:),&
                      dimU(i)*(Nx+1),MPI_DOUBLE_PRECISION,i,&
                    13*num_cores+i,PETSC_COMM_WORLD,mpi_status,ierr)
            CALL MPI_Recv(V_full(rowStartV(i):rowEndV(i),:),&
                      dimV(i)*(Nx+2),MPI_DOUBLE_PRECISION,i,&
                    17*num_cores+i,PETSC_COMM_WORLD,mpi_status,ierr)
            CALL MPI_Recv(P_full(rowStartP(i):rowEndP(i),:),&
                      dimP(i)*(Nx),MPI_DOUBLE_PRECISION,i,&
                    21*num_cores+i,PETSC_COMM_WORLD,mpi_status,ierr)
            CALL MPI_Recv(B_full(rowStartP(i):rowEndP(i),:),&
                      dimP(i)*(Nx),MPI_DOUBLE_PRECISION,i,&
                    25*num_cores+i,PETSC_COMM_WORLD,mpi_status,ierr)
            END DO
            END IF

            ! Write Solution and grid to file
            IF (save_soln > 1) THEN
            CALL write_to_file('U',Uchar*U_full,Nx,Ny,Tchar*Tf,Ny+2,&
                               111,GRID)
            CALL write_to_file('XU',Xchar*XU_full,Nx,Ny,Tchar*Tf,Ny+2,&
                               112,GRID)
            CALL write_to_file('YU',Xchar*YU_full,Nx,Ny,Tchar*Tf,Ny+2,&
                               113,GRID)
            CALL write_to_file('V',Uchar*V_full,Nx,Ny,Tchar*Tf,Ny+1,&
                               115,GRID)
            CALL write_to_file('XV',Xchar*XV_full,Nx,Ny,Tchar*Tf,Ny+1,&
                               116,GRID)
            CALL write_to_file('YV',Xchar*YV_full,Nx,Ny,Tchar*Tf,Ny+1,&
                               117,GRID)
            CALL write_to_file('P',Pchar*P_full,Nx,Ny,Tchar*Tf,Ny,&
                               119,GRID)
            CALL write_to_file('B',Bchar*B_full,Nx,Ny,Tchar*Tf,Ny,&
                               122,GRID)
            CALL write_to_file('XP',Xchar*XP0,Nx,Ny,Tchar*Tf,Ny,120,GRID)
            CALL write_to_file('YP',Xchar*YP0,Nx,Ny,Tchar*Tf,Ny,121,GRID)
            
            ! Write N to file
            fileName ='data/N_Nx'//num2str(Nx)//'Ny'//num2str(Ny)//&
                      'Grid'//num2str(GRID)//'.txt'
            OPEN(123,FILE=fileName,ACTION='write',STATUS='replace')
            WRITE(123,*) N
            CLOSE(123)

            ! Write Lengths to file
            fileName ='data/lengths_Nx'//num2str(Nx)//'Ny'//num2str(Ny)//&
                      'Grid'//num2str(GRID)//'.txt'
            OPEN(124,FILE=fileName,ACTION='write',STATUS='replace')
            WRITE(124,*) Xchar*(/ WL, Wmid, WR, HL, Hmid, HU /)
            CLOSE(124)
            END IF
            
            ! Write output to file
            fileName ='data/mio_output_'//num2str(num_cores)//&
                   '_Nx'//num2str(Nx)//'Grid'//num2str(GRID)//'.txt'
            OPEN(101,FILE=fileName,ACTION='write',STATUS='replace')
            WRITE(101,*) 'Num_cores = ', num_cores
            WRITE(101,*) 'Spatial Params: Nx = ',Nx, ', Ny = ',Ny,&
                     ', GRID = ',GRID
            WRITE(101,*) 'Temporal Params: Tf = ',Tchar*Tf,',dt = ',& 
                          Tchar*dt 
            WRITE(101,*) 'Total time: CPU time = ',t_totf-t_tot0

            WRITE(101,*) 'Time loop only: CPU time = ',t_wf-t_w0,&
                'Avg CPU/dt = ',(t_wf-t_w0)/Nt
            CLOSE(101)
            
            DEALLOCATE(U_full,XU_full,YU_full)
            DEALLOCATE(V_full,XV_full,YV_full)
            DEALLOCATE(P_full,XP0,YP0,B_full)
      
      ! Send solutions from other cores to master for saving
      ELSE IF ((my_rank > master) .AND. (save_soln > 0)) THEN
            CALL MPI_Send(U(rowStartU(my_rank):rowEndU(my_rank),:),&
                      dimU(my_rank)*(Nx+1),MPI_DOUBLE_PRECISION,&
                   master, 13*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
            CALL MPI_Send(V(rowStartV(my_rank):rowEndV(my_rank),:),&
                      dimV(my_rank)*(Nx+2),MPI_DOUBLE_PRECISION,&
                   master, 17*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
            CALL MPI_Send(P(rowStartP(my_rank):rowEndP(my_rank),:),&
                      dimP(my_rank)*(Nx),MPI_DOUBLE_PRECISION,&
                   master, 21*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
            CALL MPI_Send(B(rowStartP(my_rank):rowEndP(my_rank),:),&
                      dimP(my_rank)*(Nx),MPI_DOUBLE_PRECISION,&
                   master, 25*num_cores+my_rank,PETSC_COMM_WORLD,ierr)
      END IF

! -----------------------------------------------------------------------
! Deallocation
! -----------------------------------------------------------------------
      ! Grid Variables
      DEALLOCATE(XU,YU,XV,YV,xp,yp,XXPP,YYPP)

      ! Indexing Variables
      DEALLOCATE(dimU,rowStartU,rowEndU,&
                dimV,rowStartV,rowEndV,&
                dimP,rowStartP,rowEndP)
      DEALLOCATE(sysDimU, sysIndxU)
      DEALLOCATE(sysDimV, sysIndxV)
      DEALLOCATE(sysDimP, sysIndxP)
      
      ! Intermediate Soln Vectors
      DEALLOCATE(u_vecF90,v_vecF90,p_vecF90)
      
      ! Soln's
      DEALLOCATE(U,U0,V,V0,P,P0)
      
      ! BC's
      IF (my_rank == num_cores -1) THEN
            DEALLOCATE(V_BC)
      END IF

      ! Nonlinear Terms
      DEALLOCATE(NLU,NLU0,NLV,NLV0)
      
      ! Brinkman Terms
      DEALLOCATE(B,B0,theta)

      ! Body Forcing
      DEALLOCATE(Quo,Qu,Qvo,Qv)
      
      ! PETSc vectors, matrices and ksp contexts
      CALL VecDestroy(u_vec, ierr)
      CALL VecDestroy(fu,ierr)
      CALL KSPDestroy(ksp_u, ierr)
      Call MatDestroy(Au,ierr)
      Call MatDestroy(Bu,ierr)
      
      CALL VecDestroy(v_vec, ierr)
      CALL VecDestroy(fv,ierr)
      CALL KSPDestroy(ksp_v, ierr)
      Call MatDestroy(Av,ierr)
      Call MatDestroy(Bv,ierr)
      
      CALL VecDestroy(p_vec, ierr)
      CALL VecDestroy(fp,ierr)
      CALL KSPDestroy(ksp_p, ierr)
      Call MatDestroy(Ap,ierr)
      

! -----------------------------------------------------------------------
! End of Code
! -----------------------------------------------------------------------
      CALL PetscFinalize(ierr)
!CALL PetscFinalize(ierr)
!CALL EXIT(ierr)
      END SUBROUTINE fluids_solver
      END MODULE fluids_mod
