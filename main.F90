      PROGRAM main 
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

! Necessary .F90 Modules
      USE params_mod
      USE fluids_mod

      IMPLICIT NONE
      ! Declare Flags for Saving, etc
      PetscInt :: save_soln

      ! Declare Spatial Variables
      PetscInt ::  Nx, Ny, GRID
      
      ! Declare Temporal Variables
      PetscScalar :: Tf

      ! Declare the steady state toleracne
      PetscScalar :: tol

! --- Start Main Program ------------------------------------------------
      
      ! Write Solutions to File: 0 = no,
      !                          1 = save output,
      !                          2 = save output and solutions as txt
      save_soln = 2
      
      ! Determine tolerance for Steady State.  Must be a Double
      tol = 1D-10

      ! Number of Cells in x and y direction 
      Nx = 256 !12
      Ny = 256 !10

      ! Specify Grid Type: 0 = uniform, 1 = non-uniform
      GRID = 1

      ! Final time and temporal step size
      Tf = 1.D0

      CALL fluids_solver(Nx,Ny,GRID,Tf,tol,save_soln)

      END PROGRAM main
