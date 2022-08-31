      MODULE utilities_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

      IMPLICIT NONE

CONTAINS ! Functions often used in Matlab and helpful I/O functions

! -----------------------------------------------------------------------
! FUNCTION: Convert integer to string for I/O file names           
! Full disclosure: I looked this up on the internet 
! -----------------------------------------------------------------------
      FUNCTION num2str(i) RESULT(res)
      CHARACTER (:), ALLOCATABLE :: res
      INTEGER, INTENT(IN) :: i
      CHARACTER (range(i)+2) :: tmp
      WRITE(tmp,'(i0)') i
      res = TRIM(tmp)
      END FUNCTION

! -----------------------------------------------------------------------
! SUBROUTINE: linspace - return a vector of length n with uniform 
!                      spacing from a to b
! -----------------------------------------------------------------------
      SUBROUTINE linspace(a,b,n,x) 
      PetscScalar, INTENT(IN) :: a,b
      PetscInt, INTENT(IN) :: n
      PetscScalar, DIMENSION(:) :: x
      PetscInt :: i
      PetscScalar :: dx

      dx = (b - a) / (n - 1)
      x = (/ (a + (i-1) * dx, i = 1, n) /)

      END SUBROUTINE linspace
! -----------------------------------------------------------------------
! SUBROUTINE: Meshgrid subroutine 
! -----------------------------------------------------------------------
      SUBROUTINE meshgrid(x,y,XX,YY)
      PetscScalar, DIMENSION(:), INTENT(IN) :: x, y
      PetscScalar, DIMENSION(:,:), INTENT(OUT):: XX, YY
      PetscInt :: Nx, Ny, i, j

      Nx = SIZE(x)
      Ny = SIZE(y)
     
      ! Create Ny Copies of x
      DO i = 1, Ny
            XX(i,:) = x
      END DO

      ! Create Nx Copies of y
      DO j = 1, Nx
            YY(:,j) = y
      END DO
      END SUBROUTINE meshgrid

! -----------------------------------------------------------------------
! SUBROUTINE: fliplr subroutine 
! -----------------------------------------------------------------------
      SUBROUTINE fliplr(x,y)
      PetscScalar, DIMENSION(:), INTENT(IN) :: x
      PetscScalar, DIMENSION(:), INTENT(OUT) :: y
      PetscInt :: N, i, j

      N = SIZE(x)
      j = N
      DO i = 1,N
            y(i) = x(j)
            j = j - 1
      END DO
      END SUBROUTINE fliplr

! -----------------------------------------------------------------------
! Function: Heaviside Function
! -----------------------------------------------------------------------
      FUNCTION heaviside(x) RESULT(h)
      PetscScalar, DIMENSION(:,:) :: x
      PetscScalar, DIMENSION(SIZE(x,1),SIZE(x,2)) :: h

      h = 0.5D0 * (SIGN(1.D0,x) + 1)

      END FUNCTION heaviside

! -----------------------------------------------------------------------
! SUBROUTINE:  Writing arrays to file
! fileName = data/'var_name'_Nx'int'Ny'int'Tf'int'.txt
! -----------------------------------------------------------------------
      SUBROUTINE write_to_file(var_name,var,Nx,Ny,Tf,rows,io,GRID)
      CHARACTER (LEN=*) :: var_name
      PetscScalar, DIMENSION(:,:) :: var
      INTEGER :: Nx, Ny
      PetscScalar :: Tf
      INTEGER :: rows, io, GRID, i
      CHARACTER (LEN = 40) :: fileName
      
      fileName = 'data/'//TRIM(var_name)//'_Nx'//num2str(Nx)//&
                 'Ny'//num2str(Ny)//'Tf'//num2str(NINT(Tf))//&
                 'Grid'//num2str(GRID)//'.txt'
      
      OPEN(io, FILE = fileName, ACTION = 'write',STATUS = 'replace')
      DO i = 1, rows
            WRITE(io,*) var(i,:)
      END DO
      CLOSE(io)
      END SUBROUTINE write_to_file
      

      END MODULE utilities_mod
