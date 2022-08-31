      MODULE soln_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

! Necessary F90 Modules
      USE params_mod
      
      IMPLICIT NONE

CONTAINS
! --- Exact solution for x-comp of velocity: u -------------------------
      SUBROUTINE u_exact(x,y,t,w,u) 
      IMPLICIT NONE
      PetscScalar, DIMENSION(:,:) :: x, y
      PetscScalar :: t, w
      PetscScalar, DIMENSION(:,:) :: u
      
      u = SIN(x) * COS(y) * COS(w*t)/Uchar
      
      END SUBROUTINE u_exact

! --- Exact solution for y-comp of velocity: v -------------------------
      SUBROUTINE v_exact(x,y,t,w,v) 
      IMPLICIT NONE
      PetscScalar, DIMENSION(:,:) :: x, y
      PetscScalar :: t, w
      PetscScalar, DIMENSION(:,:) :: v

      v = -COS(x) * SIN(y) * COS(w*t)/Uchar
      
      END SUBROUTINE v_exact      

! --- Exact solution for pressure: p -----------------------------------
      SUBROUTINE p_exact(x,y,t,w,p) 
      IMPLICIT NONE
      PetscScalar, DIMENSION(:,:) :: x, y
      PetscScalar :: t, w
      PetscScalar, DIMENSION(:,:) :: p

      p = SIN(x) * SIN(y) * COS(w*t)/Pchar
      
      END SUBROUTINE p_exact      

! --- Non-Dimensional Body Forcing of U-Momentum: Q_u ------------------
      SUBROUTINE u_bodyForce(x,y,t,w,f)
      IMPLICIT NONE
      PetscScalar, DIMENSION(:,:) :: x, y
      PetscScalar :: t, w
      PetscScalar, DIMENSION(:,:) :: f
      
      f = 1.D0/Fchar * ( -rho*w*sin(x)*cos(y)*sin(w*t) &
         + rho*(cos(y)*sin(x)*cos(w*t)*cos(y)*cos(x)*cos(w*t) &
         + cos(x)*sin(y)*cos(w*t)*sin(y)*sin(x)*cos(w*t)) &
         + cos(x)*sin(y)*cos(w*t) + 2.D0*mu*cos(y)*sin(x)*cos(w*t))
      
      END SUBROUTINE u_bodyForce

! --- Non-Dimensional Body Forcing of V-Momentum: Q_v ------------------
      SUBROUTINE v_bodyForce(x,y,t,w,f)
      IMPLICIT NONE
      PetscScalar, DIMENSION(:,:) :: x, y
      PetscScalar :: t, w
      PetscScalar, DIMENSION(:,:) :: f

      f = 1.D0/Fchar*( rho*w*cos(x)*sin(y)*sin(w*t) &
          +rho*(cos(y)*sin(x)*cos(w*t)*sin(x)*sin(y)*cos(w*t) &
          + cos(x)*sin(y)*cos(w*t)*cos(x)*cos(y)*cos(w*t)) &
          + sin(x)*cos(y)*cos(w*t) - 2.D0*mu*cos(x)*sin(y)*cos(w*t))
      
      END SUBROUTINE v_bodyForce

! --- Non-Dimensional du/dy for BC -------------------------------------
      SUBROUTINE u_y(x,y,t,w,BC)
      IMPLICIT NONE
      PetscScalar, DIMENSION(:) :: x, y
      PetscScalar :: t, w
      PetscScalar, DIMENSION(:) :: BC
     
      BC = 1.D0/Uchar * ( -sin(y)*sin(x)*cos(w*t) )
      END SUBROUTINE u_y

! --- Non-Dimensional dv/dy for BC -------------------------------------
      SUBROUTINE v_y(x,y,t,w,BC)
      IMPLICIT NONE
      PetscScalar, DIMENSION(:) :: x, y
      PetscScalar :: t, w
      PetscScalar, DIMENSION(:) :: BC
     
      BC = 1.D0/Uchar * ( -cos(x)*cos(y)*cos(w*t) )
      END SUBROUTINE v_y
      
      END MODULE soln_mod
