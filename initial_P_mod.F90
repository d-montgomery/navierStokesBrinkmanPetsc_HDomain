      MODULE initial_P_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

! Necessary F90 Modules
      USE params_mod
      USE utilities_mod      
      IMPLICIT NONE

CONTAINS
! -----------------------------------------------------------------------
! Function: Initial Pressure in Wash Channel (Scalar Version) 
! -----------------------------------------------------------------------
      FUNCTION p0_w(y) RESULT(p) 
      IMPLICIT NONE
      PetscScalar, INTENT(IN) :: y ! Dimensionless
      PetscScalar :: p ! Pa
      
      p = ( Pchar*(Pw_in - Pw_out)/(Xchar*(HL + Hmid + HU)) &
             * Xchar*y ) + Pchar*Pw_out; 
     
      END FUNCTION p0_w

! -----------------------------------------------------------------------
! Subroutine: Initial Pressure in Wash Channel (Array Version) 
! -----------------------------------------------------------------------
      SUBROUTINE p2_w(y,p)
      IMPLICIT NONE
      PetscScalar, DIMENSION(:,:), INTENT(IN) :: y ! Dimensionless
      PetscScalar, DIMENSION(:,:), INTENT(OUT) :: p ! Pa

      p = ( Pchar*(Pw_in - Pw_out)/(Xchar*(HL + Hmid + HU)) &
             * Xchar*y ) + Pchar*Pw_out; 
     
      END SUBROUTINE p2_w

! -----------------------------------------------------------------------
! Function: Initial Pressure in Blood Channel (Scalar Version) 
! -----------------------------------------------------------------------
      FUNCTION p0_b(y) RESULT(p) 
      IMPLICIT NONE
      PetscScalar, INTENT(IN) :: y ! Dimensionless
      PetscScalar :: p ! Pa
      
      p = ( Pchar*(Pb_in - Pb_out)/(Xchar*(HL + Hmid + HU)) &
             * Xchar*y ) + Pchar*Pb_out     
      END FUNCTION p0_b

! -----------------------------------------------------------------------
! Subroutine: Initial Pressure in Blood Channel (Array Version) 
! -----------------------------------------------------------------------
      SUBROUTINE p2_b(y,p) 
      IMPLICIT NONE
      PetscScalar, DIMENSION(:,:), INTENT(IN) :: y ! Dimensionless
      PetscScalar, DIMENSION(:,:), INTENT(OUT) :: p ! Pa
      
      p = ( Pchar*(Pb_in - Pb_out)/(Xchar*(HL + Hmid + HU)) &
             * Xchar*y ) + Pchar*Pb_out     
      END SUBROUTINE p2_b

! -----------------------------------------------------------------------
! Subroutine: Initial Pressure in Blood and Wash Channels (Array) 
! -----------------------------------------------------------------------
      SUBROUTINE Pbw(x,y,p)
      PetscScalar, DIMENSION(:,:), INTENT(IN) :: x, y ! Dimensionless
      PetscScalar, DIMENSION(:,:), INTENT(OUT):: p ! Pa
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: pw, pb

      ALLOCATE(pw(SIZE(x,1),SIZE(x,2)),pb(SIZE(x,1),SIZE(x,2)))

      CALL p2_w(y,pw)
      CALL p2_b(y,pb)

      p = heaviside(Xchar*WL-Xchar*x)*pw &
               + heaviside(Xchar*x-Xchar*(WL+Wmid))*pb
      
      DEALLOCATE(pw,pb)

      END SUBROUTINE

! -----------------------------------------------------------------------
! Subroutine: Initial Pressure in the Injury Channel (Array) 
! -----------------------------------------------------------------------
      SUBROUTINE P_in(x,y,c,p)
      PetscScalar, DIMENSION(:,:), INTENT(IN) :: x, y ! Dimensionless
      PetscScalar, DIMENSION(:), INTENT(IN) :: c
      PetscScalar, DIMENSION(:,:), INTENT(OUT):: p ! Pa
      
      p = c(1) + c(2)*Xchar*x + c(3)*Xchar*y
      END SUBROUTINE


! -----------------------------------------------------------------------
! Subroutine: Putting it all together, Initial Pressure in H Domain 
! -----------------------------------------------------------------------
      SUBROUTINE initial_P(N,XP,YP,i_locStart,i_locEnd,P0) 
      IMPLICIT NONE
      ! Inputs
      PetscInt, DIMENSION(:), INTENT(IN) :: N
      PetscInt :: i_locStart, i_locEnd
      PetscScalar, DIMENSION(i_locStart-2:,:), INTENT(IN) :: XP, YP
      PetscScalar, DIMENSION(i_locStart-2:,:), INTENT(OUT) :: P0
      ! Variables
      PetscInt :: NxL, NxM, Nx, NyL, NyM, Ny
      PetscScalar, ALLOCATABLE, DIMENSION(:,:) :: A,AA
      PetscScalar, ALLOCATABLE, DIMENSION(:) :: c, bb, b
      PetscScalar :: denom

      ! Unpack Various N values for H-domain
      NxL = N(1); NxM = N(2); Nx = N(3);
      NyL = N(4); NyM = N(5); Ny = N(6);
      
      ! Get the initial pressure in blood and wash channels
      CALL Pbw(XP,YP,P0)
      
      ! Use Least Squares to calculate the pressure in the injury channel
      ! This task is only completed by nodes with cells in the injury chnl.
      IF ((i_locEND + 2 .GE. NyL+1) .AND. (i_locStart - 2 .LE. NyM)) THEN
            ALLOCATE(AA(4,3),A(3,3),bb(4),b(3),c(3))

            ! Build Matrix AA and vector bb for the normal equations
            ! AA^T * AA * c = AA^T * b

            AA(1,:) = Xchar*(/ 1/Xchar, WL       , HL /)
            AA(2,:) = Xchar*(/ 1/Xchar, WL       , HL + Hmid /)
            AA(3,:) = Xchar*(/ 1/Xchar, WL + Wmid, HL /)
            AA(4,:) = Xchar*(/ 1/Xchar, WL + Wmid, HL + Hmid /)

            ! RHS Vector for AA * c = bb
            bb = (/ p0_w(HL), p0_w(HL + Hmid), p0_b(HL), p0_b(HL + Hmid) /)

            ! Setup Normal equations, A*c = b, with A = AA^T * AA, and 
            ! vector b = AA^T * bb
            A = MATMUL(TRANSPOSE(AA),AA)
            b = MATMUL(TRANSPOSE(AA),bb)
      
            DEALLOCATE(AA,bb)

            ! Solve the normal equations exactly for c
            denom = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
                   - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
                   + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

            c(1) =   A(1,2)*A(2,3)*b(3) - A(1,3)*A(2,2)*b(3) &
                   - A(1,2)*A(3,3)*b(2) + A(1,3)*A(3,2)*b(2) &
                   + A(2,2)*A(3,3)*b(1) - A(2,3)*A(3,2)*b(1)
            c(2) = - A(1,1)*A(2,3)*b(3) + A(1,3)*A(2,1)*b(3) &
                   + A(1,1)*A(3,3)*b(2) - A(1,3)*A(3,1)*b(2) &
                   - A(2,1)*A(3,3)*b(1) + A(2,3)*A(3,1)*b(1)
            c(3) =   A(1,1)*A(2,2)*b(3) - A(1,2)*A(2,1)*b(3) &
                   - A(1,1)*A(3,2)*b(2) + A(1,2)*A(3,1)*b(2) &
                   + A(2,1)*A(3,2)*b(1) - A(2,2)*A(3,1)*b(1)
            c = c / denom 

            ! Use c from normal eqts to create initial pressure in inj. chnl
            CALL P_in(XP(MAX(NyL+1,i_locStart-2):MIN(NyM,i_locEnd+2),&
                         NxL+1:NxM),&
                      YP(MAX(NyL+1,i_locStart-2):MIN(NyM,i_locEnd+2),&
                         NxL+1:NxM),&
                      c,&
                      P0(MAX(NyL+1,i_locStart-2):MIN(NyM,i_locEnd+2),&
                         NxL+1:NxM))
            
            DEALLOCATE(A,b,c)

      END IF  ! End Least Squares in Injury Channel
      
      ! All cores non-dimensionalize
      P0 = P0 / Pchar       
      
      END SUBROUTINE initial_P

      
      END MODULE initial_P_mod
