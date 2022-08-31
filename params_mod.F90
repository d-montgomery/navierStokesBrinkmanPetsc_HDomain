      MODULE params_mod
! Necessary PETSc Libraries
#include <petsc/finclude/petscsys.h>  
      USE petscsys

      IMPLICIT NONE

! --- Mathematical Constants ---------------------------------------------
      PetscScalar, PARAMETER :: pi = ACOS(-1.D0)

! --- Physical Parameters ------------------------------------------------
      PetscScalar :: mu, rho                                          

! --- Lengths for the H Domain -------------------------------------------
      PetscScalar :: WL, Wmid, WR, HL, Hmid, HU    

! --- Declare Flow Rates and Pressure at Inlet/Outlets
      PetscScalar :: Q1, QM, Q3, Pw_in, Pw_out, Pb_in, Pb_out

! --- Characteristic Scales for Non-Dimensionalaization ------------------
      PetscScalar :: Uchar, Xchar, Pchar, Tchar, Fchar, Bchar, Re    

CONTAINS

! -----------------------------------------------------------------------
! Subroutine: Set parameter values and physical constants as GLOBAL vars
! -----------------------------------------------------------------------
      SUBROUTINE setparams

      ! Lengths of the H Domain (all lengths in cm, must be doubles)
      WL = 1.D-2; Wmid = 1.5D-2; WR = 1.D-2;
      HL = 0.75D-2; Hmid = 0.2D-2; HU = 0.75D-2;
      
      ! Physical Parameters (fluid)
      mu = 2.62507D-2                    !g/cm/s
      rho = 1.D0                         !g/cm^3

      ! Volumetric Flow Rates at Upper Left and Upper Right Boundaries
      Q1 = 2.77778D-2 ! Left Channel cm^2 / s
      Q3 = 1.55556D-2 ! Right Channel cm^2 /s
      QM = 4.25D-3    ! m^2 / s

      ! Initial inlet/outlet pressure from HCA g/cm/s^2 = 10 Pa
      Pw_in = 10.3471D1  ! g/cm/s^2 : Pressure Inlet (wash)
      Pw_out = 0.D0      ! g/cm/s^2 : Pressure Outlet (wash)
      Pb_in = 294.076D1  ! g/cm/s^2 :Pressure Inlet (blood)
      Pb_out = 284.213D1 ! g/cm/s^2 : Pressure Outlet (blood)
      
      ! Characteristic Scales
      Xchar = 3.D-4                           !cm
      Uchar = Q3 / WR                         !cm/s
      Tchar = Xchar / Uchar                   !s
      Pchar = rho * Uchar**2                  !g/cm/s^2
      Bchar = rho * Uchar / (mu * Xchar)      !1/cm^2
      Re = rho * Uchar * Xchar / mu           !1
      
      ! Non Dimensionalization (DON'T EDIT THIS!)
      WL = WL / Xchar; Wmid = Wmid / Xchar; WR = WR / Xchar;
      HL = HL / Xchar; Hmid = Hmid / Xchar; HU = HU / Xchar;
      
      Pw_in = Pw_in / Pchar; Pw_out = Pw_out / Pchar;
      Pb_in = Pb_in / Pchar; Pb_out = Pb_out / Pchar;
      END SUBROUTINE setparams

      END MODULE params_mod
