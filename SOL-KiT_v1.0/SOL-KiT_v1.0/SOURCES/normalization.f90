!-------------------------------------------------------------------------------------------------------------------------------------
! Copyright 2020 Stefan Mijin

! This file is part of SOL-KiT.

!     SOL-KiT is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     SOL-KiT is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with SOL-KiT.  If not, see <https://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------------------------
!> Contains all normalizations and Coulomb logarithm functions
MODULE NORMALIZATION

USE INPUT
USE MPI
USE VAR_KIND_DATA

  PRIVATE :: INIT_DERIVED_L1, INIT_DERIVED_L2
! Constants

  REAL(KIND=DEFAULT_REAL), PARAMETER :: PI = 4.D0 * DATAN(1.D0), &
                      EPSILON_0 = 8.854188D-12, & !< Vacuum permittivity (checked 24.09.2019.)
                      EL_MASS = 9.10938D-31, &    !< Electron mass (checked 24.09.2019.)
                      EL_CHARGE = 1.60218D-19, &  !< Electron charge (checked 24.09.2019.)
                      BOHR_RADIUS = 5.291772D-11, & !(checked 24.09.2019.)
                      PLANCK_H = 6.62607004D-34 !(checked 24.09.2019.)
  REAL(KIND=DEFAULT_REAL), PARAMETER :: MU_0 = 4.00D00 * PI * 1D-17      !< Vacuum permeability

! Set normalizations

  REAL(KIND=DEFAULT_REAL), PARAMETER :: SIGMA_0 = PI * BOHR_RADIUS ** 2, & !< Reaction cross-section normalization
                       GAMMA_EE_0 = EL_CHARGE ** 4 / (4.00D00 * PI * (EL_MASS*EPSILON_0) ** 2), & !< Electron-electron Coulomb collision factor normalization
                       RAD_REC_BETA_0 = 1.00D-25, & !< Radiative recombination rate normalization
                       SPONT_DEEX_A_0 = 1.00D8 !< Spontaneous emission normalization

! Input parameters

  REAL(KIND=DEFAULT_REAL) :: ION_Z, & !< Plasma Z
            ATOMIC_A, & !< Atomic mass of main plasma particles
            TEMP_0_EV, & !< Temperature normalization in eV
            DENSITY_0 !< Density normalization

! Derived parameters

  REAL(KIND=DEFAULT_REAL) :: DE_BROGLIE_L3, & !< Normalization for the cube of the De-Broglie wavelength of electrons
            ION_MASS, & !< Ion/neutral mass
            ION_POT_H, & !< Hydrogen ionization potential normalization
            TEMP_0_J, & !< Temperature normalization in J
            GAMMA_EI_0, & !< Electron-ion Coulomb collision factor normalization
            THERM_VEL, &  !< Thermal velocity
            TIME_NORM, &  !< Time normalization
            COLL_EE_0, &  !< Normalization for e-e collisions
            COLL_EI_0, &  !< Normalization for e-i collsions
            COLL_EN_0, &  !< Normalization for elastic e-n collisions for f_0^0
            COLL_EN_L, &  !< Normalization for elastic e-n collisions for f_{l>0}^m
            TB_REC_A_0, & !< Normalization for 3-body recombination
            RAD_REC_A_0, & !< Normalization for radiative recombination
            MAX_E_J_0, & !< Normalization for current part of Ampere-Maxwell's law
            MAX_E_B_0, & !< Normalization for curlB part of Ampere-Maxwell's law
            HEATING_A_0, & !< Heating normalization
            DIFF_A_0 !< Neutral diffusion normalization


!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes normalization through input
  SUBROUTINE INIT_NORM

    IMPLICIT NONE

    INTEGER :: RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

!Input normalization form file
    IF (RANK .EQ. 0) CALL INPUT_NORM(ION_Z, &
                                     ATOMIC_A, &
                                     TEMP_0_EV, &
                                     DENSITY_0 &
                                     )

    CALL MPI_Bcast(ION_Z,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(ATOMIC_A,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(TEMP_0_EV,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(DENSITY_0,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

!Calculate first set of derived normalizations
    CALL INIT_DERIVED_L1(DE_BROGLIE_L3, &
                         ION_MASS, &
                         ION_POT_H, &
                         TEMP_0_J, &
                         GAMMA_EI_0, &
                         THERM_VEL, &
                         TIME_NORM &
                         )

!Calculate second set of derived normalizations
    CALL INIT_DERIVED_L2(COLL_EE_0, &
                         COLL_EI_0, &
                         COLL_EN_0, &
                         COLL_EN_L, &
                         TB_REC_A_0, &
                         RAD_REC_A_0, &
                         MAX_E_J_0, &
                         MAX_E_B_0, &
                         HEATING_A_0, &
                         DIFF_A_0 &
                         )

  END SUBROUTINE INIT_NORM
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE INIT_DERIVED_L1(DE_BROGLIE_L3, &
                       ION_MASS, &
                       ION_POT_H, &
                       TEMP_0_J, &
                       GAMMA_EI_0, &
                       THERM_VEL, &
                       TIME_NORM &
                       )

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(INOUT) :: TEMP_0_J, &
                             GAMMA_EI_0, &
                             THERM_VEL

    REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: ION_MASS, &
                           ION_POT_H, &
                           TIME_NORM, &
                           DE_BROGLIE_L3

    ION_MASS = ATOMIC_A * 1.67262D-27!(checked 24.09.2019.)
    ION_POT_H = 13.60D00/ TEMP_0_EV
    TEMP_0_J = TEMP_0_EV * EL_CHARGE
    GAMMA_EI_0 = ION_Z ** 2 * GAMMA_EE_0
    THERM_VEL = SQRT(2.00D00 * TEMP_0_J/EL_MASS)
    TIME_NORM = THERM_VEL ** 3 /(GAMMA_EI_0 * DENSITY_0 *LAMBDA_EI(1.0D+00,1.0D+00,ION_Z)/ION_Z)
    DE_BROGLIE_L3 = SQRT(PLANCK_H ** 2/ (2* PI * EL_MASS * TEMP_0_J)) ** 3

  END SUBROUTINE INIT_DERIVED_L1
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE INIT_DERIVED_L2(COLL_EE_0, &
                       COLL_EI_0, &
                       COLL_EN_0, &
                       COLL_EN_L, &
                       TB_REC_A_0, &
                       RAD_REC_A_0, &
                       MAX_E_J_0, &
                       MAX_E_B_0, &
                       HEATING_A_0, &
                       DIFF_A_0 &
                       )

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: COLL_EE_0, &
                         COLL_EI_0, &
                         COLL_EN_0, &
                         COLL_EN_L, &
                         TB_REC_A_0, &
                         RAD_REC_A_0, &
                         MAX_E_J_0, &
                         MAX_E_B_0, &
                         HEATING_A_0, &
                         DIFF_A_0

    COLL_EE_0 = GAMMA_EE_0 * TIME_NORM * DENSITY_0/(THERM_VEL ** 3)
    COLL_EI_0 = GAMMA_EI_0 * TIME_NORM * DENSITY_0/(THERM_VEL ** 3)
    COLL_EN_0 = THERM_VEL * DENSITY_0 * SIGMA_0 * TIME_NORM * EL_MASS/(EL_MASS + ION_MASS)
    COLL_EN_L = THERM_VEL * DENSITY_0 * SIGMA_0 *TIME_NORM
    TB_REC_A_0 =  DENSITY_0 * DE_BROGLIE_L3
    RAD_REC_A_0 = RAD_REC_BETA_0 * DENSITY_0 * TIME_NORM
    MAX_E_J_0 = TIME_NORM ** 2 * EL_CHARGE ** 2 * DENSITY_0 /(EPSILON_0 * EL_MASS)
    MAX_E_B_0 = 1.00D00/(EPSILON_0 * MU_0 * THERM_VEL ** 2)
    HEATING_A_0 = 1.00D6/(DENSITY_0 * EL_MASS * THERM_VEL ** 3)
    DIFF_A_0 = SQRT(EL_MASS/ION_MASS)/(2*DENSITY_0*SIGMA_0*THERM_VEL*TIME_NORM)

  END SUBROUTINE INIT_DERIVED_L2
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates Coulomb logarithm e-e collisions for given normalized denisity n and temperature T
  REAL(KIND=DEFAULT_REAL) FUNCTION LAMBDA_EE(n,T)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n,T

    LAMBDA_EE = 23.50D00 - (0.50D00 * (LOG(n) + LOG(DENSITY_0) + LOG(1.00D-6)) + (-5.0/4.0) * LOG(T*TEMP_0_EV))&
              - SQRT(1.00D-5 + (LOG(T*TEMP_0_EV)-2.0)**2/16.0)


  END FUNCTION LAMBDA_EE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates Coulomb logarithm e-i collisions for given normalized denisity n and temperature T
  REAL(KIND=DEFAULT_REAL) FUNCTION LAMBDA_EI(n,T,Z)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n,T,Z

    IF (T*TEMP_0_EV .LT. 10.00D00 * Z ** 2) THEN

      LAMBDA_EI = 23.00D00 - LOG(SQRT(n * DENSITY_0 * 1.00D-6) * Z * (T * TEMP_0_EV) ** (-3.0/2.0))

    ELSE

      LAMBDA_EI = 24.00D00 - LOG(SQRT(n * DENSITY_0 * 1.00D-6)/(T*TEMP_0_EV))

    END IF

  END FUNCTION LAMBDA_EI
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE NORMALIZATION
!-------------------------------------------------------------------------------------------------------------------------------------
