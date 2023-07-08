!-------------------------------------------------------------------------------------------------------------------------------------
!> Contains post-processing test routines
MODULE TESTS

  USE GRID
  USE SWITCHES
  USE NORMALIZATION
  USE NEUT_AND_HEAT
  USE MPI
  USE COLL_CS_RATES
  USE VAR_KIND_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
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
!> Test heat flow q against Spitzer-Harm heat flow for given temperature profile and returns ratio
  FUNCTION SH_TEST(q,T,n)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SH_TEST(NUM_X) !< Returned ratio

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: q, & !< Input heat flow
                                            T, &!< Input temperature
                                            n !< Input density

    REAL(KIND=DEFAULT_REAL) :: KAPPA, log_ratio(NUM_X)
    INTEGER :: I

    SH_TEST = 0

    DO I = 1, NUM_X

      log_ratio(I) = LAMBDA_EI(1.00D00,1.00D00,ION_Z)/LAMBDA_EI(n(I),T(I),Z_PROF(I))

    END DO

    IF (.NOT. COLL_EE_L_SWITCH) THEN

      KAPPA = 13.581 * 3.00D00 * SQRT(PI)/16.00D00                                          !Lorentz case

    ELSE

      KAPPA = 3.20 * 3.00D00 * SQRT(PI)/16.00D00                                            !SH case with Epperlein coeff

    END IF

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      SH_TEST(1) = q(1)/(- KAPPA*ION_Z/Z_PROF(1) *log_ratio(1)* T(1) ** (5.00D00/2.00D00) * (T(2) - T(NUM_X))/dxc(1))
      SH_TEST(NUM_X) = q(NUM_X)/(- KAPPA*ION_Z/Z_PROF(NUM_X) *log_ratio(NUM_X)* T(NUM_X) &
      ** (5.00D00/2.00D00) * (T(1) - T(NUM_X - 1))/dxc(NUM_X))

    END IF

    DO I = 2, NUM_X - 1

      SH_TEST(I) = q(I)/(- KAPPA*ION_Z/Z_PROF(I) * log_ratio(I) * T(I) ** (5.00D00/2.00D00)*(T(I + 1) - T(I - 1))/dxc(I))

    END DO

  END FUNCTION SH_TEST
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates total energy contained in atomic structure (for hydrogen) and in electron thermal motion
  FUNCTION ATOMIC_EN_TEST(T, n, n_n)

    IMPLICIT NONE

    REAL(KIND=HIGH_PREC_REAL) :: ATOMIC_EN_TEST(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: T, &!< Input temperature
                                                  n !< Input electron density

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n(NUM_NEUTRALS,NUM_X) !< Input neutral state densities

    REAL(KIND=HIGH_PREC_REAL) :: E(NUM_NEUTRALS), TEST(NUM_X)

    INTEGER :: I,J

    !Calculate energies

    DO I = 1, NUM_NEUTRALS

      E(I) = ION_POT_H * (1.00D00 - 1.00D00 / REAL(I ** 2,KIND=DEFAULT_REAL))

    END DO

    TEST = 0

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        TEST(I) = TEST(I) + n_n(J,I) * E(J)

      END DO

      TEST(I) = TEST(I) + n(I) * ION_POT_H
      TEST(I) = TEST(I) + 3.00D00 * n(I)*T(I) /2.00D00

    END DO

    ATOMIC_EN_TEST = TEST

  END FUNCTION ATOMIC_EN_TEST
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates relative difference between two spatial vectors
  FUNCTION DIFF_TEST(a,b)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: DIFF_TEST(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: a, &
                                                  b

    REAL(KIND=DEFAULT_REAL) :: TEST(NUM_X)

    INTEGER :: I

    DO I = 1, NUM_X

      TEST(I) = (a(I) - b(I))

    END DO

    DIFF_TEST = TEST

  END FUNCTION DIFF_TEST
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates current total ionization rate
  FUNCTION S_ION_SK(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: S_ION_SK(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    REAL(KIND=DEFAULT_REAL) :: I_RATES(NUM_X,NUM_NEUTRALS)

    REAL(KIND=DEFAULT_REAL) :: F0(NUM_V)

    INTEGER :: I,J,K

    S_ION_SK = 0

    DO K = 1, NUM_X

      F0 = f(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        I_RATES(K,I) = COLL_RATE(F0, 10000 + 100 * I)

      END DO

    END DO

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        S_ION_SK(I) = S_ION_SK(I) + I_RATES(I,J) * f(X_POS(I) + NUM_V*NUM_H +  J)

      END DO

    END DO

  END FUNCTION S_ION_SK
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates Maxwellian ionization rate
  FUNCTION S_ION_M(T,n,f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: S_ION_M(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: T, &
                                                  n

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    REAL(KIND=DEFAULT_REAL) :: I_RATES(NUM_X,NUM_NEUTRALS)

    REAL(KIND=DEFAULT_REAL) :: F0(NUM_V)

    INTEGER :: I,J,K

    S_ION_M = 0

    DO K = 1, NUM_X

      DO J = 1, NUM_V

        F0(J) = n(K) * ((PI * T(K)) ** (- 3.00D00/2.00D00)) * EXP( - (V_GRID(J) ** 2)/ T(K))

      END DO

      DO I = 1,NUM_NEUTRALS

        I_RATES(K,I) = COLL_RATE(F0, 10000 + 100 * I)

      END DO

    END DO

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        S_ION_M(I) = S_ION_M(I) + I_RATES(I,J) * f(X_POS(I) + NUM_V*NUM_H +  J)

      END DO

    END DO

  END FUNCTION S_ION_M
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates total energy loss rate due to ionization
  FUNCTION ION_E_RATE(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: ION_E_RATE(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    REAL(KIND=DEFAULT_REAL) :: I_RATES(NUM_X,NUM_NEUTRALS)

    REAL(KIND=DEFAULT_REAL) :: F0(NUM_V)

    INTEGER :: I,J,K

    ION_E_RATE = 0

    DO K = 1, NUM_X

      F0 = f(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        I_RATES(K,I) = COLL_RATE(F0, 10000 + 100 * I)

      END DO

    END DO

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        ION_E_RATE(I) = ION_E_RATE(I) + I_RATES(I,J) * f(X_POS(I) + NUM_V*NUM_H +  J) *&
         ION_POT_H / REAL(J ** 2,KIND=HIGH_PREC_REAL)

      END DO

    END DO

  END FUNCTION ION_E_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates total energy loss rate due to excitation
  FUNCTION EX_E_RATE(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: EX_E_RATE(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    REAL(KIND=DEFAULT_REAL) :: EXC_RATES(NUM_X,NUM_NEUTRALS,NUM_NEUTRALS)

    REAL(KIND=DEFAULT_REAL) :: F0(NUM_V)

    INTEGER :: I,J,K

    EX_E_RATE = 0

    DO K = 1, NUM_X

      F0 = f(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        DO J = I + 1, NUM_NEUTRALS

          EXC_RATES(K,I,J) = COLL_RATE(F0, 10000 + I * 100 + J)

        END DO

      END DO

    END DO

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        DO K = J+1, NUM_NEUTRALS

          EX_E_RATE(I) = EX_E_RATE(I) + EXC_RATES(I,J,K) * f(X_POS(I) + NUM_V*NUM_H +  J) *&
          ION_POT_H * (1.00D00 / REAL(J ** 2,KIND=HIGH_PREC_REAL) &
          - 1.00D00 / REAL(K ** 2,KIND=HIGH_PREC_REAL))

        END DO

      END DO

    END DO

  END FUNCTION EX_E_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates total energy gain rate due to 3b recombination
  FUNCTION REC_3B_E_RATE(f,n_e,n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: REC_3B_E_RATE(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: n_e,n_i

    REAL(KIND=DEFAULT_REAL) :: REC_3B_RATES(NUM_X,NUM_NEUTRALS)

    REAL(KIND=DEFAULT_REAL) :: F0(NUM_V)

    INTEGER :: I,J,K

    REC_3B_E_RATE = 0

    DO K = 1, NUM_X

      F0 = f(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        REC_3B_RATES(K,I) = COLL_RATE(F0, 10000 + I,K)

      END DO

    END DO

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        REC_3B_E_RATE(I) = REC_3B_E_RATE(I) + TB_REC_A_0*REC_3B_RATES(I,J) * n_e(I)*n_i(I) *&
         (ION_POT_H / REAL(J ** 2,KIND=HIGH_PREC_REAL))

      END DO

    END DO

  END FUNCTION REC_3B_E_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates total energy gain rate due to deexcitation
  FUNCTION DEEX_E_RATE(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: DEEX_E_RATE(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    REAL(KIND=DEFAULT_REAL) :: DEEXC_RATES(NUM_X,NUM_NEUTRALS,NUM_NEUTRALS)

    REAL(KIND=DEFAULT_REAL) :: F0(NUM_V)

    INTEGER :: I,J,K

    DEEX_E_RATE = 0

    DO K = 1, NUM_X

      F0 = f(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        DO J = 1, I-1

          DEEXC_RATES(K,I,J) = COLL_RATE(F0, 10000 + I * 100 + J, K)

        END DO

      END DO

    END DO

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        DO K = 1, J-1

          DEEX_E_RATE(I) = DEEX_E_RATE(I) + DEEXC_RATES(I,J,K) * f(X_POS(I) + NUM_V*NUM_H +  J) *&
          ION_POT_H * (1.00D00 / REAL(K ** 2,KIND=HIGH_PREC_REAL) &
          - 1.00D00 / REAL(J ** 2,KIND=HIGH_PREC_REAL))

        END DO

      END DO

    END DO

  END FUNCTION DEEX_E_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates total energy radiation rate in deexcitation
  FUNCTION RAD_DEEX_E_RATE(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: RAD_DEEX_E_RATE(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    INTEGER :: I,J,K

    RAD_DEEX_E_RATE = 0

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        DO K = 1, J-1

          RAD_DEEX_E_RATE(I) = RAD_DEEX_E_RATE(I) + SPONT_DEEX_A_0 * TIME_NORM * DIPOLE_TRANSITION(J,K)  &
          * f(X_POS(I) + NUM_V*NUM_H +  J) * ION_POT_H * (1.00D00 / REAL(K ** 2,KIND=HIGH_PREC_REAL) &
          - 1.00D00 / REAL(J ** 2,KIND=HIGH_PREC_REAL))

        END DO

      END DO

    END DO

  END FUNCTION RAD_DEEX_E_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates total energy radiation rate in recombination
  FUNCTION RAD_REC_E_RATE(T_e,n_e,n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: RAD_REC_E_RATE(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: n_e,n_i,T_e

    INTEGER :: I,J

    RAD_REC_E_RATE = 0

    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        RAD_REC_E_RATE(I) = RAD_REC_E_RATE(I) + RAD_REC_A_0* RECOMB_RATE(J, T_e(I))* n_i(I) * n_e(I)*&
        (ION_POT_H / REAL(J ** 2,KIND=HIGH_PREC_REAL))

      END DO

    END DO

  END FUNCTION RAD_REC_E_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates numerical heating from jE term
  FUNCTION NUM_DV_HEATING(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: NUM_DV_HEATING(NUM_X)

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f

    REAL(KIND=DEFAULT_REAL) :: E

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_V) :: F1

    INTEGER :: J,K

    NUM_DV_HEATING = 0

    DO K = 1, NUM_X

      F1 = f(X_POS(K) + (NUM_H - 2) * NUM_V + 1: X_POS(K) + (NUM_H - 1) * NUM_V)
      E = f(X_POS(K) + NUM_H*NUM_V + NUM_NEUTRALS + 1)

      NUM_DV_HEATING(K) = NUM_DV_HEATING(K) + (4.0D0*PI/3.0D0)*E*V_GRID(1) ** 2 * ((V_CELL_BOUNDARY(1)) ** 2 * &
            ((1.00D00-v_interp(1))*F1(2) + v_interp(1)*F1(1))) /2.0D00

      DO J = 2, NUM_V-1

        NUM_DV_HEATING(K) = NUM_DV_HEATING(K) + (4.0D0*PI/3.0D0)*V_GRID(J)**2*E*((V_CELL_BOUNDARY(J)) ** (2) * &
              ((1.00D00-v_interp(J))*F1(J + 1) + v_interp(J)*F1(J)) - (V_CELL_BOUNDARY(J - 1)) ** (2) * &
              ((1.00D00 - v_interp(J-1))*F1(J) + v_interp(J-1)*F1(J - 1))) / 2.0D00

      END DO

      NUM_DV_HEATING(K) = NUM_DV_HEATING(K) + V_GRID(NUM_V) ** (2) * (- (V_CELL_BOUNDARY(NUM_V-1)) ** (2) * &
            ((1.00D00 - v_interp(NUM_V-1))*F1(NUM_V) + v_interp(NUM_V-1)*F1(NUM_V - 1))) / 2.00D00

    END DO

  END FUNCTION NUM_DV_HEATING
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE TESTS
!-------------------------------------------------------------------------------------------------------------------------------------
