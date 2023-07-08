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
!> Contains all collisional-radiative submatrix builders
MODULE BUILD_CRM

  USE GRID
  USE NORMALIZATION
  USE NEUT_AND_HEAT
  USE COLL_CS_RATES
  USE MATRIX_DATA
  USE INEL_GRID
  USE MPI
  USE VAR_KIND_DATA
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Builds collisional-radiative submatrix for exciataion and spontaneous de-excitation
  SUBROUTINE FILL_CRM_EX(P)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: P !< Current spatial cell
    REAL(KIND=HIGH_PREC_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    INTEGER :: I, J, ROW, COLUMN , K


!Fill diagonal elements according to adequate marker
    DO I = 1, NUM_NEUTRALS

      ROW = I
      COLUMN = I

      VAL = 0.0

      DO J = 1, I - 1

        VAL = VAL + SPONT_DEEX_A_0 * TIME_NORM * DIPOLE_TRANSITION(I,J)         !Loss due to spontaneous emission

      END DO

      DO J = I + 1, NUM_NEUTRALS

        VAL = VAL + EX_RATES(I,J,P)                                             !Collisional loss

      END DO

      DO J = 1, I - 1

        VAL = VAL + DEEX_RATES(I,J,P)                                           !Collisional loss

      END DO

      LOCAL_M%VALUE(MARKER_CRM_EX(P) + I) = LOCAL_M%VALUE(MARKER_CRM_EX(P) + I) - VAL  !Send element to global matrix
    END DO

    K = NUM_NEUTRALS + 1                                                        !Advance nonzero count

!Fill lower triangular elements
    DO I = 1, NUM_NEUTRALS

      DO J = 1, I - 1

        ROW = I
        COLUMN = J

        VAL = EX_RATES(J,I,P)                                                   !Gain due to collisions

        LOCAL_M%VALUE(MARKER_CRM_EX(P) + K) = LOCAL_M%VALUE(MARKER_CRM_EX(P) + K) + VAL !Send element to global matrix

        K = K + 1                                                               !Advance nonzero count

      END DO

    END DO

!Fill upper triangular elements
    DO I = 1, NUM_NEUTRALS

      DO J = I + 1, NUM_NEUTRALS

        ROW = I
        COLUMN = J

        VAL =  SPONT_DEEX_A_0 * TIME_NORM * DIPOLE_TRANSITION(J, I)              !Gain due to spontaneous emission

        VAL = VAL + DEEX_RATES(J,I,P)                                           !Gain due to collisions

        LOCAL_M%VALUE(MARKER_CRM_EX(P) + K) = LOCAL_M%VALUE(MARKER_CRM_EX(P) + K) + VAL !Send element to global matrix

        K = K + 1                                                               !Advance nonzero count

      END DO

    END DO

  END SUBROUTINE FILL_CRM_EX
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Builds collisional-radiative submatrix for ionization
  SUBROUTINE FILL_CRM_ION(P)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: P !< Current spatial cell

    INTEGER :: I
    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix


!Fill diagonal elements
    DO I = 1, NUM_NEUTRALS

      VAL = - ION_RATES(I,P)                                                    !Ionization loss

      LOCAL_M%VALUE(MARKER_CRM_ION(P) + I) = LOCAL_M%VALUE(MARKER_CRM_ION(P) + I) + VAL !Send element to global matrix

    END DO

  END SUBROUTINE FILL_CRM_ION
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Builds recombination CRM submatrix, including the loss of cold electrons
  SUBROUTINE FILL_CRM_RECOMB(T_e, n_i, P)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::  T_e, & !< Lagged electron temperature
                                 n_i !< Lagged ion density

    INTEGER, INTENT(IN) :: P !< Current spatial cell

    REAL(KIND=HIGH_PREC_REAL) :: VAL,& !< Current value of matrix element to be passed to global matrix
                    RAD_RATE(NUM_NEUTRALS), TB_RATE(NUM_NEUTRALS),&
                    R_TOT, TB_TOT

    INTEGER :: I,J, ROW, COLUMN, K

    R_TOT = 0
    TB_TOT = 0


!Initialize recombination rates
    DO I = 1, NUM_NEUTRALS

      RAD_RATE(I) = RECOMB_RATE(I, T_e)
      TB_RATE(I) = TB_RECOMB_RATES(I,P)

      R_TOT = R_TOT + RAD_RATE(I)
      TB_TOT = TB_TOT + TB_RATE(I)


    END DO

    VAL = (TB_REC_A_0 * TB_TOT  &
              + RAD_REC_A_0 * R_TOT) * n_i

    DO I = 1, NUM_V                                                             !Calculate cold electron loss rate due to recombination

      ROW = 1
      COLUMN = I

      LOCAL_M%VALUE(MARKER_RECOMB_CRM(P,I)) = &
                              LOCAL_M%VALUE(MARKER_RECOMB_CRM(P,I))  &
                              - VAL * V_GRID(I) ** 2 * V_GRID_WIDTH(I)/ (V_GRID(1) ** 2 * V_GRID_WIDTH(1)) !Send element to global matrix

    END DO

    K = NUM_V + 1                                                               !Advance nonzero count


    DO I = 1, NUM_NEUTRALS                                                      !Calculate CR model effect of recombination

      DO J = 1, NUM_V

        ROW = I + NUM_V
        COLUMN = J

        VAL = (TB_REC_A_0 * TB_RATE(I) * n_i &
              + RAD_REC_A_0 * RAD_RATE(I) * n_i) * &
              4.00D00 * PI * V_GRID(J) ** 2 * V_GRID_WIDTH(J)

        LOCAL_M%VALUE(MARKER_RECOMB_CRM(P,K)) = LOCAL_M%VALUE(MARKER_RECOMB_CRM(P,K)) + VAL

        K = K + 1                                                               !Advance nonzero count

      END DO

    END DO

  END SUBROUTINE FILL_CRM_RECOMB
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Builds recombination CRM matrix for full fluid mode
  SUBROUTINE FILL_CRM_RECOMB_FF(T_e, n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::  T_e(NUM_X), & !< Lagged electron temperature
                                 n_i(NUM_X) !< Lagged ion density

    REAL(KIND=HIGH_PREC_REAL) :: VAL,& !< Current value of matrix element to be passed to global matrix
                    RAD_RATE(MIN_X:MAX_X,NUM_NEUTRALS), TB_RATE(MIN_X:MAX_X,NUM_NEUTRALS)


    INTEGER :: I,J, K, P

!Initialize recombination rates

    DO P = MIN_X,MAX_X

      DO I = 1, NUM_NEUTRALS

        RAD_RATE(P,I) = RECOMB_RATE(I, T_e(P))
        TB_RATE(P,I) = TB_RECOMB_RATES(I,P)

      END DO

    END DO

    K = 1

    DO P = MIN_X,MAX_X

      IF (MOD(P,2) .EQ. 1) THEN

        DO J = 1, NUM_NEUTRALS

          VAL = TB_REC_A_0 * TB_RATE(P,J) * n_i(P) &
                    + RAD_REC_A_0 * RAD_RATE(P,J) * n_i(P)


          LOCAL_M%VALUE(MARKER_RECOMB_CRM_FF(K)) = &
                                    LOCAL_M%VALUE(MARKER_RECOMB_CRM_FF(K))  + VAL  !Send element to global matrix

          K = K + 1

        END DO

      END IF

    END DO

  END SUBROUTINE FILL_CRM_RECOMB_FF
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_CRM
!-------------------------------------------------------------------------------------------------------------------------------------
