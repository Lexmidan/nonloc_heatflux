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
!> Contains build routine for neutral diffusion
MODULE BUILD_NEUTRAL_DIFF

  USE GRID
  USE NEUT_AND_HEAT
  USE MATRIX_DATA
  USE SWITCHES
  USE NORMALIZATION
  USE MPI
  USE VAR_KIND_DATA

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: SIGMA_CX_DIFF(:)

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills neutral diffusion submatrix
  SUBROUTINE FILL_DIFF_N(f,n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(DIM_F) !< Lagged vector
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i(NUM_X) !< Lagged ion density

    REAL(KIND=DEFAULT_REAL) :: n_n(NUM_X,NUM_NEUTRALS), &
              VAL, D_plus(MIN_X:MAX_X,NUM_NEUTRALS), D_minus(MIN_X:MAX_X,NUM_NEUTRALS)

    INTEGER :: I,POS,J,P



    IF (.NOT. ALLOCATED(SIGMA_CX_DIFF)) THEN

      ALLOCATE(SIGMA_CX_DIFF(NUM_NEUTRALS))

      SIGMA_CX_DIFF = 0

      IF (SIMPLE_CX_SWITCH) THEN

        DO P = MIN_X, MAX_X

          IF (MOD(P,2) .EQ. 0) THEN

            SIGMA_CX_DIFF(1) =  3.00D-19/SIGMA_0

            IF (NUM_NEUTRALS .GT. 1) &
            SIGMA_CX_DIFF(2) = 2**4 * 1.00D-19/SIGMA_0
            IF (NUM_NEUTRALS .GT. 2) &
            SIGMA_CX_DIFF(3) = 3**4 * 7.00D-20/SIGMA_0

            DO I = 4, NUM_NEUTRALS

              SIGMA_CX_DIFF(4) = I ** 4 * 6.00D-20/SIGMA_0

            END DO

          END IF

        END DO

      END IF

    END IF

    !Get lagged neutral density
    DO I = 1, NUM_X

      DO J = 1, NUM_NEUTRALS

        n_n(I,J) = f(X_POS(I) + NUM_V * NUM_H + J)

      END DO

    END DO

    DO J = 1, NUM_NEUTRALS
      !Calculate diffusion coefficients
      DO I = MAX(2,MIN_X), MIN(MAX_X,NUM_X - 1)

        D_plus(I,J) = SQRT(NEUTRAL_TEMP) &
        / (SIGMA_EL * J**4* (n_n(I+1,1)+ n_i(I+1)) + SIGMA_CX_DIFF(J) * n_i(I+1))
        D_minus(I,J) = SQRT(NEUTRAL_TEMP) &
        / (SIGMA_EL * J**4* (n_n(I-1,1)+ n_i(I-1)) + SIGMA_CX_DIFF(J) * n_i(I-1))

      END DO

  !Handle first and last cell
      IF (MIN_X .EQ. 1) THEN

        D_minus(1,J) = 0
        D_plus(1,J) = SQRT(NEUTRAL_TEMP) / (SIGMA_EL * J**4 * (n_n(2,1)+ n_i(2))+ SIGMA_CX_DIFF(J) * n_i(2))

      END IF

      IF (MAX_X .EQ. NUM_X) THEN

        D_plus(NUM_X,J) = 0
        D_minus(NUM_X,J) = SQRT(NEUTRAL_TEMP)&
         / (SIGMA_EL * J**4 * (n_n(NUM_X-1,1)+ n_i(NUM_X-1))+ SIGMA_CX_DIFF(J) * n_i(NUM_X-1))

      END IF

    END DO

    DO I = 1, NEUTRAL_DIFF_SP%N_NZ

      POS  = (NEUTRAL_DIFF_SP%ROW(I) / NUM_0D) + 1 - NEUTRAL_DIFF_N(I) / NUM_0D

      VAL = DIFF_A_0 * &         !Extract constants and position
            ((NEUTRAL_DIFF_SP%COL(I) / NEUTRAL_DIFF_SP%ROW(I)) * (NEUTRAL_DIFF_SP%ROW(I) / NEUTRAL_DIFF_SP%COL(I)) * &          !Diagonal elements
            (-D_plus(POS,NEUTRAL_DIFF_N(I))/dxp(POS) -D_minus(POS,NEUTRAL_DIFF_N(I))/dxm(POS) + &
            (1 / POS) * D_minus(POS,NEUTRAL_DIFF_N(I))/dxm(POS) + &
            (POS / NUM_X) * D_plus(POS,NEUTRAL_DIFF_N(I))/dxp(POS)) &
            + (NEUTRAL_DIFF_SP%COL(I) / (NEUTRAL_DIFF_SP%ROW(I) + 2 * NUM_0D)) * D_plus(POS,NEUTRAL_DIFF_N(I))/ dxp(POS) &            !Upper off-diagonal elements
            + (NEUTRAL_DIFF_SP%ROW(I) / (NEUTRAL_DIFF_SP%COL(I) + 2 * NUM_0D)) * D_minus(POS,NEUTRAL_DIFF_N(I))/ dxm(POS)  &           !Lower off-diagonal elements
            ) / dxc(POS)

      LOCAL_M%VALUE(MARKER_NEUT_DIFF(I)) = LOCAL_M%VALUE(MARKER_NEUT_DIFF(I)) + VAL

    END DO

  END SUBROUTINE FILL_DIFF_N
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_NEUTRAL_DIFF
