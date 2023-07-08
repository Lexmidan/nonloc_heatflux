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
!> Contains heating term submatrix builder
MODULE BUILD_HEATING

  USE GRID
  USE NEUT_AND_HEAT
  USE NORMALIZATION
  USE SWITCHES
  USE MATRIX_DATA
  USE MPI
  USE VAR_KIND_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills diffusive  heating submatrix
  SUBROUTINE FILL_HEATING(t, n_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X) !< Lagged density vector
    INTEGER, INTENT(IN) :: t !< Timestep number

    INTEGER :: I, K, OFFSET_UP_H, MIN_C, MAX_C, RANK, IERR

    REAL(KIND=DEFAULT_REAL) :: n_h, &
              VAL

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    IF (t .LE. T_HEATING) THEN

      !Set fixed boundary condition offsets
      OFFSET_UP_H = 0

      IF (FIXED_BOUNDARY_UP_SWITCH) OFFSET_UP_H = 2

      IF (RANK .EQ. 0) THEN

        MIN_C = 1 + OFFSET_UP_H

      ELSE

        MIN_C = MIN_X + MOD(MIN_X + 1, 2)

      END IF

      MAX_C = MAX_X - MOD(MAX_X+1,2)

      IF ((MAX_C + 1)/2 .GT. N_HEATING + OFFSET_UP_H/2)  MAX_C = 2 * N_HEATING - 1 + OFFSET_UP_H

      !Calculate heating length
      n_h = 0
      DO I = 1 + OFFSET_UP_H, 2 * N_HEATING - 1 + OFFSET_UP_H, 2

          n_h = n_h + dxc(I)

      END DO

      !Fill matrix
      DO K = MIN_C , MAX_C , 2

        DO I = 1, TRI_DIAG_SP%N_NZ

          VAL = HEATING_A_0 * HEAT_POWER * &
                ((TRI_DIAG_SP%COL(I) / TRI_DIAG_SP%ROW(I)) * (TRI_DIAG_SP%ROW(I) / TRI_DIAG_SP%COL(I)) * &          !Diagonal elements
                ((-1.0D0 + (1 / TRI_DIAG_SP%ROW(I))) &
                * V_CELL_BOUNDARY(TRI_DIAG_SP%ROW(I) - 1) ** 2/dvm(TRI_DIAG_SP%ROW(I)) + &
                (-1.0D0 + (TRI_DIAG_SP%ROW(I)) / NUM_V) &
                * V_CELL_BOUNDARY(TRI_DIAG_SP%ROW(I)) ** 2/dvp(TRI_DIAG_SP%ROW(I)))  &
                + (TRI_DIAG_SP%COL(I) / (TRI_DIAG_SP%ROW(I) + 1) &
                * V_CELL_BOUNDARY(TRI_DIAG_SP%ROW(I)) ** 2 / dvp(TRI_DIAG_SP%ROW(I))) &                                                 !Upper off-diagonal elements
                + (TRI_DIAG_SP%ROW(I) / (TRI_DIAG_SP%COL(I) + 1) &
                * V_CELL_BOUNDARY(TRI_DIAG_SP%ROW(I)-1) ** 2/dvm(TRI_DIAG_SP%ROW(I)))  &                                                 !Lower off-diagonal elements
                ) / (3.0D0 * n_e(K) * n_h * V_GRID(TRI_DIAG_SP%ROW(I)) ** 2 * V_GRID_WIDTH(TRI_DIAG_SP%ROW(I)))

          LOCAL_M%VALUE(MARKER_HEATING(K) + I) = LOCAL_M%VALUE(MARKER_HEATING(K) + I) + VAL

        END DO

      END DO

    END IF


  END SUBROUTINE FILL_HEATING
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron heating term if full fluid mode
  SUBROUTINE FILL_EL_HEATING(n_e,T_e,TIMESTEP)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T_e(NUM_X), & !< Lagged electron temperature
                                           n_e(NUM_X)

    INTEGER, INTENT(IN) :: TIMESTEP

    INTEGER :: J, P

    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    IF (TIMESTEP .LE. T_HEATING) THEN

      J = 1

      DO P = MIN_X, MIN(MAX_X,2 * N_HEATING - 1)

        IF (MOD(P,2) .EQ. 1) THEN

            VAL = 4.00D00 *HEATING_A_0 * HEAT_POWER &
            / (3.00D00* n_e(P) * T_e(P) * (X_GRID(2*N_HEATING-1)+dxc(2*N_HEATING-1)/2.00D00+dxc(1)/2.00D00))

            LOCAL_M%VALUE(MARKER_EL_HEATING(J)) = LOCAL_M%VALUE(MARKER_EL_HEATING(J)) + VAL

            J = J + 1

        END IF

    END DO

    END IF

  END SUBROUTINE FILL_EL_HEATING
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_HEATING
!-------------------------------------------------------------------------------------------------------------------------------------
