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
!> Contains velocity space advection submatrix routines due to E_x field
MODULE BUILD_E_ADV

  USE SWITCHES
  USE GRID
  USE MATRIX_DATA
  USE MPI
  USE VAR_KIND_DATA

  PRIVATE :: G_F, H_F

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills velocity space advection submatrix due to E_x field
  SUBROUTINE FILL_E_ADV(f_lagged)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(DIM_F) !< Lagged variable vector

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

!Fill matrix according to adequate sparsity pattern, marker, and properties list
    DO I = 1, E_ADV_SP%N_NZ

      VAL = PROP_E_ADV%M(I) * ((PROP_E_ADV%H(I) - 1) / (2.00D00 * (PROP_E_ADV%H(I) - 1) - 1.00D00) &
            * G_F(PROP_E_ADV%V(I), f_lagged, PROP_E_ADV%H(I) - 2, PROP_E_ADV%P(I)) &
            + (PROP_E_ADV%H(I) - 1 + 1.00D00) / (2.00D00 * (PROP_E_ADV%H(I) - 1) + 3.00D00) &
            * H_F(PROP_E_ADV%V(I), f_lagged, PROP_E_ADV%H(I), PROP_E_ADV%P(I)))

      LOCAL_M%VALUE(MARKER_E_ADV + I) = LOCAL_M%VALUE(MARKER_E_ADV + I) + VAL

    END DO

  END SUBROUTINE FILL_E_ADV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates the G_l^m derivative function for given lagged vector
  REAL(KIND=DEFAULT_REAL) FUNCTION G_F(J,f_l,L,P)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: J, & !< Velocity cell number at which function is evaluated
                           L, & !< Harmonic L-number of funcion
                           P !< Spatial cell number at which function is evaluated
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_l(DIM_F) !< Lagged vector

    REAL(KIND=DEFAULT_REAL) :: F(NUM_V)

    G_F = 0

    IF ((L .LT. L_MAX + 1) .AND. (L .GT. -1)) THEN

      F = f_l(X_POS(P) + NUM_V * H_POS(L) + 1: X_POS(P) + NUM_V * (H_POS(L) + 1))

      IF (J .EQ. 1) THEN

        IF (L .EQ. 0) G_F =(1.00D00 + V_GRID(1) ** 2 / V_GRID(2) ** 2) * &
        ((1.00D00 - v_interp(1))*v_interp(1)*F(2) - v_interp(1)*F(1)) / dv

      ELSE IF (L .GT. - 1) THEN

        IF (J .EQ. NUM_V) THEN

          G_F =  V_GRID(J) ** L * ((V_CELL_BOUNDARY(J)) ** (-L) * v_interp(J)* F(J) &
          - (V_CELL_BOUNDARY(J - 1)) ** (-L) * &
          ((1.00D00 - v_interp(J-1))*F(J) + v_interp(J-1)*F(J - 1))) / V_GRID_WIDTH(J)

        ELSE

          G_F =  V_GRID(J) ** L * ((V_CELL_BOUNDARY(J)) ** (-L) * &
          ((1.00D00-v_interp(J))*F(J + 1) + v_interp(J)*F(J)) &
          - (V_CELL_BOUNDARY(J - 1)) ** (-L) * &
          ((1.00D00 - v_interp(J-1))*F(J) + v_interp(J-1)*F(J - 1))) / V_GRID_WIDTH(J)

        END IF

      END IF

    END IF

  END FUNCTION G_F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates the H_l^m derivative function for given lagged vector
  REAL(KIND=DEFAULT_REAL) FUNCTION H_F(J,f_l,L,P)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: J, & !< Velocity cell number at which function is evaluated
                           L, & !< Harmonic L-number of funcion
                           P !< Spatial cell number at which function is evaluated
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_l(DIM_F) !< Lagged vector

    REAL(KIND=DEFAULT_REAL) :: F(NUM_V)

    H_F = 0

    IF (L .LT. L_MAX + 1) THEN

      F = f_l(X_POS(P) + NUM_V * H_POS(L) + 1: X_POS(P) + NUM_V * (H_POS(L) + 1))

      IF (J .EQ. NUM_V) THEN

        H_F = V_GRID(J) ** (-L - 1) * (- (V_CELL_BOUNDARY(J-1)) ** (L + 1) * &
              ((1.00D00 - v_interp(J-1))*F(J) + v_interp(J-1)*F(J - 1))) / V_GRID_WIDTH(J)

      ELSE IF (J .EQ. 1) THEN

        H_F = V_GRID(J) ** (-L - 1) * ((V_CELL_BOUNDARY(J)) ** (L + 1) * &
              ((1.00D00-v_interp(J))*F(J + 1) + v_interp(J)*F(J))) / V_GRID_WIDTH(J)

      ELSE

        H_F =  V_GRID(J) ** (-L - 1) * ((V_CELL_BOUNDARY(J)) ** (L + 1) * &
              ((1.00D00-v_interp(J))*F(J + 1) + v_interp(J)*F(J)) - (V_CELL_BOUNDARY(J - 1)) ** (L + 1) * &
              ((1.00D00 - v_interp(J-1))*F(J) + v_interp(J-1)*F(J - 1))) / V_GRID_WIDTH(J)

      END IF

    END IF

  END FUNCTION H_F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_E_ADV
!-------------------------------------------------------------------------------------------------------------------------------------
