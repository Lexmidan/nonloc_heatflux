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
!> Contains the spatial advection submatrix and the routine used to pre-build it
MODULE BUILD_X_ADV

  USE GRID
  USE SWITCHES
  USE SPARSE
  USE MPI
  USE VAR_KIND_DATA

  IMPLICIT NONE

  TYPE(SPARSE_MAT) :: M_X_ADV      !< Spatial x-advection submatrix

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills the spatial advection submatrix
  SUBROUTINE FILL_X_ADV

    IMPLICIT NONE

    INTEGER :: I, J, K, P, N_NZ, L, c, MIN_X_XADV, MAX_X_XADV

    REAL(KIND=DEFAULT_REAL) :: VAL

    INTEGER :: RANK, SIZE, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0
!Calculate number of nonzeroes
    DO I = MIN_X, MAX_X

      c = 2

      IF (((I .EQ. 1) .AND. NO_FLOW_BOUNDARY_UP_SWITCH) .OR. ((I .EQ. NUM_X) .AND. NO_FLOW_BOUNDARY_DIV_SWITCH)) c = 1

      DO L = 0, L_MAX

        IF (((MOD(L,2) .EQ. 0) .AND. (MOD(I,2) .EQ. 1)) .OR. ((MOD(L,2) .EQ. 1) .AND. (MOD(I,2) .EQ. 0))) THEN

          IF (L - 1 .GE. 0) N_NZ = N_NZ + C * NUM_V

          IF (L + 1 .LE. L_MAX) N_NZ = N_NZ + C * NUM_V

        END IF

      END DO

    END DO

    CALL ALLOCATE_SPARSE(M_X_ADV, LOC_ROWS, DIM_F, N_NZ)

    MIN_X_XADV = MIN_X
    MAX_X_XADV = MAX_X

    IF (.NOT. FIXED_BOUNDARY_DIV_SWITCH .AND. (RANK .EQ. SIZE - 1)) MAX_X_XADV = MAX_X - 1

    K = 1
!Fill auxillary sparse matrix

    IF (RANK .EQ. 0) THEN

      IF (.NOT. FIXED_BOUNDARY_UP_SWITCH) THEN                                    !Fill first cell advection if not fixed boundary

        MIN_X_XADV = MIN_X + 1
        DO I = 1, NUM_H

          IF (MOD(I - 1,2) .EQ. 0) THEN                                        !Evolve only even harmonics

            IF (0 .LE. I - 2) THEN                                  !First set of potential nonzeroes

              DO J = 1, NUM_V

                M_X_ADV%ROW(K) = NUM_V * (NUM_H - I) + J
                M_X_ADV%COLUMN(K) = X_POS(2) + NUM_V * H_POS(I - 2) + J

                VAL = - (I - 1) / (2.00D00 * (I - 1) - 1.00D00) * V_GRID(J) / dxc(1)
                M_X_ADV%VALUE(K) = VAL

                K = K + 1                                                         !Advance nonzero count

                IF (PERIODIC_BOUNDARY_SWITCH) THEN                                !If periodic boundary include advection from last cell

                  M_X_ADV%ROW(K) = NUM_V * (NUM_H - I) + J
                  M_X_ADV%COLUMN(K) = X_POS(NUM_X) + NUM_V * H_POS(I - 2) + J

                  VAL = (I - 1) / (2.00D00 * (I - 1) - 1.00D00) * V_GRID(J) / dxc(1)

                  M_X_ADV%VALUE(K) = VAL

                  K = K + 1                                                       !Advance nonzero count

                END IF

              END DO

            END IF

            IF (I .LE. L_MAX) THEN                                     !Second set of potential nonzeroes

              DO J = 1,NUM_V

                M_X_ADV%ROW(K) = NUM_V * (NUM_H - I) + J
                M_X_ADV%COLUMN(K) = X_POS(2) + NUM_V * H_POS(I) + J

                VAL = - I / (2.00D00 * (I - 1) + 3.00D00) * V_GRID(J) / dxc(1)

                M_X_ADV%VALUE(K) = VAL

                K = K + 1                                                         !Advance nonzero count

                IF (PERIODIC_BOUNDARY_SWITCH) THEN                                !If periodic boundary include advection from last cell

                  M_X_ADV%ROW(K) = NUM_V * (NUM_H - I) + J
                  M_X_ADV%COLUMN(K) = X_POS(NUM_X) + NUM_V * H_POS(I) + J

                  VAL = I / (2.00D00 * (I - 1) + 3.00D00) * V_GRID(J) / dxc(1)

                  M_X_ADV%VALUE(K) = VAL

                  K = K + 1                                                       !Advance nonzero count

                END IF

              END DO

            END IF

          END IF

        END DO

      END IF

    END IF

!Fill bulk spatial cells
    DO P = MIN_X_XADV, MAX_X_XADV

      DO I = 1, NUM_H

        IF (((MOD(I - 1,2) .EQ. 0) .AND. (MOD(P,2) .EQ. 1)) .OR. ((MOD(I - 1,2) .EQ. 1) .AND. (MOD(P,2) .EQ. 0))) THEN  !Evolve only even harmonics on odd spatial cells (centres) and vice-versa

          IF (0 .LE. I - 2) THEN                                  !First set of potential nonzeroes

            DO J = 1, NUM_V

              M_X_ADV%ROW(K) = X_POS(P) + NUM_V * (NUM_H - I) + J
              M_X_ADV%COLUMN(K) = X_POS(P + 1) + NUM_V * H_POS(I - 2) + J

              VAL = - (I - 1) / (2.00D00 * (I - 1) - 1.00D00) * V_GRID(J) / dxc(P)

              M_X_ADV%VALUE(K) = VAL

              K = K + 1                                                         !Advance nonzero count

              M_X_ADV%ROW(K) = X_POS(P) + NUM_V * (NUM_H - I) + J
              M_X_ADV%COLUMN(K) = X_POS(P - 1) + NUM_V * H_POS(I - 2) + J

              VAL = (I - 1) / (2.00D00 * (I - 1) - 1.00D00) * V_GRID(J) / dxc(P)

              M_X_ADV%VALUE(K) = VAL

              K = K + 1                                                         !Advance nonzero count

            END DO

          END IF

          IF (I .LE. L_MAX) THEN                                     !Second set of potential nonzeroes

            DO J = 1,NUM_V

              M_X_ADV%ROW(K) = X_POS(P) + NUM_V * (NUM_H - I) + J
              M_X_ADV%COLUMN(K) = X_POS(P + 1) + NUM_V * H_POS(I) + J

              VAL = - I / (2.00D00 * (I - 1) + 3.00D00) * V_GRID(J) / dxc(P)

              M_X_ADV%VALUE(K) = VAL

              K = K + 1                                                         !Advance nonzero count

              M_X_ADV%ROW(K) = X_POS(P) + NUM_V * (NUM_H - I) + J
              M_X_ADV%COLUMN(K) = X_POS(P - 1) + NUM_V * H_POS(I) + J

              VAL = I / (2.00D00 * (I - 1) + 3.00D00) * V_GRID(J) / dxc(P)

              M_X_ADV%VALUE(K) = VAL

              K = K + 1                                                         !Advance nonzero count

            END DO

          END IF

        END IF

      END DO

    END DO

    IF (RANK .EQ. SIZE - 1) THEN

      IF (.NOT. FIXED_BOUNDARY_DIV_SWITCH) THEN                                   !Fill last cell advection if not fixed boundary

        DO  I = 1, NUM_H

          IF (((MOD(I - 1,2) .EQ. 0) .AND. (MOD(NUM_X,2) .EQ. 1)) .OR. ((MOD(I - 1,2) .EQ. 1) .AND. (MOD(NUM_X,2) .EQ. 0))) THEN !Evolve only even harmonics on odd spatial cells (centres) and vice-versa

            IF (0 .LE. I - 2) THEN                                  !First set of potential nonzeroes

              DO J = 1, NUM_V

                M_X_ADV%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                M_X_ADV%COLUMN(K) = X_POS(NUM_X - 1) + NUM_V * H_POS(I - 2) + J

                VAL =  (I - 1) / (2.00D00 * (I - 1) - 1.00D00) * V_GRID(J) / dxc(NUM_X)

                M_X_ADV%VALUE(K) = VAL

                K = K + 1                                                         !Advance nonzero count

                IF (PERIODIC_BOUNDARY_SWITCH) THEN                                !If periodic boundary add advection to first cell

                  M_X_ADV%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                  M_X_ADV%COLUMN(K) =  NUM_V * H_POS(I - 2) + J

                  VAL = - (I - 1) / (2.00D00 * (I - 1) - 1.00D00) * V_GRID(J) / dxc(NUM_X)

                  M_X_ADV%VALUE(K) = VAL

                  K = K + 1                                                       !Advance nonzero count

                END IF

              END DO

            END IF

            IF (I .LE. L_MAX) THEN                                     !Second set of potential nonzeroes

              DO J = 1,NUM_V

                M_X_ADV%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                M_X_ADV%COLUMN(K) = X_POS(NUM_X - 1) + NUM_V * H_POS(I) + J

                VAL = I / (2.00D00 * (I - 1) + 3.00D00) * V_GRID(J) / dxc(NUM_X)

                M_X_ADV%VALUE(K) = VAL

                K = K + 1                                                         !Advance nonzero count

                IF (PERIODIC_BOUNDARY_SWITCH) THEN                                !If periodic boundary add advection to first cell

                  M_X_ADV%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                  M_X_ADV%COLUMN(K) = NUM_V * H_POS(I) + J

                  VAL = - I / (2.00D00 * (I - 1) + 3.00D00) * V_GRID(J) / dxc(NUM_X)

                  M_X_ADV%VALUE(K) = VAL

                  K = K + 1                                                       !Advance nonzero count

                END IF

              END DO

            END IF

          END IF

        END DO

      END IF

    END IF


  END SUBROUTINE FILL_X_ADV
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_X_ADV
!-------------------------------------------------------------------------------------------------------------------------------------
