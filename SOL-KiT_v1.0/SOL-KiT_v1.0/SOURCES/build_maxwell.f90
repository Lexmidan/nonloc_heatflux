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
!> Contains submatrix for Maxwell-Ampere's law and the subroutine that fills it
MODULE BUILD_MAXWELL

  USE GRID
  USE SPARSE
  USE NORMALIZATION
  USE SWITCHES
  USE MPI
  USE VAR_KIND_DATA
  
  IMPLICIT NONE

  TYPE (SPARSE_MAT) :: M_MAXWELL !< Maxwell equation submatrix

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills Maxwell equation submatrix
  SUBROUTINE FILL_MAXWELL

    IMPLICIT NONE

    INTEGER :: I, J, K, N_NZ, RANK, IERR, SIZE, OFFSET

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    N_NZ = NUM_V * (nd_loc - (LOC_NUM_C + OFFSET))                                              !Number of nonzeroes

    CALL ALLOCATE_SPARSE(M_MAXWELL, LOC_ROWS, DIM_F, N_NZ)

    K = 1                                                                       !Advance nonzero count

!Fill auxillary sparse matrix
    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        DO J = 1, NUM_V

          M_MAXWELL%ROW(K) = X_POS(I) + NUM_H*NUM_V + NUM_NEUTRALS + 1
          M_MAXWELL%COLUMN(K) = X_POS(I) + NUM_V * H_POS(1) + J

          VAL = MAX_E_J_0 * 4.00D00 * PI / 3.00D00 * V_GRID(J) ** 3 * V_GRID_WIDTH(J)

          M_MAXWELL%VALUE(K) = VAL

          K = K + 1                                                             !Advance nonzero count

        END DO

      END IF

    END DO

  END SUBROUTINE FILL_MAXWELL
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_MAXWELL
!-------------------------------------------------------------------------------------------------------------------------------------
