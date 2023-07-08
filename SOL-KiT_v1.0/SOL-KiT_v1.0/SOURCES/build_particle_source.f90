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
!> Contains particle source builder
MODULE BUILD_PARTICLE_SOURCE

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
!> Fills electron source submatrix
  SUBROUTINE FILL_PART_SOURCE(t, n_e,T_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), & !< Lagged density vector
                                T_e(NUM_X)    !< Lagged temperature vector
    INTEGER, INTENT(IN) :: t !< Timestep number

    INTEGER :: I, OFFSET_UP_H, MIN_C, MAX_C, RANK, IERR

    REAL(KIND=DEFAULT_REAL) :: F_R, &
                    VAL, &
                    T_S(NUM_X)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)


    IF (t .LE. T_PART_SOURCE) THEN

!Calculate source flux constant

!Set source temp
      T_S = TEMP_PART_SOURCE

      IF (PART_SOURCE_BACKGROUND_TEMP_SWITCH) T_S = T_e

!Determine offsets and cell centre range
      OFFSET_UP_H = 0

      IF (FIXED_BOUNDARY_UP_SWITCH) OFFSET_UP_H = 2

      IF (RANK .EQ. 0) THEN

        MIN_C = 1 + OFFSET_UP_H

      ELSE

        MIN_C = MIN_X + MOD(MIN_X + 1, 2)

      END IF

      MAX_C = MAX_X - MOD(MAX_X+1,2)

      F_R = P_FLUX_IN / (X_GRID(2 * N_PART_SOURCE - 1 + OFFSET_UP_H)+dxc(2 * N_PART_SOURCE - 1 + OFFSET_UP_H)/2&
      -(X_GRID(1 + OFFSET_UP_H)-dxc(1 + OFFSET_UP_H)/2))

      IF ((MAX_C + 1)/2 .GT. N_PART_SOURCE + OFFSET_UP_H/2)  MAX_C = 2 * N_PART_SOURCE - 1 + OFFSET_UP_H
!Fill matrix
      DO I = 1, PART_SOURCE_SP%N_NZ

        VAL = F_R * 4.00D00 * PI * V_GRID(SOURCE_DATA(I)%K) ** 2 * V_GRID_WIDTH(SOURCE_DATA(I)%K) / n_e(SOURCE_DATA(I)%POS) * &
              (PI * T_S(SOURCE_DATA(I)%POS)) **(-3.00/2.00) * EXP(- V_GRID(SOURCE_DATA(I)%J)**2/T_S(SOURCE_DATA(I)%POS))

        LOCAL_M%VALUE(MARKER_PART_SOURCE(I)) = LOCAL_M%VALUE(MARKER_PART_SOURCE(I)) + VAL


      END DO

    END IF


  END SUBROUTINE FILL_PART_SOURCE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills ion source submatrix
  SUBROUTINE FILL_ION_PART_SOURCE(t, n_e,T_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), & !< Lagged electron density vector
                                           T_e(NUM_X)    !< Lagged temperature vector
    INTEGER, INTENT(IN) :: t !< Timestep number

    INTEGER :: I, J, OFFSET_UP_H, MIN_C, MAX_C, RANK, IERR

    REAL(KIND=DEFAULT_REAL) :: F_R, &
                    VAL, &
                    T_S(NUM_X), &
                    NUM_PART_RATE(MIN_X:MAX_X) !< Numerical particle rate - consistent with electron source


    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)


    IF (t .LE. T_PART_SOURCE) THEN

      NUM_PART_RATE = 0

!Calculate source fluxconstant

      T_S = TEMP_PART_SOURCE

      IF (PART_SOURCE_BACKGROUND_TEMP_SWITCH) T_S = T_e

!Determine offsets and cell centre range
      OFFSET_UP_H = 0

      IF (FIXED_BOUNDARY_UP_SWITCH) OFFSET_UP_H = 2

      IF (RANK .EQ. 0) THEN

        MIN_C = 1 + OFFSET_UP_H

      ELSE

        MIN_C = MIN_X + MOD(MIN_X + 1, 2)

      END IF

      MAX_C = MAX_X - MOD(MAX_X+1,2)

      F_R = P_FLUX_IN / (X_GRID(2 * N_PART_SOURCE - 1 + OFFSET_UP_H)+dxc(2 * N_PART_SOURCE - 1 + OFFSET_UP_H)/2&
      -(X_GRID(1 + OFFSET_UP_H)-dxc(1 + OFFSET_UP_H)/2))

      IF ((MAX_C + 1)/2 .GT. N_PART_SOURCE + OFFSET_UP_H/2)  MAX_C = 2 * N_PART_SOURCE - 1 + OFFSET_UP_H

      DO I = MIN_C, MAX_C,2

        DO J = 1, NUM_V

          NUM_PART_RATE(I)=NUM_PART_RATE(I)+ F_R * 4.00D00*PI*V_GRID(J)**2 * V_GRID_WIDTH(J) &
          * (PI * T_S(I)) **(-3.00/2.00) * EXP(- V_GRID(J)**2/T_S(I))

        END DO

      END DO
!Fill matrix
      DO I = 1, ION_PART_SOURCE_SP%N_NZ

        VAL = NUM_PART_RATE(ION_SOURCE_DATA(I)%POS)*4.00D00 * PI * &
        V_GRID(ION_SOURCE_DATA(I)%K) ** 2 * V_GRID_WIDTH(ION_SOURCE_DATA(I)%K) / n_e(ION_SOURCE_DATA(I)%POS)

        LOCAL_M%VALUE(MARKER_ION_PART_SOURCE(I)) = LOCAL_M%VALUE(MARKER_ION_PART_SOURCE(I)) + VAL

      END DO

    END IF

  END SUBROUTINE FILL_ION_PART_SOURCE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_PARTICLE_SOURCE
!-------------------------------------------------------------------------------------------------------------------------------------
