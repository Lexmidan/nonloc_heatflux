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
!> Contains neutral and heating parameters
MODULE NEUT_AND_HEAT

 USE SWITCHES
 USE INPUT
 USE MPI
 USE VAR_KIND_DATA

! Neutral parameters

  INTEGER :: NUM_NEUTRALS, & !< Number of tracked neutral state
             NUM_CS_ANGLES !< Number of resolved angles in differential cross-section

  REAL(KIND=DEFAULT_REAL) :: NEUTRAL_DENSITY, & !< Reference neutral density
            NEUTRAL_TEMP, & !< Neutral temperature
            SIGMA_EL !< Electron-neutral elastic cross-section in units of SIGMA_0

! Heating parameters

  INTEGER :: N_HEATING, & !< Number of heated upstream cell centres
             T_HEATING !< Last heating timestep

  REAL(KIND=DEFAULT_REAL) :: HEAT_POWER !< Heating power in MW/m^2

! Recycling and divertor Mach number

  REAL(KIND=DEFAULT_REAL) :: MACH_N_DIV, & !< Divertor Mach number
                             REC_R !< Recycling ratio

! Particle source parameters

  INTEGER :: N_PART_SOURCE, & !< Number of source cells upstream
             T_PART_SOURCE !< Last source timestep

  REAL(KIND=DEFAULT_REAL) :: P_FLUX_IN, & !< Effective input flux
                  TEMP_PART_SOURCE !< Source temperature if not using background temperature
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes neutral and heating/particle source parameters from input
  SUBROUTINE INIT_HEAT_NEUT

    IMPLICIT NONE

    INTEGER :: RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    NUM_NEUTRALS = 0

    IF (NEUTRAL_TRACK_SWITCH) THEN

      IF (RANK .EQ. 0) CALL INPUT_NEUTRALS(NUM_NEUTRALS, &                                       !Input neutral data
                                           NUM_CS_ANGLES,&
                                           NEUTRAL_DENSITY, &
                                           NEUTRAL_TEMP, &
                                           SIGMA_EL&
                                           )

      CALL MPI_Bcast(NUM_NEUTRALS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(NUM_CS_ANGLES,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(NEUTRAL_DENSITY,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(NEUTRAL_TEMP,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(SIGMA_EL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)



    END IF

    IF (RECYCLING_SWITCH .OR. PLASMA_SINK_SWITCH) THEN

      IF (RANK .EQ. 0) CALL INPUT_REC_R(MACH_N_DIV, &
                                        REC_R)                             !Input recycling rate and Mach number

      CALL MPI_Bcast(MACH_N_DIV,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(REC_R,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    END IF

    IF (HEATING_SWITCH) THEN

      IF (RANK .EQ. 0) CALL INPUT_HEATING(N_HEATING, &                         !Input heating data
                                          T_HEATING, &
                                          HEAT_POWER)

      CALL MPI_Bcast(N_HEATING,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(T_HEATING,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(HEAT_POWER,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    END IF

    IF (PART_SOURCE_SWITCH) THEN

      IF (RANK .EQ. 0) CALL INPUT_PART_SOURCE(N_PART_SOURCE, &                         !Input particle source data
                                              T_PART_SOURCE, &
                                              P_FLUX_IN, &
                                              TEMP_PART_SOURCE)

      CALL MPI_Bcast(N_PART_SOURCE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(T_PART_SOURCE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(P_FLUX_IN,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Bcast(TEMP_PART_SOURCE,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    END IF

  END SUBROUTINE INIT_HEAT_NEUT
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE NEUT_AND_HEAT
!-------------------------------------------------------------------------------------------------------------------------------------
