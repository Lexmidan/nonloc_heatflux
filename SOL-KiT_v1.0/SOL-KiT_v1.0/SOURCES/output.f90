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
!> Contains all output routines
MODULE OUTPUT

  USE GRID
  USE NEUT_AND_HEAT
  USE NORMALIZATION
  USE SWITCHES
  USE SOLVER_PARAMS
  USE F_INIT
  USE MPI
  USE VAR_KIND_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output spatial vector vec at timestep t for quantity s
  SUBROUTINE OUTPUT_VECTOR_X(vec,s,t)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: vec(NUM_X) !< Output vector
    CHARACTER (LEN = *), INTENT(IN) :: s !< Output name
    INTEGER, INTENT(IN) :: t !< Output timestep

    CHARACTER (LEN = 60) :: FILENAME, &
                            x1

    INTEGER :: IR , I

    WRITE(x1,'(I5.5)') t

    FILENAME = 'OUTPUT/'// s // '/' // s // '_' // TRIM(x1) // '.txt'

    OPEN(10, FILE = FILENAME, ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    DO I = 1, NUM_X

      WRITE(10,'(ES22.14)') vec(I)

    END DO

    CLOSE(10, IOSTAT = IR)

  END SUBROUTINE OUTPUT_VECTOR_X
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output spatial and velocity grid
  SUBROUTINE OUTPUT_GRIDS

    IMPLICIT NONE

    INTEGER :: IR , I

    OPEN(10, FILE = 'OUTPUT/GRIDS/X_GRID.txt', ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    DO I = 1, NUM_X

      WRITE(10, '(G22.14)') X_GRID(I)

    END DO

    CLOSE(10, IOSTAT = IR)

    OPEN(10, FILE = 'OUTPUT/GRIDS/V_GRID.txt', ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    DO I = 1, NUM_V

      WRITE(10, '(G22.14)') V_GRID(I)

    END DO

    CLOSE(10, IOSTAT = IR)

    OPEN(10, FILE = 'OUTPUT/GRIDS/V_GRID_WIDTH.txt', ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    DO I = 1, NUM_V

      WRITE(10, '(G22.14)') V_GRID_WIDTH(I)

    END DO

    CLOSE(10, IOSTAT = IR)

  END SUBROUTINE OUTPUT_GRIDS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output neutral state densities n_n at time t
  SUBROUTINE OUTPUT_NEUTRALS(n_n,t)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n(NUM_NEUTRALS,NUM_X) !< Matrix containing neutral state density vectors

    INTEGER, INTENT(IN) :: t !< Output time

    CHARACTER (LEN  = 60) :: FILENAME, x1, ROWFMT

    INTEGER :: I, IR, J

    WRITE(x1,'(I5.5)') t

    WRITE(ROWFMT, '(A,I3,A)') '(',NUM_NEUTRALS,'(1X,ES22.14E3))'

    FILENAME = 'OUTPUT/NEUTRAL_DENS/NEUTRAL_DENS_' // trim(x1) // '.txt'

    OPEN(10,FILE = FILENAME, ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    DO I = 1, NUM_X

      WRITE(10, ROWFMT) (n_n(J,I) , J = 1, NUM_NEUTRALS)

    END DO

    CLOSE(10, IOSTAT = IR)

  END SUBROUTINE OUTPUT_NEUTRALS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output distribution functions at time t
  SUBROUTINE OUTPUT_F(f,t)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(0:L_MAX, NUM_V, NUM_X) !< Chopped-up distribution function into harmonics
    INTEGER, INTENT(IN) :: t !< Output time

    CHARACTER (LEN = 60) :: ROWFMT, FILENAME, x1, x2

    INTEGER :: I, K, IR, P

    WRITE(ROWFMT,'(A,I4,A)') '(',NUM_X,'(1X,ES22.14E3))'

    WRITE(x2,'(I5.5)') t

    DO I = 0, L_MAX

      WRITE(x1,'(I3)') I

          FILENAME = 'OUTPUT/DIST_F/F_L'// TRIM(x1) // '_' // TRIM(x2) // '.txt'

          OPEN(10, FILE = FILENAME, ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

          DO K = 1, NUM_V

            WRITE(10,ROWFMT) (f(I,K,P), P = 1, NUM_X)

          END DO

          CLOSE(10, IOSTAT = IR)

        x1 = ''

    END DO

  END SUBROUTINE OUTPUT_F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output full F vector into restart folder
  SUBROUTINE OUTPUT_RESTART_F(vec)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: vec(DIM_F) !< Output vector

    INTEGER :: IR , I

    OPEN(10, FILE = 'INPUT/RESTART/VAR_VEC_INPUT.txt', ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    DO I = 1, DIM_F

      WRITE(10,'(ES22.14E3)') vec(I)

    END DO

    CLOSE(10, IOSTAT = IR)

  END SUBROUTINE OUTPUT_RESTART_F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output restart data
  SUBROUTINE OUTPUT_RESTART_INFO(TIME,TIMESTEP)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: TIME !< Elapse total time in run
    INTEGER, INTENT(IN) :: TIMESTEP !< Current timestep

    CHARACTER(LEN = *), PARAMETER :: FORMAT_F = '(1X, A40," = ", G19.12, 5X, A)' , & !< Float format specifier
                                     FORMAT_I = '(1X, A40," = ", I14, 5X, A)' , & !< Integer format specifier
                                     FORMAT_L = '(1X, A40," = ", L14, 5X, A)',& !< Logical format specifier
                                     FORMAT_C = '(1X, A40," = ", A14, 5X, A)' !< Character format specifier

    INTEGER :: IR

    OPEN(11, FILE = 'INPUT/RESTART/RESTART_INFO.txt', ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')


    WRITE(11,FORMAT_F) 'TIME',TIME,'Elapsed real time'
    WRITE(11,FORMAT_I) 'TIMESTEP',TIMESTEP,'Current timestep'
    WRITE(11,FORMAT_F) 'TOTAL_L',TOTAL_L,'Simulation length in meters'
    WRITE(11, *)

    WRITE(11, FORMAT_I) 'NUM_NEUTRALS', NUM_NEUTRALS, 'Number of tracked neutral state'
    WRITE(11, FORMAT_I) 'L_MAX', L_MAX, 'L-number of maximum resolved harmonic'
    WRITE(11, FORMAT_I) 'NUM_V', NUM_V, 'Number of cells in velocity space'
    WRITE(11, FORMAT_F) 'dv', dv, 'Size of velocity cells'
    WRITE(11, FORMAT_L) 'FULL_FLUID_MODE', FULL_FLUID_MODE, 'Full fluid mode on'
    WRITE(11, *)

    CLOSE(11, IOSTAT = IR)

  END SUBROUTINE OUTPUT_RESTART_INFO
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output sheath data
  SUBROUTINE OUTPUT_SHEATH_DATA(gamma, phi, v_c, flow_lim, bv, TIMESTEP )

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: gamma, & !< Sheath transmission coefficient
                                phi, & !< Sheath potential normalized to boundary temperature
                                v_c, & !< Cut-off Velocity
                                flow_lim !< Extrapolation flow limiter

    INTEGER, INTENT(IN) :: TIMESTEP

    LOGICAL, INTENT(IN) :: bv !< True if using Bohm value at boundary

    CHARACTER (LEN = 60) :: FILENAME, &
                            x1

    INTEGER :: IR

    WRITE(x1,'(I5.5)') TIMESTEP

    FILENAME = 'OUTPUT/SHEATH_DATA/SHEATH_DATA_' // TRIM(x1) // '.txt'

    OPEN(10, FILE = FILENAME, ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    WRITE(10,'(ES22.14)') gamma
    WRITE(10,'(ES22.14)') phi
    WRITE(10,'(ES22.14)') v_c
    WRITE(10,'(ES22.14)') flow_lim
    IF (bv) THEN

      WRITE(10,'(ES22.14)') 1.0D00

    ELSE

      WRITE(10,'(ES22.14)') 0.0D00

    END IF

    CLOSE(10, IOSTAT = IR)

  END SUBROUTINE OUTPUT_SHEATH_DATA
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output timestep data
  SUBROUTINE OUTPUT_TIMESTEP_DATA(dt, TIMESTEP,TIME)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: dt, & !< Current timestep length
                                TIME !< Current elapsed time

    INTEGER, INTENT(IN) :: TIMESTEP


    CHARACTER (LEN = 60) :: FILENAME, &
                            x1

    INTEGER :: IR

    WRITE(x1,'(I5.5)') TIMESTEP

    FILENAME = 'OUTPUT/TIMESTEP_DATA/TIMESTEP_DATA_' // TRIM(x1) // '.txt'

    OPEN(10, FILE = FILENAME, ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    WRITE(10,'(ES19.12)') dt
    WRITE(10,'(ES19.12)') TIME

    CLOSE(10, IOSTAT = IR)

  END SUBROUTINE OUTPUT_TIMESTEP_DATA
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Output total density data
  SUBROUTINE OUTPUT_TOT_DENS_DATA(TIMESTEP,n_tot,n_n_tot,old_tot_dens,old_n_tot_dens)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_tot, & !< Current total density
                                n_n_tot, &!< Current neutral total density
                                old_tot_dens , &!< Old total density
                                old_n_tot_dens !< Old neutral total density

    INTEGER, INTENT(IN) :: TIMESTEP


    CHARACTER (LEN = 60) :: FILENAME, &
                            x1

    INTEGER :: IR

    WRITE(x1,'(I5.5)') TIMESTEP

    FILENAME = 'OUTPUT/TOT_DENS_DATA/TOT_DENS_DATA_' // TRIM(x1) // '.txt'

    OPEN(10, FILE = FILENAME, ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    WRITE(10,'(ES22.14)') n_tot
    WRITE(10,'(ES22.14)') n_n_tot
    WRITE(10,'(ES22.14)') ABS(old_tot_dens - n_tot)/n_tot
    IF(n_n_tot .GT. 0.0D00) WRITE(10,'(ES22.14)') ABS(old_n_tot_dens - n_n_tot)/n_n_tot

    CLOSE(10, IOSTAT = IR)

  END SUBROUTINE OUTPUT_TOT_DENS_DATA
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE OUTPUT
!-------------------------------------------------------------------------------------------------------------------------------------
