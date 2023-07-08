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
!> Contains all file reading routines
MODULE INPUT

  USE MPI
  USE VAR_KIND_DATA

  IMPLICIT NONE

  PRIVATE :: IR, FORMAT_F, FORMAT_I, FORMAT_L

  INTEGER :: IR !< IOSTAT holder

  CHARACTER(LEN = *), PARAMETER :: FORMAT_F = '(44X, G19.12E3)' , & !< Float format specifier
                                   FORMAT_I = '(44X, I14)' , & !< Integer format specifier
                                   FORMAT_L = '(44X, L14)', & !< Logical format specifier
                                   FORMAT_C = '(44X, A14)' !< Character format specifier
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads normalization input
SUBROUTINE INPUT_NORM(ION_Z, &
                        ATOMIC_A, &
                        TEMP_0_EV, &
                        DENSITY_0 &
                        )
  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: ION_Z, & !< Plasma Z
                         ATOMIC_A, & !< Atomic mass of main plasma particles
                         TEMP_0_EV, & !< Temperature normalization in eV
                         DENSITY_0 !< Density normalization

  OPEN(10, FILE = "INPUT/NORMALIZATION_INPUT.txt", ACTION = 'READ', IOSTAT = IR, STATUS = 'OLD')

  READ(10, FORMAT_F) ION_Z
  READ(10, FORMAT_F) ATOMIC_A
  READ(10, FORMAT_F) TEMP_0_EV
  READ(10, FORMAT_F) DENSITY_0

  CLOSE(10, IOSTAT = IR)

END SUBROUTINE INPUT_NORM
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads switches
SUBROUTINE INPUT_SWITCHES(MAXWELL_SWITCH, &
                          X_ADV_SWITCH, &
                          E_ADV_SWITCH,&
                          COLL_EE_0_SWITCH, &
                          COLL_EE_L_SWITCH, &
                          DIAG_EE_L_SWITCH, &
                          COLL_EI_L_SWITCH, &
                          COLL_EN_EL_0_SWITCH, &
                          COLL_EN_EL_L_SWITCH, &
                          COLL_EN_EX, &
                          COLL_EN_ION, &
                          COLL_RECOMB, &
                          LOGARITHMIC_GRID_SWITCH, &
                          PERIODIC_BOUNDARY_SWITCH, &
                          FIXED_BOUNDARY_UP_SWITCH, &
                          FIXED_BOUNDARY_DIV_SWITCH, &
                          NO_FLOW_BOUNDARY_UP_SWITCH, &
                          NO_FLOW_BOUNDARY_DIV_SWITCH, &
                          NEUTRAL_TRACK_SWITCH, &
                          NEUTRAL_DIFFUSION_SWITCH, &
                          RECYCLING_SWITCH, &
                          FAST_DETAILED_BALANCE_SWITCH, &
                          HEATING_SWITCH, &
                          PLASMA_SINK_SWITCH, &
                          OUTPUT_DENSITY, &
                          OUTPUT_TEMP, &
                          OUTPUT_FLOW_VEL, &
                          OUTPUT_HEAT_FLOW, &
                          OUTPUT_E_FIELD, &
                          OUTPUT_SH_TEST, &
                          OUTPUT_RATE_DATA, &
                          OUTPUT_NEUTRAL_DATA, &
                          OUTPUT_ATOMIC_EN_TEST_SWITCH,&
                          OUTPUT_QN_TEST, &
                          OUTPUT_CURRENT_TEST, &
                          LINEAR_INIT, &
                          DROP_INIT, &
                          TWO_POINT_M_INIT, &
                          UNIFORM_NEUTRAL_INIT, &
                          NEUTRAL_CLOUD_INIT, &
                          DENSITY_FROM_FILE_INIT, &
                          TEMPERATURE_FROM_FILE_INIT, &
                          NEUTRAL_GROUND_DENS_FROM_FILE_INIT, &
                          ION_VEL_FROM_FILE_INIT, &
                          LOCAL_INIT_SWITCH,  &
                          LOCAL_SAHA_BOLTZMANN_INIT_SWITCH, &
                          RESTART_SWITCH, &
                          SAVE_RESTART_SWITCH, &
                          CONTINUE_RUN_SWITCH, &
                          ADAPTIVE_RESTART_SWITCH, &
                          COLD_ION_FLUID_SWITCH, &
                          ION_CONT_OFF_SWITCH, &
                          ION_CONV_UPWINDING_SWITCH, &
                          SIMPLE_CX_SWITCH, &
                          ION_EL_TEMP_SWITCH, &
                          NO_EXTRAPOLATION_SWITCH, &
                          SONIC_OUTFLOW_DIV_SWITCH, &
                          PART_SOURCE_SWITCH, &
                          PART_SOURCE_BACKGROUND_TEMP_SWITCH, &
                          ADAPTIVE_TIMESTEP_SWITCH, &
                          FULL_FLUID_MODE, &
                          Z_PROFILE_FROM_FILE, &
                          X_GRID_FROM_FILE &
                          )

  IMPLICIT NONE

  LOGICAL, INTENT(OUT) :: MAXWELL_SWITCH, &
                          X_ADV_SWITCH, &
                          E_ADV_SWITCH,&
                          COLL_EE_0_SWITCH, &
                          COLL_EE_L_SWITCH, &
                          DIAG_EE_L_SWITCH, &
                          COLL_EI_L_SWITCH, &
                          COLL_EN_EL_0_SWITCH, &
                          COLL_EN_EL_L_SWITCH, &
                          COLL_EN_EX, &
                          COLL_EN_ION, &
                          COLL_RECOMB, &
                          LOGARITHMIC_GRID_SWITCH, &
                          PERIODIC_BOUNDARY_SWITCH, &
                          FIXED_BOUNDARY_UP_SWITCH, &
                          FIXED_BOUNDARY_DIV_SWITCH, &
                          NO_FLOW_BOUNDARY_UP_SWITCH, &
                          NO_FLOW_BOUNDARY_DIV_SWITCH, &
                          NEUTRAL_TRACK_SWITCH, &
                          NEUTRAL_DIFFUSION_SWITCH, &
                          RECYCLING_SWITCH, &
                          FAST_DETAILED_BALANCE_SWITCH, &
                          HEATING_SWITCH, &
                          PLASMA_SINK_SWITCH, &
                          OUTPUT_DENSITY, &
                          OUTPUT_TEMP, &
                          OUTPUT_FLOW_VEL, &
                          OUTPUT_HEAT_FLOW, &
                          OUTPUT_E_FIELD, &
                          OUTPUT_SH_TEST, &
                          OUTPUT_RATE_DATA, &
                          OUTPUT_NEUTRAL_DATA, &
                          OUTPUT_ATOMIC_EN_TEST_SWITCH, &
                          OUTPUT_QN_TEST, &
                          OUTPUT_CURRENT_TEST, &
                          LINEAR_INIT, &
                          DROP_INIT, &
                          TWO_POINT_M_INIT, &
                          UNIFORM_NEUTRAL_INIT, &
                          NEUTRAL_CLOUD_INIT, &
                          DENSITY_FROM_FILE_INIT, &
                          TEMPERATURE_FROM_FILE_INIT, &
                          NEUTRAL_GROUND_DENS_FROM_FILE_INIT, &
                          ION_VEL_FROM_FILE_INIT, &
                          LOCAL_INIT_SWITCH, &
                          LOCAL_SAHA_BOLTZMANN_INIT_SWITCH, &
                          RESTART_SWITCH, &
                          SAVE_RESTART_SWITCH, &
                          CONTINUE_RUN_SWITCH, &
                          ADAPTIVE_RESTART_SWITCH, &
                          COLD_ION_FLUID_SWITCH, &
                          ION_CONT_OFF_SWITCH, &
                          ION_CONV_UPWINDING_SWITCH, &
                          SIMPLE_CX_SWITCH, &
                          ION_EL_TEMP_SWITCH, &
                          NO_EXTRAPOLATION_SWITCH, &
                          SONIC_OUTFLOW_DIV_SWITCH, &
                          PART_SOURCE_SWITCH, &
                          PART_SOURCE_BACKGROUND_TEMP_SWITCH, &
                          ADAPTIVE_TIMESTEP_SWITCH, &
                          FULL_FLUID_MODE, &
                          Z_PROFILE_FROM_FILE, &
                          X_GRID_FROM_FILE


  OPEN(11, FILE = 'INPUT/SWITCHES_INPUT.txt', ACTION = 'READ', IOSTAT = IR, STATUS = 'OLD')

  READ(11,*)
  READ(11, FORMAT_L) MAXWELL_SWITCH
  READ(11, FORMAT_L) X_ADV_SWITCH
  READ(11, FORMAT_L) E_ADV_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) COLL_EE_0_SWITCH
  READ(11, FORMAT_L) COLL_EE_L_SWITCH
  READ(11, FORMAT_L) DIAG_EE_L_SWITCH
  READ(11, FORMAT_L) COLL_EI_L_SWITCH
  READ(11, FORMAT_L) COLL_EN_EL_0_SWITCH
  READ(11, FORMAT_L) COLL_EN_EL_L_SWITCH
  READ(11, FORMAT_L) COLL_EN_EX
  READ(11, FORMAT_L) COLL_EN_ION
  READ(11, FORMAT_L) COLL_RECOMB
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) LOGARITHMIC_GRID_SWITCH
  READ(11, FORMAT_L) PERIODIC_BOUNDARY_SWITCH
  READ(11, FORMAT_L) FIXED_BOUNDARY_UP_SWITCH
  READ(11, FORMAT_L) FIXED_BOUNDARY_DIV_SWITCH
  READ(11, FORMAT_L) NO_FLOW_BOUNDARY_UP_SWITCH
  READ(11, FORMAT_L) NO_FLOW_BOUNDARY_DIV_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) NEUTRAL_TRACK_SWITCH
  READ(11, FORMAT_L) NEUTRAL_DIFFUSION_SWITCH
  READ(11, FORMAT_L) RECYCLING_SWITCH
  READ(11, FORMAT_L) FAST_DETAILED_BALANCE_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) HEATING_SWITCH
  READ(11, FORMAT_L) PLASMA_SINK_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) OUTPUT_DENSITY
  READ(11, FORMAT_L) OUTPUT_TEMP
  READ(11, FORMAT_L) OUTPUT_FLOW_VEL
  READ(11, FORMAT_L) OUTPUT_HEAT_FLOW
  READ(11, FORMAT_L) OUTPUT_E_FIELD
  READ(11, FORMAT_L) OUTPUT_SH_TEST
  READ(11, FORMAT_L) OUTPUT_RATE_DATA
  READ(11, FORMAT_L) OUTPUT_NEUTRAL_DATA
  READ(11, FORMAT_L) OUTPUT_ATOMIC_EN_TEST_SWITCH
  READ(11, FORMAT_L) OUTPUT_QN_TEST
  READ(11, FORMAT_L) OUTPUT_CURRENT_TEST
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) LINEAR_INIT
  READ(11, FORMAT_L) DROP_INIT
  READ(11, FORMAT_L) TWO_POINT_M_INIT
  READ(11, FORMAT_L) UNIFORM_NEUTRAL_INIT
  READ(11, FORMAT_L) NEUTRAL_CLOUD_INIT
  READ(11, FORMAT_L) DENSITY_FROM_FILE_INIT
  READ(11, FORMAT_L) TEMPERATURE_FROM_FILE_INIT
  READ(11, FORMAT_L) NEUTRAL_GROUND_DENS_FROM_FILE_INIT
  READ(11, FORMAT_L) ION_VEL_FROM_FILE_INIT
  READ(11, FORMAT_L) LOCAL_INIT_SWITCH
  READ(11, FORMAT_L) LOCAL_SAHA_BOLTZMANN_INIT_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) RESTART_SWITCH
  READ(11, FORMAT_L) SAVE_RESTART_SWITCH
  READ(11, FORMAT_L) CONTINUE_RUN_SWITCH
  READ(11, FORMAT_L) ADAPTIVE_RESTART_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) COLD_ION_FLUID_SWITCH
  READ(11, FORMAT_L) ION_CONT_OFF_SWITCH
  READ(11, FORMAT_L) ION_CONV_UPWINDING_SWITCH
  READ(11, FORMAT_L) SIMPLE_CX_SWITCH
  READ(11, FORMAT_L) ION_EL_TEMP_SWITCH
  READ(11, FORMAT_L) NO_EXTRAPOLATION_SWITCH
  READ(11, FORMAT_L) SONIC_OUTFLOW_DIV_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) PART_SOURCE_SWITCH
  READ(11, FORMAT_L) PART_SOURCE_BACKGROUND_TEMP_SWITCH
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_L) ADAPTIVE_TIMESTEP_SWITCH
  READ(11, FORMAT_L) FULL_FLUID_MODE
  READ(11, FORMAT_L) Z_PROFILE_FROM_FILE
  READ(11, FORMAT_L) X_GRID_FROM_FILE

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_SWITCHES
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads neutral parameters
SUBROUTINE INPUT_NEUTRALS(NUM_NEUTRALS, &
                          NUM_CS_ANGLES,&
                          NEUTRAL_DENSITY, &
                          NEUTRAL_TEMP, &
                          SIGMA_EL&
                          )

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: NUM_NEUTRALS, & !< Number of tracked neutral state
                          NUM_CS_ANGLES !< Number of resolved angles in differential cross-section

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: NEUTRAL_DENSITY, & !< Reference neutral density
                         NEUTRAL_TEMP, & !< Neutral temperature
                         SIGMA_EL !< Electron-neutral elastic cross-section in units of SIGMA_0

  OPEN(11, FILE = 'INPUT/NEUT_AND_HEAT_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  READ(11,*)
  READ(11, FORMAT_I)  NUM_NEUTRALS
  READ(11, FORMAT_I)  NUM_CS_ANGLES
  READ(11, FORMAT_F)  NEUTRAL_DENSITY
  READ(11, FORMAT_F)  NEUTRAL_TEMP
  READ(11, FORMAT_F)  SIGMA_EL
  READ(11, *)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_NEUTRALS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads heating parameters
SUBROUTINE INPUT_HEATING(N_HEATING, &
                         T_HEATING, &
                         HEAT_POWER)

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: N_HEATING, & !< Number of heated upstream cell centres
                          T_HEATING !< Last heating timestep

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: HEAT_POWER !< Heating power in MW/m^2

  INTEGER :: I

  OPEN(11, FILE = 'INPUT/NEUT_AND_HEAT_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  DO I = 1,7

    READ(11,*)

  END DO

  READ(11,*)
  READ(11, FORMAT_I)  N_HEATING
  READ(11, FORMAT_I)  T_HEATING
  READ(11, FORMAT_F)  HEAT_POWER
  READ(11,*)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_HEATING
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads particle source parameters
SUBROUTINE INPUT_PART_SOURCE(N_PART_SOURCE, &
                            T_PART_SOURCE, &
                            P_FLUX_IN, &
                            TEMP_PART_SOURCE)
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: N_PART_SOURCE, &     !< Total number of source cells
                          T_PART_SOURCE        !< Last source timestep


  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: P_FLUX_IN, &    !< Input particle flux
                               TEMP_PART_SOURCE !< Input particle temperature (if not background temperature)
  INTEGER :: I

  OPEN(11, FILE = 'INPUT/NEUT_AND_HEAT_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  DO I = 1,16

    READ(11,*)

  END DO

  READ(11,*)
  READ(11, FORMAT_I)  N_PART_SOURCE
  READ(11, FORMAT_I)  T_PART_SOURCE
  READ(11, FORMAT_F)  P_FLUX_IN
  READ(11, FORMAT_F)  TEMP_PART_SOURCE
  READ(11,*)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_PART_SOURCE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads recycling rate and Mach number
SUBROUTINE INPUT_REC_R(MACH_N_DIV,REC_R)

  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: MACH_N_DIV, & !< Mach number at divertor boundary
                                          REC_R !< Recycling ratio

  INTEGER :: I

  OPEN(11, FILE = 'INPUT/NEUT_AND_HEAT_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  DO I = 1,12

    READ(11,*)

  END DO

  READ(11,*)
  READ(11, FORMAT_F)  MACH_N_DIV
  READ(11, FORMAT_F)  REC_R
  READ(11,*)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_REC_R
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads grid parameters
SUBROUTINE INPUT_GRID(TIMESTEP_NUM, &
                      PRETIMESTEP_NUM, &
                      dt, &
                      pre_dt, &
                      T_SAVE, &
                      L_MAX, &
                      NUM_V, &
                      dv, &
                      V_GRID_MULT, &
                      NUM_C, &
                      dx, &
                      SMALL_dx &
                      )

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: TIMESTEP_NUM , & !< Number of full length timesteps
                          PRETIMESTEP_NUM !< Number of small pre-timesteps

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: dt , & !< Timestep size
                         pre_dt !< Pre-timestep size

  INTEGER, INTENT(OUT) :: T_SAVE !< Save interval (save every T_SAVE timesteps)

  INTEGER, INTENT(OUT) :: L_MAX !< L-number of maximum resolved harmonic

  INTEGER, INTENT(OUT) :: NUM_V !< Number of cells in velocity space

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: dv, & !< Size of velocity cells
                               V_GRID_MULT !< Velocity cell width common ratio

  INTEGER, INTENT(OUT) :: NUM_C !< Number of spatial cells

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: dx , & !< Size of spatial cells
                         SMALL_dx  !< Size of divertor cell

  OPEN(11, FILE = 'INPUT/GRID_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  READ(11, *)
  READ(11, FORMAT_I) TIMESTEP_NUM
  READ(11, FORMAT_I)  PRETIMESTEP_NUM
  READ(11, FORMAT_F)  dt
  READ(11, FORMAT_F)  pre_dt
  READ(11, FORMAT_I)  T_SAVE
  READ(11, *)

  READ(11, *)
  READ(11, FORMAT_I) L_MAX
  READ(11, *)

  READ(11, *)
  READ(11, FORMAT_I)  NUM_V
  READ(11, FORMAT_F)  dv
  READ(11, FORMAT_F)  V_GRID_MULT
  READ(11, *)

  READ(11, *)
  READ(11, FORMAT_I) NUM_C
  READ(11, FORMAT_F) dx
  READ(11, FORMAT_F) SMALL_dx
  READ(11, *)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_GRID
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads parameters for periodic initialization
SUBROUTINE INPUT_PERIODIC_INIT(T_AVG, &
                               T_AMP, &
                               T_PHASE, &
                               T_FREQ, &
                               DENS_AVG, &
                               DENS_AMP, &
                               DENS_PHASE, &
                               DENS_FREQ &
                               )

  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: T_AVG, & !< Average temperature
                         T_AMP, & !< Temperature perturbation amplitude
                         T_PHASE !< Temperature perturbation phase

  INTEGER, INTENT(OUT) :: T_FREQ !< Temperture perturbation frequency

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: DENS_AVG, & !< Average density
                         DENS_AMP, & !< Density perturbation amplitude
                         DENS_PHASE !< Density perturbation phase

  INTEGER, INTENT(OUT) :: DENS_FREQ !< Density perturbation frequency

  OPEN(11, FILE = 'INPUT/F_INIT_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  READ(11,*)
  READ(11, FORMAT_F) T_AVG
  READ(11, FORMAT_F) T_AMP
  READ(11, FORMAT_F) T_PHASE
  READ(11, FORMAT_I) T_FREQ
  READ(11, FORMAT_F) DENS_AVG
  READ(11, FORMAT_F) DENS_AMP
  READ(11, FORMAT_F) DENS_PHASE
  READ(11, FORMAT_I) DENS_FREQ
  READ(11,*)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_PERIODIC_INIT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads monotonic initialization parameters
SUBROUTINE INPUT_MONOTONIC_INIT(T_UP, &
                                T_DIV, &
                                DENS_UP, &
                                DENS_DIV, &
                                NUM_DROP, &
                                PLASMA_RAMP_WIDTH &
                                )

  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: T_UP, & !< Upstream temperature
                         T_DIV, & !< Divertor temperature
                         DENS_UP, & !< Upstream density
                         DENS_DIV !< Divertor density

  INTEGER, INTENT(OUT) :: NUM_DROP, & !< X_GRID point at which drop/jump happens
                          PLASMA_RAMP_WIDTH !< Width of exponential ramp in cells
  INTEGER :: I

  OPEN(11, FILE = 'INPUT/F_INIT_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  DO I = 1,10

    READ(11,*)

  END DO


  READ(11,*)
  READ(11, FORMAT_F) T_UP
  READ(11, FORMAT_F) T_DIV
  READ(11, FORMAT_F) DENS_UP
  READ(11, FORMAT_F) DENS_DIV
  READ(11, FORMAT_I) NUM_DROP
  READ(11, FORMAT_I) PLASMA_RAMP_WIDTH
  READ(11,*)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_MONOTONIC_INIT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Reads position of neutral cloud
SUBROUTINE INPUT_NEUTRAL_CLOUD(NUM_CLOUD)

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: NUM_CLOUD

  INTEGER :: I

  OPEN(11, FILE = 'INPUT/F_INIT_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  DO I = 1,18

    READ(11,*)

  END DO


  READ(11,*)
  READ(11, FORMAT_I) NUM_CLOUD
  READ(11,*)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_NEUTRAL_CLOUD
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!>Reads input vector x of size DIM from file INPUT/str_INPUT.txt
SUBROUTINE INPUT_INIT(str,x,DIM)

  IMPLICIT NONE

  CHARACTER(LEN = *), INTENT(IN) :: str
  INTEGER, INTENT(IN) :: DIM
  REAL(KIND=DEFAULT_REAL), INTENT(OUT), DIMENSION(:) :: x
  INTEGER :: I

  OPEN(11, FILE = 'INPUT/'// str //'_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  READ(11,*) (x(I), I = 1, DIM)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_INIT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!>Reads solver parameters for Lis and for nonlinear iterations
SUBROUTINE INPUT_SOLVER_PARAMS(SOLVER_TOL, &
                               SOLVER_MAX_ITER,&
                               MAX_NONLIN, &
                               NONLIN_TOL, &
                               BIS_TOL)


  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) ::  SOLVER_TOL!< Solver tolerance


  INTEGER, INTENT(OUT) :: MAX_NONLIN, & !< Maximum number of iterations
                          SOLVER_MAX_ITER !< Solver max iteration
  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: NONLIN_TOL !< Nonlinear tolerance

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: BIS_TOL

  OPEN(11, FILE = 'INPUT/SOLVER_PARAMS_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  READ(11,*)
  READ(11, FORMAT_F) SOLVER_TOL
  READ(11, FORMAT_I) SOLVER_MAX_ITER
  READ(11, *)

  READ(11,*)
  READ(11, FORMAT_I) MAX_NONLIN
  READ(11, FORMAT_F) NONLIN_TOL
  READ(11,*)

  READ(11,*)
  READ(11, FORMAT_F) BIS_TOL
  READ(11,*)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_SOLVER_PARAMS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Read dipole transition probabilites from file into matrix M
  SUBROUTINE INPUT_DIPOL_TRANS(M,max_n)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: M(:,:)

    INTEGER, INTENT(IN) :: max_n

    CHARACTER (LEN=*), PARAMETER :: FORMAT = '(G14.7E2,47X,I8,8X,I8)'

    REAL(KIND=DEFAULT_REAL) :: VAL

    INTEGER :: I, J, K

    M = 0

    OPEN(11, FILE = 'INPUT/H_DIPOLE_TRANS_PROB_INPUT.txt', ACTION = 'READ', IOSTAT = IR, STATUS = 'OLD')

    DO K = 1, 20 * 19 / 2

      READ(11, FORMAT) VAL, J, I

      IF (I .LE. max_n) M(I,J) = VAL

    END DO

    CLOSE(11,IOSTAT = IR)

  END SUBROUTINE INPUT_DIPOL_TRANS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!>Reads input vector x of size DIM from file INPUT/RESTART/str_INPUT.txt
SUBROUTINE INPUT_RESTART_INIT(str,x,DIM)

  IMPLICIT NONE

  CHARACTER(LEN = *), INTENT(IN) :: str
  CHARACTER(*), PARAMETER :: fmt = "(ES22.14E3)"
  INTEGER, INTENT(IN) :: DIM
  REAL(KIND=DEFAULT_REAL), INTENT(OUT), DIMENSION(:) :: x
  INTEGER :: I

  OPEN(11, FILE = 'INPUT/RESTART/'// str //'_INPUT.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  READ(11,fmt) (x(I), I = 1, DIM)

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_RESTART_INIT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!>Reads time and current timestep from restart info file
SUBROUTINE INPUT_RESTART(TIME,TIMESTEP)

  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: TIME
  INTEGER, INTENT(OUT) :: TIMESTEP

  OPEN(11, FILE = 'INPUT/RESTART/RESTART_INFO.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  READ(11,FORMAT_F) TIME
  READ(11,FORMAT_I) TIMESTEP

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_RESTART
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!>Inputs data needed for adaptive restarts
SUBROUTINE INPUT_ADAPTIVE_RESTART(nn, nL, nv, ndv,old_ff)

  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: ndv
  INTEGER, INTENT(OUT) :: nn, nL, nv

  LOGICAL, INTENT(OUT) :: old_ff

  INTEGER :: I

  OPEN(11, FILE = 'INPUT/RESTART/RESTART_INFO.txt', ACTION = "READ", IOSTAT = IR, STATUS = "OLD")

  DO I = 1, 4

    READ(11,*)

  END DO

  READ(11,FORMAT_I) nn
  READ(11,FORMAT_I) nL
  READ(11, FORMAT_I) nv
  READ(11, FORMAT_F) ndv
  READ(11, FORMAT_L) old_ff

  CLOSE(11, IOSTAT = IR)

END SUBROUTINE INPUT_ADAPTIVE_RESTART
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE INPUT
!-------------------------------------------------------------------------------------------------------------------------------------
