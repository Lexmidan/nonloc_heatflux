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
!> Contains main time loop and calls to all initializations
PROGRAM MAIN

#include "petsc/finclude/petscksp.h"

  USE OUTPUT
  USE POST_PROCESSING
  USE NORMALIZATION
  USE F_INIT
  USE GRID
  USE SWITCHES
  USE NEUT_AND_HEAT
  USE EVOLVE
  USE SPARSE
  USE COLL_CS_RATES
  USE INEL_GRID
  USE PRINTING
  USE SOLVER_PARAMS
  USE MATRIX_DATA
  USE BUILD_COULOMB_COLL
  USE EXT_SOLVER
  USE MPI
  USE petscksp
  USE VAR_KIND_DATA

  IMPLICIT NONE

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: f_old(:), & !< Variable vector from previous timestep
                               f_new(:) !< Current variable vector
  REAL(KIND=DEFAULT_REAL) :: TIME !< Real elapsed time
  INTEGER :: TIMESTEP, & !< Timestep counter
             INIT_TIMESTEP !< Initial timestep if restarting run

  INTEGER :: RANK, IERR, IR

!Initialize PETSc and MPI
  CALL PetscInitialize(PETSC_NULL_CHARACTER,IERR)
  CALL MPI_Comm_rank(PETSC_COMM_WORLD,RANK,IERR)

  IF (RANK .EQ. 0) THEN

    CALL PRINT_START                                                              !Print logo

  END IF

  TIME = 0                                                                      !Initialize time
  INIT_TIMESTEP = 0                                                             !Default initial timestep

  CALL INITIALIZE(f_old,f_new,TIME,INIT_TIMESTEP)                               !Call initializations of all quantities (see below)

  IF (RANK .EQ. 0) THEN

    CALL PRINT_ECHO('-----------------------------------------------------------------------------------------&
        &------------------------')
    CALL PRINT_ECHO('Starting pre-timesteps')
    CALL PRINT_ECHO('-----------------------------------------------------------------------------------------&
        &------------------------')

  END IF

!Set FIRST_PASS module variable in ext_solver to true - build PETSc objects
  CALL SET_FIRST_PASS

  IF (.NOT. RESTART_SWITCH) THEN

    DO TIMESTEP = 1, PRETIMESTEP_NUM                                              !Take short timesteps

      IF (RANK .EQ. 0) CALL PRINT_START_EVOLVE(TIMESTEP,'pre-timestep')
      CALL EVOLVE_F(pre_dt, f_old, f_new,TIMESTEP)                                !Push to next pretimestep

      TIME = TIME + ADAPTIVE_dt * TIME_NORM                                                !Update elapsed time

      IF (RANK .EQ. 0) CALL PRINT_TIMESTEP(TIMESTEP, TIME, 'pre-timestep')
      IF (RANK .EQ. 0) CALL PRINT_ECHO('---------------------------------------------------------&
      &--------------------------------------------------------')

      f_old = f_new                                                               !Update old variable vector
    END DO

    IF (RANK .EQ. 0) CALL PRINT_ECHO('---------------------------------------------------------&
    &--------------------------------------------------------')
    IF (RANK .EQ. 0) CALL PRINT_ECHO('Finished last pre-timestep')

  END IF

  IF (RANK .EQ. 0) CALL PRINT_ECHO('Starting full length timesteps')
  IF (RANK .EQ. 0) CALL PRINT_ECHO('---------------------------------------------------------&
  &--------------------------------------------------------')

  DO TIMESTEP = INIT_TIMESTEP + 1, TIMESTEP_NUM                                                 !Main evolution loop

    IF (RANK .EQ. 0) CALL PRINT_START_EVOLVE(TIMESTEP,'timestep')
    CALL EVOLVE_F(dt, f_old, f_new,TIMESTEP)                                    !Push to next timestep

    TIME = TIME + ADAPTIVE_dt * TIME_NORM                                            !Update elapsed time

    IF (RANK .EQ. 0) CALL PRINT_TIMESTEP(TIMESTEP, TIME, 'timestep')
    IF (RANK .EQ. 0) CALL PRINT_ECHO('---------------------------------------------------------&
    &--------------------------------------------------------')

    IF (RANK .EQ. 0) THEN

      IF (MOD(TIMESTEP,T_SAVE) .EQ. 0) CALL POST_PROC(f_new, TIMESTEP,TIME)          !Call post-processing and output

      IF (SAVE_RESTART_SWITCH) THEN

        CALL PRINT_ECHO('Outputting restart data')
        CALL PRINT_ECHO('-----------------------------------------------------------------------------------------------&
        &------------------')
        CALL OUTPUT_RESTART_F(f_new)
        CALL OUTPUT_RESTART_INFO(TIME,TIMESTEP)

      END IF

    END IF

    f_old = f_new                                                               !Update old variable vector

  END DO

!Deallocate all
  DEALLOCATE(f_new)
  DEALLOCATE(f_old)
  CALL DEALLOCATE_SPARSE(M_MAXWELL)
  CALL DEALLOCATE_SPARSE(M_X_ADV)

!Destroy PETSc objects
  CALL VecDestroy(sol,IERR)
  CALL VecDestroy(rhs,IERR)
  CALL MatDestroy(PETSc_mat,IERR)
  CALL KSPDestroy(solver,IERR)

!Finish run
  IF (RANK .EQ. 0) THEN

    CALL PRINT_ECHO('Run Complete')

    OPEN(11, FILE = 'OUTPUT/STATUS.txt', ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

    WRITE(11,*) 'Normal run finish'

  END IF

  CALL PetscFinalize(IERR)
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calls all initialization routines from linked modules
  SUBROUTINE INITIALIZE(f_old, f_new,TIME,INIT_TIMESTEP)

    REAL(KIND=DEFAULT_REAL), ALLOCATABLE, INTENT(INOUT) :: f_old(:), f_new(:)
    REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: TIME
    INTEGER, INTENT(OUT) :: INIT_TIMESTEP
    INTEGER :: RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('normalization')
    CALL INIT_NORM                                                              !Initialize normalizations

    IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('switches')
    CALL INIT_SWITCHES                                                          !Load switches

    IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('neutral and heating data')
    CALL INIT_HEAT_NEUT                                                         !Load neutral and heating data

    IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('grids')
    CALL INIT_GRID                                                              !Initialize grids

    IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('solver parameters')
    CALL INIT_SOLVER_PARAMS                                                     !Initialize solver parameters for both external solver and internal nonlinear iterations

    IF (RANK .EQ. 0) CALL PRINT_ECHO('Outputting grid vectors')
    IF (RANK .EQ. 0) CALL OUTPUT_GRIDS                                          !Output grid data

    IF (NEUTRAL_TRACK_SWITCH) THEN

      IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('inelastic collision grid mappings')
      CALL INIT_INEL_GRID                                                       !Initialize inelastic collision mappings

    END IF

    IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('global matrix and sparsity patterns')

    CALL INIT_MATRIX_DATA                                                       !Initialize all matrix data, including the global sparse matrix, as well as all sparsity patterns and element markers
    CALL MPI_Barrier(PETSC_COMM_WORLD,IERR)

    IF (NEUTRAL_TRACK_SWITCH) THEN

      IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('neutral collisional and radiative data')
      CALL INIT_SIGMA                                                           !Initializes cross-section data
      CALL INIT_DIPOLE_TRANS_PROB                                               !Load dipole transition probablity data

    END IF

    IF (COLL_EE_L_SWITCH) CALL INIT_ROSENBLUTH                                  !Initialize arrays arising from rosenbluth potential integral I and J (cf. Shkarofsky)

!Allocate variable vectors
    ALLOCATE(f_old(DIM_F))
    ALLOCATE(f_new(DIM_F))

    IF (RESTART_SWITCH) THEN                                                    !Load initial conditions from restart folder

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Restart on - loading initial conditions from restart folder')
      IF (RANK .EQ. 0) THEN

        IF (ADAPTIVE_RESTART_SWITCH) THEN

          CALL INIT_F_ADAPTIVE(f_old)

        ELSE

         CALL INPUT_RESTART_INIT('VAR_VEC',f_old,DIM_F)

        END IF

      END IF

      CALL MPI_Bcast(f_old,DIM_F,MPI_REAL8,0,PETSC_COMM_WORLD,IERR)

      IF (CONTINUE_RUN_SWITCH) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('Continue run on - setting time and timestep to those recorded in restart folder')
        IF (RANK .EQ. 0) CALL INPUT_RESTART(TIME,INIT_TIMESTEP)

        CALL MPI_Bcast(TIME,1,MPI_REAL8,0,PETSC_COMM_WORLD,IERR)
        CALL MPI_Bcast(INIT_TIMESTEP,1,MPI_INTEGER,0,PETSC_COMM_WORLD,IERR)

      ELSE

        IF (RANK .EQ. 0) CALL POST_PROC(f_old,0,0.00D00)                                                     !Perform post-processing and output on initial conditions

      END IF

    ELSE

      IF (RANK .EQ. 0) CALL PRINT_INITIALIZING('variable vector')
      IF (RANK .EQ. 0) CALL INIT_F(f_old)                                                          !Set f_old to initial conditions

      CALL MPI_Bcast(f_old,DIM_F,MPI_REAL8,0,PETSC_COMM_WORLD,IERR)

      IF (RANK .EQ. 0) CALL POST_PROC(f_old,0,0.00D00)                                                     !Perform post-processing and output on initial conditions

    END IF

    IF (RANK .EQ. 0) CALL PRINT_ECHO('Initializations complete')
    IF (RANK .EQ. 0) CALL PRINT_ECHO('--------------------------------------------------------------------------&
         &---------------------------------------')

  END SUBROUTINE INITIALIZE
!-------------------------------------------------------------------------------------------------------------------------------------
END PROGRAM MAIN
!-------------------------------------------------------------------------------------------------------------------------------------
