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
!>Contains main evolution routines and calls to matrix builder
MODULE EVOLVE

  USE GRID
  USE NEUT_AND_HEAT
  USE MOMENTS
  USE SWITCHES
  USE SPARSE
  USE EXT_SOLVER
  USE SOLVER_PARAMS
  USE MAT_BUILD
  USE PRINTING
  USE MATRIX_DATA
  USE COLL_CS_RATES
  USE BUILD_DIV_BOUNDARY
  USE MPI
  USE VAR_KIND_DATA

  REAL(KIND=DEFAULT_REAL) :: ADAPTIVE_dt !< Current adaptive timestep

  INTEGER :: CURRENT_TIMESTEP_NUM !< Number of timesteps performed in current coupling mode

  LOGICAL :: FIRST_STEP !< True if this is first step of simulation
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Performs a timestep of length dt
  SUBROUTINE EVOLVE_F(dt, f_old, f_new, TIMESTEP)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: dt !< Timestep size
    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(INOUT) :: f_old, & !< Old vector
                                                     f_new !< New/unknown vector
    INTEGER, INTENT(IN) :: TIMESTEP !< Current timestep number

    TYPE (SPARSE_MAT) :: M_SOLVER !< Solver matrix

    INTEGER :: I, J, &
               N_NONLIN!< Number of nonlinear iterations

    REAL(KIND=DEFAULT_REAL) :: DELTA_F
    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X) :: n_old, &      !< Old densities
                                            u_old, &      !< Old electron velocity
                                            T_old, &      !< Old temperatures
                                            n_lagged, &   !< Lagged densities
                                            T_lagged , &    !< Lagged temperature
                                            u_lagged !< Lagged electron velocity

    REAL(KIND=DEFAULT_REAL) :: f_lagged(DIM_F), &  !< Lagged vector
                    n_i_lagged(NUM_X), & !< Lagged ion densities
                    u_i_lagged(NUM_X), & !< Lagged ion velocities
                    n_i_old(NUM_X), & !< Old ion densities
                    u_i_old(NUM_X), & !< Old ion velocities
                    n_TOT(NUM_X), & !< Total old density vector
                    timestep_mult(NUM_X)!< Timestep density multiplier vector

    REAL(KIND=DEFAULT_REAL) :: n_new(NUM_X) !< New electron density

    INTEGER :: RANK, IERR, SIZE

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    IF (RANK .EQ. 0) CALL PRINT_ECHO('Initializing old density and temperature...')

    IF (FULL_FLUID_MODE) THEN

      DO I = 1, NUM_X

        n_old(I) = f_old(X_POS(I) + NUM_0D - 4)
        u_old(I) = f_old(X_POS(I) + NUM_0D - 3)
        T_old(I) = f_old(X_POS(I) + NUM_0D - 2)

      END DO

    ELSE

      u_old = 0

      DO I = 1, NUM_X                                                        !Calculate old density

        n_old(I) = DENSITY_MOMENT(f_old(X_POS(I) + (NUM_H - 1) * NUM_V + 1 : X_POS(I) + NUM_H * NUM_V))

      END DO

      IF (L_MAX .GT. 0) THEN                                                      !Calculate old temperature

          DO I = 1, NUM_X

            T_old(I) = TEMPERATURE_MOMENT(f_0 = f_old(X_POS(I) + (NUM_H - 1) * NUM_V + 1 : X_POS(I) + NUM_H * NUM_V),&
                                        f_10 = f_old(X_POS(I) + H_POS(1)* NUM_V + 1 : X_POS(I) + (1 + H_POS(1))* NUM_V) &
                                        )

          END DO

      ELSE

        DO I = 1, NUM_X

          T_old(I) = TEMPERATURE_MOMENT(f_0 = f_old(X_POS(I) + (NUM_H - 1) * NUM_V + 1 : X_POS(I) + NUM_H * NUM_V))

        END DO

      END IF

    END IF

!Update dexcitation and recombination cross sections
      IF ((COLL_EN_EX .OR. COLL_RECOMB) .AND. (FAST_DETAILED_BALANCE_SWITCH)) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('Updating detailed balance cross-sections')

        IF (COLL_EN_EX) CALL UPDATE_SIGMA_L_DEEX(T_old)

        IF (COLL_RECOMB)   CALL UPDATE_SIGMA_L_RECOMB(T_old)

      END IF

!Modify timestep if collisionally adaptive

    IF (ADAPTIVE_TIMESTEP_SWITCH) THEN

      n_TOT = 0.00D00

      DO I = 1, NUM_X

        n_TOT(I) = n_TOT(I) + n_old(I)

        DO J = 1, NUM_NEUTRALS

          n_TOT(I) = n_TOT(I) + f_old(X_POS(I) + NUM_H * NUM_V + J)

        END DO

        timestep_mult(I) = SQRT(T_old(I))**3/n_TOT(I)

      END DO

      ADAPTIVE_dt = dt * MINVAL(timestep_mult)

    ELSE

      ADAPTIVE_dt = dt

    END IF

!Extract old ion density and flow velocity

    IF (COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) THEN

      DO I = 1, NUM_X

        n_i_old(I) = f_old(X_POS(I) + NUM_0D - 1)
        u_i_old(I) = f_old(X_POS(I) + NUM_0D)

      END DO

      IF (ION_CONT_OFF_SWITCH) n_i_old = n_old/Z_PROF

    ELSE

      n_i_old = 0
      u_i_old = 0

    END IF

    DELTA_F = 1.0D00
    N_NONLIN = 0

    IF (RANK .EQ. 0) CALL PRINT_ECHO('Starting nonlinear loop')
    IF (RANK .EQ. 0) CALL PRINT_ECHO('---------------------------------------------------------------&
    &--------------------------------------------------')

    DO WHILE ((DELTA_F .GT. NONLIN_TOL) .AND. (N_NONLIN .LT. MAX_NONLIN))       !Begin nonlinear iterations

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Starting nonlinear iteration')
      IF (RANK .EQ. 0) CALL PRINT_ECHO('Calculating lagged quantities...')

!Calculate lagged quantities (variable vector, density, temperature)
      IF (N_NONLIN .GT. 0) THEN

        f_lagged = f_new

        IF (FULL_FLUID_MODE) THEN

          DO I = 1, NUM_X

            n_lagged(I) = f_lagged(X_POS(I) + NUM_0D - 4)
            u_lagged(I) = f_lagged(X_POS(I) + NUM_0D - 3)
            T_lagged(I) = f_lagged(X_POS(I) + NUM_0D - 2)

          END DO

        ELSE

          u_lagged = 0

          DO I = 1, NUM_X

            n_lagged(I) = DENSITY_MOMENT(f_lagged(X_POS(I) + (NUM_H - 1) * NUM_V + 1 : X_POS(I) + NUM_H * NUM_V))

          END DO

          IF (L_MAX .GT. 0) THEN

            DO I = 1, NUM_X

                T_lagged(I) = TEMPERATURE_MOMENT(f_0 = f_lagged(X_POS(I) + (NUM_H - 1) * NUM_V + 1 : X_POS(I) + NUM_H * NUM_V),&
                                        f_10 = f_lagged(X_POS(I) + H_POS(1)* NUM_V + 1 : X_POS(I) + (1 + H_POS(1))* NUM_V) &
                                        )

            END DO

          ELSE

            DO I = 1, NUM_X

              T_lagged(I) = TEMPERATURE_MOMENT(f_0 = f_lagged(X_POS(I) + (NUM_H - 1) * NUM_V + 1 : X_POS(I) + NUM_H * NUM_V))

            END DO

          END IF

        END IF

        IF (COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) THEN

          DO I = 1, NUM_X

            n_i_lagged(I) = f_lagged(X_POS(I) + NUM_0D - 1)
            u_i_lagged(I) = f_lagged(X_POS(I) + NUM_0D)

          END DO

          IF (ION_CONT_OFF_SWITCH) n_i_lagged = n_lagged/Z_PROF


        ELSE

          n_i_lagged = 0
          u_i_lagged = 0

        END IF

      ELSE

        f_lagged = f_old
        n_lagged = n_old
        T_lagged = T_old
        u_lagged = u_old
        n_i_lagged = n_i_old
        u_i_lagged = u_i_old

      END IF

!Update dexcitation and recombination cross sections
      IF ((COLL_EN_EX .OR. COLL_RECOMB) .AND. (.NOT. FAST_DETAILED_BALANCE_SWITCH)) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('Updating detailed balance cross-sections')

        IF (COLL_EN_EX) CALL UPDATE_SIGMA_L_DEEX(T_lagged)

        IF (COLL_RECOMB)   CALL UPDATE_SIGMA_L_RECOMB(T_lagged)

      END IF

      IF (COLL_EN_EX .OR. COLL_RECOMB .OR. COLL_EN_ION) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('Updating collision rates')
        IF (COLL_EN_EX) CALL UPDATE_EX_DEEX_RATES(f_lagged)
        IF (COLL_RECOMB) CALL UPDATE_TB_RECOMB_RATES(f_lagged)
        IF (COLL_EN_ION) CALL UPDATE_ION_RATES(f_lagged)

      END IF

!Allocate boundary distribution vector
      IF (PLASMA_SINK_SWITCH) CALL ALLOCATE_F_BOUNDARY

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Starting matrix construction')

      CALL FILL_M(f_lagged, n_old, T_old, n_lagged, T_lagged, n_i_lagged, u_i_lagged,u_lagged,&
                   TIMESTEP, N_NONLIN) !Call matrix filling routines

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Preparing sparse matrix for solver')

      CALL MAT_PREP(LOCAL_M,M_SOLVER,ADAPTIVE_dt)                                         !Prepare matrix for solver (M_SOLVER = I - dt*M)

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Flushing previous entries')

      CALL MAT_FLUSH                                                            !Set all changing matrix elements to 0

!Wait for all processes to prepare their respective matrix band
      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Calling solver')
      f_new = f_lagged                                                          !Set unknown vector to lagged vector in case solver set to not assume zeroes
      CALL SOLVE_M(M_SOLVER, f_old, f_new,TIMESTEP)                                      !Call external solver

      IF (FULL_FLUID_MODE) THEN

        DO I = 1, NUM_X

          n_new(I) = f_new(X_POS(I) + NUM_0D - 4)

        END DO

      ELSE

        DO I = 1, NUM_X

          n_new(I) = DENSITY_MOMENT(f_new(X_POS(I) + (NUM_H - 1) * NUM_V + 1 : X_POS(I) + NUM_H * NUM_V))

        END DO

      END IF


      DO I = 1, NUM_X - 1

        IF (MOD(I,2) .EQ. 0) THEN

          n_new(I) = ((X_GRID(I) - X_GRID(I - 1)) * n_new(I+1) &
          + (X_GRID(I + 1) - X_GRID(I)) * n_new(I-1))/dxc(I)

        END IF

      END DO

      IF (PERIODIC_BOUNDARY_SWITCH) n_new(NUM_X) = 0.5D00 * (n_new(NUM_X-1)+n_new(1))

      IF ((COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) .AND. ION_CONT_OFF_SWITCH) THEN

        DO I = 1, NUM_X

          f_new(X_POS(I) + NUM_0D - 1) = n_new(I)/Z_PROF(I)

        END DO

      END IF

      N_NONLIN = N_NONLIN + 1                                                   !Advance nonlinear iteration counter

!Broadcasting boundary vector and temperature
      IF (PLASMA_SINK_SWITCH .AND. (.NOT. FULL_FLUID_MODE)) THEN

        DO I = 1, NUM_H

          CALL MPI_Bcast(f_boundary(I,:),NUM_V,MPI_REAL8,SIZE-1,PETSC_COMM_WORLD,IERR)

        END DO

        IF (SONIC_OUTFLOW_DIV_SWITCH .OR. (.NOT.(COLD_ION_FLUID_SWITCH))) THEN

          CALL MPI_Bcast(gamma_e,1,MPI_REAL8,SIZE-1,PETSC_COMM_WORLD,IERR)
          CALL MPI_Bcast(pot_drop,1,MPI_REAL8,SIZE-1,PETSC_COMM_WORLD,IERR)

        END IF

      END IF

      IF (RANK .EQ. 0) THEN

        CALL PRINT_ECHO('Interpolating vector')
        CALL INTERPOLATE(f_new)                                                   !Interpolate quantities not evolved
        CALL PRINT_ECHO('Interpolation complete')

        IF (FULL_FLUID_MODE) CALL UPDATE_MAXWELLIAN(f_new)   !Update Maxwellian components

        CALL PRINT_ECHO('Calculating nonlinear residual')
        DELTA_F = NONLIN_ERR(f_new,f_lagged)                                      !Calculate nonlinear residual

        CALL PRINT_END_NONLIN(N_NONLIN,DELTA_F)
        CALL PRINT_ECHO('-----------------------------------------------------------------------------------------------&
        &------------------------------------')

      END IF

      CALL MPI_Bcast(DELTA_F,1,MPI_REAL8,0,PETSC_COMM_WORLD,IERR)
      CALL MPI_Bcast(f_new,DIM_F,MPI_REAL8,0,PETSC_COMM_WORLD,IERR)

      CALL DEALLOCATE_SPARSE(M_SOLVER)

    END DO

    IF (RANK .EQ. 0) CALL PRINT_END_EVOLVE(N_NONLIN, DELTA_F)
    IF (RANK .EQ. 0) CALL PRINT_ECHO('--------------------------------------------------&
    &---------------------------------------------------------------')

  END SUBROUTINE EVOLVE_F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Interpolates those quantities not evolved
  SUBROUTINE INTERPOLATE(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(INOUT) :: f(DIM_F)

    INTEGER :: I, J, K

    REAL(KIND=DEFAULT_REAL) :: a_minus, a_plus

    DO I = 2, NUM_X - 1

      a_minus = (X_GRID(I + 1) - X_GRID(I)) / dxc(I)
      a_plus = (X_GRID(I) - X_GRID(I - 1)) / dxc(I)

      IF (.NOT. FULL_FLUID_MODE) THEN

!Interpolate distribution functions
        DO J = 1, NUM_H

          IF (((MOD(I,2) .EQ. 0) .AND. (MOD(J - 1,2) .EQ. 0)) .OR. ((MOD(I,2) .EQ. 1) .AND. (MOD(J - 1,2) .EQ. 1))) THEN

            DO K = 1, NUM_V

              f(X_POS(I) + NUM_V*(NUM_H - J) + K) = &
                    a_minus * f(X_POS(I-1) + NUM_V*(NUM_H - J) + K) + a_plus * f(X_POS(I+1) + NUM_V*(NUM_H - J) + K)

            END DO

          END IF

        END DO

      END IF

      a_minus = (X_GRID(I + 1) - X_GRID(I)) / dxc(I)
      a_plus = (X_GRID(I) - X_GRID(I - 1)) / dxc(I)

      IF (MOD(I,2) .EQ. 1) THEN

!Interpolate fields
        IF (MAXWELL_SWITCH) THEN

          f(X_POS(I) + NUM_H*NUM_V + NUM_NEUTRALS + 1 ) = &
                a_minus * f(X_POS(I-1) + NUM_H*NUM_V + NUM_NEUTRALS + 1)  + a_plus * f(X_POS(I+1) + NUM_H*NUM_V + NUM_NEUTRALS + 1)

        END IF

!Interpolate ion flow

        IF (COLD_ION_FLUID_SWITCH) THEN

          f(X_POS(I) + NUM_0D) = &
                a_minus * f(X_POS(I-1) + NUM_0D)  + a_plus * f(X_POS(I+1) + NUM_0D)

        END IF

!Interpolate electron flow
        IF (FULL_FLUID_MODE) THEN

          f(X_POS(I) + NUM_0D-3) = &
                a_minus * f(X_POS(I-1) + NUM_0D - 3)  + a_plus * f(X_POS(I+1) + NUM_0D - 3)

        END IF

      ELSE

!Interpolate neutral density
        IF (NEUTRAL_TRACK_SWITCH) THEN

          DO J = 1, NUM_NEUTRALS

            f(X_POS(I) + NUM_H * NUM_V + J) = &
                a_minus* f(X_POS(I - 1) + NUM_H * NUM_V + J) + a_plus * f(X_POS(I + 1) + NUM_H * NUM_V + J)

          END DO

        END IF

!Interpolate ion density

        IF (COLD_ION_FLUID_SWITCH) THEN

          f(X_POS(I) + NUM_0D - 1) = &
                  a_minus * f(X_POS(I-1) + NUM_0D - 1)  + a_plus * f(X_POS(I+1) + NUM_0D - 1)

        END IF

!Interpolate electron density and temperature
        IF (FULL_FLUID_MODE) THEN

          f(X_POS(I) + NUM_0D - 4) = &
                a_minus * f(X_POS(I-1) + NUM_0D - 4)  + a_plus * f(X_POS(I+1) + NUM_0D - 4)

          f(X_POS(I) + NUM_0D - 2) = &
                a_minus * f(X_POS(I-1) + NUM_0D - 2)  + a_plus * f(X_POS(I+1) + NUM_0D - 2)

        END IF

      END IF

    END DO



!Interpolate quantities in first and last cell if periodic boundary
    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      IF (.NOT. FULL_FLUID_MODE) THEN

        DO J = 1, NUM_H

          IF (MOD(J - 1,2) .EQ. 1) THEN

            DO K = 1, NUM_V

              f(NUM_V * (NUM_H - J) + K) = &
                                          0.5D00 * (f(X_POS(2) + NUM_V * (NUM_H - J) + K) &
                                          + f(X_POS(NUM_X) + NUM_V * (NUM_H - J) + K))

             END DO

          ELSE

            DO K = 1, NUM_V

              f(X_POS(NUM_X) + NUM_V * (NUM_H - J) + K) = &
                                          0.5D00 * (f(NUM_V * (NUM_H - J) + K) + f(X_POS(NUM_X - 1) + NUM_V * (NUM_H - J) + K))

            END DO

          END IF

        END DO

      END IF

      IF (MAXWELL_SWITCH) THEN

        f(NUM_H*NUM_V + NUM_NEUTRALS + 1) = &
                                    0.5D00 * (f(X_POS(2) + NUM_H*NUM_V + NUM_NEUTRALS + 1) &
                                     + f(X_POS(NUM_X) + NUM_H*NUM_V + NUM_NEUTRALS + 1))

      END IF

      IF (COLD_ION_FLUID_SWITCH) THEN

        f(X_POS(NUM_X) + NUM_0D - 1) = &
              0.5D00 * (f(X_POS(NUM_X-1) + NUM_0D - 1)  +  f(NUM_0D - 1))

        f(NUM_0D) = &
              0.5D00 * (f(X_POS(NUM_X) + NUM_0D)  +  f(X_POS(2) + NUM_0D))

      END IF

      IF (NEUTRAL_TRACK_SWITCH) THEN

        DO J = 1, NUM_NEUTRALS

          f(X_POS(NUM_X) + NUM_H * NUM_V + J) = 0.5D00 * (f(X_POS(1) + NUM_H * NUM_V + J) &
          + f(X_POS(NUM_X - 1) + NUM_H * NUM_V + J))

        END DO

      END IF

      IF (FULL_FLUID_MODE) THEN

        f(X_POS(NUM_X) + NUM_0D - 4) = &
              0.5D00 * (f(X_POS(NUM_X-1) + NUM_0D - 4)  +  f(NUM_0D - 4))

        f(X_POS(NUM_X) + NUM_0D - 2) = &
              0.5D00 * (f(X_POS(NUM_X-1) + NUM_0D - 2)  +  f(NUM_0D - 2))

        f(NUM_0D - 3) = &
              0.5D00 * (f(X_POS(NUM_X) + NUM_0D - 3)  +  f(X_POS(2) + NUM_0D - 3))

      END IF

    END IF

!Set interpolated quantities if no flow boundaries
    IF (NO_FLOW_BOUNDARY_UP_SWITCH) THEN

      IF (.NOT. FULL_FLUID_MODE) THEN

        DO J = 1, NUM_H

          IF (MOD(J - 1,2) .EQ. 1) THEN

            DO K = 1, NUM_V

              f(NUM_V*(NUM_H - J) + K) = 0.5D00 * f(X_POS(2) + NUM_V*(NUM_H - J) + K)

            END DO

          END IF

        END DO

      END IF

      IF (MAXWELL_SWITCH) THEN

        f(NUM_H*NUM_V + NUM_NEUTRALS + 1) = 0.5D00 * f(X_POS(2) + NUM_H*NUM_V + NUM_NEUTRALS + 1)

      END IF

      IF (COLD_ION_FLUID_SWITCH) THEN

        f(NUM_0D) = &
              0.5D00 * f(X_POS(2) + NUM_0D)

      END IF

      IF (FULL_FLUID_MODE) THEN

        f(NUM_0D-3) = &
              0.5D00 * f(X_POS(2) + NUM_0D-3)

      END IF

    END IF

    IF (NO_FLOW_BOUNDARY_DIV_SWITCH) THEN

      IF (PLASMA_SINK_SWITCH .OR. RECYCLING_SWITCH) THEN

        IF (.NOT. FULL_FLUID_MODE) THEN

          DO J = 1, NUM_H

            IF (MOD(J - 1,2) .EQ. 1) THEN

                DO K = 1, NUM_V

                  f(X_POS(NUM_X) + NUM_V*(NUM_H - J) + K) = 0.50D00*(f(X_POS(NUM_X - 1) + NUM_V*(NUM_H - J) + K) &
                  +  f_boundary(J,K))

                END DO

            END IF

          END DO

        END IF

        IF (MAXWELL_SWITCH) THEN

          f(X_POS(NUM_X) + NUM_H*NUM_V + NUM_NEUTRALS + 1) =  f(X_POS(NUM_X - 1) + NUM_H*NUM_V + NUM_NEUTRALS + 1) + &
                  (f(X_POS(NUM_X - 1) + NUM_H*NUM_V + NUM_NEUTRALS + 1) - f(X_POS(NUM_X - 3) + NUM_H*NUM_V + NUM_NEUTRALS + 1))/&
                  (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)) * (X_GRID(NUM_X) - X_GRID(NUM_X - 1))

          IF (f(X_POS(NUM_X) + NUM_H*NUM_V + NUM_NEUTRALS + 1) .LT. 0) f(X_POS(NUM_X) + NUM_H*NUM_V + NUM_NEUTRALS + 1) = 0.0D00

        END IF

        IF (COLD_ION_FLUID_SWITCH) THEN

          IF (ION_EL_TEMP_SWITCH) THEN

            IF (BOHM_VALUE) THEN

              f(X_POS(NUM_X) + NUM_0D) = 0.5D00*(f(X_POS(NUM_X - 1) + NUM_0D) &
               + MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(ION_MASS)))

            ELSE

              f(X_POS(NUM_X) + NUM_0D) = FLOW_LIMITER/n_e_boundary&
                             *(f(X_POS(NUM_X - 1) + NUM_0D-1)*f(X_POS(NUM_X - 1) + NUM_0D) + &
                             (f(X_POS(NUM_X - 1) + NUM_0D-1)*f(X_POS(NUM_X - 1) + NUM_0D) -&
                             f(X_POS(NUM_X - 3) + NUM_0D-1)*f(X_POS(NUM_X - 3) + NUM_0D))/&
                             (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)) * (X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

            END IF

          ELSE

            IF (BOHM_VALUE) THEN

              f(X_POS(NUM_X) + NUM_0D) = 0.5D00*(f(X_POS(NUM_X - 1) + NUM_0D) &
              + MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(2*ION_MASS)))

            ELSE
              f(X_POS(NUM_X) + NUM_0D) = FLOW_LIMITER/n_e_boundary&
                             *(f(X_POS(NUM_X - 1) + NUM_0D-1)*f(X_POS(NUM_X - 1) + NUM_0D) + &
                            (f(X_POS(NUM_X - 1) + NUM_0D-1)*f(X_POS(NUM_X - 1) + NUM_0D) -&
                            f(X_POS(NUM_X - 3) + NUM_0D-1)*f(X_POS(NUM_X - 3) + NUM_0D))/&
                            (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)) * (X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

            END IF

          END IF

        END IF


        IF (FULL_FLUID_MODE) THEN

          IF (ION_EL_TEMP_SWITCH) THEN

            IF (BOHM_VALUE) THEN

              f(X_POS(NUM_X) + NUM_0D-3) = 0.5D00*(f(X_POS(NUM_X - 1) + NUM_0D-3) &
              + MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(ION_MASS)))

            ELSE

              f(X_POS(NUM_X) + NUM_0D-3) = f(X_POS(NUM_X) + NUM_0D)

            END IF

          ELSE

            IF (BOHM_VALUE) THEN

              f(X_POS(NUM_X) + NUM_0D-3) = 0.5D00*(f(X_POS(NUM_X - 1) + NUM_0D-3) &
              + MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(2*ION_MASS)))

            ELSE

              f(X_POS(NUM_X) + NUM_0D-3) = f(X_POS(NUM_X) + NUM_0D)

            END IF

          END IF

        END IF

      ELSE

        IF (.NOT. FULL_FLUID_MODE) THEN

          DO J = 1, NUM_H

            IF (MOD(J - 1,2) .EQ. 1) THEN

              DO K = 1, NUM_V

                f(X_POS(NUM_X) + NUM_V*(NUM_H - J) + K) =  0.5D00 * f(X_POS(NUM_X - 1) + NUM_V*(NUM_H - J) + K)

              END DO

            END IF

          END DO

        END IF

        IF (MAXWELL_SWITCH) THEN

          f(X_POS(NUM_X) + NUM_H*NUM_V + NUM_NEUTRALS + 1) =  0.5D00 * f(X_POS(NUM_X - 1) + NUM_H*NUM_V + NUM_NEUTRALS + 1)

        END IF

        IF (COLD_ION_FLUID_SWITCH) THEN

          f(X_POS(NUM_X) + NUM_0D) = 0.5D00 * f(X_POS(NUM_X - 1) + NUM_0D)

        END IF

        IF (FULL_FLUID_MODE) THEN

          f(X_POS(NUM_X) + NUM_0D-3) = 0.5D00 * f(X_POS(NUM_X - 1) + NUM_0D-3)

        END IF

      END IF

    END IF

  END SUBROUTINE INTERPOLATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates nonlinear iteration error between f_new and f_lagged
  REAL(KIND=DEFAULT_REAL) FUNCTION NONLIN_ERR(f_new, f_lagged)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(IN) :: f_new, &
                                                f_lagged

    INTEGER :: I, J, K
    REAL(KIND=DEFAULT_REAL) :: DDF(NUM_X), &
              er

    DDF = 0

    DO I = 1, NUM_X

      er = 0

!Skip harmonics and electric field if not kinetic timestep
      IF (.NOT. FULL_FLUID_MODE) THEN

!Add residual from distribution functions
        DO J = 1, NUM_H

          IF (((MOD(J - 1,2) .EQ. 0) .AND. (MOD(I,2) .EQ. 1)) .OR. ((MOD(J - 1,2) .EQ. 1) .AND. (MOD(I,2) .EQ. 0))) THEN

            DO K = 1, NUM_V

              er = (f_new(X_POS(I) + (NUM_H - J) * NUM_V + K) - f_lagged(X_POS(I) + (NUM_H - J) * NUM_V + K)) ** 2

                DDF(I) = DDF(I) + er

            END DO

          END IF

        END DO

!Add residual from fields
        IF (MAXWELL_SWITCH) THEN

          IF (MOD(I,2) .EQ. 0) THEN

            er =  (f_new(X_POS(I) + NUM_V*NUM_H + NUM_NEUTRALS + 1) &
                - f_lagged(X_POS(I) + NUM_V*NUM_H + NUM_NEUTRALS + 1)) ** 2

            DDF(I) = DDF(I) + er

          END IF

        END IF


      END IF

!Add residual from neutral densities
      IF (NEUTRAL_TRACK_SWITCH) THEN

        IF (MOD(I,2) .EQ. 1) THEN

          DO J = 1, NUM_NEUTRALS

            er = (f_new(X_POS(I) + NUM_H * NUM_V + J) - f_lagged(X_POS(I) + NUM_H * NUM_V + J)) ** 2


            DDF(I) = DDF(I) + er

          END DO

        END IF

      END IF

!Add residual from ion quantities
      IF (COLD_ION_FLUID_SWITCH) THEN

        IF (MOD(I,2) .EQ. 1) THEN

          er = (f_new(X_POS(I) + NUM_0D - 1) - f_lagged(X_POS(I) + NUM_0D - 1)) ** 2

          DDF(I) = DDF(I) + er

        END IF

        IF (MOD(I,2) .EQ. 0) THEN

          er = (f_new(X_POS(I) + NUM_0D) - f_lagged(X_POS(I) + NUM_0D)) ** 2

          DDF(I) = DDF(I) + er

        END IF

      END IF

      IF (FULL_FLUID_MODE) THEN

        IF (MOD(I,2) .EQ. 1) THEN

          er = (f_new(X_POS(I) + NUM_0D - 4) - f_lagged(X_POS(I) + NUM_0D - 4)) ** 2

          DDF(I) = DDF(I) + er

          er = (f_new(X_POS(I) + NUM_0D - 2) - f_lagged(X_POS(I) + NUM_0D - 2)) ** 2

          DDF(I) = DDF(I) + er

        END IF

        IF (MOD(I,2) .EQ. 0) THEN

          er = (f_new(X_POS(I) + NUM_0D - 3) - f_lagged(X_POS(I) + NUM_0D - 3)) ** 2

          DDF(I) = DDF(I) + er

        END IF

      END IF

    END DO

    DDF = SQRT(DDF)

    NONLIN_ERR = MAXVAL(DDF(:))

  END FUNCTION NONLIN_ERR
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update first two harmonics to Maxwellian values for full fluid mode
  SUBROUTINE UPDATE_MAXWELLIAN(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), DIMENSION(DIM_F), INTENT(INOUT) :: f

    INTEGER :: I,J

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X) :: n, T, u, n_num

    REAL(KIND=DEFAULT_REAL) :: f0(NUM_X,NUM_V)

    n_num = 0

    DO I = 1, NUM_X

      n(I) = f(X_POS(I) + NUM_0D - 4)
      u(I) = f(X_POS(I) + NUM_0D - 3)
      T(I) = f(X_POS(I) + NUM_0D - 2)

    END DO

    DO I = 1, NUM_X

      DO J = 1, NUM_V

        f0(I,J) = n(I) * ((PI * T(I)) ** (- 3.00D00/2.00D00)) * EXP( - (V_GRID(J) ** 2)/ T(I))    !Maxwellian initialization
        n_num(I) = n_num(I) + 4.00D00 * PI * V_GRID(J) ** 2 * V_GRID_WIDTH(J) * f0(I,J)

      END DO

      DO J = 1, NUM_V

        f(X_POS(I) + NUM_V * (NUM_H - 1) + J) = f0(I,J) * n(I)/n_num(I)
        f(X_POS(I) + NUM_V * (NUM_H - 2) + J) = 2.00D00 * u(I) * V_GRID(J) * f(X_POS(I) + NUM_V * (NUM_H - 1) + J) / T(I)

      END DO

    END DO

  END SUBROUTINE UPDATE_MAXWELLIAN
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE EVOLVE
!-------------------------------------------------------------------------------------------------------------------------------------
