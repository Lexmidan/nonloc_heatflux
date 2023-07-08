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
!> Contains all vector initialization
MODULE F_INIT

  USE GRID
  USE SWITCHES
  USE INPUT
  USE NEUT_AND_HEAT
  USE NORMALIZATION
  USE MPI
  USE VAR_KIND_DATA

! Periodic initialization parameters
  REAL(KIND=DEFAULT_REAL) :: T_AVG, & !< Average temperature
            T_AMP, & !< Temperature perturbation amplitude
            T_PHASE !< Temperature perturbation phase

  INTEGER :: T_FREQ !< Temperture perturbation frequency

  REAL(KIND=DEFAULT_REAL) :: DENS_AVG, & !< Average density
            DENS_AMP, & !< Density perturbation amplitude
            DENS_PHASE !< Density perturbation phase

  INTEGER :: DENS_FREQ !< Density perturbation frequency

! Monotonic initialization parameters
  REAL(KIND=DEFAULT_REAL) :: T_UP, & !< Upstream temperature
            T_DIV, & !< Divertor temperature
            DENS_UP, & !< Upstream density
            DENS_DIV !< Divertor density

  INTEGER :: NUM_DROP, & !< X_GRID point at which drop/jump happens
             PLASMA_RAMP_WIDTH !< Width of plasma exponential ramp in cells

! Neutral initialization
  INTEGER :: NUM_CLOUD !< X_GRID point at which neutral cloud starts
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes f vector
  SUBROUTINE INIT_F(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: f(DIM_F) !< Vector to be initialized
    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X) :: T, & !< Temperature vector
                              n, & !< Density vector
                              n_n, & !< Neutral density vector
                              u_i !< Ion velocity
    REAL(KIND=DEFAULT_REAL) :: f0(NUM_X,NUM_V), &
                    E(NUM_X), &
                    A(NUM_X), &
                    f1, &
                    GRAD_F_X(NUM_X,NUM_V), &
                    GRAD_F_V(NUM_X,NUM_V), &
                    GRAD_T(NUM_X), &
                    GRAD_n(NUM_X), &
                    n_TOT(NUM_X), &
                    ION_SUM, &
                    EX_SUM,&
                    b, &
                    n_num(NUM_X)

    REAL(KIND=DEFAULT_REAL) , ALLOCATABLE, DIMENSION(:,:) :: n_neut

    INTEGER:: I, &
              J, &
              F_0_POS, &
              f_1_POS

    f = 0

    IF (PERIODIC_BOUNDARY_SWITCH) THEN                                          !Periodic initialization

      CALL INPUT_PERIODIC_INIT(T_AVG, &
                               T_AMP, &
                               T_PHASE, &
                               T_FREQ, &
                               DENS_AVG, &
                               DENS_AMP, &
                               DENS_PHASE, &
                               DENS_FREQ &
                               )

      DO I = 1, NUM_X

        T(I) = T_AVG + T_AMP * SIN(2.00D00 * PI * T_FREQ * X_GRID(I)/(NUM_C * dx) + T_PHASE)
        n_TOT(I) = DENS_AVG + DENS_AMP * SIN(2.00D00 * PI * DENS_FREQ * X_GRID(I)/(NUM_C * dx) + DENS_PHASE)

      END DO

    ELSE IF (LINEAR_INIT .OR. DROP_INIT .OR. TWO_POINT_M_INIT) THEN             !Monotonic initialization

      CALL INPUT_MONOTONIC_INIT(T_UP, &
                                T_DIV, &
                                DENS_UP, &
                                DENS_DIV, &
                                NUM_DROP, &
                                PLASMA_RAMP_WIDTH &
                                )

      IF (LINEAR_INIT) THEN                                                     !Linear initialization

        DO I = 1, NUM_X

          T(I) = T_UP - (T_UP - T_DIV) * X_GRID(I)/X_GRID(NUM_X)
          n_TOT(I) = DENS_UP - (DENS_UP - DENS_DIV) * X_GRID(I)/X_GRID(NUM_X)

        END DO

      ELSE IF (DROP_INIT) THEN                                                  !Drop/jump initialization

        DO I = 1, NUM_DROP

            T(I) = T_UP
            n_TOT(I) = DENS_UP

        END DO

        IF (PLASMA_RAMP_WIDTH .GT. 0) THEN

          DO I = NUM_DROP + 1, NUM_DROP + 2 * PLASMA_RAMP_WIDTH

            T(I) = T_UP * EXP(-REAL(I - NUM_DROP,KIND=DEFAULT_REAL)&
            /REAL(2 * PLASMA_RAMP_WIDTH,KIND=DEFAULT_REAL) * LOG(T_UP/T_DIV))
            n_TOT(I) = DENS_UP * EXP(-REAL(I - NUM_DROP,KIND=DEFAULT_REAL)&
            /REAL(2 * PLASMA_RAMP_WIDTH,KIND=DEFAULT_REAL) * LOG(DENS_UP/DENS_DIV))

          END DO

        END IF

        DO I = NUM_DROP + 2* PLASMA_RAMP_WIDTH + 1, NUM_X

          T(I) = T_DIV
          n_TOT(I) = DENS_DIV

        END DO

      ELSE IF (TWO_POINT_M_INIT) THEN                                           !Two-Point model initialization

        DO I = 1, NUM_X

          T(I) = (T_UP ** (7.00D00/2.00D00) + X_GRID(I)/X_GRID(NUM_X) * &
          (T_DIV ** (7.00D00/2.00D00) - T_UP ** (7.00D00/2.00D00))) ** (2.00D00/7.00D00)
          n_TOT(I) = DENS_UP * T_UP/T(I)

        END DO

      END IF

    END IF

    IF (DENSITY_FROM_FILE_INIT) CALL INPUT_INIT('DENS', n_TOT, NUM_X)                 !Initialize from file
    IF (TEMPERATURE_FROM_FILE_INIT) CALL INPUT_INIT('TEMPERATURE', T, NUM_X)

!Local Saha-Boltzmann initialization
    IF (NEUTRAL_TRACK_SWITCH .AND. LOCAL_SAHA_BOLTZMANN_INIT_SWITCH) THEN

      ALLOCATE(n_neut(NUM_X, NUM_NEUTRALS))

      n_neut = 0

      DO J = 1, NUM_X

        ION_SUM = 0
        EX_SUM = 0

!Calculate ionization and excitation sums
        DO I = 1, NUM_NEUTRALS

          ION_SUM = ION_SUM + I ** 2 * EXP((ION_POT_H/REAL(I**2,KIND=DEFAULT_REAL))/T(J))
          EX_SUM = EX_SUM + I ** 2 * EXP(-(ION_POT_H * (1.00D00 - 1.00D00 / REAL(I ** 2,KIND=DEFAULT_REAL)))/T(J))

        END DO

        b = SQRT(T(J)) ** 3 / (DE_BROGLIE_L3 * DENSITY_0 * ION_SUM)
        n(J) =  0.50D00 * b * (SQRT(1.00D00 + 4.000D00 * n_TOT(J) / b) - 1.00D00)

        DO I = 1, NUM_NEUTRALS

          n_neut(J, I) = (n_TOT(J) - n(J)) * I ** 2 &
          * EXP(-(ION_POT_H * (1.00D00 - 1.00D00 / REAL(I ** 2,KIND=DEFAULT_REAL)))/T(J)) / EX_SUM

        END DO

      END DO

    ELSE

      n = n_TOT

    END IF

!Initialize f_0
    DO I = 1, NUM_X

      n_num(I) = 0

      F_0_POS = X_POS(I) + NUM_V * (NUM_H - 1)

      DO J = 1, NUM_V

        f0(I,J) = n(I) * ((PI * T(I)) ** (- 3.00D00/2.00D00)) * EXP( - (V_GRID(J) ** 2)/ T(I))    !Maxwellian initialization
        n_num(I) = n_num(I) + 4.00D00 * PI * V_GRID(J) ** 2 * V_GRID_WIDTH(J) * f0(I,J)


      END DO
!Correct for numerical integration density error
      DO J = 1, NUM_V

        f(F_0_POS + J) = f0(I,J) * n(I)/n_num(I)

      END DO

    END DO

!Initialize f_1 and E based on local collisional values
    IF (LOCAL_INIT_SWITCH) THEN

!Calculate density and temperature gradients
      GRAD_T = 0
      GRAD_n = 0

      DO I = 2, NUM_X - 1

        GRAD_T(I) = (T(I+1) - T(I-1))/(X_GRID(I+1)-X_GRID(I-1))
        GRAD_n(I) = (n(I+1) - n(I-1))/(X_GRID(I+1)-X_GRID(I-1))

      END DO

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        GRAD_T(1) = (T(2) - T(NUM_X))/dx
        GRAD_n(1) = (n(2) - n(NUM_X))/dx

        GRAD_T(NUM_X) = (T(1) - T(NUM_X-1))/dx
        GRAD_n(NUM_X) = (n(1) - n(NUM_X-1))/dx

      ELSE

        IF (NO_FLOW_BOUNDARY_UP_SWITCH) THEN

          GRAD_T(1) = (T(2) - T(1))/dx
          GRAD_n(1) = (n(2) - n(1))/dx

        END IF

        IF (NO_FLOW_BOUNDARY_DIV_SWITCH) THEN

          GRAD_T(NUM_X) = (T(NUM_X) - T(NUM_X-1))/(2*(X_GRID(NUM_X)-X_GRID(NUM_X-1)))
          GRAD_n(NUM_X) = (n(NUM_X) - n(NUM_X-1))/(2*(X_GRID(NUM_X)-X_GRID(NUM_X-1)))

        END IF

      END IF

!Calculate E-field
      DO I = 1, NUM_X

        E(I) = - (GRAD_n(I) * T(I)/n(I) + 5.00D00 * GRAD_T(I) / 2.00D00) / 2.00D00

        f(X_POS(I) + NUM_H * NUM_V + NUM_NEUTRALS + 1) = E(I)

      END DO
!Calculate maxwellian gradients
      DO I = 1, NUM_X

        DO J = 1, NUM_V

          GRAD_F_V(I,J) = - 2 * V_GRID(J) * f0(I,J)/T(I)

          GRAD_F_X(I,J) = f0(I,J) * GRAD_n(I) / n(I) + (V_GRID(J) ** 2/T(I) - 3.00D00/2.00D00) * f0(I,J) * GRAD_T(I) / T(I)

        END DO

      END DO
!Set local normalization of f_1 so that heat flux is Braginskii
      DO I = 1, NUM_X

        A(I) = 32.00D00 * 4.00D00 * n(I) /(PI *3.00D00 *3.20D00)

      END DO

!Initialize f_1
      DO I = 1, NUM_X

          f_1_POS = X_POS(I) + NUM_V * H_POS(1)

          DO J = 1, NUM_V

            f1 = - (V_GRID(J) ** 4 * GRAD_F_X(I,J)  - V_GRID(J) ** 3 * GRAD_F_V(I,J) * E(I)) / A(I)

            f(f_1_POS + J) = f1

          END DO

        END DO

      END IF

      IF (NEUTRAL_TRACK_SWITCH) THEN                                            !Neutral initialization

        IF (.NOT. LOCAL_SAHA_BOLTZMANN_INIT_SWITCH) THEN

          n_n = 0

          IF (UNIFORM_NEUTRAL_INIT) n_n = NEUTRAL_DENSITY                         !Uniform neutrals

          IF (NEUTRAL_CLOUD_INIT) THEN                                            !Neutral cloud

            CALL INPUT_NEUTRAL_CLOUD(NUM_CLOUD)

            IF ((DROP_INIT) .AND. (PLASMA_RAMP_WIDTH .GT. 0)) THEN                !Initialize with constant total density

              DO I = 1, NUM_X

                n_n(I) = DENS_UP - n_TOT(I)

              END DO

            ELSE

              DO I = NUM_CLOUD + 1, NUM_X

                n_n(I) = NEUTRAL_DENSITY

              END DO

            END IF

          END IF

          IF (NEUTRAL_GROUND_DENS_FROM_FILE_INIT) CALL INPUT_INIT('NEUTRAL_DENS', n_n, NUM_X) !Initialize neutrals from file

!Initialize neutral state densities if not Saha-Boltzmann
          DO I = 1, NUM_X

            F_0_POS = X_POS(I) + NUM_V * NUM_H + 1

            f(F_0_POS) = n_n(I)

          END DO

        ELSE
!Initialize neutral state densities if Saha-Boltzmann
          DO I = 1, NUM_X

            DO J = 1, NUM_NEUTRALS

              F_0_POS = X_POS(I) + NUM_V * NUM_H + J

              f(F_0_POS) = n_neut(I,J)

            END DO

          END DO

        END IF

      END IF

!Cold ion initializations
      IF (COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) THEN

        u_i = 0

        IF (ION_VEL_FROM_FILE_INIT) CALL INPUT_INIT('ION_VEL',u_i,NUM_X) !Initialize ion velocity from file

!Initialize f_1 to local slow drifting Maxwellian
        DO I = 1, NUM_X

          DO J = 1, NUM_V

            F_0_POS = X_POS(I) + NUM_V * (NUM_H - 1) + J

            f(F_0_POS - NUM_V) = 2.00D00 * u_i(I) * V_GRID(J) * f(F_0_POS) / T(I)

          END DO

        END DO

!Initialize ion quanitites
        DO I = 1, NUM_X

          F_0_POS = X_POS(I) + NUM_0D

          f(F_0_POS) = u_i(I)                                              !Ion velocity
          f(F_0_POS - 1) = n(I)

        END DO

      END IF
!Initialize fluid electron quanitities
      IF (FULL_FLUID_MODE) THEN

        DO I = 1, NUM_X

          F_0_POS = X_POS(I) + NUM_0D

          f(F_0_POS - 2) = T(I)
          f(F_0_POS - 3) = u_i(I)
          f(F_0_POS - 4) = n(I)

        END DO

      END IF

    END SUBROUTINE INIT_F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes f_old using adaptive restarting
  SUBROUTINE INIT_F_ADAPTIVE(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(OUT) :: f(DIM_F) !< Vector to be initialized

    REAL(KIND=DEFAULT_REAL) :: old_dv

    REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: old_f(:)

    INTEGER :: old_num_neutral, old_num_v, old_num_l, old_num_0D, old_dim_f

    INTEGER :: I, J, K

    LOGICAL :: OLD_FF_MODE

    CALL INPUT_ADAPTIVE_RESTART(old_num_neutral,old_num_l, old_num_v, old_dv,OLD_FF_MODE)

    f = 0

    old_num_0D = NUM_0D - NUM_V * NUM_H - NUM_NEUTRALS + old_num_neutral + old_num_v * (1 + old_num_l)

    IF (OLD_FF_MODE .AND. (.NOT. FULL_FLUID_MODE)) old_num_0D = old_num_0D + 3

    old_dim_f = NUM_X * old_num_0D

    ALLOCATE(old_f(old_dim_f))

    CALL INPUT_RESTART_INIT('VAR_VEC',old_f,old_dim_f)

    DO I = 1, NUM_X

!Input distribution function as they are
      DO J = 1, old_num_l + 1

        DO K = 1, old_num_v

          f(X_POS(I) + (NUM_H - J) * NUM_V + K) = old_f((I-1) * old_num_0D + old_num_v * (old_num_l + 1 - J) + K)

        END DO

      END DO

      DO J = 1, old_num_neutral

        f(X_POS(I) + NUM_V * NUM_H + J) = old_f((I-1) * old_num_0D + old_num_v * (old_num_l + 1) + J)

      END DO

!Input electric field
      f(X_POS(I) + NUM_V*NUM_H + NUM_NEUTRALS + 1) = old_f((I-1) * old_num_0D + old_num_v * (old_num_l + 1) + old_num_neutral + 1)

!Input fluid quantities

      IF (FULL_FLUID_MODE) THEN

        f(X_POS(I) + NUM_0D - 2) = old_f(I * old_num_0D - 2)
        f(X_POS(I) + NUM_0D - 3) = old_f(I * old_num_0D - 3)
        f(X_POS(I) + NUM_0D - 4) = old_f(I * old_num_0D - 4)

      END IF

      f(X_POS(I) + NUM_0D - 1) = old_f(I * old_num_0D - 1)

      f(X_POS(I) + NUM_0D ) = old_f(I * old_num_0D)

    END DO

  END SUBROUTINE INIT_F_ADAPTIVE
!-------------------------------------------------------------------------------------------------------------------------------------
END  MODULE F_INIT
!-------------------------------------------------------------------------------------------------------------------------------------
