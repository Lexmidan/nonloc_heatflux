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
!> Contains post-processing routine and output calls
MODULE POST_PROCESSING

  USE GRID
  USE SWITCHES
  USE NEUT_AND_HEAT
  USE MOMENTS
  USE OUTPUT
  USE TESTS
  USE PRINTING
  USE NORMALIZATION
  USE MPI
  USE BUILD_DIV_BOUNDARY
  USE EVOLVE
  USE VAR_KIND_DATA

  REAL(KIND=DEFAULT_REAL),ALLOCATABLE :: OLD_TOT_DENS(:), OLD_N_TOT_DENS(:)

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Contains all function post-processing and output calls
  SUBROUTINE POST_PROC(f,TIMESTEP,TIME)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(DIM_F), & !< Input vector
                                TIME !< Current elapsed time
    INTEGER, INTENT(IN) :: TIMESTEP !< Input timestep

    REAL(KIND=DEFAULT_REAL) :: f_chopped(0:L_MAX,NUM_V,NUM_X) !< Chopped-up vector into local harmonics

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X) :: n, &                                            !Moments and output
                                T, &
                                u_x, &
                                q_x, &
                                q_x_ratio, &
                                E_x, &
                                a_test, &
                                u_i, &
                                n_i, &
                                qn_test, &
                                c_test, &
                                S_i_SK, &
                                S_i_M, &
                                q_tot, &
                                ION_E, &
                                EX_E, &
                                REC_3B_E, &
                                DEEX_E, &
                                RAD_DEEX_E, &
                                RAD_REC_E, &
                                num_dv_err


    REAL(KIND=DEFAULT_REAL) :: n_n(NUM_NEUTRALS,NUM_X), &                                        !Neutrals
              q(3), &
              n_tot, &
              n_n_tot

    INTEGER :: I, J, K

    CALL PRINT_ECHO('Starting post-processing and output')
    n = 0 ; T = 0 ; u_x = 0 ; q_x = 0 ; q_x_ratio = 0 ; E_x = 0 ; n_n = 0; q_tot = 0

    f_chopped = 0

    DO I = 1, NUM_X                                                             !Chop f into local harmonic
      DO J = 1, NUM_H
        DO K = 1, NUM_V

          f_chopped(J-1,K,I) = f(X_POS(I) + NUM_V * (NUM_H - J) + K)

        END DO
      END DO
    END DO

    CALL OUTPUT_F(f_chopped, TIMESTEP)                                          !Output harmonics

    IF (OUTPUT_DENSITY) THEN

      IF (FULL_FLUID_MODE) THEN

        DO I = 1, NUM_X

          n(I) = f(X_POS(I) + NUM_0D - 4)

        END DO

      ELSE

        DO I = 1, NUM_X

          n(I) = DENSITY_MOMENT(f_chopped(0,:,I))

        END DO

      END IF

      CALL PRINT_ECHO('Outputting density vector...')
      CALL OUTPUT_VECTOR_X(n,'DENSITY',TIMESTEP)                                !Output density

    END IF

    IF (OUTPUT_TEMP) THEN

      IF (FULL_FLUID_MODE) THEN

        DO I = 1, NUM_X

          T(I) = f(X_POS(I) + NUM_0D - 2)

        END DO

      ELSE

        IF (L_MAX .EQ. 0) THEN

          DO I = 1, NUM_X

            T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I))

          END DO

        ELSE

          DO I = 1, NUM_X

            T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I), &
                                      f_10 = f_chopped(1,:,I))

          END DO

        END IF

      END IF

      CALL PRINT_ECHO('Outputting temperature vector...')
      CALL OUTPUT_VECTOR_X(T,'TEMPERATURE',TIMESTEP)                            !Output temperature

    END IF

    IF (OUTPUT_FLOW_VEL) THEN

      IF (FULL_FLUID_MODE) THEN

        DO I = 1, NUM_X

          u_x(I) = f(X_POS(I) + NUM_0D - 3)

        END DO

      ELSE

        IF (L_MAX .GT. 0) THEN

          DO I = 1, NUM_X

            u_x(I) = FLOW_V_X_MOMENT(f_1 = f_chopped(1,:,I), &
                                     f_0 = f_chopped(0,:,I))

          END DO

        END IF

      END IF

      CALL PRINT_ECHO('Outputting flow velocity vector x-component...')
      CALL OUTPUT_VECTOR_X(u_x,'FLOW_VEL_X',TIMESTEP)                         !Output parallel flow velocity

    END IF

    IF (OUTPUT_HEAT_FLOW) THEN                                                  !Output heatflows

      IF (L_MAX .GT. 1) THEN

        DO I = 1, NUM_X

          q = HEAT_FLOW_MOMENT(f_0 = f_chopped(0,:,I), &
                               f_10 = f_chopped(1,:,I), &
                               f_20 = f_chopped(2,:,I) &
                               )
          q_x(I) = q(1)

        END DO

        CALL PRINT_ECHO('Outputting heat_flow vector...')
        CALL OUTPUT_VECTOR_X(q_x,'HEAT_FLOW_X',TIMESTEP)

      ELSE IF (L_MAX .EQ. 1) THEN

        DO I = 1, NUM_X

          q = HEAT_FLOW_MOMENT(f_0 = f_chopped(0,:,I), &
                               f_10 = f_chopped(1,:,I) &
                               )
          q_x(I) = q(1)


        END DO

        CALL PRINT_ECHO('Outputting heat_flow vector...')
        CALL OUTPUT_VECTOR_X(q_x,'HEAT_FLOW_X',TIMESTEP)




      END IF

        DO I = 1, NUM_X

          q_tot(I)=MOMENT_V(f_chopped(1,:,I),3)/6.00D00

        END DO

        CALL PRINT_ECHO('Outputting total heat_flow vector...')
        CALL OUTPUT_VECTOR_X(q_tot,'HEAT_FLOW_X_TOT',TIMESTEP)

    END IF

  IF (OUTPUT_E_FIELD) THEN                                                      !Output E-field

    IF ((MAXWELL_SWITCH) .OR. (E_ADV_SWITCH)) THEN

      DO I = 1, NUM_X

        E_x(I) = f(X_POS(I) + NUM_V*NUM_H + NUM_NEUTRALS + 1)

      END DO
      num_dv_err = NUM_DV_HEATING(f)

      CALL PRINT_ECHO('Outputting E-field...')
      CALL OUTPUT_VECTOR_X(E_x,'E_FIELD_X',TIMESTEP)
      CALL OUTPUT_VECTOR_X(num_dv_err,'NUM_DV_HEATING',TIMESTEP)

    END IF

  END IF

  IF (OUTPUT_NEUTRAL_DATA) THEN

    IF (NEUTRAL_TRACK_SWITCH) THEN

      DO I = 1, NUM_X

        DO J = 1, NUM_NEUTRALS

          n_n(J,I) = f(X_POS(I) + NUM_V*NUM_H +J)

        END DO

      END DO

      CALL PRINT_ECHO('Outputting neutral densities...')
      CALL OUTPUT_NEUTRALS(n_n, TIMESTEP)                                       !Output neutral density

    END IF

  END IF

  IF (OUTPUT_SH_TEST) THEN                                                      !Output ratio of heatflow to SH heatflow

    IF (L_MAX .GT. 1) THEN

      DO I = 1, NUM_X

        q = HEAT_FLOW_MOMENT(f_0 = f_chopped(0,:,I), &
                             f_10 = f_chopped(1,:,I), &
                             f_20 = f_chopped(2,:,I) &
                             )
        q_x(I) = q(1)

        T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I), &
                                  f_10 = f_chopped(1,:,I) &
                                  )

      END DO

      q_x_ratio = SH_TEST(q_x, T,n)

      CALL PRINT_ECHO('Outputting SH ratios...')
      CALL OUTPUT_VECTOR_X(q_x_ratio,'SH_q_ratio',TIMESTEP)

    ELSE IF (L_MAX .EQ. 1) THEN

      DO I = 1, NUM_X

        q = HEAT_FLOW_MOMENT(f_0 = f_chopped(0,:,I), &
                               f_10 = f_chopped(1,:,I) &
                               )
        q_x(I) = q(1)

        T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I), &
                                    f_10 = f_chopped(1,:,I) &
                                    )

      END DO

      q_x_ratio = SH_TEST(q_x, T,n)

      CALL PRINT_ECHO('Outputting SH ratios...')
      CALL OUTPUT_VECTOR_X(q_x_ratio,'SH_q_ratio',TIMESTEP)


    END IF

  END IF

  IF ((OUTPUT_ATOMIC_EN_TEST_SWITCH) .AND. (NEUTRAL_TRACK_SWITCH)) THEN

    IF (FULL_FLUID_MODE) THEN

      DO I = 1, NUM_X

        T(I) = f(X_POS(I) + NUM_0D - 2)

      END DO

      DO I = 1, NUM_X

        n(I) = f(X_POS(I) + NUM_0D - 4)

      END DO

    ELSE

      IF (L_MAX .GT. 0) THEN

        DO I = 1, NUM_X

          T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I), &
                                    f_10 = f_chopped(1,:,I) &
                                    )

        END DO

      ELSE

        DO I = 1, NUM_X

          T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I))

        END DO

      END IF

      DO I = 1, NUM_X

        n(I) = DENSITY_MOMENT(f_chopped(0,:,I))

      END DO

    END IF

    IF (.NOT. OUTPUT_NEUTRAL_DATA) THEN

      DO I = 1, NUM_X

        DO J = 1, NUM_NEUTRALS

          n_n(J,I) = f(X_POS(I) + NUM_V*NUM_H +J)

        END DO

      END DO

    END IF

    a_test = ATOMIC_EN_TEST(T,n,n_n)

    CALL PRINT_ECHO('Outputting atomic energy test data...')
    CALL OUTPUT_VECTOR_X(a_test,'ATOMIC_EN_TEST',TIMESTEP)

  END IF

  IF (COLD_ION_FLUID_SWITCH) THEN

    IF (OUTPUT_DENSITY) THEN

      DO I = 1, NUM_X

        n_i(I) = f(X_POS(I) + NUM_0D - 1)

      END DO

      CALL PRINT_ECHO('Outputting ion density...')
      CALL OUTPUT_VECTOR_X(n_i,'ION_DENS',TIMESTEP)

    END IF

    IF (OUTPUT_FLOW_VEL) THEN

      DO I = 1, NUM_X

        u_i(I) = f(X_POS(I) + NUM_0D)

      END DO

      CALL PRINT_ECHO('Outputting ion velocity...')
      CALL OUTPUT_VECTOR_X(u_i,'ION_VEL',TIMESTEP)

    END IF

    IF (OUTPUT_QN_TEST) THEN

      qn_test = DIFF_TEST(n,n_i)

      CALL PRINT_ECHO('Outputting QN test...')
      CALL OUTPUT_VECTOR_X(qn_test,'QN_TEST',TIMESTEP)

    END IF

    IF (OUTPUT_CURRENT_TEST) THEN

      c_test = DIFF_TEST(u_x, Z_PROF * u_i)

      CALL PRINT_ECHO('Outputting current test...')
      CALL OUTPUT_VECTOR_X(c_test,'CURRENT_TEST',TIMESTEP)

    END IF

  END IF

  IF ((OUTPUT_RATE_DATA) .AND. (NEUTRAL_TRACK_SWITCH)) THEN

    IF (FULL_FLUID_MODE) THEN

      DO I = 1, NUM_X

        T(I) = f(X_POS(I) + NUM_0D - 2)

      END DO

    ELSE

      IF (L_MAX .EQ. 0) THEN

        DO I = 1, NUM_X

          T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I))

        END DO

      ELSE

        DO I = 1, NUM_X

          T(I) = TEMPERATURE_MOMENT(f_0 = f_chopped(0,:,I), &
                                    f_10 = f_chopped(1,:,I))

        END DO

      END IF

    END IF

    IF (FULL_FLUID_MODE) THEN

      DO I = 1, NUM_X

        n(I) = f(X_POS(I) + NUM_0D - 4)

      END DO

    ELSE


      DO I = 1, NUM_X

        n(I) = DENSITY_MOMENT(f_chopped(0,:,I))

      END DO

    END IF

    IF (COLD_ION_FLUID_SWITCH) THEN

      DO I = 1, NUM_X

        n_i(I) = f(X_POS(I) + NUM_0D - 1)

      END DO

    ELSE

      n_i = n

    END IF

    IF (COLL_EN_ION) THEN

      S_i_SK = S_ION_SK(f)
      S_i_M = S_ION_M(T,n,f)

      ION_E = ION_E_RATE(f)

      CALL PRINT_ECHO('Outputting ionization rate data...')
      CALL OUTPUT_VECTOR_X(S_i_SK,'S_ION_SK',TIMESTEP)
      CALL OUTPUT_VECTOR_X(S_i_M,'S_ION_M',TIMESTEP)
      CALL PRINT_ECHO('Outputting energy loss rate due to ionization...')
      CALL OUTPUT_VECTOR_X(ION_E,'ION_E_RATE',TIMESTEP)

    END IF

    IF (COLL_EN_EX) THEN

      EX_E = EX_E_RATE(f)
      DEEX_E = DEEX_E_RATE(f)
      RAD_DEEX_E = RAD_DEEX_E_RATE(f)

      CALL UPDATE_SIGMA_L_DEEX(T,.TRUE.)

      CALL PRINT_ECHO('Outputting excitation and deexcitation energy exchange rates...')
      CALL OUTPUT_VECTOR_X(EX_E,'EX_E_RATE',TIMESTEP)
      CALL OUTPUT_VECTOR_X(DEEX_E,'DEEX_E_RATE',TIMESTEP)
      CALL OUTPUT_VECTOR_X(RAD_DEEX_E,'RAD_DEEX_E_RATE',TIMESTEP)

    END IF

    IF (COLL_RECOMB) THEN

      REC_3B_E = REC_3B_E_RATE(f,n,n_i)
      RAD_REC_E = RAD_REC_E_RATE(T,n,n_i)

      CALL UPDATE_SIGMA_L_RECOMB(T,.TRUE.)

      CALL PRINT_ECHO('Outputting recombination energy exchange rates...')
      CALL OUTPUT_VECTOR_X(REC_3B_E,'REC_3B_E_RATE',TIMESTEP)
      CALL OUTPUT_VECTOR_X(RAD_REC_E,'RAD_REC_E_RATE',TIMESTEP)

    END IF

  END IF

  IF (PLASMA_SINK_SWITCH) THEN

    IF (SONIC_OUTFLOW_DIV_SWITCH .OR. (.NOT.(COLD_ION_FLUID_SWITCH))) THEN

      CALL PRINT_ECHO('Outputting sheath data...')
      CALL OUTPUT_SHEATH_DATA(gamma_e,pot_drop,cut_off_v,FLOW_LIMITER,BOHM_VALUE,TIMESTEP)

    END IF

  END IF

  IF (ADAPTIVE_TIMESTEP_SWITCH) THEN

    CALL PRINT_ECHO('Outputting timestep data...')
    CALL OUTPUT_TIMESTEP_DATA(ADAPTIVE_dt,TIMESTEP,TIME)

  END IF

  IF (OUTPUT_DENSITY) THEN

    n_tot = 0
    n_n_tot = 0

    IF (COLD_ION_FLUID_SWITCH .AND. (.NOT. ION_CONT_OFF_SWITCH)) THEN

      DO I = 1, NUM_X

        IF (MOD(I,2) .EQ. 1) n_tot = n_tot + n_i(I)*dxc(I)/(X_GRID(NUM_X)+dxc(NUM_X)/2+dxc(1)/2)

      END DO

    ELSE

      DO I = 1, NUM_X

        IF (MOD(I,2) .EQ. 1) n_tot = n_tot + n(I)*dxc(I)/(X_GRID(NUM_X)+dxc(NUM_X)/2+dxc(1)/2)

      END DO

    END IF

    IF (NEUTRAL_TRACK_SWITCH) THEN

      DO I = 1,NUM_X

        IF (MOD(I,2) .EQ. 1) THEN

          DO J = 1, NUM_NEUTRALS

            n_n_tot = n_n_tot + n_n(J,I)*dxc(I)/(X_GRID(NUM_X)+dxc(NUM_X)/2+dxc(1)/2)

          END DO

        END IF

      END DO

    END IF

    n_tot = n_tot + n_n_tot

    IF (.NOT. ALLOCATED(OLD_TOT_DENS)) THEN

      ALLOCATE(OLD_TOT_DENS(1))
      ALLOCATE(OLD_N_TOT_DENS(1))

      OLD_TOT_DENS = n_tot
      OLD_N_TOT_DENS = n_n_tot

    END IF

    CALL PRINT_ECHO('Outputting total density data...')
    CALL OUTPUT_TOT_DENS_DATA(TIMESTEP,n_tot,n_n_tot,OLD_TOT_DENS(1),OLD_N_TOT_DENS(1))

    OLD_TOT_DENS = n_tot
    OLD_N_TOT_DENS = n_n_tot

  END IF


  CALL PRINT_ECHO('Post-processing complete')
  CALL PRINT_ECHO('---------------------------------------------------------------------------------------------------------------&
                  &--')

END SUBROUTINE POST_PROC
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE POST_PROCESSING
