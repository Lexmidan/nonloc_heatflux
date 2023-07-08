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
!>Contains calls to submatrix builders and adds submatrices into main matrix
MODULE MAT_BUILD

  USE GRID
  USE SWITCHES
  USE BUILD_0D
  USE BUILD_VLASOV
  USE BUILD_MAXWELL
  USE BUILD_HEATING
  USE BUILD_NEUTRAL_DIFF
  USE PRINTING
  USE BUILD_DIV_BOUNDARY
  USE BUILD_ION_FLUID
  USE BUILD_FLUID_EL
  USE BUILD_PARTICLE_SOURCE
  USE VAR_KIND_DATA
  USE MPI

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills evolution matrix
  SUBROUTINE FILL_M(f_lagged, n_old, T_old, n_lagged, T_lagged, n_i_lagged, u_i_lagged, u_lagged,&
                   TIMESTEP, N_NONLIN)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(DIM_F) !< Lagged vector
    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: n_old, & !< Old density
                                             T_old, & !< Old temperature
                                             n_lagged, & !< Lagged density
                                             u_lagged, & !< Lagged electron
                                             T_lagged !< Lagged temperature
    INTEGER, INTENT(IN) :: TIMESTEP, & !< Current timestep
                           N_NONLIN !< Current nonlinear iteration

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i_lagged(NUM_X), &  !< Lagged ion density
                                  u_i_lagged(NUM_X) !< Lagged ion velocity

    INTEGER :: RANK, IERR, SIZE

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    IF (RANK .EQ. 0) CALL PRINT_ECHO('Building collision integral and CRM submatrices')
    CALL FILL_0D(f_lagged, n_old, T_old, n_lagged, T_lagged, n_i_lagged, u_i_lagged, N_NONLIN)          !Call 0D submatrix builders

    IF (.NOT. FULL_FLUID_MODE) THEN

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Building Vlasov term submatrices')
      CALL FILL_VLASOV(f_lagged)                                                  !Call Vlasov term submatrix builders

      IF (HEATING_ON) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('Building heating submatrix')
        CALL FILL_HEATING(TIMESTEP,n_lagged)                                     !Call heating operator builder

      END IF

  !Particle source builders
      IF (PART_SOURCE_ON) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('Building particle source submatrix')
        CALL FILL_PART_SOURCE(TIMESTEP,n_lagged,T_lagged)
        IF (COLD_ION_FLUID_SWITCH .AND. (.NOT. ION_CONT_OFF_SWITCH)) &
        CALL FILL_ION_PART_SOURCE(TIMESTEP,n_lagged,T_lagged)

      END IF

    END IF

    IF (NEUTRAL_DIFF_ON) THEN

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Building neutral diffusion submatrix')

      IF (COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) THEN

        CALL FILL_DIFF_N(f_lagged, n_i_lagged)                                      !Call neutral diffusion operator builder

      ELSE

        CALL FILL_DIFF_N(f_lagged, n_lagged)                                      !Call neutral diffusion operator builder

      END IF

    END IF

    IF ((REC_ON .OR. PLASMA_SINK_ON) .AND. (RANK .EQ. SIZE - 1)) THEN           !Call divertor boundary submatrix builder if last processor

      CALL PRINT_ECHO('Building divertor boundary submatrix')

      IF (COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) THEN

        IF (SONIC_OUTFLOW_DIV_SWITCH) THEN

          CALL FILL_DIV_BOUNDARY(n_lagged,T_lagged,&
                                f_lagged, &
                                u_i_lagged, n_i_lagged)


          IF (ION_CONT_ON) CALL FILL_ION_SINK(n_i_lagged)


        ELSE

          CALL FILL_DIV_BOUNDARY(n_lagged,T_lagged,&
                                f_lagged, &
                                u_i_lagged, n_i_lagged)


          IF (ion_flux .GT. 0 .AND. ION_CONT_ON) CALL FILL_ION_SINK(n_i_lagged)


        END IF

      ELSE

        CALL FILL_DIV_BOUNDARY(n_lagged,T_lagged,&
                            f_lagged)
      END IF

    END IF

    IF (PLASMA_SINK_SWITCH) THEN

        CALL MPI_Bcast(T_e_boundary,1,MPI_REAL8,SIZE-1,MPI_COMM_WORLD,IERR)
        CALL MPI_Bcast(ion_flux,1,MPI_REAL8,SIZE-1,MPI_COMM_WORLD,IERR)
        CALL MPI_Bcast(cut_off_v,1,MPI_REAL8,SIZE-1,MPI_COMM_WORLD,IERR)
        CALL MPI_Bcast(n_e_boundary,1,MPI_REAL8,SIZE-1,MPI_COMM_WORLD,IERR)
        CALL MPI_Bcast(BOHM_VALUE,1,MPI_LOGICAL,SIZE-1,MPI_COMM_WORLD,IERR)
        CALL MPI_Bcast(FLOW_LIMITER,1,MPI_REAL8,SIZE-1,MPI_COMM_WORLD,IERR)

    END IF

    IF (COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) THEN                                             !Call cold fluid operator builders

      IF (X_ADV_SWITCH) THEN
        IF (RANK .EQ. 0) CALL PRINT_ECHO('Building fluid ion submatrices')
        IF (ION_FLOW_ON) CALL FILL_ION_CONV(u_i_lagged)
        IF (ION_FLOW_ON .AND. ION_EL_TEMP_SWITCH) CALL FILL_ION_PRESSURE(n_i_lagged,T_lagged)
        IF (ION_CONT_ON) CALL FILL_ION_CONT(n_i_lagged)

      END IF

      IF (MAXWELL_SWITCH .AND. ION_FLOW_ON) CALL FILL_MAXWELL_ION(n_i_lagged)

      IF (E_ADV_SWITCH .AND. ION_FLOW_ON) CALL FILL_ION_LORENTZ

      IF (ION_CONT_ON) THEN

        IF (COLL_EN_ION) THEN

            CALL FILL_ION_GAIN

        END IF

        IF (COLL_RECOMB) THEN

           IF (FULL_FLUID_MODE) THEN

               CALL FILL_ION_LOSS_FF(n_i_lagged,T_lagged)

           ELSE

            CALL FILL_ION_LOSS(n_i_lagged,T_lagged)

           END IF

        END IF

      END IF

      IF (ION_FLOW_ON .AND. (COLL_EN_ION .OR. COLL_RECOMB)) THEN

          CALL FILL_ION_MOM_SOURCE(n_i_lagged,n_lagged,T_lagged,f_lagged)

      END IF

      IF (ION_FLOW_ON .AND. SIMPLE_CX_SWITCH .AND. NEUTRAL_TRACK_SWITCH) &
       CALL FILL_SIMPLE_CX(f_lagged,u_i_lagged)

      IF (COLL_EI_L_SWITCH .AND. ION_FLOW_ON) THEN

        IF (FULL_FLUID_MODE) THEN

          CALL FILL_ION_T_DRAG(n_lagged,n_i_lagged)

          CALL FILL_ION_U_DRAG(n_old,T_old,n_lagged,T_lagged,n_i_lagged)

        ELSE

          CALL FILL_ION_DRAG(n_i_lagged)

        END IF

      END IF

      IF (ION_FLOW_ON .AND. ((REC_ON .OR. PLASMA_SINK_ON) &
          .AND. ((MIN_X .LE. NUM_X - 1) .AND. (MAX_X .GE. NUM_X - 1)))) THEN

        IF (SONIC_OUTFLOW_DIV_SWITCH) THEN

          CALL FILL_ION_CONV_DIV(u_i_lagged,n_i_lagged)

        ELSE

          IF (ion_flux .GT. 0) CALL FILL_ION_CONV_DIV(u_i_lagged,n_i_lagged)

        END IF

      END IF

    END IF

    IF (FULL_FLUID_MODE) THEN

      IF (RANK .EQ. 0) CALL PRINT_ECHO('Building fluid electron submatrices')

      IF (EL_FLOW_ON) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('Building fluid electron flow submatrices')

        IF (MAXWELL_SWITCH) CALL FILL_MAXWELL_EL(n_lagged)
        IF (E_ADV_SWITCH) CALL FILL_EL_LORENTZ

        IF ((REC_ON .OR. PLASMA_SINK_ON) &
            .AND. ((MIN_X .LE. NUM_X - 1) .AND. (MAX_X .GE. NUM_X - 1))) THEN

          IF (ion_flux .GT. 0) CALL FILL_EL_CONV_DIV(u_lagged,n_lagged)

        END IF

        IF (COLL_EN_ION .OR. COLL_RECOMB) THEN

          CALL FILL_EL_MOM_SOURCE(n_i_lagged,n_lagged,T_lagged,f_lagged)

        END IF

        IF (X_ADV_SWITCH) THEN

          CALL FILL_EL_CONV(u_lagged)
          CALL FILL_EL_PRESSURE(n_lagged,T_lagged)
          IF (COLL_EI_L_SWITCH) THEN

            CALL FILL_EL_T_DRAG

            CALL FILL_EL_U_DRAG(n_old,T_old,n_lagged,T_lagged)

          END IF

          IF (COLL_EN_EX .OR. COLL_EN_ION .OR. COLL_EN_EL_L_SWITCH .OR. COLL_RECOMB) THEN

                CALL FILL_EL_EN_DRAG(f_lagged,n_i_lagged,n_lagged)

          END IF

        END IF

      END IF

      IF (EL_CONT_ON) THEN

        IF (RANK .EQ. 0) CALL PRINT_ECHO('   Building fluid electron continuity submatrices')

        IF (X_ADV_SWITCH) CALL FILL_EL_CONT(n_lagged)
        IF (COLL_EN_ION) THEN

          CALL FILL_EL_GAIN

        END IF
        IF (COLL_RECOMB) THEN

          CALL FILL_EL_LOSS(n_i_lagged,T_lagged)

        END IF
        IF ((REC_ON .OR. PLASMA_SINK_ON) .AND. (RANK .EQ. SIZE - 1) .AND. (ion_flux .GT. 0)) &
        CALL FILL_EL_SINK(n_lagged)
        IF (X_ADV_SWITCH) CALL FILL_EL_TEMP_CONV(u_lagged)
        IF (X_ADV_SWITCH) CALL FILL_EL_TEMP_PDV(T_lagged)
        IF ((REC_ON .OR. PLASMA_SINK_ON) .AND. (RANK .EQ. SIZE - 1) .AND. (ion_flux .GT. 0) .AND. X_ADV_SWITCH) &
        CALL FILL_EL_TEMP_PDV_DIV(u_lagged,n_lagged)
        IF (COLL_EI_L_SWITCH) THEN

          CALL FILL_EL_TEMP_Q_T(T_lagged,n_lagged)

          IF (X_ADV_SWITCH) CALL FILL_EL_TEMP_Q_U(T_lagged,n_lagged)
          IF (X_ADV_SWITCH) THEN

            CALL FILL_EL_T_DRAG_TEMP(u_lagged)

          END IF
          IF (X_ADV_SWITCH) &
          CALL FILL_EL_U_DRAG_TEMP(n_old,T_old,u_lagged,T_lagged,u_i_lagged)

        END IF
        IF ((REC_ON .OR. PLASMA_SINK_ON) .AND. (RANK .EQ. SIZE - 1) .AND. (ion_flux .GT. 0)) THEN

          CALL FILL_EL_TEMP_Q_DIV(n_lagged)

        END IF
        IF (COLL_EN_ION) THEN

            CALL FILL_EL_GAIN_TEMP(n_lagged,u_lagged,T_old)

        END IF
        IF (COLL_RECOMB) THEN

          CALL FILL_EL_LOSS_TEMP(n_i_lagged, T_old, n_lagged, u_lagged)

        END IF

        IF (X_ADV_SWITCH) THEN

          IF (COLL_EN_EX .OR. COLL_EN_ION .OR. COLL_EN_EL_L_SWITCH .OR. COLL_RECOMB) THEN

            CALL FILL_EL_EN_DRAG_TEMP(n_lagged,u_lagged)

          END IF

        END IF

        IF (HEATING_ON) CALL FILL_EL_HEATING(n_lagged,T_lagged,TIMESTEP)

        IF (COLL_EN_EX .OR. COLL_EN_ION .OR. COLL_RECOMB) THEN

            CALL FILL_EL_EN_INCOLL_TEMP(n_lagged,T_lagged,n_i_lagged)

        END IF

        IF (COLL_EN_EL_L_SWITCH) CALL FILL_EL_EN_ELCOLL_TEMP(T_lagged)

        IF (COLL_RECOMB) THEN

            CALL FILL_CRM_RECOMB_FF(T_lagged, n_i_lagged)

        END IF

      END IF

    END IF

  IF (RANK .EQ. 0) CALL PRINT_ECHO('Matrix building complete')

END SUBROUTINE FILL_M
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE MAT_BUILD
!-------------------------------------------------------------------------------------------------------------------------------------
