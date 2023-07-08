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
!> Contains submatrix building routines for fluid ion operators
MODULE BUILD_ION_FLUID

  USE GRID
  USE NORMALIZATION
  USE SWITCHES
  USE MATRIX_DATA
  USE COLL_CS_RATES
  USE NEUT_AND_HEAT
  USE BUILD_COULOMB_COLL
  USE BUILD_DIV_BOUNDARY
  USE MPI
  USE VAR_KIND_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build Maxwell-Ampere submatrix elements due to ion motion
  SUBROUTINE FILL_MAXWELL_ION(n_i_lagged)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i_lagged(NUM_X)

    INTEGER :: I, K

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    K = 1

  !Fill auxillary sparse matrix
    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        VAL = - MAX_E_J_0 * Z_PROF(I) * n_i_lagged(I)

        LOCAL_M%VALUE(MARKER_MAXWELL_ION(K)) = LOCAL_M%VALUE(MARKER_MAXWELL_ION(K)) + VAL

        K = K + 1

      END IF

    END DO

  END SUBROUTINE FILL_MAXWELL_ION
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build ion Lorentz force terms
  SUBROUTINE FILL_ION_LORENTZ

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    DO I = 1, ION_LORENTZ_SP%N_NZ

      VAL = Z_PROF(ION_LORENTZ_POS(I)) * EL_MASS / ION_MASS

      LOCAL_M%VALUE(MARKER_ION_LORENTZ(I)) = LOCAL_M%VALUE(MARKER_ION_LORENTZ(I)) + VAL

    END DO

  END SUBROUTINE FILL_ION_LORENTZ
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build ion convective term
  SUBROUTINE FILL_ION_CONV(u_i_lagged)

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_i_lagged(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL, &  !< Current value of matrix element to be passed to global matrix
                               loc_dx(MIN_X:MAX_X), &
                               loc_mult(-1:1,MIN_X:MAX_X)

    IF (ION_CONV_UPWINDING_SWITCH) THEN

      loc_mult = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0 ) THEN

          IF (u_i_lagged(I) .GE. 0) THEN

            loc_mult(1,I) = -1.00D00
            loc_mult(-1,I) = 0
            loc_mult(0,I) = 1.00D00
            loc_dx(I) = dxm(I)

          ELSE

            loc_mult(0,I) = - 1.00D00
            loc_mult(-1,I) = 1.00D00
            loc_mult(1,I) = 0
            loc_dx(I) = dxp(I)

          END IF

        END IF

      END DO

      DO I = 1, ION_CONV_SP%N_NZ

        VAL = - loc_mult(ION_CONV_PROP%SIGN(I),ION_CONV_PROP%POS(I))&
        * u_i_lagged(ION_CONV_PROP%POS(I))/loc_dx(ION_CONV_PROP%POS(I))

        LOCAL_M%VALUE(MARKER_ION_CONV(I)) = LOCAL_M%VALUE(MARKER_ION_CONV(I)) + VAL

      END DO

    ELSE

      DO I = 1, ION_CONV_SP%N_NZ

        VAL = ION_CONV_PROP%SIGN(I) * u_i_lagged(ION_CONV_PROP%POS(I)) / &
        (dxp(ION_CONV_PROP%POS(I)) + dxm(ION_CONV_PROP%POS(I)))

        LOCAL_M%VALUE(MARKER_ION_CONV(I)) = LOCAL_M%VALUE(MARKER_ION_CONV(I)) + VAL

      END DO

    END IF

  END SUBROUTINE FILL_ION_CONV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build ion convective term
  SUBROUTINE FILL_ION_CONT(n_i_lagged)

    IMPLICIT NONE

    INTEGER :: I, COL_X

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i_lagged(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    DO I = 1, ION_CONT_SP%N_NZ

      COL_X = (ION_CONT_SP%COL(I) - NUM_0D)/NUM_0D + 1

      VAL = ION_CONT_PROP%SIGN(I) * &
            n_i_lagged(COL_X) / dxc(ION_CONT_PROP%POS(I))

      LOCAL_M%VALUE(MARKER_ION_CONT(I)) = LOCAL_M%VALUE(MARKER_ION_CONT(I)) + VAL

    END DO

  END SUBROUTINE FILL_ION_CONT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build ion source term - ionization
  SUBROUTINE FILL_ION_GAIN

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=HIGH_PREC_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    DO I = 1, ION_GAIN_SP%N_NZ

      VAL = ION_RATES(ION_GAIN_N(I),ION_GAIN_POS(I))

      LOCAL_M%VALUE(MARKER_ION_GAIN(I)) = LOCAL_M%VALUE(MARKER_ION_GAIN(I)) + VAL

    END DO

  END SUBROUTINE FILL_ION_GAIN
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build ion loss term - volume recombination
  SUBROUTINE FILL_ION_LOSS(n_i, T_e)

    IMPLICIT NONE

    INTEGER :: I, K

    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::  n_i(NUM_X), & !< Lagged ion density
                                  T_e(NUM_X) !< Lagged electron temperature

    REAL(KIND=HIGH_PREC_REAL) :: VAL,& !< Current value of matrix element to be passed to global matrix
                    RAD_RATE(NUM_NEUTRALS,MIN_X:MAX_X), TB_RATE(NUM_NEUTRALS,MIN_X:MAX_X), &
                    R_TOT(MIN_X:MAX_X), TB_TOT(MIN_X:MAX_X)


    R_TOT = 0
    TB_TOT = 0
!Calculate recombination rates
    DO K = MIN_X,MAX_X

      DO I = 1, NUM_NEUTRALS

        RAD_RATE(I, K) = RECOMB_RATE(I, T_e(K))
        TB_RATE(I, K) = TB_RECOMB_RATES(I,K)

        R_TOT(K) = R_TOT(K) + RAD_RATE(I, K)
        TB_TOT(K) = TB_TOT(K) + TB_RATE(I, K)

      END DO

    END DO

    DO I = 1, ION_LOSS_SP%N_NZ

      VAL = (TB_REC_A_0 * TB_TOT(ION_RECOMB_POS(I))  + &
            RAD_REC_A_0 * R_TOT(ION_RECOMB_POS(I))) * n_i(ION_RECOMB_POS(I)) * &
            4.00D00 * PI * V_GRID(ION_RECOMB_V(I)) ** 2 * V_GRID_WIDTH(ION_RECOMB_V(I))

      LOCAL_M%VALUE(MARKER_ION_LOSS(I)) = LOCAL_M%VALUE(MARKER_ION_LOSS(I)) - VAL

    END DO

  END SUBROUTINE FILL_ION_LOSS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build ion loss term - volume recombination - full fluid case
  SUBROUTINE FILL_ION_LOSS_FF(n_i, T_e)

    IMPLICIT NONE

    INTEGER :: I, K

    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::  n_i(NUM_X), & !< Lagged ion density
                                  T_e(NUM_X) !< Lagged electron temperature

    REAL(KIND=HIGH_PREC_REAL) :: VAL,& !< Current value of matrix element to be passed to global matrix
                    RAD_RATE(NUM_NEUTRALS,MIN_X:MAX_X), TB_RATE(NUM_NEUTRALS,MIN_X:MAX_X), &
                    R_TOT(MIN_X:MAX_X), TB_TOT(MIN_X:MAX_X)


    R_TOT = 0
    TB_TOT = 0

!Calculate recombination rates
    DO K = MIN_X,MAX_X

      DO I = 1, NUM_NEUTRALS

        RAD_RATE(I, K) = RECOMB_RATE(I, T_e(K))
        TB_RATE(I, K) = TB_RECOMB_RATES(I,K)

        R_TOT(K) = R_TOT(K) + RAD_RATE(I, K)
        TB_TOT(K) = TB_TOT(K) + TB_RATE(I, K)

      END DO

    END DO

    DO I = 1, ION_LOSS_FF_SP%N_NZ

      VAL = (TB_REC_A_0 * TB_TOT(ION_RECOMB_FF_POS(I))  + &
            RAD_REC_A_0 * R_TOT(ION_RECOMB_FF_POS(I))) * n_i(ION_RECOMB_FF_POS(I))

      LOCAL_M%VALUE(MARKER_ION_LOSS_FF(I)) = LOCAL_M%VALUE(MARKER_ION_LOSS_FF(I)) - VAL

    END DO

  END SUBROUTINE FILL_ION_LOSS_FF
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build e-i drag term in ion momentum equation
  SUBROUTINE FILL_ION_DRAG(n_i)

    IMPLICIT NONE

    INTEGER :: I, K, P

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i(NUM_X)!< Lagged ion density

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    K = 1
    DO P = MIN_X,MAX_X

      IF (MOD(P,2) .EQ. 0) THEN

        DO I = 1,NUM_V+2

          VAL = - (4.00D00 * PI / 3.00D00) * EL_MASS/ION_MASS * EI_DRAG_TERMS(P)%VAL(I)/n_i(P)

          LOCAL_M%VALUE(MARKER_ION_DRAG(K)) = LOCAL_M%VALUE(MARKER_ION_DRAG(K)) + VAL

          K = K + 1

        END DO

      END IF

    END DO

  END SUBROUTINE FILL_ION_DRAG
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build divertor ion convection term
  SUBROUTINE FILL_ION_CONV_DIV(u_i,n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_i(NUM_X), & !< Lagged ion velociy
                                           n_i(NUM_X) !< Lagged ion density

    REAL(KIND=DEFAULT_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                               loc_mult, &
                               loc_dx

    IF (ION_CONV_UPWINDING_SWITCH) THEN

      IF (u_i(NUM_X-1) .GT. 0) THEN

        loc_mult = 0
        loc_dx = dxm(NUM_X-1)

      ELSE

        loc_mult = 1.00D00
        loc_dx = dxc(NUM_X)

      END IF

    ELSE

      loc_mult = 1.00D00
      loc_dx = dxm(NUM_X-1) + dxc(NUM_X)

    END IF

    IF (SONIC_OUTFLOW_DIV_SWITCH .AND. BOHM_VALUE) THEN

      IF (ION_EL_TEMP_SWITCH) THEN

        VAL = - MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(ION_MASS))

      ELSE

        VAL = - MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(2.00D00*ION_MASS))

      END IF

      LOCAL_M%VALUE(MARKER_ION_CONV_DIV(1)) = LOCAL_M%VALUE(MARKER_ION_CONV_DIV(1)) + loc_mult*VAL/loc_dx

    ELSE

      VAL = - u_i(NUM_X - 1)*n_i(NUM_X-1)/n_e_boundary* &
      (1.00D00 + (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1))) / (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)))

      LOCAL_M%VALUE(MARKER_ION_CONV_DIV(1)) = &
      LOCAL_M%VALUE(MARKER_ION_CONV_DIV(1)) + FLOW_LIMITER*loc_mult*VAL/loc_dx

      VAL =  u_i(NUM_X - 1)*n_i(NUM_X-3)/n_e_boundary* &
      (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1))) / (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3))

      LOCAL_M%VALUE(MARKER_ION_CONV_DIV(2)) = &
      LOCAL_M%VALUE(MARKER_ION_CONV_DIV(2)) + FLOW_LIMITER*loc_mult*VAL/loc_dx

    END IF

  END SUBROUTINE FILL_ION_CONV_DIV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build divertor ion sink term
  SUBROUTINE FILL_ION_SINK(n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i(NUM_X)!< Lagged ion density

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    IF (SONIC_OUTFLOW_DIV_SWITCH .AND. BOHM_VALUE) THEN

      VAL = - ion_flux / n_i(NUM_X)
      LOCAL_M%VALUE(MARKER_ION_SINK(1)) = LOCAL_M%VALUE(MARKER_ION_SINK(1)) + VAL/ dxc(NUM_X)

    ELSE

      VAL = - n_i(NUM_X - 1)* (1.00D00 + (dxc(NUM_X) / dxm(NUM_X-1)))

      LOCAL_M%VALUE(MARKER_ION_SINK(2)) = LOCAL_M%VALUE(MARKER_ION_SINK(2)) + FLOW_LIMITER*VAL/dxc(NUM_X)

      VAL =  n_i(NUM_X - 3)* (dxc(NUM_X) / dxm(NUM_X-1))

      LOCAL_M%VALUE(MARKER_ION_SINK(3)) = LOCAL_M%VALUE(MARKER_ION_SINK(3)) + FLOW_LIMITER*VAL/dxc(NUM_X)

    END IF

  END SUBROUTINE FILL_ION_SINK
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build ion source term in momentum equation
  SUBROUTINE FILL_ION_MOM_SOURCE(n_i,n_e,T_e,f0)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN), DIMENSION(NUM_X) :: n_i, &!< Lagged ion density
                                                  T_e, & !< Lagged electron temperature
                                                  n_e!< Lagged electron density

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f0(DIM_F) !< Lagged distribution function

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    REAL(KIND=DEFAULT_REAL), DIMENSION(MIN_X:MAX_X) :: TOTAL_GAIN, TOTAL_LOSS

    INTEGER :: I, J, K

    TOTAL_GAIN = 0
    TOTAL_LOSS = 0

    IF (COLL_EN_ION) THEN

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) THEN

          DO J = 1, NUM_NEUTRALS

            TOTAL_GAIN(I) = TOTAL_GAIN(I) + &
            ION_RATES(J,I) * f0(X_POS(I) + NUM_H * NUM_V + J)

          END DO

        END IF

      END DO

    END IF

    IF (COLL_RECOMB) THEN

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) THEN

          DO J = 1, NUM_NEUTRALS

            TOTAL_LOSS(I) = TOTAL_LOSS(I) + RECOMB_RATE(J, T_e(I)) * RAD_REC_A_0 * n_e(I) * n_i(I)

            TOTAL_LOSS(I) = TOTAL_LOSS(I) + &
            TB_RECOMB_RATES(J,I) * TB_REC_A_0 * n_e(I) * n_i(I)

          END DO

        END IF

      END DO

    END IF

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        VAL = - (TOTAL_GAIN(I) - TOTAL_LOSS(I)) / n_i(I)

        LOCAL_M%VALUE(MARKER_ION_MOM_SOURCE(K)) = LOCAL_M%VALUE(MARKER_ION_MOM_SOURCE(K)) + VAL

        K = K + 1

      END IF

    END DO

  END SUBROUTINE FILL_ION_MOM_SOURCE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build pressure gradient term in ion momentum equation
  SUBROUTINE FILL_ION_PRESSURE(n_i,T_e)

    IMPLICIT NONE

    INTEGER :: K, P, I

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i(NUM_X), &!< Lagged ion density
                                T_e(NUM_X)!< Lagged electron temperature

    REAL(KIND=DEFAULT_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    K = 1

    IF ((.NOT. ION_CONT_OFF_SWITCH) .OR. (FULL_FLUID_MODE)) THEN

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 0) THEN

          VAL = ION_PRESSURE_SIGN(K) * 0.50D00 * EL_MASS/ION_MASS * T_e(ION_PRESSURE_POS(K))/n_i(P)

          LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) = LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) + VAL/dxc(P)

          K = K + 1

          VAL = ION_PRESSURE_SIGN(K) * 0.50D00 * EL_MASS/ION_MASS * T_e(ION_PRESSURE_POS(K))/n_i(P)

          LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) = LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) + VAL/dxc(P)

          K = K + 1

        END IF

      END DO

    ELSE

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 0) THEN

          DO I = 1, NUM_V

            VAL = ION_PRESSURE_SIGN(K) * 2.00D00 * PI * EL_MASS/ION_MASS &
            * T_e(ION_PRESSURE_POS(K))/n_i(P) * V_GRID(I) ** 2 * V_GRID_WIDTH(I)

            LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) = LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) + VAL/dxc(P)

            K = K + 1

            VAL = ION_PRESSURE_SIGN(K) * 2.00D00 * PI * EL_MASS/ION_MASS &
            * T_e(ION_PRESSURE_POS(K))/n_i(P) * V_GRID(I) ** 2 * V_GRID_WIDTH(I)

            LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) = LOCAL_M%VALUE(MARKER_ION_PRESSURE(K)) + VAL/dxc(P)

            K = K + 1

          END DO

        END IF

      END DO

    END IF

  END SUBROUTINE FILL_ION_PRESSURE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build simple CX drag term in ion momentum equation
  SUBROUTINE FILL_SIMPLE_CX(f,u_i)

    IMPLICIT NONE

    INTEGER :: I, K, P

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(DIM_F), & !< Lagged unknown vector
                                           u_i(NUM_X) !< Lagged ion velocity

    REAL(KIND=DEFAULT_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                               SIMPLE_CX_RATE(MIN_X:MAX_X)

    SIMPLE_CX_RATE = 0

    DO P = MIN_X, MAX_X

      IF (MOD(P,2) .EQ. 0) THEN

        SIMPLE_CX_RATE(P) = SIMPLE_CX_RATE(P) + f(X_POS(P)+NUM_V*NUM_H + 1) * ABS(u_i(P)) &
                           * 3.00D-19/SIGMA_0

        IF (NUM_NEUTRALS .GT. 1)&
        SIMPLE_CX_RATE(P) = SIMPLE_CX_RATE(P) + f(X_POS(P)+NUM_V*NUM_H + 2) * ABS(u_i(P)) &
                           * 2**4 * 1.00D-19/SIGMA_0

        IF (NUM_NEUTRALS .GT. 2)&
        SIMPLE_CX_RATE(P) = SIMPLE_CX_RATE(P) + f(X_POS(P)+NUM_V*NUM_H + 3) * ABS(u_i(P)) &
                           * 3**4 * 7.00D-20/SIGMA_0

        DO I = 4, NUM_NEUTRALS

          SIMPLE_CX_RATE(P) = SIMPLE_CX_RATE(P) + f(X_POS(P)+NUM_V*NUM_H + I) * ABS(u_i(P)) &
                            * I ** 4 * 6.00D-20/SIGMA_0

        END DO

      END IF

    END DO

    K = 1

    DO P = MIN_X,MAX_X

      IF (MOD(P,2) .EQ. 0) THEN

          VAL = - SIMPLE_CX_RATE(P)

          LOCAL_M%VALUE(MARKER_SIMPLE_CX(K)) = LOCAL_M%VALUE(MARKER_SIMPLE_CX(K)) + COLL_EN_L*VAL

          K = K + 1

      END IF

    END DO

  END SUBROUTINE FILL_SIMPLE_CX
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build thermal drag force term for ions
  SUBROUTINE FILL_ION_T_DRAG(n_e,n_i)

    IMPLICIT NONE

    INTEGER :: K, P

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i(NUM_X), & !< Lagged ion density
                                           n_e(NUM_X) !< Lagged electron density

    REAL(KIND=DEFAULT_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    K = 1

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 0) THEN

          VAL = - ION_T_DRAG_SIGN(K) * 0.50D00  * 0.71D00 * EL_MASS/ION_MASS * n_e(P)/n_i(P)

          LOCAL_M%VALUE(MARKER_ION_T_DRAG(K)) = LOCAL_M%VALUE(MARKER_ION_T_DRAG(K)) + VAL/dxc(P)

          K = K + 1

          VAL = - ION_T_DRAG_SIGN(K) * 0.50D00  * 0.71D00 * EL_MASS/ION_MASS * n_e(P)/n_i(P)

          LOCAL_M%VALUE(MARKER_ION_T_DRAG(K)) = LOCAL_M%VALUE(MARKER_ION_T_DRAG(K)) + VAL/dxc(P)

          K = K + 1

        END IF

      END DO

  END SUBROUTINE FILL_ION_T_DRAG
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build e-i drag force term for ions
  SUBROUTINE FILL_ION_U_DRAG(n_o,T_o,n_e,T_e,n_i)

    IMPLICIT NONE

    INTEGER :: K, P

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), &!< Lagged electron density
                                           n_i(NUM_X), &!< Lagged ion density
                                           T_e(NUM_X) , &!< Lagged electron temperature
                                           n_o(NUM_X), &!< Old electron density
                                           T_o(NUM_X) !< Old electron temperature

    REAL(KIND=DEFAULT_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    K = 1

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 0) THEN

          VAL = - EL_MASS/ION_MASS * ION_U_DRAG_SIGN(K) &
                * DENSITY_0 * TIME_NORM/(TEMP_0_EV**(3.00D00/2.00D00)) * 1.48D-12 &
               * LAMBDA_EI(n_o(P),T_o(P),Z_PROF(P)) * Z_PROF(P)/ION_Z * n_e(P) ** 2/(n_i(P)*T_e(P)**(3.00D00/2.00D00))

          LOCAL_M%VALUE(MARKER_ION_U_DRAG(K)) = LOCAL_M%VALUE(MARKER_ION_U_DRAG(K)) + VAL

          K = K + 1

          VAL = - EL_MASS/ION_MASS * ION_U_DRAG_SIGN(K) &
              * DENSITY_0 * TIME_NORM/(TEMP_0_EV**(3.00D00/2.00D00)) * 1.48D-12 &
               * LAMBDA_EI(n_o(P),T_o(P),Z_PROF(P)) * Z_PROF(P)/ION_Z * n_e(P) ** 2/(n_i(P)*T_e(P)**(3.00D00/2.00D00))

          LOCAL_M%VALUE(MARKER_ION_U_DRAG(K)) = LOCAL_M%VALUE(MARKER_ION_U_DRAG(K)) + VAL

          K = K + 1

        END IF

      END DO

  END SUBROUTINE FILL_ION_U_DRAG
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_ION_FLUID
!-------------------------------------------------------------------------------------------------------------------------------------
