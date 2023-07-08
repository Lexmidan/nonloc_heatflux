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
!> Contains submatrix building routines for fluid electron operators
MODULE BUILD_FLUID_EL

  USE GRID
  USE NORMALIZATION
  USE SWITCHES
  USE MATRIX_DATA
  USE COLL_CS_RATES
  USE NEUT_AND_HEAT
  USE BUILD_DIV_BOUNDARY
  USE MPI
  USE VAR_KIND_DATA

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: R_ex(:,:), & !< Contains excitation drag matrix elements
                                          R_deex(:,:), & !< Contains deexcitation drag matrix elements
                                          R_ion(:,:), & !< Contains ionization drag matrix elements
                                          R_recomb(:), & !< Contains recombination drag matrix elements
                                          R_el(:,:) !< Contains elastic matrix elements

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build Maxwell-Ampere submatrix elements due to electron motion
  SUBROUTINE FILL_MAXWELL_EL(n_e_lagged)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e_lagged(NUM_X)

    INTEGER :: I, K

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    K = 1

  !Fill auxillary sparse matrix
    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        VAL =  MAX_E_J_0 * n_e_lagged(I)

        LOCAL_M%VALUE(MARKER_MAXWELL_EL(K)) = LOCAL_M%VALUE(MARKER_MAXWELL_EL(K)) + VAL

        K = K + 1

      END IF

    END DO

  END SUBROUTINE FILL_MAXWELL_EL
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron Lorentz force terms
  SUBROUTINE FILL_EL_LORENTZ

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    VAL = - 1.00D00

    DO I = 1, EL_LORENTZ_SP%N_NZ

      LOCAL_M%VALUE(MARKER_EL_LORENTZ(I)) = LOCAL_M%VALUE(MARKER_EL_LORENTZ(I)) + VAL

    END DO

  END SUBROUTINE FILL_EL_LORENTZ
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron convective term
  SUBROUTINE FILL_EL_CONV(u_e_lagged)

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_e_lagged(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL, &  !< Current value of matrix element to be passed to global matrix
                               loc_dx(MIN_X:MAX_X), &
                               loc_mult(-1:1,MIN_X:MAX_X)

    IF (ION_CONV_UPWINDING_SWITCH) THEN

      loc_mult = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0 ) THEN

          IF (u_e_lagged(I) .GE. 0) THEN

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

      DO I = 1, EL_CONV_SP%N_NZ

        VAL = - loc_mult(EL_CONV_PROP%SIGN(I),EL_CONV_PROP%POS(I))&
        * u_e_lagged(EL_CONV_PROP%POS(I))/loc_dx(EL_CONV_PROP%POS(I))

        LOCAL_M%VALUE(MARKER_EL_CONV(I)) = LOCAL_M%VALUE(MARKER_EL_CONV(I)) + VAL

      END DO

    ELSE

      DO I = 1, EL_CONV_SP%N_NZ

        VAL = EL_CONV_PROP%SIGN(I) * u_e_lagged(EL_CONV_PROP%POS(I)) / &
        (dxp(EL_CONV_PROP%POS(I)) + dxm(EL_CONV_PROP%POS(I)))

        LOCAL_M%VALUE(MARKER_EL_CONV(I)) = LOCAL_M%VALUE(MARKER_EL_CONV(I)) + VAL

      END DO

    END IF

  END SUBROUTINE FILL_EL_CONV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron convective term
  SUBROUTINE FILL_EL_CONT(n_e_lagged)

    IMPLICIT NONE

    INTEGER :: I, COL_X

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e_lagged(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    DO I = 1, EL_CONT_SP%N_NZ

      COL_X = (EL_CONT_SP%COL(I) - (NUM_0D-3))/NUM_0D + 1

      VAL = EL_CONT_PROP%SIGN(I) * &
            n_e_lagged(COL_X) / dxc(EL_CONT_PROP%POS(I))

      LOCAL_M%VALUE(MARKER_EL_CONT(I)) = LOCAL_M%VALUE(MARKER_EL_CONT(I)) + VAL

    END DO

  END SUBROUTINE FILL_EL_CONT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron source term - ionization
  SUBROUTINE FILL_EL_GAIN

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix


    DO I = 1, EL_GAIN_SP%N_NZ

      VAL = ION_RATES(EL_GAIN_N(I),EL_GAIN_POS(I))

      LOCAL_M%VALUE(MARKER_EL_GAIN(I)) = LOCAL_M%VALUE(MARKER_EL_GAIN(I)) + VAL

    END DO

  END SUBROUTINE FILL_EL_GAIN
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build el loss term - volume recombination
  SUBROUTINE FILL_EL_LOSS(n_i, T_e)

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

    DO I = 1, EL_LOSS_SP%N_NZ

      VAL = (TB_REC_A_0 * TB_TOT(EL_RECOMB_POS(I))  + &
            RAD_REC_A_0 * R_TOT(EL_RECOMB_POS(I))) * n_i(EL_RECOMB_POS(I))

      LOCAL_M%VALUE(MARKER_EL_LOSS(I)) = LOCAL_M%VALUE(MARKER_EL_LOSS(I)) - VAL

    END DO

  END SUBROUTINE FILL_EL_LOSS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build divertor electron convection term
  SUBROUTINE FILL_EL_CONV_DIV(u_e,n_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_e(NUM_X), & !< Lagged electron velociy
                                           n_e(NUM_X) !< Lagged electron density

    REAL(KIND=DEFAULT_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                               loc_mult, &
                               loc_dx

    IF (ION_CONV_UPWINDING_SWITCH) THEN

      IF (u_e(NUM_X-1) .GT. 0) THEN

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

      LOCAL_M%VALUE(MARKER_EL_CONV_DIV(1)) = LOCAL_M%VALUE(MARKER_EL_CONV_DIV(1)) + loc_mult*VAL/loc_dx

    ELSE

      VAL = - u_e(NUM_X - 1)*n_e(NUM_X-1)/n_e_boundary* &
      (1.00D00 + (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1))) / (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)))

      LOCAL_M%VALUE(MARKER_EL_CONV_DIV(1)) = &
      LOCAL_M%VALUE(MARKER_EL_CONV_DIV(1)) + FLOW_LIMITER*loc_mult*VAL/loc_dx

      VAL =  u_e(NUM_X - 1)*n_e(NUM_X-3)/n_e_boundary* &
      (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1))) / (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3))

      LOCAL_M%VALUE(MARKER_EL_CONV_DIV(2)) = &
      LOCAL_M%VALUE(MARKER_EL_CONV_DIV(2)) + FLOW_LIMITER*loc_mult*VAL/loc_dx

    END IF

  END SUBROUTINE FILL_EL_CONV_DIV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build divertor electron sink term
  SUBROUTINE FILL_EL_SINK(n_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X)!< Lagged electron density

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    IF (SONIC_OUTFLOW_DIV_SWITCH .AND. BOHM_VALUE) THEN

      VAL = - ion_flux / n_e(NUM_X)

      LOCAL_M%VALUE(MARKER_EL_SINK(1)) = LOCAL_M%VALUE(MARKER_EL_SINK(1)) + VAL/ dxc(NUM_X)

    ELSE

      VAL = - n_e(NUM_X - 1)* (1.00D00 + (dxc(NUM_X) / dxm(NUM_X-1)))

      LOCAL_M%VALUE(MARKER_EL_SINK(2)) = LOCAL_M%VALUE(MARKER_EL_SINK(2)) + FLOW_LIMITER*VAL/dxc(NUM_X)

      VAL =  n_e(NUM_X - 3)* (dxc(NUM_X) / dxm(NUM_X-1))

      LOCAL_M%VALUE(MARKER_EL_SINK(3)) = LOCAL_M%VALUE(MARKER_EL_SINK(3)) + FLOW_LIMITER*VAL/dxc(NUM_X)

    END IF

  END SUBROUTINE FILL_EL_SINK
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron source term in momentum equation
  SUBROUTINE FILL_EL_MOM_SOURCE(n_i,n_e,T_e,f0)

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

        VAL = - (TOTAL_GAIN(I) - TOTAL_LOSS(I)) / n_e(I)

        LOCAL_M%VALUE(MARKER_EL_MOM_SOURCE(K)) = LOCAL_M%VALUE(MARKER_EL_MOM_SOURCE(K)) + VAL

        K = K + 1

      END IF

    END DO

  END SUBROUTINE FILL_EL_MOM_SOURCE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build pressure gradient term in electron momentum equation
  SUBROUTINE FILL_EL_PRESSURE(n_e,T_e)

    IMPLICIT NONE

    INTEGER :: K, P

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), &!< Lagged ion density
                                T_e(NUM_X)!< Lagged electron temperature

    REAL(KIND=DEFAULT_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    K = 1

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 0) THEN

          VAL = EL_PRESSURE_SIGN(K) * 0.50D00  * T_e(EL_PRESSURE_POS(K))/n_e(P)

          LOCAL_M%VALUE(MARKER_EL_PRESSURE(K)) = LOCAL_M%VALUE(MARKER_EL_PRESSURE(K)) + VAL/dxc(P)

          K = K + 1

          VAL = EL_PRESSURE_SIGN(K) * 0.50D00  * T_e(EL_PRESSURE_POS(K))/n_e(P)

          LOCAL_M%VALUE(MARKER_EL_PRESSURE(K)) = LOCAL_M%VALUE(MARKER_EL_PRESSURE(K)) + VAL/dxc(P)

          K = K + 1

        END IF

      END DO

  END SUBROUTINE FILL_EL_PRESSURE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build thermal drag force term for electrons
  SUBROUTINE FILL_EL_T_DRAG

    IMPLICIT NONE

    INTEGER :: K, P

    REAL(KIND=DEFAULT_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    K = 1

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 0) THEN

          VAL = EL_T_DRAG_SIGN(K) * 0.50D00  * 0.71D00

          LOCAL_M%VALUE(MARKER_EL_T_DRAG(K)) = LOCAL_M%VALUE(MARKER_EL_T_DRAG(K)) + VAL/dxc(P)

          K = K + 1

          VAL = EL_T_DRAG_SIGN(K) * 0.50D00  * 0.71D00

          LOCAL_M%VALUE(MARKER_EL_T_DRAG(K)) = LOCAL_M%VALUE(MARKER_EL_T_DRAG(K)) + VAL/dxc(P)

          K = K + 1

        END IF

      END DO

  END SUBROUTINE FILL_EL_T_DRAG
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build e-i drag force term for electrons
  SUBROUTINE FILL_EL_U_DRAG(n_o,T_o,n_e,T_e)

    IMPLICIT NONE

    INTEGER :: K, P

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), &!< Lagged electron density
                                           T_e(NUM_X) , &!< Lagged electron temperature
                                           n_o(NUM_X), &!< Old electron density
                                           T_o(NUM_X) !< Old electron temperature

    REAL(KIND=DEFAULT_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    K = 1

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 0) THEN

          VAL = EL_U_DRAG_SIGN(K) * DENSITY_0 * TIME_NORM/(TEMP_0_EV**(3.00D00/2.00D00)) * 1.48D-12 &
               * LAMBDA_EI(n_o(P),T_o(P),Z_PROF(P)) * Z_PROF(P)/ION_Z * n_e(P)/T_e(P)**(3.00D00/2.00D00)

          LOCAL_M%VALUE(MARKER_EL_U_DRAG(K)) = LOCAL_M%VALUE(MARKER_EL_U_DRAG(K)) + VAL

          K = K + 1

          VAL = EL_U_DRAG_SIGN(K) * DENSITY_0 * TIME_NORM/(TEMP_0_EV**(3.00D00/2.00D00)) * 1.48D-12 &
               * LAMBDA_EI(n_o(P),T_o(P),Z_PROF(P)) * Z_PROF(P)/ION_Z * n_e(P)/T_e(P)**(3.00D00/2.00D00)

          LOCAL_M%VALUE(MARKER_EL_U_DRAG(K)) = LOCAL_M%VALUE(MARKER_EL_U_DRAG(K)) + VAL

          K = K + 1

        END IF

      END DO

  END SUBROUTINE FILL_EL_U_DRAG
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron temperature convective term
  SUBROUTINE FILL_EL_TEMP_CONV(u_e_lagged)

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_e_lagged(NUM_X)

    INTEGER :: u_sign(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL, &  !< Current value of matrix element to be passed to global matrix
                               loc_dx(MIN_X:MAX_X), &
                               loc_mult(-1:1,MIN_X:MAX_X)

      u_sign = 0
      loc_mult = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 1 ) THEN

          IF (u_e_lagged(I) .GE. 0) THEN

            u_sign(I) = -1
            loc_mult(1,I) = -1.00D00
            loc_mult(-1,I) = 0
            loc_mult(0,I) = 1.00D00
            loc_dx(I) = dxm(I)

          ELSE

            u_sign(I) = -1
            loc_mult(0,I) = - 1.00D00
            loc_mult(-1,I) = 1.00D00
            loc_mult(1,I) = 0
            loc_dx(I) = dxp(I)

          END IF

        END IF

      END DO

      IF ((u_e_lagged(1) .GE. 0) .AND. NO_FLOW_BOUNDARY_UP_SWITCH) u_sign(1) = 0
      IF ((u_e_lagged(NUM_X) .LT. 0) .AND. NO_FLOW_BOUNDARY_DIV_SWITCH) THEN

        u_sign(NUM_X) = 0

        IF (PLASMA_SINK_SWITCH) PRINT*,'WARNING: Velocity in last cell centre negative!'

      END IF

      DO I = 1, EL_CONV_TEMP_SP%N_NZ

        VAL = u_sign(EL_CONV_TEMP_PROP%POS(I))* loc_mult(EL_CONV_TEMP_PROP%SIGN(I),EL_CONV_TEMP_PROP%POS(I))&
        * u_e_lagged(EL_CONV_TEMP_PROP%POS(I))/loc_dx(EL_CONV_TEMP_PROP%POS(I))

        LOCAL_M%VALUE(MARKER_EL_CONV_TEMP(I)) = LOCAL_M%VALUE(MARKER_EL_CONV_TEMP(I)) + VAL

      END DO

  END SUBROUTINE FILL_EL_TEMP_CONV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron pdV temperature term
  SUBROUTINE FILL_EL_TEMP_PDV(T_e_lagged)

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T_e_lagged(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    DO I = 1, EL_TEMP_PDV_SP%N_NZ

      VAL = EL_TEMP_PDV_PROP%SIGN(I) * 2.00D00/3.00D00 * &
            T_e_lagged(EL_TEMP_PDV_PROP%POS(I)) / dxc(EL_TEMP_PDV_PROP%POS(I))

      LOCAL_M%VALUE(MARKER_EL_TEMP_PDV(I)) = LOCAL_M%VALUE(MARKER_EL_TEMP_PDV(I)) + VAL

    END DO

  END SUBROUTINE FILL_EL_TEMP_PDV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build divertor electron temperature pdv term
  SUBROUTINE FILL_EL_TEMP_PDV_DIV(u_e,n_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_e(NUM_X), & !< Lagged electron velociy
                                           n_e(NUM_X)!< Lagged electron density

    REAL(KIND=DEFAULT_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                               loc_dx

    loc_dx = dxc(NUM_X)


    IF (SONIC_OUTFLOW_DIV_SWITCH .AND. BOHM_VALUE) THEN

      IF (ION_EL_TEMP_SWITCH) THEN

        VAL = - MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(ION_MASS))

      ELSE

        VAL = - MACH_N_DIV*SQRT(T_e_boundary * EL_MASS/(2.00D00*ION_MASS))

      END IF

      LOCAL_M%VALUE(MARKER_EL_TEMP_PDV_DIV(1)) = LOCAL_M%VALUE(MARKER_EL_TEMP_PDV_DIV(1)) + (2.00D00/3.00D00)*VAL/loc_dx

    ELSE

      VAL = - u_e(NUM_X - 1)*n_e(NUM_X-1)/n_e_boundary* &
      (1.00D00 + (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1))) / (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)))


      VAL = VAL + u_e(NUM_X - 3)*n_e(NUM_X-3)/n_e_boundary* &
      (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1))) / (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3))

      LOCAL_M%VALUE(MARKER_EL_TEMP_PDV_DIV(1)) = &
      LOCAL_M%VALUE(MARKER_EL_TEMP_PDV_DIV(1)) + (2.00D00/3.00D00)*FLOW_LIMITER*VAL/loc_dx

    END IF

  END SUBROUTINE FILL_EL_TEMP_PDV_DIV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron temperature div q_t term
  SUBROUTINE FILL_EL_TEMP_Q_T(T_e_lagged,n_e_lagged)

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T_e_lagged(NUM_X), &
                                           n_e_lagged(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL, &  !< Current value of matrix element to be passed to global matrix
                               loc_mult(-1:1,MIN_X:MAX_X), &
                               KAPPA, log_ratio(NUM_X)

      loc_mult = 0

      DO I = 1, NUM_X

        log_ratio(I) = LAMBDA_EI(1.00D00,1.00D00,ION_Z)&
        /LAMBDA_EI(n_e_lagged(I),T_e_lagged(I),Z_PROF(I))

      END DO


      IF (.NOT. COLL_EE_L_SWITCH) THEN

        KAPPA = 13.581 * 3.00D00 * SQRT(PI)/16.00D00                                          !Lorentz case

      ELSE

        KAPPA = 3.20 * 3.00D00 * SQRT(PI)/16.00D00                                            !SH case with Epperlein coeff

      END IF

        DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 1 ) THEN

            IF (I .EQ. 1) THEN

              IF (PERIODIC_BOUNDARY_SWITCH) THEN

                loc_mult(0,I) = -T_e_lagged(I + 1)**(5.00D00/2.00D00)&
                                 *log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I) - &
                                  T_e_lagged(NUM_X)**(5.00D00/2.00D00)&
                                  *log_ratio(NUM_X)*ION_Z/Z_PROF(NUM_X)/dxm(I)
                loc_mult(-1,I) =  T_e_lagged(I + 1)**(5.00D00/2.00D00)*log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I)
                loc_mult(1,I) = T_e_lagged(NUM_X)**(5.00D00/2.00D00)*log_ratio(NUM_X)*ION_Z/Z_PROF(NUM_X)/dxm(I)

              ELSE

                loc_mult(0,I) = - T_e_lagged(I + 1)**(5.00D00/2.00D00)*log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I)
                loc_mult(-1,I) = T_e_lagged(I + 1)**(5.00D00/2.00D00)*log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I)
                loc_mult(1,I) = 0.00D00

              END IF

            ELSE IF (I .GT. NUM_X - 2 ) THEN

              IF (PERIODIC_BOUNDARY_SWITCH) THEN

                loc_mult(0,I) = - T_e_lagged(I + 1)**(5.00D00/2.00D00)*log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I) - &
                                  T_e_lagged(I - 1)**(5.00D00/2.00D00)*log_ratio(I-1)*ION_Z/Z_PROF(I-1)/dxm(I)
                loc_mult(-1,I) = T_e_lagged(I + 1)**(5.00D00/2.00D00)*log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I)
                loc_mult(1,I) = T_e_lagged(I - 1)**(5.00D00/2.00D00)*log_ratio(I-1)*ION_Z/Z_PROF(I-1)/dxm(I)

              ELSE

                loc_mult(0,I) = - T_e_lagged(I - 1)**(5.00D00/2.00D00)*log_ratio(I-1)*ION_Z/Z_PROF(I-1)/dxm(I)
                loc_mult(-1,I) = 0.00D00
                loc_mult(1,I) = T_e_lagged(I - 1)**(5.00D00/2.00D00)*log_ratio(I-1)*ION_Z/Z_PROF(I-1)/dxm(I)

              END IF

            ELSE

              loc_mult(0,I) = - T_e_lagged(I + 1)**(5.00D00/2.00D00)*log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I) - &
                              T_e_lagged(I - 1)**(5.00D00/2.00D00)*log_ratio(I-1)*ION_Z/Z_PROF(I-1)/dxm(I)
              loc_mult(-1,I) = T_e_lagged(I + 1)**(5.00D00/2.00D00)*log_ratio(I+1)*ION_Z/Z_PROF(I+1)/dxp(I)
              loc_mult(1,I) =  T_e_lagged(I - 1)**(5.00D00/2.00D00)*log_ratio(I-1)*ION_Z/Z_PROF(I-1)/dxm(I)

            END IF

          END IF

        END DO

      DO I = 1, EL_CONV_TEMP_SP%N_NZ

        VAL = 4.00D00*KAPPA&
              *loc_mult(EL_CONV_TEMP_PROP%SIGN(I),EL_CONV_TEMP_PROP%POS(I))&
             /(3.00D00 * dxc(EL_CONV_TEMP_PROP%POS(I)) * n_e_lagged(EL_CONV_TEMP_PROP%POS(I)))

        LOCAL_M%VALUE(MARKER_EL_CONV_TEMP(I)) = LOCAL_M%VALUE(MARKER_EL_CONV_TEMP(I))   + VAL

      END DO

  END SUBROUTINE FILL_EL_TEMP_Q_T
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron div q_u term
  SUBROUTINE FILL_EL_TEMP_Q_U(T_e_lagged,n_e_lagged)

    IMPLICIT NONE

    INTEGER :: I, COL_X

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T_e_lagged(NUM_X), n_e_lagged(NUM_X)

    REAL(KIND=DEFAULT_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    DO I = 1, EL_TEMP_Q_U_SP%N_NZ

      COL_X = (EL_TEMP_Q_U_SP%COL(I) - (NUM_0D-3))/NUM_0D + 1

      VAL = EL_TEMP_Q_U_PROP%SIGN(I) * 2.00D00/3.00D00 * 0.71D00 *&
            n_e_lagged(COL_X)*T_e_lagged(COL_X) / &
            (n_e_lagged(EL_TEMP_Q_U_PROP%POS(I)) * dxc(EL_TEMP_Q_U_PROP%POS(I)))

      LOCAL_M%VALUE(MARKER_EL_TEMP_Q_U(I)) = LOCAL_M%VALUE(MARKER_EL_TEMP_Q_U(I)) + VAL

    END DO

  END SUBROUTINE FILL_EL_TEMP_Q_U
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build divertor electron temperature div q term
  SUBROUTINE FILL_EL_TEMP_Q_DIV(n_e,mult)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X) !< Lagged electron density
    REAL(KIND=DEFAULT_REAL), INTENT(IN), OPTIONAL :: mult

    REAL(KIND=DEFAULT_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                               loc_dx, &
                               gamma_e, &
                               temp_mult

    loc_dx = dxc(NUM_X)
    temp_mult = 0
    IF (ION_EL_TEMP_SWITCH) temp_mult = 1.00D00

    gamma_e = 2.00D00 - 0.5D00*LOG(2*PI*(1.00D00+temp_mult)*EL_MASS/ION_MASS)

    IF (PRESENT(mult)) gamma_e = mult

    VAL = - (gamma_e - 2.5D00)/n_e(NUM_X) * ion_flux

    LOCAL_M%VALUE(MARKER_EL_TEMP_PDV_DIV(1)) = LOCAL_M%VALUE(MARKER_EL_TEMP_PDV_DIV(1)) + (2.00D00/3.00D00)*VAL/loc_dx

  END SUBROUTINE FILL_EL_TEMP_Q_DIV
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron source temp term - ionization
  SUBROUTINE FILL_EL_GAIN_TEMP(n_e,u_e,T_e)

    IMPLICIT NONE
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), & !< Lagged electron density
                                           u_e(NUM_X), & !< Lagged electron velocity
                                           T_e(NUM_X) !< Old electron temperature (using old for conservation properties)
    INTEGER :: I

    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    DO I = 1, EL_GAIN_TEMP_SP%N_NZ

      VAL = - ION_RATES(EL_GAIN_TEMP_N(I),EL_GAIN_TEMP_POS(I)) * &
            (T_e(EL_GAIN_TEMP_POS(I))-2.00D00*u_e(EL_GAIN_TEMP_POS(I))**2/3.00D00)&
            /n_e(EL_GAIN_TEMP_POS(I))

      LOCAL_M%VALUE(MARKER_EL_GAIN_TEMP(I)) = LOCAL_M%VALUE(MARKER_EL_GAIN_TEMP(I)) + VAL

    END DO

  END SUBROUTINE FILL_EL_GAIN_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build el loss temp term - volume recombination
  SUBROUTINE FILL_EL_LOSS_TEMP(n_i, T_e, n_e, u_e)

    IMPLICIT NONE

    INTEGER :: I, K
    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::  n_i(NUM_X), & !< Lagged ion density
                                            T_e(NUM_X), & !< Old electron temperature
                                            n_e(NUM_X), & !< Lagged electron density
                                            u_e(NUM_X) !< Lagged electron velocity

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

    DO I = 1, EL_LOSS_TEMP_SP%N_NZ

      VAL = - (TB_REC_A_0 * TB_TOT(EL_RECOMB_TEMP_POS(I))  + &
            RAD_REC_A_0 * R_TOT(EL_RECOMB_TEMP_POS(I))) * n_i(EL_RECOMB_TEMP_POS(I))* &
            (T_e(EL_RECOMB_TEMP_POS(I))-2.00D00*u_e(EL_RECOMB_TEMP_POS(I))**2/3.00D00)&
            /n_e(EL_RECOMB_TEMP_POS(I))

      LOCAL_M%VALUE(MARKER_EL_LOSS_TEMP(I)) = LOCAL_M%VALUE(MARKER_EL_LOSS_TEMP(I)) - VAL

    END DO

  END SUBROUTINE FILL_EL_LOSS_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build thermal drag force heating term for electrons
  SUBROUTINE FILL_EL_T_DRAG_TEMP(u_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::  u_e(NUM_X) !< Lagged electron velocity

    INTEGER :: K, P

    REAL(KIND=DEFAULT_REAL) :: VAL, &!< Current value of matrix element to be passed to global matrix
                               loc_dx(MIN_X:MAX_X)

    K = 1


    DO P = MIN_X,MAX_X

      loc_dx(P) = dxp(P) + dxm(P)

    END DO

    IF ((MIN_X .EQ. 1) .AND. (NO_FLOW_BOUNDARY_UP_SWITCH)) loc_dx(1) = dxp(1)
    IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) loc_dx = dxm(NUM_X)

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 1) THEN

          VAL = - EL_T_DRAG_TEMP_SIGN(K) * 2.00D00  * 0.71D00 * u_e(P) / 3.00D00

          LOCAL_M%VALUE(MARKER_EL_T_DRAG_TEMP(K)) = LOCAL_M%VALUE(MARKER_EL_T_DRAG_TEMP(K)) + VAL/loc_dx(P)

          K = K + 1

          VAL = - EL_T_DRAG_TEMP_SIGN(K) * 2.00D00  * 0.71D00 * u_e(P) / 3.00D00

          LOCAL_M%VALUE(MARKER_EL_T_DRAG_TEMP(K)) = LOCAL_M%VALUE(MARKER_EL_T_DRAG_TEMP(K)) + VAL/loc_dx(P)

          K = K + 1

        END IF

      END DO

  END SUBROUTINE FILL_EL_T_DRAG_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build e-i drag force term for electron temp eq
  SUBROUTINE FILL_EL_U_DRAG_TEMP(n_o,T_o,u_e,T_e,u_i)

    IMPLICIT NONE

    INTEGER :: K, P

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_e(NUM_X), &!< Lagged electron velocity
                                           T_e(NUM_X) , &!< Lagged electron temperature
                                           n_o(NUM_X), &!< Old electron density
                                           u_i(NUM_X), &!< Lagged ion velocity
                                           T_o(NUM_X) !< Old electron temperature

    REAL(KIND=DEFAULT_REAL) :: VAL!< Current value of matrix element to be passed to global matrix

    K = 1

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 1) THEN

          VAL =  DENSITY_0 * TIME_NORM/(TEMP_0_EV**(3.00D00/2.00D00)) * 1.48D-12&
               * LAMBDA_EI(n_o(P),T_o(P),Z_PROF(P)) *Z_PROF(P)/ION_Z* u_e(P)*(u_e(P)-u_i(P))/T_e(P)**(3.00D00/2.00D00)

          LOCAL_M%VALUE(MARKER_EL_U_DRAG_TEMP(K)) = LOCAL_M%VALUE(MARKER_EL_U_DRAG_TEMP(K)) + 4.00D00*VAL/3.00D00

          K = K + 1

        END IF

      END DO

  END SUBROUTINE FILL_EL_U_DRAG_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update (de)excitation drag matrix elements
  SUBROUTINE UPDATE_R_EX_DEEX(f_1)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::  f_1(NUM_V,NUM_X) !< Lagged f1 vector

    INTEGER :: I, P, nl, nu, ROW, COLUMN

    REAL(KIND=DEFAULT_REAL) :: VAL

    IF (.NOT. ALLOCATED(R_ex)) ALLOCATE(R_ex(NUM_X,NUM_NEUTRALS))
    IF (.NOT. ALLOCATED(R_deex)) ALLOCATE(R_deex(NUM_X,NUM_NEUTRALS))


    DO P = MIN_X, MAX_X

      DO nl = 1, NUM_NEUTRALS

        VAL = 0

        DO nu = nl + 1, NUM_NEUTRALS

          DO I = 1, INEL_SP(nl,nu)%N_NZ

            ROW = INEL_SP(nl,nu)%ROW(I)
            COLUMN = INEL_SP(nl,nu)%COL(I)

            VAL = VAL - COLL_EN_L * f_1(COLUMN,P) * V_GRID(ROW) ** 4 * V_GRID_WIDTH(ROW) * &
                  SIGMA_L(0, ROW , nl, nu)  * IN_DAT(nl,nu)%W(ROW,COLUMN)

            VAL = VAL + COLL_EN_L * f_1(COLUMN,P) * V_GRID(COLUMN)* V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW) * &
                  (SIGMA_L(0, COLUMN , nl, nu) - SIGMA_L(1,COLUMN, nl, nu)) * IN_DAT(nl,nu)%W(ROW,COLUMN)

          END DO

        END DO

        R_ex(P,nl) = 4.00D00 * PI * VAL / 3.00D00

      END DO

      DO nu = 1, NUM_NEUTRALS

        VAL = 0

        DO nl = 1, nu - 1

          DO I = 1, INEL_SP(nu,nl)%N_NZ

            ROW = INEL_SP(nu,nl)%ROW(I)
            COLUMN = INEL_SP(nu,nl)%COL(I)

            VAL = VAL - COLL_EN_L * f_1(COLUMN,P) * V_GRID(ROW) ** 4 * V_GRID_WIDTH(ROW) * &
                  SIGMA_L_DEEX(0, nu, nl, P, COLUMN)  * IN_DAT(nu,nl)%W(ROW,COLUMN)

            VAL = VAL + COLL_EN_L * f_1(COLUMN,P) * V_GRID(COLUMN)* V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW) * &
                  (SIGMA_L_DEEX(0, nu, nl, P, COLUMN ) - SIGMA_L_DEEX(1, nu, nl, P, COLUMN)) * IN_DAT(nu,nl)%W(ROW,COLUMN)

          END DO

        END DO

        R_deex(P,nu) = 4.00D00 * PI * VAL / 3.00D00

      END DO

    END DO


  END SUBROUTINE UPDATE_R_EX_DEEX
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update ionization drag matrix elements
  SUBROUTINE UPDATE_R_ION(f_1)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_1(NUM_V,NUM_X) !< Lagged f1 vector

    INTEGER :: I, P, nl, ROW, COLUMN

    REAL(KIND=DEFAULT_REAL) :: VAL

    IF (.NOT. ALLOCATED(R_ion)) ALLOCATE(R_ion(NUM_X,NUM_NEUTRALS))

    DO P = MIN_X, MAX_X

      DO nl = 1, NUM_NEUTRALS

        VAL = 0

          DO I = 1, INEL_SP(nl,0)%N_NZ

            ROW = INEL_SP(nl,0)%ROW(I)
            COLUMN = INEL_SP(nl,0)%ROW(I)

            VAL = VAL - COLL_EN_L * f_1(COLUMN,P) * V_GRID(ROW) ** 4 * V_GRID_WIDTH(ROW) * &
                  SIGMA_L(0, ROW , nl, 0)  * IN_DAT(nl,0)%W(ROW,COLUMN)

            VAL = VAL + COLL_EN_L * f_1(COLUMN,P) * V_GRID(COLUMN)* V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW) * &
                  (SIGMA_L(0, COLUMN , nl, 0) - SIGMA_L(1,COLUMN, nl, 0)) * IN_DAT(nl,0)%W(ROW,COLUMN)

          END DO

        R_ion(P,nl) = 4.00D00 * PI * VAL / 3.00D00

      END DO

    END DO


  END SUBROUTINE UPDATE_R_ION
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update recombination drag matrix elements
  SUBROUTINE UPDATE_R_RECOMB(f_1,n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_1(NUM_V,NUM_X), & !< Lagged f1 vector
                                           n_i(NUM_X) !< Lagged ion density

    INTEGER :: I, P, nu, ROW, COLUMN

    REAL(KIND=DEFAULT_REAL) :: VAL

    IF (.NOT. ALLOCATED(R_recomb)) ALLOCATE(R_recomb(NUM_X))

    DO P = MIN_X, MAX_X

      VAL = 0

      DO nu = 1, NUM_NEUTRALS

          DO I = 1, INEL_SP(0,nu)%N_NZ

            ROW = INEL_SP(0,nu)%ROW(I)
            COLUMN = INEL_SP(0,nu)%COL(I)

            VAL = VAL - DENSITY_0 * DE_BROGLIE_L3 * COLL_EN_L * f_1(COLUMN,P) * &
                  n_i(P) * V_GRID(ROW) ** 4 * V_GRID_WIDTH(ROW) * SIGMA_L_RECOMB(0,nu, P,ROW) * IN_DAT(0,nu)%E(I)              !Diagonal elements

            VAL = VAL + DENSITY_0 * DE_BROGLIE_L3 * COLL_EN_L * f_1(COLUMN,P) &
                   * n_i(P)  * V_GRID(COLUMN) * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW) * &
                  (SIGMA_L_RECOMB(0,nu, P,COLUMN) -  SIGMA_L_RECOMB(1,nu, P,COLUMN)) * IN_DAT(0,nu)%W(ROW,COLUMN)  !Absorbing elements

          END DO

      END DO

      R_recomb(P) = 4.00D00 * PI * VAL / 3.00D00

    END DO


  END SUBROUTINE UPDATE_R_RECOMB
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update elastic drag matrix elements
  SUBROUTINE UPDATE_R_EL(f_1)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_1(NUM_V,NUM_X) !< Lagged f1 vector

    INTEGER :: I, P, nu

    REAL(KIND=DEFAULT_REAL) :: VAL

    IF (.NOT. ALLOCATED(R_el)) ALLOCATE(R_el(NUM_X,NUM_NEUTRALS))

    DO P = MIN_X, MAX_X

     DO nu = 1, NUM_NEUTRALS

        VAL = 0

        DO I = 1, NUM_V

          VAL = VAL - COLL_EN_L * f_1(I,P) * V_GRID(I) ** 4 * V_GRID_WIDTH(I) * SIGMA_L_EL(1,I,nu)

        END DO

        R_el(P,nu) = 4.00D00 * PI * VAL / 3.00D00

      END DO

    END DO


  END SUBROUTINE UPDATE_R_EL
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build e-n drag force term for electrons
  SUBROUTINE FILL_EL_EN_DRAG(f,n_i,n_e)

    IMPLICIT NONE

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), &!< Lagged electron density
                                           n_i(NUM_X), &!< Lagged ion density
                                           f(DIM_F) !< Lagged unknown vector

    REAL(KIND=DEFAULT_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                               f_1(NUM_V,NUM_X), &
                               TOTAL_R(NUM_X,NUM_NEUTRALS)  !< Sum of matrix elements for included e-n collisions (excluding recombination)

    DO I = 1, NUM_X

        f_1(:,I) = f(X_POS(I) + NUM_V*(NUM_H-2) + 1:X_POS(I) + NUM_V*(NUM_H-1))

    END DO

    TOTAL_R = 0

    IF (COLL_EN_EL_L_SWITCH) THEN

      CALL UPDATE_R_EL(f_1)
      TOTAL_R = TOTAL_R + R_el

    END IF

    IF (COLL_EN_EX) THEN

      CALL UPDATE_R_EX_DEEX(f_1)
      TOTAL_R = TOTAL_R + R_ex
      TOTAL_R = TOTAL_R + R_deex

    END IF

    IF (COLL_EN_ION) THEN

      CALL UPDATE_R_ION(f_1)
      TOTAL_R = TOTAL_R + R_ion

    END IF

    IF (EL_EN_DRAG_N_SP%N_NZ .GT. 0) THEN

      DO I = 1, EL_EN_DRAG_N_SP%N_NZ

          VAL = EL_EN_DRAG_N_MULT(I) * TOTAL_R(EL_EN_DRAG_N_POS(I),EL_EN_DRAG_N(I)) &
                /n_e(EL_EN_DRAG_N_POS(I))

          LOCAL_M%VALUE(MARKER_EL_EN_DRAG(I)) = LOCAL_M%VALUE(MARKER_EL_EN_DRAG(I)) + VAL

      END DO

    END IF

    IF (COLL_RECOMB) THEN

      CALL UPDATE_R_RECOMB(f_1,n_i)

      DO I = 1, EL_EN_DRAG_RECOMB_SP%N_NZ

        VAL = EL_EN_DRAG_REC_MULT(I) * R_recomb(EL_EN_DRAG_REC_POS(I)) &
              /n_e(EL_EN_DRAG_REC_POS(I))

        LOCAL_M%VALUE(MARKER_EL_EN_DRAG_RECOMB(I)) = LOCAL_M%VALUE(MARKER_EL_EN_DRAG_RECOMB(I)) + VAL

      END DO

    END IF

  END SUBROUTINE FILL_EL_EN_DRAG
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build e-n drag temperature term for electrons
  SUBROUTINE FILL_EL_EN_DRAG_TEMP(n_e,u_e)

    IMPLICIT NONE

    INTEGER :: I, J, K

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_e(NUM_X), & !< Lagged electron velocity
                                           n_e(NUM_X)!< Lagged electron density
    REAL(KIND=DEFAULT_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                               TOTAL_R(NUM_X,NUM_NEUTRALS)  !< Sum of matrix elements for included e-n collisions (excluding recombination)


    TOTAL_R = 0

    IF (COLL_EN_EL_L_SWITCH) THEN

      TOTAL_R = TOTAL_R + R_el

    END IF

    IF (COLL_EN_EX) THEN

      TOTAL_R = TOTAL_R + R_ex
      TOTAL_R = TOTAL_R + R_deex

    END IF

    IF (COLL_EN_ION) THEN

      TOTAL_R = TOTAL_R + R_ion

    END IF

    K = 1

    IF (EL_EN_DRAG_N_SP%N_NZ .GT. 0) THEN

      DO I = MIN_X,MAX_X

        IF (MOD(I,2) .EQ. 1) THEN

          DO J = 1, NUM_NEUTRALS

            VAL = - u_e(I) * TOTAL_R(I,J)/n_e(I)

            LOCAL_M%VALUE(MARKER_EL_EN_DRAG_TEMP(K)) = LOCAL_M%VALUE(MARKER_EL_EN_DRAG_TEMP(K)) + 4.00D00*VAL/3.00D00

            K = K + 1

          END DO

        END IF

      END DO

    END IF

    IF (COLL_RECOMB) THEN

      DO I = MIN_X,MAX_X

        IF (MOD(I,2) .EQ. 1) THEN

          VAL = -u_e(I) * R_recomb(I)/n_e(I)

          LOCAL_M%VALUE(MARKER_EL_EN_DRAG_TEMP(K)) = LOCAL_M%VALUE(MARKER_EL_EN_DRAG_TEMP(K)) + 4.00D00*VAL/3.00D00

          K = K + 1

        END IF

      END DO

    END IF

  END SUBROUTINE FILL_EL_EN_DRAG_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron energy loss/gain term due to e-n inelastic collisions
  SUBROUTINE FILL_EL_EN_INCOLL_TEMP(n_e,T_e, n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: t_e(NUM_X), & !< Lagged electron temperature
                                           n_e(NUM_X), &!< Lagged electron density
                                           n_i(NUM_X) !< Lagged ion density

    INTEGER :: I, J, P

    REAL(KIND=HIGH_PREC_REAL) :: VAL, & !< Current value of matrix element to be passed to global matrix
                                 EN_RATE(MIN_X:MAX_X,NUM_NEUTRALS+1)

    EN_RATE = 0

    IF (COLL_EN_EX) THEN

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 1) THEN

          DO I = 1, NUM_NEUTRALS

            DO J = I + 1, NUM_NEUTRALS

              EN_RATE(P,I) = EN_RATE(P,I) - EX_RATES(I,J,P) * ION_POT_H * &
              (1.00D00 / REAL(I ** 2,KIND=HIGH_PREC_REAL) - 1.00D00 / REAL(J ** 2,KIND=HIGH_PREC_REAL))

            END DO

            DO J = 1, I - 1

              EN_RATE(P,I) = EN_RATE(P,I) - DEEX_RATES(I,J,P) * ION_POT_H * &
              (1.00D00 / REAL(I ** 2,KIND=HIGH_PREC_REAL) - 1.00D00 / REAL(J ** 2,KIND=HIGH_PREC_REAL))

            END DO

          END DO

        END IF

      END DO

    END IF

    IF (COLL_EN_ION) THEN

      DO P = MIN_X,MAX_X

        IF (MOD(P,2) .EQ. 1) THEN

          DO I = 1, NUM_NEUTRALS

            EN_RATE(P,I) = EN_RATE(P,I) - ION_RATES(I,P) * &
            ION_POT_H / REAL(I ** 2,KIND=HIGH_PREC_REAL)

          END DO

        END IF

      END DO

    END IF

    IF (COLL_RECOMB) THEN

      DO P = MIN_X,MAX_X

        DO I = 1, NUM_NEUTRALS

          EN_RATE(P,NUM_NEUTRALS+1) = EN_RATE(P,NUM_NEUTRALS+1) + &
          (RAD_REC_A_0*RECOMB_RATE(I, T_e(P)) + TB_REC_A_0 *TB_RECOMB_RATES(I,P))*n_i(P) *&
          ION_POT_H / REAL(I ** 2,KIND=HIGH_PREC_REAL)

        END DO

      END DO

    END IF

    J = 1

    DO P = MIN_X, MAX_X

      IF (MOD(P,2) .EQ. 1) THEN

        DO I = 1, NUM_NEUTRALS + 1

          VAL = 2.00D00*EN_RATE(P,I)/(3.00D00*n_e(P))

          LOCAL_M%VALUE(MARKER_EL_EN_INCOLL_TEMP(J)) = LOCAL_M%VALUE(MARKER_EL_EN_INCOLL_TEMP(J)) + VAL

          J = J + 1

        END DO

      END IF

    END DO

  END SUBROUTINE FILL_EL_EN_INCOLL_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Build electron energy loss/gain term due to e-n elastic collisions (assuming constant cross-section)
  SUBROUTINE FILL_EL_EN_ELCOLL_TEMP(T_e)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T_e(NUM_X) !< Lagged electron temperature

    INTEGER :: I, J, P

    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    J = 1

    DO P = MIN_X, MAX_X

      IF (MOD(P,2) .EQ. 1) THEN

        DO I = 1, NUM_NEUTRALS

          VAL = - 2.00D00**4/(3.00D00*SQRT(PI))* COLL_EN_0 * SQRT(T_e(P))*(T_e(P)-NEUTRAL_TEMP) * SIGMA_L_EL(1,1,I)

          LOCAL_M%VALUE(MARKER_EL_EN_ELCOLL_TEMP(J)) = LOCAL_M%VALUE(MARKER_EL_EN_ELCOLL_TEMP(J)) + VAL

          J = J + 1

        END DO

      END IF

    END DO

  END SUBROUTINE FILL_EL_EN_ELCOLL_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_FLUID_EL
!-------------------------------------------------------------------------------------------------------------------------------------
