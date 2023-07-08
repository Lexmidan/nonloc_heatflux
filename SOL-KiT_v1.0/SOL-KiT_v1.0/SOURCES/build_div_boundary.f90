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
!> Contains plasma sink and recycling submatrix building routines
MODULE BUILD_DIV_BOUNDARY

  USE GRID
  USE NORMALIZATION
  USE NEUT_AND_HEAT
  USE SOLVER_PARAMS
  USE SWITCHES
  USE MATRIX_DATA
  USE MPI
  USE COLL_CS_RATES
  USE VAR_KIND_DATA

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: P_FIXED(:,:) !< Fixed part of Legendre tensor

  REAL(KIND=DEFAULT_REAL) :: n_e_boundary, & !< Cut-off distribution density
                  T_e_boundary, & !< Cut-off distribution temperature
                  ion_flux, & !< Ion flux
                  gamma_e, & !< Electron sheath heat transmission coefficient
                  pot_drop !< Sheath potential drop normalized to temperature


  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: f_boundary(:,:), & !< Boundary distribution function harmonics to be used in interpolation
                               shift_mult(:) !< Boundary distribution harmonic rescaling factors

  LOGICAL :: BOHM_VALUE

  REAL(KIND=DEFAULT_REAL) :: FLOW_LIMITER, & !< Factor limiting ion flow when extrapolating
                             cut_off_v
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills divertor boundary submatrix, including sink and recycling terms
  SUBROUTINE FILL_DIV_BOUNDARY(n_e,T_e,f,u_i,n_i)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e(NUM_X), & !< Lagged electron density
                                T_e(NUM_X), & !< Lagged electron temperature
                                f(DIM_F) !< lagged f

    REAL(KIND=DEFAULT_REAL), INTENT(IN), OPTIONAL :: u_i(NUM_X), & !< Lagged ion velocity
                                          n_i(NUM_X)    !< Lagged ion density

    REAL(KIND=DEFAULT_REAL) :: VAL, &
              v_f, & !< Cut-off velocity
              a1, & !< Right interpolation constant
              a2, & !< Left interpolation constant
              el_flux, & !< Calculated electron flux
              bis_err, & !< False position error
              v_p1, &  !< Right interpolated velocity cell centre
              v_p2, & !< Left interpolated velocity cell centre
              v_f1, v_f2, y1, y2, & ! False position method variables
              f_ext(NUM_V,NUM_H), & !< Extrapolated electron distribution function on last cell boundary
              P_LL(NUM_H,NUM_H,NUM_V), & !< Legendre decomposition tensor for boundary cut-off
              dv1, & !< Right interpolated cell width
              dv2, & !< Left interpolated cell width
              ion_flux1, & ! False position method variables
              ion_flux2, &  ! False position method variables
              T_i, & ! Ion temperature
              ion_flux_bohm, &
              shift_mult_even, &
              shift_mult_odd


    INTEGER :: I, K, N, J, L, P

    LOGICAL :: POSITIVE_FLUX

    BOHM_VALUE = .FALSE.

    cut_off_v = 0.00D00

    ion_flux = 0.0D00

    FLOW_LIMITER = 1.00D00

    T_i = 0.00D00

    n_e_boundary = n_e(NUM_X)*n_e(NUM_X)/n_e(NUM_X-1)

    T_e_boundary = T_e(NUM_X)

    IF (ION_EL_TEMP_SWITCH) T_i = T_e_boundary

    ion_flux_bohm = n_e_boundary * MACH_N_DIV*SQRT((T_e_boundary + T_i) * EL_MASS / ( 2.00D00 * ION_MASS)) !Ion Bohm flux
    !Perform flux initialization
    IF (PRESENT(u_i)) THEN

      ion_flux = u_i(NUM_X - 1)*n_i(NUM_X - 1) + &
                (u_i(NUM_X-1)*n_i(NUM_X -1) - u_i(NUM_X-3)*n_i(NUM_X - 3))/&
                (dxm(NUM_X-1)) * dxc(NUM_X)

      POSITIVE_FLUX = .TRUE.

      IF (SONIC_OUTFLOW_DIV_SWITCH) THEN

        IF ((ion_flux .LT. ion_flux_bohm) .OR. NO_EXTRAPOLATION_SWITCH) THEN

          ion_flux = ion_flux_bohm

          BOHM_VALUE = .TRUE.

        END IF

      ELSE

        IF (ion_flux .LE. 0) POSITIVE_FLUX = .FALSE.

      END IF


    ELSE

      BOHM_VALUE = .TRUE.
      ion_flux = ion_flux_bohm

      POSITIVE_FLUX = .TRUE.

    END IF

!Add sink
    IF (PLASMA_SINK_ON .AND. (.NOT. FULL_FLUID_MODE)) THEN

      IF (.NOT. ALLOCATED(P_FIXED)) CALL P_FIXED_INIT

      IF (.NOT. ALLOCATED(shift_mult)) ALLOCATE(shift_mult(SINK_SP%N_NZ))

      !Exponential density extrapolation multipliers
      shift_mult_even = n_e(NUM_X)/n_e(NUM_X-1)
      shift_mult_odd = (n_e(NUM_X)/n_e(NUM_X-1))**2
      !Set extrapolation multipliers (to avoid f_0 going negative)
      shift_mult = shift_mult_even + (shift_mult_odd &
      - shift_mult_even)*REAL(SINK_SHIFT,KIND=DEFAULT_REAL)

      IF (POSITIVE_FLUX) THEN

        !Extrapolate to boundary

        f_ext = 0

        DO J = 1, NUM_H

          IF (MOD(J,2) .EQ. 1) THEN

            DO I = 1, NUM_V

              f_ext(I,J) = shift_mult_even*f(X_POS(NUM_X)+NUM_V*(NUM_H-J)+I)

            END DO

          ELSE

            DO I = 1, NUM_V

              f_ext(I,J) = shift_mult_odd*f(X_POS(NUM_X-1)+NUM_V*(NUM_H-J)+I)

            END DO

          END IF

        END DO

        el_flux = 0

!Determine closest spatial cell to cut-off
        DO P = NUM_V - 1, 2, -1

          v_f = V_CELL_BOUNDARY(P-1)

          el_flux = FLUX(f_ext, v_f, P)

          K = P

          IF (el_flux .GT. ion_flux) THEN

            EXIT

          END IF

        END DO

        IF (K .EQ. 2) THEN

          PRINT*, 'WARNING: Logical boundary condition cut-off close to 0, ion and electron fluxes potentially unequal'
          PRINT*, 'Adjusting ion flow'

          FLOW_LIMITER = el_flux/ion_flux
          ion_flux = el_flux

        END IF

!Initialize false position method variables
        bis_err = 1.00D00

        v_f1 = V_CELL_BOUNDARY(K-1)

        v_f2 = V_CELL_BOUNDARY(K)

        ion_flux1 = ion_flux
        ion_flux2 = ion_flux

        y1 = FLUX(f_ext,v_f1,K) - ion_flux1
        y2 = FLUX(f_ext,v_f2,K) - ion_flux2

!Calculate precise cut-off
        DO WHILE (bis_err .GT. BIS_TOL)

           v_f = (v_f1 * y2 - v_f2 * y1)/ (y2 - y1)

           el_flux = FLUX(f_ext, v_f, K)

          IF (el_flux .GT. ion_flux) THEN

             v_f1 = v_f

             y1 = el_flux - ion_flux

          ELSE

             v_f2 = v_f

             y2 = el_flux - ion_flux

          END IF

          bis_err = ABS(el_flux - ion_flux) / ion_flux

        END DO

        cut_off_v = v_f

!Calculate interpolated cell data
        CALL INTERP_CELLS(v_f,K,v_p1,v_p2,dv1,dv2,a1,a2)

!Calculate Legendre tensor components
        DO I = 1, NUM_H

          P_LL(I,:,:) = P_ROW(I-1,v_f,K)

        END DO

        N = 1

!Fill matrix
!Lower grid cells
        DO I = 1, K-2

          DO J = 1, NUM_H

            VAL = (-1.00D00/3.00D00 * V_GRID(I) /dxc(NUM_X)) *  P_LL(2,J,I)

            LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

            N = N + 1

          END DO

          DO J = 2,NUM_H

            IF (MOD(J,2) .EQ. 1) THEN

              DO L = 1, NUM_H

                VAL = (-REAL(J-1,KIND=DEFAULT_REAL)/REAL(2*J-3,KIND=DEFAULT_REAL) &
                                      * V_GRID(I)/dxc(NUM_X) * P_LL(J-1,L,I)&
                                      -REAL(J,KIND=DEFAULT_REAL)/REAL(2*J+1,KIND=DEFAULT_REAL) &
                                      * V_GRID(I)/dxc(NUM_X)* P_LL(J+1,L,I))

                LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

                N = N + 1

              END DO

            END IF

          END DO

        END DO

! K-1
        DO J = 1, NUM_H

          VAL = (-1.00D00/3.00D00 * v_p2**3 /dxc(NUM_X)) * P_LL(2,J,K-1) * a2 * dv2 /(V_GRID_WIDTH(K-1)*V_GRID(K-1)**2)

          LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

          N = N + 1

        END DO

        DO J = 2,NUM_H

          IF (MOD(J,2) .EQ. 1) THEN

            DO L = 1, NUM_H

              VAL = (-REAL(J-1,KIND=DEFAULT_REAL)/REAL(2*J-3,KIND=DEFAULT_REAL) * v_p2**3/dxc(NUM_X)* P_LL(J-1,L,K-1)&
                      -REAL(J,KIND=DEFAULT_REAL)/REAL(2*J+1,KIND=DEFAULT_REAL) * v_p2**3/dxc(NUM_X)* P_LL(J+1,L,K-1)) &
                      *a2 * dv2 /(V_GRID_WIDTH(K-1)*V_GRID(K-1)**2)

              LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

              N = N + 1

            END DO

          END IF

        END DO

! K
        DO J = 1, NUM_H

          VAL = (-1.00D00/3.00D00 * v_p2**3 /dxc(NUM_X)) *  P_LL(2,J,K-1) * (1.00D00-a2) * dv2 /(V_GRID_WIDTH(K)*V_GRID(K)**2)

          VAL = VAL + (-1.00D00/3.00D00 * v_p1**3 /dxc(NUM_X)) *  P_LL(2,J,K) * a1 * dv1 /(V_GRID_WIDTH(K)*V_GRID(K)**2)


          LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

          N = N + 1

        END DO

        DO J = 2,NUM_H

          IF (MOD(J,2) .EQ. 1) THEN

            DO L = 1, NUM_H

              VAL = (-REAL(J-1,KIND=DEFAULT_REAL)/REAL(2*J-3,KIND=DEFAULT_REAL) &
                                    * v_p2**3/dxc(NUM_X)* P_LL(J-1,L,K-1)&
                                    -REAL(J,KIND=DEFAULT_REAL)/REAL(2*J+1,KIND=DEFAULT_REAL) &
                                    * v_p2**3/dxc(NUM_X)* P_LL(J+1,L,K-1)) &
                                    * (1.00D00-a2) * dv2 /(V_GRID_WIDTH(K)*V_GRID(K)**2)

              VAL = VAL + (-REAL(J-1,KIND=DEFAULT_REAL)/REAL(2*J-3,KIND=DEFAULT_REAL) &
                                          * v_p1**3/dxc(NUM_X) * P_LL(J-1,L,K)&
                                          -REAL(J,KIND=DEFAULT_REAL)/REAL(2*J+1,KIND=DEFAULT_REAL) &
                                          * v_p1**3/dxc(NUM_X)* P_LL(J+1,L,K)) &
                                          * a1 * dv1 /(V_GRID_WIDTH(K)*V_GRID(K)**2)

              LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

              N = N + 1

            END DO

          END IF

        END DO

! K + 1
        DO J = 1, NUM_H

          VAL =   (-1.00D00/3.00D00 * v_p1**3 /dxc(NUM_X)) *  P_LL(2,J,K+1) * (1.00D00-a1) * dv1 /(V_GRID_WIDTH(K+1)*V_GRID(K+1)**2)

          VAL = VAL + (-1.00D00/3.00D00 * V_GRID(K+1) /dxc(NUM_X)) *  P_LL(2,J,K+1)


          LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

          N = N + 1

        END DO

        DO J = 2,NUM_H

          IF (MOD(J,2) .EQ. 1) THEN

            DO L = 1, NUM_H


              VAL = (-REAL(J-1,KIND=DEFAULT_REAL)/REAL(2*J-3,KIND=DEFAULT_REAL) &
                                    * v_p1**3/dxc(NUM_X) * P_LL(J-1,L,K+1)&
                                    -REAL(J,KIND=DEFAULT_REAL)/REAL(2*J+1,KIND=DEFAULT_REAL) &
                                    * v_p1**3/dxc(NUM_X)* P_LL(J+1,L,K+1)) &
                                    * (1.00D00-a1) * dv1 /(V_GRID_WIDTH(K+1)*V_GRID(K+1)**2)

              VAL = VAL + (-REAL(J-1,KIND=DEFAULT_REAL)/REAL(2*J-3,KIND=DEFAULT_REAL) &
                                    * V_GRID(K+1)/dxc(NUM_X) * P_LL(J-1,L,K+1)&
                                    -REAL(J,KIND=DEFAULT_REAL)/REAL(2*J+1,KIND=DEFAULT_REAL) &
                                    * V_GRID(K+1)/dxc(NUM_X)* P_LL(J+1,L,K+1))

              LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

              N = N + 1

            END DO

          END IF

        END DO

!Upper grid cells
        DO I = K+2,NUM_V

          DO J = 1, NUM_H

            VAL = (-1.00D00/3.00D00 * V_GRID(I) /dxc(NUM_X)) *  P_LL(2,J,I)


            LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

            N = N + 1

          END DO

          DO J = 2,NUM_H

            IF (MOD(J,2) .EQ. 1) THEN

              DO L = 1, NUM_H

                VAL = (-REAL(J-1,KIND=DEFAULT_REAL)/REAL(2*J-3,KIND=DEFAULT_REAL) &
                                      * V_GRID(I)/dxc(NUM_X) * P_LL(J-1,L,I)&
                                      -REAL(J,KIND=DEFAULT_REAL)/REAL(2*J+1,KIND=DEFAULT_REAL) &
                                      * V_GRID(I)/dxc(NUM_X)* P_LL(J+1,L,I))

                LOCAL_M%VALUE(MARKER_SINK(N)) = LOCAL_M%VALUE(MARKER_SINK(N)) + shift_mult(N)*VAL

                N = N + 1

              END DO

            END IF

          END DO

        END DO

      END IF

!Calculate sheath heat transmission coefficient if Bohm switch on

      IF (SONIC_OUTFLOW_DIV_SWITCH) THEN

        VAL = INTERP_MOM(1,3,v_f,f_ext,K)/6.00D00 !Heat flux
        gamma_e = 2*VAL/(ion_flux*T_e_boundary)
        pot_drop = v_f**2/T_e_boundary

      END IF

      DO I = 1, NUM_H

        DO J = 1, NUM_V

          f_boundary(I,J) = 0

          DO K = 1, NUM_H

            f_boundary(I,J) = f_boundary(I,J) + P_LL(I,K,J)*f_ext(J,K)

          END DO

        END DO

      END DO

    END IF



!Add recycling
    IF (REC_ON .AND. (POSITIVE_FLUX)) THEN

      IF (FULL_FLUID_MODE) THEN

        VAL = REC_R * (ion_flux / n_e(NUM_X)) / &
              (2.00D00 * (X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

        LOCAL_M%VALUE(MARKER_REC(1)) = LOCAL_M%VALUE(MARKER_REC(1)) + VAL

      ELSE

        IF (PRESENT(u_i) .AND. (.NOT. BOHM_VALUE)) THEN


          VAL =    n_i(NUM_X - 1) + &
                      n_i(NUM_X -1) /&
                      (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)) * (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

          LOCAL_M%VALUE(MARKER_CI_REC(1)) = LOCAL_M%VALUE(MARKER_CI_REC(1)) &
                                + FLOW_LIMITER * REC_R * VAL / (2.00D00 * (X_GRID(NUM_X) - X_GRID(NUM_X - 1)))


          VAL =  - n_i(NUM_X - 3)/&
          (X_GRID(NUM_X - 1) - X_GRID(NUM_X - 3)) * (2.00D00*(X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

          LOCAL_M%VALUE(MARKER_CI_REC(2)) = LOCAL_M%VALUE(MARKER_CI_REC(2)) &
                                + FLOW_LIMITER * REC_R * VAL / (2.00D00 * (X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

        ELSE

          IF ((SONIC_OUTFLOW_DIV_SWITCH .AND. COLD_ION_FLUID_SWITCH) .AND. (.NOT. ION_CONT_OFF_SWITCH)) THEN

              VAL = REC_R * (ion_flux / n_i(NUM_X)) / &
                    (2.00D00 * (X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

              LOCAL_M%VALUE(MARKER_REC(1)) = LOCAL_M%VALUE(MARKER_REC(1)) + VAL

          ELSE

            DO I = 1, NUM_V

              VAL = REC_R * (n_e_boundary / n_e(NUM_X)) * MACH_N_DIV*SQRT((T_e_boundary+T_i)* EL_MASS / ( 2.00D00 * ION_MASS))&
               * 4.00D00 * PI * V_GRID(I) ** 2 * V_GRID_WIDTH(I) / &
                    (2.00D00 * (X_GRID(NUM_X) - X_GRID(NUM_X - 1)))

              LOCAL_M%VALUE(MARKER_REC(I)) = LOCAL_M%VALUE(MARKER_REC(I)) + VAL

            END DO

          END IF

        END IF

      END IF

    END IF

  END SUBROUTINE FILL_DIV_BOUNDARY
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates electron flux including cut-off
  REAL(KIND=DEFAULT_REAL) FUNCTION FLUX(f_ext, v_f, K)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_ext(NUM_V,NUM_H), & !< extrapolated harmonics on boundary
                          v_f !< Cut-off velocity

    INTEGER, INTENT(IN) :: K !< Velocity cell closest to cutoff

    REAL(KIND=DEFAULT_REAL) :: fl

    fl = INTERP_MOM(1,1,v_f,f_ext,K)

    FLUX = fl / 3.00D00

  END FUNCTION FLUX
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates Lth row of Legendre decomposition tensor
   FUNCTION P_ROW(L, v_f,K_c)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) ::   v_f !< Cut-off velocity

    INTEGER, INTENT(IN) :: L, &!< Row L number
                           K_c

    REAL(KIND=DEFAULT_REAL) :: P_R(NUM_H,NUM_V), P_ROW(NUM_H,NUM_V)


    INTEGER :: I,K,P

    REAL(KIND=DEFAULT_REAL) :: Pl(NUM_V,0:NUM_H), &
                    z_min(NUM_V), &
                    v_p1

    P_R = 0.0D00
    v_p1 = (V_CELL_BOUNDARY(K_c) + v_f) / 2.00D00            !Interpolated right velocity cell in integral
    z_min = - v_f/V_GRID
    z_min(K_c) = - v_f/v_p1
    CALL p_polynomial_value(NUM_V,NUM_H,z_min,Pl)

    IF (L .EQ. 0) THEN

      DO I = 1, NUM_H

        IF (I - 1 .EQ. 0) THEN

          DO K = K_c, NUM_V

            P_R(I,K) = 0.50D00 - 0.50D00 * z_min(K)

            P_R(I,K) = (-1.00D00)**(I-1) * P_R(I,K) + P_FIXED(L+1,I)

          END DO

        ELSE

          IF (I .NE. 1) THEN

            DO K = K_c, NUM_V

              P_R(I,K) = - z_min(K) * Pl(K,L)*Pl(K,I-1)/REAL(I+L,KIND=DEFAULT_REAL)
              P_R(I,K) = P_R(I,K) - REAL(I-1,KIND=DEFAULT_REAL)*Pl(K,L)*Pl(K,I-2)/REAL((L-I+1)*(L+I),KIND=DEFAULT_REAL)
              P_R(I,K) = 0.50D00*REAL(2*L + 1,KIND=DEFAULT_REAL)*P_R(I,K)

              P_R(I,K) = (-1.00D00)**(I-1) * P_R(I,K) + P_FIXED(L+1,I)

            END DO

          ELSE

            DO K = K_c, NUM_V

              P_R(I,K) = - z_min(K) * Pl(K,L)*Pl(K,I-1)/REAL(I+L,KIND=DEFAULT_REAL)
              P_R(I,K) = 0.50D00*REAL(2*L + 1,KIND=DEFAULT_REAL)*P_R(I,K)

              P_R(I,K) = (-1.00D00)**(I-1) * P_R(I,K) + P_FIXED(L+1,I)

            END DO

        END IF

      END IF

      END DO

    ELSE

      DO I = 1, NUM_H

        IF (I - 1 .EQ. L) THEN

          DO K = K_c, NUM_V

            P_R(I,K) = 0.50D00 - 0.50D00 * z_min(K) * Pl(K,I-1) ** 2

            DO P = 1, L - 1

              P_R(I,K) = P_R(I,K) - Pl(K,P)*(z_min(K)*Pl(K,P) - Pl(K,P+1))

            END DO

            P_R(I,K) = (-1)**(I-1) * P_R(I,K) + P_FIXED(L+1,I)

          END DO

        ELSE

          IF (I .NE. 1) THEN

            DO K = K_c, NUM_V

              P_R(I,K) = - z_min(K) * Pl(K,L)*Pl(K,I-1)/REAL(I+L,KIND=DEFAULT_REAL)
              P_R(I,K) = P_R(I,K) + (REAL(L,KIND=DEFAULT_REAL)*Pl(K,L-1)*Pl(K,I-1) -&
               REAL(I-1,KIND=DEFAULT_REAL)*Pl(K,L)*Pl(K,I-2))/REAL((L-I+1)*(L+I),KIND=DEFAULT_REAL)
              P_R(I,K) = 0.50D00*REAL(2*L + 1,KIND=DEFAULT_REAL)*P_R(I,K)

              P_R(I,K) = (-1.00D00)**(I-1) * P_R(I,K) + P_FIXED(L+1,I)

            END DO

          ELSE

            DO K = K_c, NUM_V

              P_R(I,K) = - z_min(K) * Pl(K,L)*Pl(K,I-1)/REAL(I+L,KIND=DEFAULT_REAL)
              P_R(I,K) = P_R(I,K) + (REAL(L,KIND=DEFAULT_REAL)*Pl(K,L-1)*Pl(K,I-1))/REAL((L-I+1)*(L+I),KIND=DEFAULT_REAL)
              P_R(I,K) = 0.50D00*REAL(2*L + 1,KIND=DEFAULT_REAL)*P_R(I,K)

              P_R(I,K) = (-1.00D00)**(I-1) * P_R(I,K) + P_FIXED(L+1,I)

            END DO

          END IF

        END IF

      END DO

    END IF

    DO I = 1, NUM_H

      IF (I-1 .EQ. L) THEN

        DO K = 1, K_c - 1

          P_R(I,K) = (-1.00D00) ** (I-1) + P_FIXED(L+1,I)

        END DO

      ELSE

        DO K = 1, K_c - 1

          P_R(I,K) = P_FIXED(L+1,I)

        END DO

      END IF

    END DO

    P_ROW = P_R

  END FUNCTION P_ROW
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fixed part of Legendre tensor
   SUBROUTINE P_FIXED_INIT

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: P_F(0:NUM_H,0:NUM_H) !< Placeholder

    INTEGER :: I,J

    REAL(KIND=DEFAULT_REAL) :: Pl0(1,0:NUM_H)

    REAL(KIND=DEFAULT_REAL) :: VALS(1)

    P_F = 0.0D00
    VALS = 0.00D00
    CALL p_polynomial_value(1,NUM_H,VALS,Pl0)

    ALLOCATE(P_FIXED(NUM_H,NUM_H))

    DO I = 1,NUM_H

      DO J = 1, NUM_H

        IF (MOD(J-1,2) .EQ. 1) THEN

          IF (I .EQ. J) THEN

            P_F(I,J) = 1.00D00

          ELSE IF (MOD(I-1,2) .EQ. 0) THEN

            IF (I .EQ. 1) THEN

              P_F(I,J) = - REAL(J-1,KIND=DEFAULT_REAL)*Pl0(1,J-2)

            ELSE IF (J .EQ. 1) THEN

              P_F(I,J) = REAL(I-1,KIND=DEFAULT_REAL)*Pl0(1,I-2)

            ELSE

              P_F(I,J) = REAL(I-1,KIND=DEFAULT_REAL)*Pl0(1,I-2)*Pl0(1,J-1) - REAL(J-1,KIND=DEFAULT_REAL)*Pl0(1,J-2)*Pl0(1,I-1)

            END IF

            P_F(I,J) = P_F(I,J) * REAL(2*(I-1) + 1,KIND=DEFAULT_REAL)/REAL((I - J) * (I + J -1),KIND=DEFAULT_REAL)

          END IF

        END IF

      END DO

    END DO

    P_FIXED = P_F

  END SUBROUTINE P_FIXED_INIT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates electron/ion density at cut-off
  REAL(KIND=DEFAULT_REAL) FUNCTION DIV_DENS(f_ext, v_f, K)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_ext(NUM_V,NUM_H), & !< extrapolated harmonics on boundary
                          v_f !< Cut-off velocity

    INTEGER, INTENT(IN) :: K !< Velocity cell closest to cutoff

    REAL(KIND=DEFAULT_REAL) :: fl

    fl = INTERP_MOM(0,0,v_f,f_ext,K)

    DIV_DENS = fl

  END FUNCTION DIV_DENS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates electron temperature at cut-off
  REAL(KIND=DEFAULT_REAL) FUNCTION DIV_TEMP(f_ext, v_f, K)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_ext(NUM_V,NUM_H), & !< extrapolated harmonics on boundary
                          v_f !< Cut-off velocity

    INTEGER, INTENT(IN) :: K !< Velocity cell closest to cutoff

    INTEGER :: T_i_factor

    REAL(KIND=DEFAULT_REAL) :: fl

    fl = INTERP_MOM(0,2,v_f,f_ext,K)

    T_i_factor = 0
    IF (ION_EL_TEMP_SWITCH) T_i_factor = 1
    DIV_TEMP = 2.00D00*fl/(3.00D00*n_e_boundary*(1.00D00+(1.00D00+T_i_factor)*EL_MASS/(3.00D00*ION_MASS)))

  END FUNCTION DIV_TEMP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates interpolation data around cut-off
  SUBROUTINE INTERP_CELLS(v_f,K,v1,v2,dv1,dv2,a1,a2)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: v_f !< Cut-off velocity
    INTEGER, INTENT(IN) :: K !< Cut-off cell

    REAL(KIND=DEFAULT_REAL), INTENT(INOUT) :: v1, & !< Right interpolated cell centre
                                   v2, & !< Left interpolated cell centre
                                   dv1, & !< Right interpolated cell width
                                   dv2, & !< Left interpolated cell width
                                   a1, & !< Right cell interpolation constant
                                   a2 !< Left cell interpolation constant


    v1 = (V_CELL_BOUNDARY(K) + v_f) / 2.00D00
    v2 = (V_CELL_BOUNDARY(K-2) + v_f) / 2.00D00

    a1 = 1.00D00 - (v1-V_GRID(K))/V_GRID_WIDTH(K)
    a2 = 1.00D00 - (v2-V_GRID(K-1))/V_GRID_WIDTH(K-1)            !This will be nonnegative iff V_GRID_MULT < 1.5

    dv1 = V_CELL_BOUNDARY(K) - v_f
    dv2 = v_f - V_CELL_BOUNDARY(K-2)

  END SUBROUTINE INTERP_CELLS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates interpolated M-th moment of L-th harmonic
  REAL(KIND=DEFAULT_REAL) FUNCTION INTERP_MOM(L,M,v_f,f_ext,K)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_ext(NUM_V,NUM_H), & !< extrapolated harmonics on boundary
                          v_f !< Cut-off velocity

    INTEGER, INTENT(IN) :: L,&!< Harmonic to use
                           M, & !< Moment number
                           K !< Velocity cell closest to cutoff

    INTEGER :: I, J

    REAL(KIND=DEFAULT_REAL) :: a1, &
              v_p1, &
              a2, &
              v_p2, &
              fl, &
              P_1L(NUM_H,NUM_V), &
              loc_f1, &
              dv1, dv2

    !Calculate interpolation data
    CALL INTERP_CELLS(v_f,K,v_p1,v_p2,dv1,dv2,a1,a2)

    P_1L = P_ROW(L,v_f,K)                            !Grab expansion coefficients
    fl = 0
!Add all lower contributions
    DO I = 1, K - 2

      loc_f1 = 0

      DO J = 1, NUM_H

        loc_f1 = loc_f1 + P_1L(J,I) * f_ext(I,J)

      END DO

      fl = fl + 4.00D00 * PI * loc_f1 * V_GRID(I) ** (2 + M) * V_GRID_WIDTH(I)

    END DO

!Add left interpolated cell

    loc_f1 = 0

    DO J = 1, NUM_H

      loc_f1 = loc_f1 + P_1L(J,K-1) * f_ext(K-1,J)

    END DO

    fl = fl + 4.00D00 * PI * v_p2 ** (2 + M) * dv2 * a2*loc_f1
    loc_f1 = 0

    DO J = 1, NUM_H

      loc_f1 = loc_f1  + P_1L(J,K-1) * f_ext(K,J)

    END DO
    fl = fl + 4.00D00 * PI * v_p2 ** (2 + M) * dv2 * (1.00D00-a2) * loc_f1

!Add right interpolated cell
    loc_f1 = 0

    DO J = 1, NUM_H

      loc_f1 = loc_f1  + P_1L(J,K) * f_ext(K,J)

    END DO

    fl = fl + 4.00D00 * PI * v_p1**(2 + M) * dv1 * a1 * loc_f1
    loc_f1 = 0

    DO J = 1, NUM_H

      loc_f1 = loc_f1 + P_1L(J,K+1) * f_ext(K+1,J)

    END DO

    fl = fl + 4.00D00 * PI * v_p1**(2 + M) * dv1 * (1.00D00-a1) * loc_f1

!Add all upper contributions

    DO I = K + 1, NUM_V

      loc_f1 = 0

      DO J = 1, NUM_H

        loc_f1 = loc_f1 + P_1L(J,I) * f_ext(I,J)

      END DO

      fl = fl + 4.00D00 * PI * loc_f1 * V_GRID(I) ** (2 + M) * V_GRID_WIDTH(I)

    END DO

    INTERP_MOM = fl

  END FUNCTION INTERP_MOM
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Allocates f_boundary
SUBROUTINE ALLOCATE_F_BOUNDARY

  IMPLICIT NONE

  IF (.NOT. ALLOCATED(f_boundary)) ALLOCATE(f_boundary(NUM_H,NUM_V))

END SUBROUTINE ALLOCATE_F_BOUNDARY
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_DIV_BOUNDARY
!-------------------------------------------------------------------------------------------------------------------------------------
