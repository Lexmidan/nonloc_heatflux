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
!> Contains all matrix builders and related routines for Coulomb collision integrals
MODULE BUILD_COULOMB_COLL

  USE GRID
  USE NORMALIZATION
  USE SWITCHES
  USE MATRIX_DATA
  USE MPI
  USE VAR_KIND_DATA

  IMPLICIT NONE

!> Contains row contracted drag terms for ion flow (and e-e collision momentum differncing error)
  TYPE DRAG_TERMS

    REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: VAL(:) !< Value

  END TYPE DRAG_TERMS

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE, DIMENSION(:,:) :: OLD_CC_WEIGHTS, & !< Old Chang-Cooper weights
                                         NEW_CC_WEIGHTS !< New Chang-Cooper weights

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE, DIMENSION(:,:,:) :: I_L_G, & !< ROSENBLUTH_I(L)
                                           I_L2_G, & !< ROSENBLUTH_I(L+2)
                                           J_L1_G, & !< ROSENBLUTH_J(1-L)
                                           J_L_MINUS_G !< ROSENBLUTH_(L)

  PRIVATE :: OLD_CC_WEIGHTS, NEW_CC_WEIGHTS, UPDATE_WEIGHTS, C00_DRAG, C00_DIFF, &
             ROSENBLUTH_I, ROSENBLUTH_J, D2F, R_INT_I, R_INT_J

  TYPE (DRAG_TERMS), ALLOCATABLE, DIMENSION(:) :: EI_DRAG_TERMS, EE_FD_ERR

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Updates Chang-Cooper weights
  SUBROUTINE UPDATE_WEIGHTS(N_NONLIN, f_lagged,P)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N_NONLIN, & !< Current nonlinear iteration
                           P !< Current spatial cell

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(NUM_V) !< Lagged f_0^0

    REAL(KIND=DEFAULT_REAL) :: W

    INTEGER :: I

    IF (N_NONLIN .EQ. 0) THEN                                                   !Initialize weights if first nonlinear iteration

      OLD_CC_WEIGHTS(:,P) = 0.5D00
      OLD_CC_WEIGHTS(NUM_V,P) = 0.0D00

    ELSE

      OLD_CC_WEIGHTS(:,P) = NEW_CC_WEIGHTS(:,P)

    END IF

!Calculate weights
    DO I = 1, NUM_V - 1

      W = dvp(I) * C00_DRAG(I, f_lagged) / C00_DIFF(I, f_lagged,P)

      NEW_CC_WEIGHTS(I,P) = 1.00D00/W - 1.00D00/(EXP(W) - 1.00D00)

    END DO

    NEW_CC_WEIGHTS(0,P) = 0.50D00
    NEW_CC_WEIGHTS(NUM_V,P) = 0.0D00

  END SUBROUTINE UPDATE_WEIGHTS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!>Allocate CC weights
  SUBROUTINE ALLOCATE_WEIGHTS

    IMPLICIT NONE

    IF (.NOT. ALLOCATED(OLD_CC_WEIGHTS)) THEN

      ALLOCATE(OLD_CC_WEIGHTS(0:NUM_V,MIN_X:MAX_X))
      OLD_CC_WEIGHTS = 0.0D00

    END IF

    IF (.NOT. ALLOCATED(NEW_CC_WEIGHTS)) THEN

      ALLOCATE(NEW_CC_WEIGHTS(0:NUM_V,MIN_X:MAX_X))
      NEW_CC_WEIGHTS = 0.0D00

    END IF

  END SUBROUTINE ALLOCATE_WEIGHTS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculate Chang-Cooper drag coefficient
  REAL(KIND=DEFAULT_REAL) FUNCTION C00_DRAG(I, f_lagged)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: I !< Velocity cell index at which coefficient is evaluated
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(NUM_V) !< Lagged f_0^0

    INTEGER :: J
    REAL(KIND=DEFAULT_REAL) :: C

    C = 0

    DO J = 1, I

      C = C + 4.00D00 * PI * V_GRID(J) ** 2 * f_lagged(J) * V_GRID_WIDTH(J)

    END DO

    C00_DRAG = C

  END FUNCTION C00_DRAG
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculate Chang-Cooper diffusion coefficient
  REAL(KIND=DEFAULT_REAL) FUNCTION C00_DIFF(I, f_lagged,P)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: I,& !< Velocity cell index at which coefficient is evaluated
                           P   !< Current spatial cell
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(NUM_V) !< Lagged f_0^0

    INTEGER :: J, K

    REAL(KIND=DEFAULT_REAL) :: D, V_INV, V1

    D = 0

    IF ((I .GT. 0) .AND. (I .LT. NUM_V)) THEN

      V_INV = 1.00D00/ (V_GRID(I + 1) + V_GRID(I))

      DO J = 1, I

        V1 = 4.00D00 * PI * V_GRID(J) ** 2 * V_GRID_WIDTH(J) * V_INV

        DO K = J, NUM_V - 1

          D = D +  ((1.00D00 - OLD_CC_WEIGHTS(K,P)) * f_lagged(K + 1) + &
                    OLD_CC_WEIGHTS(K,P) * f_lagged(K)) * (V_GRID(K + 1) + V_GRID(K)) * V1 * dvp(K)

        END DO

      END DO

    END IF

    C00_DIFF = D

  END FUNCTION C00_DIFF
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills e-e collision matrix for f_0^0
  SUBROUTINE FILL_EE_00(f, n_o, T_o, N_NONLIN,P)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(NUM_V), & !< Lagged f_0^0
                          n_o, & !< Old density
                          T_o !< Old temperature
    INTEGER, INTENT(IN) :: N_NONLIN, &!< Nonlinear iteration number
                           P!< Current spatial cell

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: LAMBDA, & !< Coulomb logarithm
              VAL, & !< Current value of matrix element to be passed to global matrix
              C(0:NUM_V), & !< CC drag coefficients
              D(0:NUM_V) !< CC diffusion coefficients

    LAMBDA = LAMBDA_EE(n_o, T_o)                                                !Initialiaze Coulomb logarithm
    CALL UPDATE_WEIGHTS(N_NONLIN, f,P)                                          !Update CC weights

!Initialize CC drag and diffusion coefficients
    DO I = 0, NUM_V

      C(I) = C00_DRAG(I,f)
      D(I) = C00_DIFF(I,f,P)

    END DO

!Fill matrix according to tri-diagonal sparsity pattern and adequate markers
    DO I = 1, TRI_DIAG_SP%N_NZ

      VAL = COLL_EE_0 * LAMBDA * &
            ((TRI_DIAG_SP%COL(I) / TRI_DIAG_SP%ROW(I)) * (TRI_DIAG_SP%ROW(I) / TRI_DIAG_SP%COL(I)) * &            !Diagonal elements
            (C(TRI_DIAG_SP%ROW(I)) * NEW_CC_WEIGHTS(TRI_DIAG_SP%ROW(I),P) &
            - D(TRI_DIAG_SP%ROW(I))/dvp(TRI_DIAG_SP%ROW(I)) &
            - C(TRI_DIAG_SP%ROW(I) - 1) * (1.00D00 - NEW_CC_WEIGHTS(TRI_DIAG_SP%ROW(I) - 1,P)) &
            - D(TRI_DIAG_SP%ROW(I) - 1)/dvm(TRI_DIAG_SP%ROW(I))) &
            + (TRI_DIAG_SP%COL(I) / (TRI_DIAG_SP%ROW(I) + 1)) * &                                                 !Upper off-diagonal elements
            (C(TRI_DIAG_SP%ROW(I)) * (1.00D00 - NEW_CC_WEIGHTS(TRI_DIAG_SP%ROW(I),P)) &
            + D(TRI_DIAG_SP%ROW(I))/dvp(TRI_DIAG_SP%ROW(I))) &
            + (TRI_DIAG_SP%ROW(I) / (TRI_DIAG_SP%COL(I) +1)) * &                                                  !Lower off-diagonal elements
            (- C(TRI_DIAG_SP%ROW(I) - 1) * NEW_CC_WEIGHTS(TRI_DIAG_SP%ROW(I) - 1,P) &
            + D(TRI_DIAG_SP%ROW(I) - 1)/dvm(TRI_DIAG_SP%ROW(I)))) &
            / (V_GRID(TRI_DIAG_SP%ROW(I)) ** 2 * V_GRID_WIDTH(TRI_DIAG_SP%ROW(I)))

      LOCAL_M%VALUE(MARKER_EE_0(P) + I) = &                                      !Send value to global matrix
                          LOCAL_M%VALUE(MARKER_EE_0(P) + I) + VAL

    END DO

  END SUBROUTINE FILL_EE_00
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills e-i collision integral submatrix for given L-number
  SUBROUTINE FILL_EI_L(n_i, n_o, T_o, L,P,H)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i, & !< Lagged ion density
                          n_o, & !< Old density
                          T_o !< Old temperature

    INTEGER, INTENT(IN) :: L, & !< Harmonic L-number
                           P, & !< Current position (cell boundary)
                           H !< Harmonic index of current harmonic
    INTEGER :: I, LL, L1

    REAL(KIND=DEFAULT_REAL) :: LAMBDA, & !< Coulomb logarithm
              VAL !< Current value of matrix element to be passed to global matrix

    LAMBDA = LAMBDA_EI(n_o, T_o,Z_PROF(P))                                                !Initialiaze Coulomb logarithm

    LL = - (L * (L + 1))/ 2

    IF (.NOT. ALLOCATED(EI_DRAG_TERMS)) THEN

      ALLOCATE(EI_DRAG_TERMS(MIN_X:MAX_X))

        DO I = MIN_X, MAX_X

          ALLOCATE(EI_DRAG_TERMS(I)%VAL(NUM_V + 2))

        END DO

    END IF

    IF (L .EQ. 1) THEN

      EI_DRAG_TERMS(P)%VAL = 0

    END IF

    L1 = ((L+1)/2)*(2/(L+1))

!Fill matrix according to simple diagonal pattern and adequate markers
    DO I = 1, NUM_V

      VAL = COLL_EI_0 * LAMBDA * LL * n_i * (Z_PROF(P) ** 2/ION_Z**2) / (V_GRID(I) ** 3)

      LOCAL_M%VALUE(MARKER_EI_L(P,H) + I) = LOCAL_M%VALUE(MARKER_EI_L(P,H) + I) + VAL !Send value to global matrix

      EI_DRAG_TERMS(P)%VAL(I) =  EI_DRAG_TERMS(P)%VAL(I) + L1 * VAL * V_GRID(I) ** 3 * V_GRID_WIDTH(I)

    END DO

  END SUBROUTINE FILL_EI_L
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills e-i collision integral submatrix for given L-number - cold ion case
  SUBROUTINE FILL_CI_EI_L(n_i, n_o, T_o, L,P,H,u_i,f_0)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i, & !< Lagged ion density
                          n_o, & !< Old density
                          T_o, & !< Old temperature
                          u_i, & !< Lagged ion velocity
                          f_0(NUM_V)

    INTEGER, INTENT(IN) :: L, & !< Harmonic L-number
                           P, & !< Current position (cell boundary)
                           H !< Harmonic index of current harmonic

    INTEGER :: I, LL, K, ROW, COL, L1

    REAL(KIND=DEFAULT_REAL) :: LAMBDA, & !< Coulomb logarithm
              VAL, & !< Current value of matrix element to be passed to global matrix
              C(6), &
              a_plus, &
              a_minus

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_V) :: I_0, &
                                      I_2, &
                                      J_m, &
                                      I_2L, &
                                      J_ml, &
                                      J_1ml,&
                                      Il, &
                                      A, &
                                      B1, &
                                      B2, &
                                      A_l, &
                                      B_l, &
                                      C_l, &
                                      C_0, &
                                      df0, &
                                      ddf0, &
                                      C_n, &
                                      B

    LAMBDA = LAMBDA_EI(n_o, T_o,Z_PROF(P)) * Z_PROF(P)**2/ION_Z**2                                          !Initialiaze Coulomb logarithm

    LL = - (L * (L + 1))/ 2

    DO I = 1, NUM_V

      IF (V_GRID(I) .GT. u_i) THEN

        K = I

        EXIT

      END IF

    END DO

!Initialize interpolation constants
    IF ((.NOT. PERIODIC_BOUNDARY_SWITCH) .AND. (MOD(L,2) .EQ. 1)) THEN

      a_minus = (X_GRID(P + 1) - X_GRID(P)) / dxc(P)
      a_plus = (X_GRID(P) - X_GRID(P - 1)) / dxc(P)

    ELSE

      a_minus = 0.50D00
      a_plus = 0.50D00

    END IF

!Initialize drag terms for ions
    IF (.NOT. ALLOCATED(EI_DRAG_TERMS)) THEN

      ALLOCATE(EI_DRAG_TERMS(MIN_X:MAX_X))

        DO I = MIN_X, MAX_X

          ALLOCATE(EI_DRAG_TERMS(I)%VAL(NUM_V + 2))

        END DO

    END IF

    IF (L .EQ. 1) THEN

      EI_DRAG_TERMS(P)%VAL = 0

    END IF

    L1 = ((L+1)/2)*(2/(L+1))

!Initialize f_0 derivatives
    DO I = 1, NUM_V

      df0(I) = D1F(I,f_0)
      ddf0(I) = D2F(I,f_0)

    END DO

!Initialize constants
    C(1) = (L + 1.0D00) * (L + 2.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(2) = - (L - 1.0D00) * L / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(3) = - (L * (L + 1.0D00) / 2.0D00 + L + 1.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(4) = (- L * (L + 1.0D00) / 2.0D00 + L + 2.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(5) = (L * (L + 1.0D00) / 2.0D00 + L - 1.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L - 1.0D00))
    C(6) = (L * (L + 1.0D00) / 2.0D00 - L) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L - 1.0D00))

!Calculate needed Rosenbluth coeffs

    I_0 = COLD_ION_I(u_i,n_i,0,0,K)
    I_2 = COLD_ION_I(u_i,n_i,0,2,K)
    Il = COLD_ION_I(u_i,n_i,L,L,K)
    I_2L = COLD_ION_I(u_i,n_i,L,L+ 2,K)
    J_m =  COLD_ION_J(u_i,n_i,0,-1,K)
    J_ml = COLD_ION_J(u_i,n_i,L, -L - 1,K)
    J_1ml = COLD_ION_J(u_i,n_i,L, 1 - L,K)

!Introduce shorthands
    DO I = 1, NUM_V

      A(I) = (I_2(I) + J_m(I))/(3 * V_GRID(I)*V_GRID_WIDTH(I))

      B1(I) = (-I_2(I) + 2* J_m(I) + 3*I_0(I))/(3*V_GRID(I) ** 2)

      B2(I) = B1(I) * LL / V_GRID(I)

      A_l(I) = (C(1)*I_2L(I) + C(1) * J_ml(I) + C(2) * Il(I) + C(2) * J_1ml(I)) / (2 * V_GRID(I))

      B_l(I) = (C(3)*I_2L(I) + C(4) * J_ml(I) + C(5) * Il(I) + C(6) * J_1ml(I)) / (V_GRID(I) ** 2)

      C_l(I) = -  ((L+1.00D00)/(2*L+1.00D00) * Il(I) - L * J_ml(I)/(2*L+1.00D00))/(V_GRID(I) ** 2)

      C_0(I) = - 1/(V_GRID(I) ** 2) * I_0(I)

      C_n(I) = A_l(I) * ddf0(I) + (C_l(I) + B_l(I)) * df0(I)

      C_n(I) = C_n(I)/n_i

    END DO

    DO I = 1, NUM_V

      B(I) = (C_0(I) + B1(I))/(dvp(I)+dvm(I))

    END DO

!Handle v=0 exception
    B(1) = B(1) * (dvp(1)+dvm(1)) / V_GRID(2)

!Fill in main tri-diagonal submatrix - f_l
    DO I = 1, TRI_DIAG_SP%N_NZ

      ROW = TRI_DIAG_SP%ROW(I)
      COL = TRI_DIAG_SP%COL(I)

      VAL = COLL_EI_0 * LAMBDA * (&
            ((ROW/COL) * (COL/ROW))*(- A(ROW) * (1.00D00/dvp(ROW)+1.00D00/dvm(ROW) &
            - ((ROW+1)/2)*(2/(ROW+1))*(1.00D00/dvm(1) - L1/V_GRID(1)))  + B2(ROW)) &              !Diagonal elements
            + (COL / (ROW + 1)) * &                                                 !Upper off-diagonal elements
            ( A(ROW)/dvp(ROW) + B(ROW)) &
            + (ROW / (COL + 1)) * &                                                  !Lower off-diagonal elements
            (A(ROW)/dvm(ROW) - B(ROW)))

      LOCAL_M%VALUE(MARKER_EI_L(P,H) + I) = &                                      !Send value to global matrix
                          LOCAL_M%VALUE(MARKER_EI_L(P,H) + I) + VAL

        EI_DRAG_TERMS(P)%VAL(COL) =  EI_DRAG_TERMS(P)%VAL(COL) + L1 * VAL * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW)

    END DO

!Fill off-diagonal f_0 terms for even L
    IF (MOD(L,2) .EQ. 0) THEN

      DO I = 1,NUM_V

        VAL = COLL_EI_0 * LAMBDA * C_n(I)

        LOCAL_M%VALUE(MARKER_CI_EI_R(P,H,I)) = &                                      !Send value to global matrix
                                 LOCAL_M%VALUE(MARKER_CI_EI_R(P,H,I)) + VAL

      END DO

!Fill last cell f_0 terms for odd L
    ELSE IF (P  .EQ. NUM_X) THEN

      DO I = 1,NUM_V

        ROW = I
        COL = NUM_V

        VAL = 0.50D00 * COLL_EI_0 * LAMBDA * C_n(I)

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

            LOCAL_M%VALUE(MARKER_CI_EI_R(P,H,I)) = &         !Send value to global matrix
                                 LOCAL_M%VALUE(MARKER_CI_EI_R(P,H,I)) + VAL


            EI_DRAG_TERMS(P)%VAL(COL + 1) =  EI_DRAG_TERMS(P)%VAL(COL + 1) + L1*VAL * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW)


        END IF

             LOCAL_M%VALUE(MARKER_CI_EI_L(P,H,I)) = &                               !Send value to global matrix
                                 LOCAL_M%VALUE(MARKER_CI_EI_L(P,H,I)) + VAL

             EI_DRAG_TERMS(P)%VAL(COL + 2) = EI_DRAG_TERMS(P)%VAL(COL + 2) + L1 *  VAL * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW)

      END DO

!Fill f_0 terms for odd L in bulk
    ELSE

      DO I = 1,NUM_V

        ROW = I
        COL = NUM_V

        VAL = COLL_EI_0 * LAMBDA * C_n(I)


        LOCAL_M%VALUE(MARKER_CI_EI_R(P,H,I)) = &         !Send value to global matrix
                                 LOCAL_M%VALUE(MARKER_CI_EI_R(P,H,I)) + a_plus * VAL

        EI_DRAG_TERMS(P)%VAL(COL + 1) =  EI_DRAG_TERMS(P)%VAL(COL + 1) + a_plus * L1 * VAL * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW)

        LOCAL_M%VALUE(MARKER_CI_EI_L(P,H,I)) = &                               !Send value to global matrix
                                 LOCAL_M%VALUE(MARKER_CI_EI_L(P,H,I)) + a_minus * VAL

        EI_DRAG_TERMS(P)%VAL(COL + 2) = EI_DRAG_TERMS(P)%VAL(COL + 2) + a_minus * L1 * VAL * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW)

      END DO

    END IF

  END SUBROUTINE FILL_CI_EI_L
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates matrix elements arising from Rosenbluth integral I (cf. Shkarofsky)
  FUNCTION ROSENBLUTH_I(J)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_V, NUM_V) :: ROSENBLUTH_I, I_J
    INTEGER, INTENT(IN) :: J
    INTEGER :: K, L

    I_J = 0

    DO K = 1, NUM_V

      DO L = 1, K

        IF (L .EQ. K) THEN

          I_J(K,L) = 2.00D00 * PI * V_GRID(L) ** 2 * V_GRID_WIDTH(L)

        ELSE

          I_J(K,L) = 4.00D00 * PI * (V_GRID(L)/V_GRID(K)) ** J * V_GRID(L) ** 2 * V_GRID_WIDTH(L)

        END IF

      END DO

    END DO

    ROSENBLUTH_I = I_J

  END FUNCTION ROSENBLUTH_I
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates matrix elements arising from Rose integral J (cf. Shkarofsky)
  FUNCTION ROSENBLUTH_J(J)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_V, NUM_V) :: ROSENBLUTH_J, J_J
    INTEGER, INTENT(IN) :: J
    INTEGER :: K, L

    J_J = 0

    DO K = 1, NUM_V

      DO L = K, NUM_V

        IF (L .EQ. K) THEN

          J_J(K,L) = 2.00D00 * PI * V_GRID(K) ** 2 * V_GRID_WIDTH(L)

        ELSE

          J_J(K,L) = 4.00D00 * PI * (V_GRID(L)/V_GRID(K)) ** J * V_GRID(L) ** 2 * V_GRID_WIDTH(L)

        END IF

      END DO

    END DO

    ROSENBLUTH_J = J_J

  END FUNCTION ROSENBLUTH_J
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates Rosenbluth integral I of lagged f at velocity cell N
  REAL(KIND=DEFAULT_REAL) FUNCTION R_INT_I(N, J, f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(NUM_V)

    INTEGER, INTENT(IN) :: N, &
                           J

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: INT, &
              RI(NUM_V, NUM_V)

    INT = 0

    IF (J .EQ. 2) THEN

      RI = I_L2_G(0,:,:)

    ELSE IF (J .EQ. 0) THEN

      RI = I_L_G(0,:,:)

    END IF

    DO I = 1, NUM_V

      INT = INT + RI(N,I) * f(I)

    END DO

    R_INT_I = INT

  END FUNCTION R_INT_I
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates Rosenbluth integral J of lagged f at velocity cell N
  REAL(KIND=DEFAULT_REAL) FUNCTION R_INT_J(N, J, f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(NUM_V)
    INTEGER, INTENT(IN) :: N, &
                         J

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: INT, &
              RJ(NUM_V, NUM_V)

    INT = 0

    IF (J .EQ. -1) THEN

      RJ = J_L_MINUS_G(0,:,:)

    END IF

    DO I = 1, NUM_V

      INT = INT + RJ(N,I) * f(I)

    END DO

    R_INT_J = INT

  END FUNCTION R_INT_J
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates first derivative of f_0^0 at cell N
  REAL(KIND=DEFAULT_REAL) FUNCTION D1F(N,f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(NUM_V)
    INTEGER, INTENT(IN) :: N

    REAL(KIND=DEFAULT_REAL) :: d

    IF (N .EQ. 1) THEN                                                          !Extrapolate f_0^0 at v = 0

      d = (f(2) - (f(1) - f(2) * V_GRID(1) ** 2 / V_GRID(2) ** 2) / &
          (1.00D00 - V_GRID(1) ** 2 / V_GRID(2) ** 2)) / V_GRID(2)

    ELSE IF (N .EQ. NUM_V) THEN

      d = - f(NUM_V - 1) / (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1))

    ELSE

      d = (f(N + 1) - f(N - 1)) / (V_GRID(N + 1) - V_GRID(N - 1))

    END IF

    D1F = d

  END FUNCTION D1F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates second derivative of f_0^0 at cell N
  REAL(KIND=DEFAULT_REAL) FUNCTION D2F(N,f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(NUM_V)
    INTEGER, INTENT(IN) :: N

    REAL(KIND=DEFAULT_REAL) :: d

    IF (N .EQ. 1) THEN

      d = (f(2) - f(1)) / (dv * dvp(1))

    ELSE IF (N .EQ. NUM_V) THEN

      d = - (f(NUM_V) - f(NUM_V - 1)) / (V_GRID_WIDTH(NUM_V)*dvm(NUM_V))

    ELSE

      d = ((f(N + 1) - f(N))/dvp(N) + (f(N - 1)-f(N))/dvm(N)) / V_GRID_WIDTH(N)

    END IF

    D2F = d

  END FUNCTION D2F
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills e-e collision integral for L>0 contribution to matrix (cf. Tzoufras 2011)
  SUBROUTINE FILL_EE_L(n_o, T_o, f, L,P,H)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_o, & !< Old density
                          T_o, & !< Old temperature
                          f(NUM_V) !< Lagged f_0^0
    INTEGER, INTENT(IN) :: L, & !< Harmonic L-number
                           P, & !< Current position (cell boundary)
                           H !< Harmonic index of current harmonic

    REAL(KIND=DEFAULT_REAL) :: LAMBDA, & !< Coulomb logarith
              VAL, & !< Current value of matrix element to be passed to global matrix
              !Various calculation constants and shorthands
              A(NUM_V), &
              B(NUM_V), &
              I_L(NUM_V,NUM_V), &
              I_L2(NUM_V,NUM_V), &
              J_L1(NUM_V, NUM_V), &
              J_L_MINUS(NUM_V,NUM_V), &
              I2, &
              J_MINUS, &
              I0, &
              C(6), &
              D1(NUM_V), &
              D2(NUM_V), &
              AL(NUM_V), &
              A_OD(NUM_V), &
              BL_M(NUM_V), &
              BL_P(NUM_V)

    INTEGER :: LL, I, ROW, COLUMN, L1

    LL = - (L * (L + 1))/2

    LAMBDA = LAMBDA_EE(n_o,T_o)                                                 !Initialize Coulomb logarithm

    L1 = ((L+1)/2)*(2/(L+1))

!Allocate differencing error terms
    IF (.NOT. ALLOCATED(EE_FD_ERR)) THEN

      ALLOCATE(EE_FD_ERR(MIN_X:MAX_X))

      DO I = MIN_X, MAX_X

        ALLOCATE(EE_FD_ERR(I)%VAL(NUM_V + 2))

      END DO

    END IF

    IF (L .EQ. 1) THEN

      EE_FD_ERR(P)%VAL = 0

    END IF

!Initialize constants
    C(1) = (L + 1.0D00) * (L + 2.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(2) = - (L - 1.0D00) * L / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(3) = - (L * (L + 1.0D00) / 2.0D00 + L + 1.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(4) = (- L * (L + 1.0D00) / 2.0D00 + L + 2.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L + 3.0D00))
    C(5) = (L * (L + 1.0D00) / 2.0D00 + L - 1.0D00) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L - 1.0D00))
    C(6) = (L * (L + 1.0D00) / 2.0D00 - L) / ((2.0D00 * L + 1.0D00) * (2.0D00 * L - 1.0D00))

!Calculate Rosenbluth integral terms
    IF (DIAG_EE_L_SWITCH) THEN

      I_L = 0
      I_L2 = 0
      J_L1 = 0
      J_L_MINUS = 0
      D1 = 0
      D2 = 0

    ELSE

      I_L = I_L_G(L,:,:)
      I_L2 = I_L2_G(L,:,:)
      J_L1 = J_L1_G(L,:,:)
      J_L_MINUS = J_L_MINUS_G(L,:,:)

      DO I = 1, NUM_V

        D1(I) = D1F(I, f)
        D2(I) = D2F(I, f)

      END DO

    END IF

!Calculate shorthands
    DO I = 1, NUM_V

      I2 = R_INT_I(I,2,f)
      J_MINUS = R_INT_J(I,-1,f)
      I0 = R_INT_I(I,0,f)

      A(I) = (- I2 + 2.0D00 * J_MINUS + 3.0D00 * I0) / 3.0D00
      B(I) = (I2 + J_MINUS) / 3.0D00

      AL(I) = LL * A(I) / V_GRID(I) ** 3

      BL_M(I) = - B(I) / (V_GRID(I) * V_GRID_WIDTH(I))
      BL_P(I) = B(I) / (V_GRID(I) * V_GRID_WIDTH(I))

    END DO

    DO I = 2, NUM_V

      A_OD(I) = A(I) / (V_GRID(I) ** 2 * (dvp(I)+dvm(I)))

    END DO

!Hande v=0 exception
    A_OD(1) = A(1) / (V_GRID(1) ** 2 * V_GRID(2))

!Fill tridiagonal bit of submatrix according to adequate sparsity pattern and markers
    DO I = 1, TRI_DIAG_SP%N_NZ

      ROW = TRI_DIAG_SP%ROW(I)
      COLUMN = TRI_DIAG_SP%COL(I)

!Tridiagonal component
      VAL = COLL_EE_0 * LAMBDA * &
            ((COLUMN/ ROW) * (ROW / COLUMN) * &                                 !Diagonal elements
            (AL(ROW) + 8.00D00 * PI * f(ROW) + BL_M(ROW)*(1.00D00/dvp(ROW) + 1.0D00/dvm(ROW))) &
            + (COLUMN / (ROW + 1)) * &                                          !Upper off-diagonal elements
            (A_OD(ROW) + BL_P(ROW)/dvp(ROW)) &
            + (ROW / (COLUMN +1)) * &                                           !Lower off-diagonal elements
            ( - A_OD(ROW) + BL_P(ROW)/dvm(ROW)))

!Triangular matrix contributions to tridiagonal components
      VAL = VAL + COLL_EE_0 * LAMBDA * &
            (1.00D00/(2.0D00 * V_GRID(ROW)) * D2(ROW) * &
            (C(1) * I_L2(ROW, COLUMN) &
            + C(2) * I_L(ROW, COLUMN)) &
            + 1.00D00/(V_GRID(ROW) ** 2) * D1(ROW) * &
            (C(3) * I_L2(ROW, COLUMN) &
            + C(5) * I_L(ROW, COLUMN)))

      VAL = VAL + COLL_EE_0 * LAMBDA * &
            (1.00D00/(2.0D00 * V_GRID(ROW)) * D2(ROW) * &
            (C(1) * J_L_MINUS(ROW, COLUMN) &
            + C(2) * J_L1(ROW, COLUMN)) &
            + 1.00D00/(V_GRID(ROW) ** 2) * D1(ROW) * &
            (C(4) * J_L_MINUS(ROW, COLUMN) &
            + C(6) * J_L1(ROW, COLUMN)))

      LOCAL_M%VALUE(MARKER_EE_L(P,H) + I) = LOCAL_M%VALUE(MARKER_EE_L(P,H) + I) + VAL !Send value to global matrix


      EE_FD_ERR(P)%VAL(COLUMN) =  EE_FD_ERR(P)%VAL(COLUMN) + L1*VAL * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW)

    END DO

!Add triangular submatrices according to dense sparsity pattern
    IF (.NOT. DIAG_EE_L_SWITCH) THEN

      DO I = TRI_DIAG_SP%N_NZ + 1, DENSE_SP%N_NZ

        ROW = DENSE_SP%ROW(I)
        COLUMN = DENSE_SP%COL(I)

        VAL = COLL_EE_0 * LAMBDA * &
              (1.00D00/(2.0D00 * V_GRID(ROW)) * D2(ROW) * &
              (C(1) * I_L2(ROW, COLUMN) &
              + C(2) * I_L(ROW, COLUMN)) &
              + 1.00D00/(V_GRID(ROW) ** 2) * D1(ROW) * &
              (C(3) * I_L2(ROW, COLUMN) &
              + C(5) * I_L(ROW, COLUMN)))

        VAL = VAL + COLL_EE_0 * LAMBDA * &
              (1.00D00/(2.0D00 * V_GRID(ROW)) * D2(ROW) * &
              (C(1) * J_L_MINUS(ROW, COLUMN) &
              + C(2) * J_L1(ROW, COLUMN)) &
              + 1.00D00/(V_GRID(ROW) ** 2) * D1(ROW) * &
              (C(4) * J_L_MINUS(ROW, COLUMN) &
              + C(6) * J_L1(ROW, COLUMN)))

        LOCAL_M%VALUE(MARKER_EE_L(P,H) + I) = LOCAL_M%VALUE(MARKER_EE_L(P,H) + I) + VAL !Send value to global matrix

        EE_FD_ERR(P)%VAL(COLUMN) =  EE_FD_ERR(P)%VAL(COLUMN) + L1*VAL * V_GRID(ROW) ** 3 * V_GRID_WIDTH(ROW)

      END DO

    END IF

    IF (COLD_ION_FLUID_SWITCH .AND. COLL_EI_L_SWITCH) EI_DRAG_TERMS(P)%VAL = EI_DRAG_TERMS(P)%VAL + EE_FD_ERR(P)%VAL   !Add differencing error from e-e collisions as ion momentum

  END SUBROUTINE FILL_EE_L
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculate cold ion I_j(F_L)
  FUNCTION COLD_ION_I(u_i,n_i, L,J,K)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_i, n_i
    INTEGER, INTENT(IN) :: L,K,J

    REAL(KIND=DEFAULT_REAL):: COLD_ION_I(NUM_V)

    REAL(KIND=DEFAULT_REAL) :: CI_I(NUM_V), LL

    INTEGER :: I

    LL =  (2 * L + 1) * n_i

    CI_I = 0

    DO I = K, NUM_V

      CI_I(I) = LL * u_i ** J/ V_GRID(I) ** J

    END DO

    COLD_ION_I = CI_I

  END FUNCTION COLD_ION_I
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculate cold ion J_j(F_L)
  FUNCTION COLD_ION_J(u_i,n_i, L,J,K)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: u_i, n_i
    INTEGER, INTENT(IN) :: L,K,J

    REAL(KIND=DEFAULT_REAL) :: COLD_ION_J(NUM_V)

    REAL(KIND=DEFAULT_REAL) :: CI_J(NUM_V), LL

    INTEGER :: I

    LL =  (2 * L + 1) * n_i

    CI_J = 0

    DO I = 1, K - 1

      CI_J(I) = LL * u_i ** J/ V_GRID(I) ** J

    END DO

    COLD_ION_J = CI_J

  END FUNCTION COLD_ION_J
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes precalculated I and J integrals
  SUBROUTINE INIT_ROSENBLUTH

    IMPLICIT NONE

    INTEGER :: I

    ALLOCATE(I_L_G(0:L_MAX,NUM_V,NUM_V))
    ALLOCATE(I_L2_G(0:L_MAX,NUM_V,NUM_V))
    ALLOCATE(J_L1_G(0:L_MAX,NUM_V,NUM_V))
    ALLOCATE(J_L_MINUS_G(0:L_MAX,NUM_V,NUM_V))

    DO I = 0, L_MAX

      I_L_G(I,:,:) = ROSENBLUTH_I(I)
      I_L2_G(I,:,:) = ROSENBLUTH_I(I + 2)
      J_L1_G(I,:,:) = ROSENBLUTH_J(1 - I)
      J_L_MINUS_G(I,:,:) =  ROSENBLUTH_J(-1 - I)

    END DO

  END SUBROUTINE INIT_ROSENBLUTH
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_COULOMB_COLL
!-------------------------------------------------------------------------------------------------------------------------------------
