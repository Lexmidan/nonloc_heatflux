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
!> Contains all e-n collision integral matrix build routines
MODULE BUILD_NEUTRAL_COLL

  USE GRID
  USE NORMALIZATION
  USE NEUT_AND_HEAT
  USE COLL_CS_RATES
  USE INEL_GRID
  USE MATRIX_DATA
  USE MPI
  USE VAR_KIND_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Builds elastic e-n collision integral submatrix for f_0^0
  SUBROUTINE FILL_EN_EL_00(n_n,P,N)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n !< Neutral state density

    INTEGER, INTENT(IN) :: P, & !< Current position
                           N !< Atomic state number
    INTEGER :: I

    REAL(KIND=HIGH_PREC_REAL) :: SIGMA(NUM_V + 1), &
              C(0:NUM_V), &
              D(0:NUM_V), &
              WEIGHT(0:NUM_V), &
              VAL !< Current value of matrix element to be passed to global matrix

    SIGMA = SIGMA_L_EL(1,:,N)

!Define Chang-Cooper-like drag, diffusion, and weights
    C = 0
    D = 0
    WEIGHT = 0

    DO I = 1, NUM_V

      WEIGHT(I) = v_interp(I)

      C(I) = (V_CELL_BOUNDARY(I)) ** 4 * ((1.00D00 - v_interp(I))*SIGMA(I + 1) + v_interp(I)*SIGMA(I))

      D(I) = (V_CELL_BOUNDARY(I)) ** 3 * ((1.00D00 - v_interp(I))*SIGMA(I + 1) + v_interp(I)*SIGMA(I))  * NEUTRAL_TEMP / 2.00D00

    END DO

    WEIGHT(NUM_V) = 0

!Fill matrix according to tri-diagonal sparsity pattern and adequate markers
    DO I = 1, TRI_DIAG_SP%N_NZ

      VAL = COLL_EN_0 * n_n * &
            ((TRI_DIAG_SP%COL(I) / TRI_DIAG_SP%ROW(I)) * (TRI_DIAG_SP%ROW(I) / TRI_DIAG_SP%COL(I)) * &            !Diagonal elements
            (C(TRI_DIAG_SP%ROW(I)) * WEIGHT(TRI_DIAG_SP%ROW(I)) &
            - D(TRI_DIAG_SP%ROW(I))/dvp(TRI_DIAG_SP%ROW(I)) &
            - C(TRI_DIAG_SP%ROW(I) - 1) * (1.00D00 - WEIGHT(TRI_DIAG_SP%ROW(I) - 1)) &
            - D(TRI_DIAG_SP%ROW(I) - 1)/dvm(TRI_DIAG_SP%ROW(I))) &
            + (TRI_DIAG_SP%COL(I) / (TRI_DIAG_SP%ROW(I) + 1)) * &                                                 !Upper off-diagonal elements
            (C(TRI_DIAG_SP%ROW(I)) * (1.00D00 - WEIGHT(TRI_DIAG_SP%ROW(I))) &
            + D(TRI_DIAG_SP%ROW(I))/dvp(TRI_DIAG_SP%ROW(I))) &
            + (TRI_DIAG_SP%ROW(I) / (TRI_DIAG_SP%COL(I) +1)) * &                                                  !Lower off-diagonal elements
            (- C(TRI_DIAG_SP%ROW(I) - 1) * WEIGHT(TRI_DIAG_SP%ROW(I) - 1) &
            + D(TRI_DIAG_SP%ROW(I) - 1)/dvm(TRI_DIAG_SP%ROW(I)))) /&
             (V_GRID(TRI_DIAG_SP%ROW(I)) ** 2 * V_GRID_WIDTH(TRI_DIAG_SP%ROW(I)))

      LOCAL_M%VALUE(MARKER_EN_EL_0(P) + I) = &                                      !Send value to global matrix
                          LOCAL_M%VALUE(MARKER_EN_EL_0(P) + I) + VAL

    END DO

  END SUBROUTINE FILL_EN_EL_00
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Builds e-n collision integral for higher harmonics
  SUBROUTINE FILL_EN_EL_L(n_n, L,P,H,N)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n !< Neutral state density

    INTEGER, INTENT(IN) :: L !< Harmonic L-number
    INTEGER, INTENT(IN) :: P, & !< Current position
                           H, & !< Harmonic index of current harmonic
                           N  !< Atomic state number
    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: SIGMA(NUM_V + 1), &
              VAL !< Current value of matrix element to be passed to global matrix

    SIGMA = SIGMA_L_EL(L,:,N)

!Fill matrix according to simple diagonal pattern and adequate markers
    DO I = 1, NUM_V

      VAL = - COLL_EN_L * n_n * V_GRID(I) * SIGMA(I)

      LOCAL_M%VALUE(MARKER_EN_EL_L(P,H) + I) = LOCAL_M%VALUE(MARKER_EN_EL_L(P,H) + I) + VAL !Send value to global matrix

    END DO

  END SUBROUTINE FILL_EN_EL_L
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calls individual electron impact and de-excitation excitation collision integral matrix builders and adds submatrices to get total excitation submatrix
  SUBROUTINE FILL_EN_EX(n_n, L,P, H)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n(NUM_NEUTRALS) !< Vector of neutral state densities

    INTEGER, INTENT(IN) :: L !< Harmonic L-number
    INTEGER, INTENT(IN) :: P, & !< Current position
                           H !< Harmonic index of current harmonic
    INTEGER :: I, J

!Call individual transition (de)excitation submatrix builders
    DO I = 1, NUM_NEUTRALS

      DO J = I + 1, NUM_NEUTRALS

        CALL FILL_EN_EX_NM(n_n(I), I, J, L,P, H)

      END DO

      DO J = 1, I - 1

        CALL FILL_EN_DEEX_NM(n_n(I), I, J, L,P, H)

      END DO

    END DO

  END SUBROUTINE FILL_EN_EX
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills individual e-n excitation reaction submatrices
  SUBROUTINE FILL_EN_EX_NM(n_n, nl, nu, L, P, H)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n !< Density of lower state
    INTEGER, INTENT(IN) :: nl, & !< Lower state number
                           nu, & !< Upper/excited state number
                           L !< Harmonic L-number

    INTEGER, INTENT(IN) :: P, & !< Current position
                           H !< Harmonic index of current harmonic

    INTEGER :: I, ROW, COLUMN
    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    VAL = 0


    IF (L .GT. 0) THEN

!Fill matrix according to adequate sparsity pattern and markers
      DO I = 1, INEL_SP(nl,nu)%N_NZ

          ROW = INEL_SP(nl,nu)%ROW(I)
          COLUMN = INEL_SP(nl,nu)%COL(I)

          VAL = - COLL_EN_L * n_n * V_GRID(ROW) * SIGMA_L(0, ROW , nl, nu) * IN_DAT(nl,nu)%E(I)              !Diagonal elements
          VAL = VAL  + COLL_EN_L * n_n * V_GRID(COLUMN) * &
                (SIGMA_L(0, COLUMN , nl, nu) - SIGMA_L(L,COLUMN, nl, nu)) * IN_DAT(nl,nu)%W(ROW,COLUMN)  !Absorbing elements

          LOCAL_M%VALUE(INEL_M(P,H,nl,nu)%LIST(I)) = LOCAL_M%VALUE(INEL_M(P,H,nl,nu)%LIST(I)) + VAL !Send value to global matrix

      END DO

    ELSE

!Fill matrix according to adequate sparsity pattern and markers
      DO I = 1, INEL_SP(nl,nu)%N_NZ

          ROW = INEL_SP(nl,nu)%ROW(I)
          COLUMN = INEL_SP(nl,nu)%COL(I)

          VAL = - COLL_EN_L * n_n * V_GRID(ROW) * SIGMA_L(0, ROW , nl, nu) * IN_DAT(nl,nu)%E(I)              !Diagonal elements
          VAL = VAL + COLL_EN_L * n_n * V_GRID(COLUMN) * &
                SIGMA_L(0, COLUMN , nl, nu) * IN_DAT(nl,nu)%W(ROW,COLUMN)  !Absorbing elements

          LOCAL_M%VALUE(INEL_M(P,H,nl,nu)%LIST(I)) = LOCAL_M%VALUE(INEL_M(P,H,nl,nu)%LIST(I)) + VAL !Send value to global matrix

      END DO

    END IF

  END SUBROUTINE FILL_EN_EX_NM
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills individual e-n de-excitation reaction submatrices
  SUBROUTINE FILL_EN_DEEX_NM(n_n, nu, nl, L, P, H)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n !< Density of state nu

    INTEGER, INTENT(IN) :: nl, & !< Lower state number
                           nu, & !< Upper/excited state number
                           L !< Harmonic L-number

    INTEGER, INTENT(IN) :: P, & !< Current position
                           H !< Harmonic index of current harmonic

    INTEGER :: I, ROW, COLUMN
    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    VAL = 0

    IF (L .GT. 0) THEN

!Fill matrix according to adequate sparsity pattern and markers
      DO I = 1, INEL_SP(nu,nl)%N_NZ

          ROW = INEL_SP(nu,nl)%ROW(I)
          COLUMN = INEL_SP(nu,nl)%COL(I)

          VAL = - COLL_EN_L * n_n * V_GRID(ROW) * SIGMA_L_DEEX(0,nu,nl, P,ROW) * IN_DAT(nu,nl)%E(I)            !Diagonal elements
          VAL = VAL + COLL_EN_L * n_n * V_GRID(COLUMN) * &
                (SIGMA_L_DEEX(0,nu,nl, P,COLUMN) - SIGMA_L_DEEX(L,nu,nl, P,COLUMN)) *  IN_DAT(nu,nl)%W(ROW,COLUMN)  !Absorbing elements

          LOCAL_M%VALUE(INEL_M(P,H,nu,nl)%LIST(I)) = LOCAL_M%VALUE(INEL_M(P,H,nu,nl)%LIST(I)) + VAL  !Send value to global matrix

      END DO

    ELSE

      DO I = 1, INEL_SP(nu,nl)%N_NZ

          ROW = INEL_SP(nu,nl)%ROW(I)
          COLUMN = INEL_SP(nu,nl)%COL(I)

          VAL = - COLL_EN_L  *n_n* V_GRID(ROW) * SIGMA_L_DEEX(0,nu,nl, P,ROW) * IN_DAT(nu,nl)%E(I)              !Diagonal elements
          VAL = VAL + COLL_EN_L *n_n * V_GRID(COLUMN) * &
                SIGMA_L_DEEX(0,nu,nl, P,COLUMN) *  IN_DAT(nu,nl)%W(ROW,COLUMN)  !Absorbing elements

          LOCAL_M%VALUE(INEL_M(P,H,nu,nl)%LIST(I)) = LOCAL_M%VALUE(INEL_M(P,H,nu,nl)%LIST(I)) + VAL  !Send value to global matrix

      END DO

    END IF

  END SUBROUTINE FILL_EN_DEEX_NM
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills individual e-n recombination reaction submatrices
  SUBROUTINE FILL_EN_RECOMB(n_e,n_i, L, P, H)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_e, & !< Lagged electron density
                                n_i !< Lagged ion density

    INTEGER, INTENT(IN) :: L !< Harmonic L-number

    INTEGER, INTENT(IN) :: P, & !< Current position
                           H !< Harmonic index of current harmonic

    REAL(KIND=HIGH_PREC_REAL), DIMENSION(NUM_V + 1) :: SIGMA_0, & !< Total collision cross-section
                                    SIGMA_L1 !< L moment of collision cross-section

!Particular inelastic weights
    REAL(KIND=HIGH_PREC_REAL):: W(NUM_V,NUM_V)

    INTEGER :: I, J, ROW, COLUMN
    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix

    DO J = 1, NUM_NEUTRALS

      SIGMA_0 = SIGMA_L_RECOMB(0,J, P,:)                            !Get total collision cross-section

      IF (L .GT. 0) THEN

        SIGMA_L1 = SIGMA_L_RECOMB(L, J, P,:)                         !Get L moment of collision cross-section if L>0

      ELSE

        SIGMA_L1 = 0

      END IF

!Grab inelastic grid data
      W = IN_DAT(0,J)%W

!Fill matrix according to adequate sparsity pattern and markers
      DO I = 1, INEL_SP(0,J)%N_NZ

          ROW = INEL_SP(0,J)%ROW(I)
          COLUMN = INEL_SP(0,J)%COL(I)

          VAL = - DENSITY_0 * DE_BROGLIE_L3 * COLL_EN_L * n_e * n_i * V_GRID(ROW) * SIGMA_0(ROW) * IN_DAT(0,J)%E(I)              !Diagonal elements
          VAL = VAL + DENSITY_0 * DE_BROGLIE_L3 * COLL_EN_L * n_e * n_i  * V_GRID(COLUMN) * &
                (SIGMA_0(COLUMN) - SIGMA_L1(COLUMN)) * W(ROW,COLUMN)  !Absorbing elements

          LOCAL_M%VALUE(INEL_M(P,H,0,J)%LIST(I)) = LOCAL_M%VALUE(INEL_M(P,H,0,J)%LIST(I)) + VAL !Send value to global matrix

      END DO

    END DO

  END SUBROUTINE FILL_EN_RECOMB
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Builds e-n impact ionization submatrix
  SUBROUTINE FILL_EN_ION (n_n, L, P, H)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_n(NUM_NEUTRALS)!< Density vector of exicted states
    INTEGER, INTENT(IN) :: P, & !< Current position
                           H !< Harmonic index of current harmonic

    INTEGER, INTENT(IN) :: L

    INTEGER :: I
    REAL(KIND=HIGH_PREC_REAL) :: VAL !< Current value of matrix element to be passed to global matrix


!Fill excitation-like elements
    DO I = 1, NUM_NEUTRALS                                                      !Add individual excitation-like ionization submatrices

      CALL FILL_EN_EX_NM(n_n(I),I,0,L, P,H)

    END DO

    IF (L .EQ. 0) THEN                                                          !Add the secondary slow electron production matrix

      DO I = 1, NUM_NEUTRALS

        VAL = ION_RATES(I,P) / (4.00D00 * PI * V_GRID(1) ** 2 *dv)

        LOCAL_M%VALUE(MARKER_SEC_EL(P,I)) = LOCAL_M%VALUE(MARKER_SEC_EL(P,I)) + VAL !Send value to global matrix

      END DO

    END IF

  END SUBROUTINE FILL_EN_ION
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_NEUTRAL_COLL
!-------------------------------------------------------------------------------------------------------------------------------------
