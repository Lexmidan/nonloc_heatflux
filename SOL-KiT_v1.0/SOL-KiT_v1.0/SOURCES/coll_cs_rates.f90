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
! The subroutine p_polynomial_value is distributed under the GNU LGPL license. See below.
!-------------------------------------------------------------------------------------------------------------------------------------
! The subroutine e1xa is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged. See below.
!-------------------------------------------------------------------------------------------------------------------------------------
!> Contains all cross-section and rate data and routines
MODULE COLL_CS_RATES

  USE GRID
  USE NORMALIZATION
  USE NEUT_AND_HEAT
  USE INPUT
  USE INEL_GRID
  USE MATRIX_DATA
  USE MPI
  USE VAR_KIND_DATA

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: SIGMA_DAT(:,:,:), & !< Raw integral cross-section data
                         SIGMA_L(:,:,:,:), &  !< Legendre moments of cross-section data
                         SIGMA_L_EL(:,:,:), & !< Legendre moments of elastic cross-section data
                         SIGMA_L_DEEX(:,:,:,:,:), & !< Current deexcitation cross-section
                         SIGMA_L_RECOMB(:,:,:,:) !< Current recombination cross-section

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: DIPOLE_TRANSITION(:,:), & !< Dipole transition probability data
                         LEG_P(:,:), &  !< Legendre polynomials
                         ION_RATES(:,:), & !< Current ionization rates
                         TB_RECOMB_RATES(:,:) , & !< Current recombination rates
                         EX_RATES(:,:,:), & !< Current excitation rates
                         DEEX_RATES(:,:,:) !< Current de-excitaiton rates
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes cross-section data
  SUBROUTINE INIT_SIGMA

    IMPLICIT NONE

    INTEGER :: I, J

    ALLOCATE(SIGMA_DAT(NUM_V + 1, NUM_NEUTRALS, 0:NUM_NEUTRALS))

    CALL INIT_LEG_P                                                             !Initialize Legendre polynomials

    SIGMA_DAT = 0

!Loop over all excitation processes
    DO I = 1, NUM_NEUTRALS

      DO J = I + 1, NUM_NEUTRALS

        IF ((I .EQ. 1) .AND. (J .EQ. 2)) THEN

          SIGMA_DAT(:,I,J) = SIGMA_EX_1_TO_2()

        ELSE IF (I .EQ. 1) THEN

          IF (J .LT. 6) THEN

            SIGMA_DAT(:,I,J) = SIGMA_EX_1_TO_3_4_5(J)

          ELSE

            SIGMA_DAT(:,I,J) = SIGMA_EX_1_TO_6PLUS(J)

          END IF

        ELSE IF ((I .EQ. 2) .AND. (J .EQ. 3)) THEN

          SIGMA_DAT(:,I,J) = SIGMA_EX_2_TO_3()

        ELSE

          SIGMA_DAT(:,I,J) = SIGMA_EX_2PLUS_TO_4PLUS(I,J)

        END IF

      END DO

    END DO

!Loop over all ionization processes
    DO I = 1, NUM_NEUTRALS

      IF  (I .LT. 4) THEN

        SIGMA_DAT(:,I,0) = SIGMA_ION_1_2_3(I)

      ELSE

        SIGMA_DAT(:,I,0) = SIGMA_ION_4PLUS(I)

      END IF

    END DO

    CALL INIT_SIGMA_L                                                           !Initialize Legendre moments of cross-section data

  END SUBROUTINE INIT_SIGMA
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates the deexcitation cross-section based on detailed balance
  FUNCTION SIGMA_DEEX(L,I,J,T)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_DEEX(NUM_V + 1)

    INTEGER, INTENT(IN) :: I, & !< Higher excited level
                           J, & !< Lower level
                           L !< Harmonic l number

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T !< Lagged temperature

    INTEGER :: M, N

    REAL(KIND=DEFAULT_REAL) :: SIG(NUM_V+1) , E

    SIG = 0.00D00

    E = ION_POT_H * (1.00D00/REAL(J**2,KIND=DEFAULT_REAL) - 1.00D00/REAL(I**2,KIND=DEFAULT_REAL))

    DO M = 1, NUM_V

      DO N = 1, INEL_EM(M,J,I)%NE

        SIG(M) = SIG(M) + REAL(J**2,KIND=DEFAULT_REAL)/REAL(I**2,KIND=DEFAULT_REAL) * EXP(E/T) * &
                IN_DAT(J,I)%W(M,INEL_EM(M,J,I)%E(N)) * SIGMA_L(L,INEL_EM(M,J,I)%E(N),J,I) * V_GRID(INEL_EM(M,J,I)%E(N)) &
                / V_GRID(M) * EXP(-(V_GRID(INEL_EM(M,J,I)%E(N)) ** 2 - V_GRID(M) ** 2)/T) &
                 * INEL_ABS(INEL_EM(M,J,I)%E(N),J,I)%EMIT_SW

      END DO

    END DO
    SIGMA_DEEX = SIG

  END FUNCTION SIGMA_DEEX
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates the 3B-recombination cross-section based on detailed balance
  FUNCTION SIGMA_RECOMB(L,J,T)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_RECOMB(NUM_V + 1)

    INTEGER, INTENT(IN) :: J, & !< Lower level
                           L !< Harmonic l number

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T !< Lagged temperature

    INTEGER :: M, N

    REAL(KIND=DEFAULT_REAL) :: SIG(NUM_V+1) , E

    SIG = 0.00D00

    E = ION_POT_H /REAL(J**2,KIND=DEFAULT_REAL)

    DO M = 1, NUM_V

      DO N = 1, INEL_EM(M,J,0)%NE

        SIG(M) = SIG(M) + REAL(J**2,KIND=DEFAULT_REAL)/(SQRT(T) ** 3) * EXP(E/T) * &
                IN_DAT(J,0)%W(M,INEL_EM(M,J,0)%E(N)) * SIGMA_L(L,INEL_EM(M,J,0)%E(N),J,0) * V_GRID(INEL_EM(M,J,0)%E(N)) &
                / V_GRID(M) * EXP(-(V_GRID(INEL_EM(M,J,0)%E(N)) ** 2 - V_GRID(M) ** 2)/T) &
                * INEL_ABS(INEL_EM(M,J,0)%E(N),J,0)%EMIT_SW

      END DO

    END DO

    SIGMA_RECOMB = SIG

  END FUNCTION SIGMA_RECOMB
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fetches cross-section with given ID and takes the L-th angular moment
  FUNCTION SIGMA_GET(L,ID)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_GET(NUM_V + 1)
    INTEGER, INTENT(IN) :: L, & !< L-number
                           ID !< Process ID

    REAL(KIND=DEFAULT_REAL) :: S(NUM_V + 1, NUM_CS_ANGLES), &  !< Angle resolved cross-section data (differential cross-sections)
              SIG(NUM_V + 1), &
              SIG_SUM(NUM_V + 1)

    INTEGER :: I, J

    IF (ID .EQ. 1) THEN

      SIG = SIGMA_EL                                                              !Model elastic scattering

    ELSE

      SIG = SIGMA_DAT(:,MOD(ID, 10000)/100, MOD(ID, 100))                       !Takes data from stored cross-section

    END IF

!Set differential cross-section using integral cross-section and assuming isotropic scattering

      DO I = 1, NUM_V + 1

        DO J = 1, NUM_CS_ANGLES

          S(I,J) = SIG(I) / (4.00D00 * PI)

        END DO

      END DO

      SIG_SUM = 0

      DO I = 1, NUM_V + 1

        DO J = 1, NUM_CS_ANGLES

          IF (J .EQ. 1) THEN

            IF (L .EQ. 0) THEN

              SIG_SUM(I) = SIG_SUM(I) +  PI * S(I,J) * SIN(REAL(J,KIND=DEFAULT_REAL) &
              * PI / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL)) * PI / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL)

            ELSE

              SIG_SUM(I) = SIG_SUM(I) +  PI * S(I,J) * SIN(REAL(J,KIND=DEFAULT_REAL) &
              * PI / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL)) * PI / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL) * (1.00D00 - LEG_P(J,L))

            END IF

          ELSE

            IF (L .EQ. 0) THEN

               SIG_SUM(I) = SIG_SUM(I) + PI * (&
                            S(I,J) * SIN(REAL(J,KIND=DEFAULT_REAL) * PI/ REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL)) &
                            + S(I,J - 1) * SIN(REAL(J - 1,KIND=DEFAULT_REAL) * PI &
                            / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL))) * PI / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL)

            ELSE

              SIG_SUM(I) = SIG_SUM(I) + PI * (&
                           S(I,J) * SIN(REAL(J,KIND=DEFAULT_REAL) * PI/ &
                           REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL)) * (1.00D00 - LEG_P(J,L)) &
                           + S(I,J - 1) * SIN(REAL(J - 1,KIND=DEFAULT_REAL) * PI &
                           / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL))&
                            * (1.00D00- LEG_P(J - 1,L))) * PI / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL)

            END IF

          END IF

        END DO

      END DO

    SIGMA_GET = SIG_SUM

  END FUNCTION SIGMA_GET
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes Legendre Polynomials
  SUBROUTINE INIT_LEG_P

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: X(NUM_CS_ANGLES)

    INTEGER :: I

!Set angle space
    DO I = 1, NUM_CS_ANGLES

      X(I) = COS(REAL(I,KIND=DEFAULT_REAL) * PI / REAL(NUM_CS_ANGLES,KIND=DEFAULT_REAL))

    END DO

    ALLOCATE(LEG_P(NUM_CS_ANGLES, 0:L_MAX))

    CALL p_polynomial_value(NUM_CS_ANGLES, L_MAX, X, LEG_P)                     !Call Legendre polynomial calculator subroutine

  END SUBROUTINE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates collisional rate for process with given ID and using the distribution function f_0
  REAL(KIND=DEFAULT_REAL) FUNCTION COLL_RATE(f_0,ID,K)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_0(NUM_V) !< f_0^0
    INTEGER, INTENT(IN) :: ID !< Process ID

    INTEGER, INTENT(IN), OPTIONAL :: K !< Position

    INTEGER :: I

    REAL(KIND=DEFAULT_REAL) :: RATE

    RATE = 0

    IF ((MOD(ID, 100) .GT. 0) .AND. (MOD(ID, 100) .LT. MOD(ID, 10000)/100)) THEN

      DO I = 1, NUM_V

        RATE = RATE + 4.00D00 * PI * COLL_EN_L * V_GRID(I) ** 3 * f_0(I) &
        * SIGMA_L_DEEX(0,MOD(ID, 10000)/100,MOD(ID, 100),K,I)  * V_GRID_WIDTH(I)

      END DO

    ELSE IF (MOD(ID, 10000)/100 .EQ. 0) THEN

      DO I = 1, NUM_V

        RATE = RATE + 4.00D00 * PI * COLL_EN_L * V_GRID(I) ** 3 * f_0(I) &
         * SIGMA_L_RECOMB(0,MOD(ID, 100),K,I) * V_GRID_WIDTH(I)

      END DO

    ELSE

      DO I = 1, NUM_V

        RATE = RATE + 4.00D00 * PI * COLL_EN_L * V_GRID(I) ** 3 * f_0(I) &
        * SIGMA_L(0,I,MOD(ID, 10000)/100, MOD(ID, 100)) * V_GRID_WIDTH(I) &
        * INEL_ABS(I,MOD(ID, 10000)/100,MOD(ID, 100))%EMIT_SW

      END DO

    END IF


    COLL_RATE = RATE

  END FUNCTION COLL_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates recombination rate for radiative recombination
  REAL(KIND=DEFAULT_REAL) FUNCTION RECOMB_RATE(n, T)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T !< Lagged electron temperature
    INTEGER, INTENT(IN) :: n !< Final recombined state


    REAL(KIND=DEFAULT_REAL), PARAMETER :: Ry = 13.6D00 !< Rydberg constant

    REAL(KIND=DEFAULT_REAL) :: T_k, BETA, E1, I

    T_k = T * TEMP_0_EV                                                         !Temperature in eV
    I = Ry / REAL(n ** 2,KIND=DEFAULT_REAL)                                                             !Ionization potential in eV
    BETA = I / T_k

!Radiative recombination

      IF (n .EQ. 1) THEN

        RECOMB_RATE = 3.92D00 * SQRT(I/Ry) * BETA ** (3.0D00/2.0D00)/(BETA + 0.35D00) * 1.00D-20 / RAD_REC_BETA_0

      ELSE IF (n .EQ. 2) THEN

        RECOMB_RATE = (2.47D00 * SQRT(I/Ry) * BETA ** (3.0D00/2.0D00)/(BETA + 0.12D00) &
                      + 6.22D00 * SQRT(I/Ry) * BETA ** (3.0D00/2.0D00)/(BETA + 0.61D00)) * 1.00D-20 / RAD_REC_BETA_0

      ELSE

        CALL E1XA(BETA, E1)                                                     !Call exponential integral calculator subroutine

        RECOMB_RATE = 5.201D00 * BETA ** (3.0D00/2.0D00) * E1 * EXP(BETA) * 1.00D-20 / RAD_REC_BETA_0

      END IF

  END FUNCTION RECOMB_RATE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!*************************************************************************************************************************************
!Cross-section routines (cf. Janev)
!*************************************************************************************************************************************
  FUNCTION SIGMA_EX_1_TO_2()

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_EX_1_TO_2(NUM_V + 1)
    REAL(KIND=DEFAULT_REAL), PARAMETER :: deltaE = 10.20D00, &
                         sigma0 = 5.984D00, &
                         a = 0.228D00, &
                         b = 0.1865D00, &
                         c = 0.5025D00, &
                         A0 = 4.4979D00, &
                         Av(1:5) = (/ 1.4182D00, -20.877D00, 49.735D00, -46.249D00, 17.442D00/)

    REAL(KIND=DEFAULT_REAL) :: E(NUM_V + 1), &
              x(NUM_V + 1), &
              s(NUM_V + 1)

    INTEGER :: I, J

    s = 0

    DO I = 1, NUM_V

      E(I) = TEMP_0_EV * V_GRID(I) ** 2

    END DO

    E(NUM_V + 1) = TEMP_0_EV * (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1)) ** 2

    x = E / deltaE

    DO I = 1, NUM_V + 1

      IF (E(I) .GT. deltaE) THEN

        IF (E(I) .LT. 11.56D00) THEN

          s(I) = a + b * (E(I) - deltaE)

        ELSE IF (E(I) .LT. 12.23D00) THEN

          s(I) = c

        ELSE

          s(I) = sigma0 * A0 * LOG(x(I)) / (deltaE * x(I))

          DO J = 1, 5

            s(I) = s(I) + sigma0 * Av(J) / (deltaE * x(I) ** J)

          END DO

        END IF

      END IF

    END DO

    SIGMA_EX_1_TO_2 = s * 1.00D-20 / SIGMA_0

  END FUNCTION SIGMA_EX_1_TO_2
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION SIGMA_EX_1_TO_3_4_5(n)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_EX_1_TO_3_4_5(NUM_V + 1)
    INTEGER, INTENT(IN) :: n

    REAL(KIND=DEFAULT_REAL), PARAMETER :: s0 = 5.984D00
    REAL(KIND=DEFAULT_REAL) :: E(NUM_V +1), &
              x(NUM_V + 1), &
              s(NUM_V + 1), &
              deltaE, &
              alpha, &
              A0, &
              A(4)

    INTEGER :: I, J
    s = 0

    DO I = 1, NUM_V

      E(I) = TEMP_0_EV * V_GRID(I) ** 2

    END DO

    E(NUM_V + 1) = TEMP_0_EV * (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1)) ** 2

    IF (n .EQ. 3) THEN

      deltaE = 12.09D00
      alpha = 0.38277D00
      A0 = 0.75448D00
      A(1) = 0.42956D00
      A(2)= -0.58288D00
      A(3) = 1.0693D00
      A(4) = 0.00D00

    ELSE IF (n .EQ. 4) THEN

      deltaE = 12.75D00
      alpha = 0.41844D00
      A0 = 0.24300D00
      A(1) = 0.24846D00
      A(2)= 0.19701D00
      A(3) = 0.00D00
      A(4) = 0.00D00

    ELSE IF (n .EQ. 5) THEN

      deltaE = 13.06D00
      alpha = 0.45929D00
      A0 = 0.11508D00
      A(1) = 0.13092D00
      A(2)= 0.23581D00
      A(3) = 0.00D00
      A(4) = 0.00D00

    END IF

    x = E / deltaE

    DO I = 1, NUM_V + 1

      IF (X(I) .GT. 1) THEN

      s(I) = s0 * ((1.00D00 - 1.00D00/x(I)) ** alpha) * A0 * LOG(x(I)) / (deltaE * x(I))

      DO J = 1, 4

        s(I) = s(I) + s0 * ((1.00D00 - 1.00D00/x(I)) ** alpha) * A(J)  / (deltaE * x(I) ** J)

      END DO

      END IF

    END DO

    SIGMA_EX_1_TO_3_4_5 = s * 1.00D-20 / SIGMA_0

  END FUNCTION SIGMA_EX_1_TO_3_4_5
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  REAL(KIND=DEFAULT_REAL) FUNCTION JOHNSON_y(n,m)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, m

    JOHNSON_y = 1.00D00 - (REAL(n, KIND=DEFAULT_REAL)/REAL(m, KIND=DEFAULT_REAL)) ** 2

  END FUNCTION JOHNSON_y
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  REAL(KIND=DEFAULT_REAL) FUNCTION g_factor(n,ynm)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: ynm

    REAL(KIND=DEFAULT_REAL) :: g0, g1, g2

    IF (n .eq. 1) THEN

      g0 = 1.1330D00
      g1 = -0.4059D00
      g2 = 0.0714D00

    ELSE IF (n .eq. 2) THEN

      g0 = 1.0785D00
      g1 = -0.2319D00
      g2 = 0.0295D00

    ELSE

      g0 = 0.9935D00 + 0.2328D00  / REAL(n,KIND=DEFAULT_REAL) - 0.1296D00  / REAL(n ** 2,KIND=DEFAULT_REAL)
      g1 = -(0.6282D00 - 0.5598D00/REAL(n,KIND=DEFAULT_REAL) + 0.5299D00/REAL(n ** 2,KIND=DEFAULT_REAL)) / REAL(n,KIND=DEFAULT_REAL)
      g2 = (0.3887D00 -1.181D00/REAL(n,KIND=DEFAULT_REAL) + 1.1470D00/REAL(n ** 2,KIND=DEFAULT_REAL))/REAL(n ** 2,KIND=DEFAULT_REAL)

    END IF

    g_factor = g0 + g1/ynm + g2/ynm ** 2

  END FUNCTION g_factor
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  REAL(KIND=DEFAULT_REAL) FUNCTION f_osc(n,m)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, m

    f_osc = 32.00D00 * REAL(n,KIND=DEFAULT_REAL) * g_factor(n,JOHNSON_y(n,m)) &
    / (JOHNSON_y(n,m) ** 3 * REAL(m,KIND=DEFAULT_REAL) ** 3 * 3.00D00 * SQRT(3.00D00) * PI)

  END FUNCTION f_osc
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION SIGMA_EX_1_TO_6PLUS(n)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_EX_1_TO_6PLUS(NUM_V + 1)

    INTEGER, INTENT(IN) :: n
    REAL(KIND=DEFAULT_REAL) :: E(NUM_V + 1), &
              x(NUM_V + 1), &
              deltaE, &
              A, &
              r, &
              B, &
              s(NUM_V + 1)

    INTEGER :: I

    s = 0

    DO I = 1, NUM_V

      E(I) = TEMP_0_EV * V_GRID(I) ** 2

    END DO

    E(NUM_V + 1) = TEMP_0_EV * (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1)) ** 2

    deltaE = 13.6D00 * (1.00D00 - 1.00D00/REAL(n ** 2,KIND=DEFAULT_REAL))
    x = E / deltaE
    A = 2.00D00 * f_osc(1, n) / JOHNSON_y(1, n)
    r = 0.45D00
    B = 4.00D00 * (1.00D00 + 4.00D00 /(3.00D00 * JOHNSON_y(1, n)) &
        - 0.603D00/ JOHNSON_y(1, n) ** 2) / (REAL(n ** 3,KIND=DEFAULT_REAL) * JOHNSON_y(1, n) ** 2)

    DO I = 1, NUM_V + 1

      IF (E(I) .GT. deltaE) THEN

        s(I) = 1.76D00 * (1.00D00 - EXP(-r * JOHNSON_y(1,n) * x(I))) * (A * (LOG(x(I)) + 1.00D00 / (2.00D00 * x(I))) + &
               (B - A * LOG(2.00D00/JOHNSON_y(1,n))) * (1.00D00 - 1.00D00/x(I))) / (JOHNSON_y(1,n) * x(I))

      END IF

    END DO

    SIGMA_EX_1_TO_6PLUS = s * 1.00D-20 / SIGMA_0

  END FUNCTION SIGMA_EX_1_TO_6PLUS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION SIGMA_EX_2_TO_3()

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_EX_2_TO_3(NUM_V + 1)

    REAL(KIND=DEFAULT_REAL), PARAMETER :: s0 = 5.984D00, &
                         alpha = 1.3196D00, &
                         A0 = 38.906D00, &
                         deltaE = 13.6 * (1.00D00/4.00D00 - 1.00D00/9.00D00), &
                         A(1:4) = (/ 5.2373D00, 119.25D00, -595.39D00, 816.71D00 /)

    INTEGER :: I, J
    REAL(KIND=DEFAULT_REAL) :: s(NUM_V + 1), &
              E(NUM_V + 1), &
              x(NUM_V + 1)

    s = 0

    DO I = 1, NUM_V

      E(I) = TEMP_0_EV * V_GRID(I) ** 2

    END DO

    E(NUM_V + 1) = TEMP_0_EV * (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1)) ** 2

    x = E / deltaE

    DO I = 1, NUM_V + 1

      IF (x(I) .GT. 1) THEN

        s(I) = s0 * ((1.00D00 - 1.00D00/x(I)) ** alpha) * A0 * LOG(x(I)) / (deltaE * x(I))

        DO J = 1, 4

          s(I) = s(I) + s0 * ((1.00D00 - 1.00D00/x(I)) ** alpha) * A(J)  / (deltaE * x(I) ** J)

        END DO

      END IF

    END DO

    SIGMA_EX_2_TO_3 = s * 1.00D-20 / SIGMA_0



  END FUNCTION SIGMA_EX_2_TO_3
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION SIGMA_EX_2PLUS_TO_4PLUS(n,m)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_EX_2PLUS_TO_4PLUS(NUM_V + 1)
    INTEGER, INTENT(IN) :: n, m

    REAL(KIND=DEFAULT_REAL) :: E(NUM_V + 1), &
              x(NUM_V + 1), &
              s(NUM_V + 1), &
              r, &
              deltaE, &
              A, &
              B, &
              bb

    INTEGER :: I

    s = 0

    DO I = 1, NUM_V

      E(I) = TEMP_0_EV * V_GRID(I) ** 2

    END DO

    E(NUM_V + 1) = TEMP_0_EV * (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1)) ** 2

    deltaE = 13.6D00 * (1/REAL(n ** 2,KIND=DEFAULT_REAL) - 1/REAL(m ** 2,KIND=DEFAULT_REAL))
    x = E / deltaE
    r = 1.94D00/REAL(n,KIND=DEFAULT_REAL) ** 1.57D00
    A = 2.00D00 * REAL(n ** 2,KIND=DEFAULT_REAL) * f_osc(n,m) / JOHNSON_y(n,m)
    bb = (4.00D00 - 18.63D00/REAL(n,KIND=DEFAULT_REAL) + 36.24/REAL(n ** 2,KIND=DEFAULT_REAL)&
     - 28.09D00/REAL(n ** 3,KIND=DEFAULT_REAL)) / REAL(n,KIND=DEFAULT_REAL)
    B = 4.00D00 * REAL(n ** 4,KIND=DEFAULT_REAL)&
     * (1.00D00 + 4.00D00/(3.00D00 * JOHNSON_y(n,m)) + bb/JOHNSON_y(n,m) ** 2) &
     / (REAL(m ** 3,KIND=DEFAULT_REAL) * JOHNSON_y(n,m) ** 2)

    DO I = 1, NUM_V + 1

      IF (E(I) .GT. deltaE) THEN

        s(i) = 1.76D00 * REAL(n ** 2,KIND=DEFAULT_REAL) * (1.00D00 - EXP(- r * JOHNSON_y(n,m) * x(I))) * &
               (A * (LOG(x(I)) + 1.00D00/(2.00D00 * x(I))) &
               + (B - A * LOG(2.00D00 * REAL(n ** 2,KIND=DEFAULT_REAL) /JOHNSON_y(n,m)))&
                * (1.00D00 - 1.00D00/x(I))) / (JOHNSON_y(n,m) * x(I))

      END IF

    END DO

    SIGMA_EX_2PLUS_TO_4PLUS = s * 1.00D-20 / SIGMA_0

  END FUNCTION SIGMA_EX_2PLUS_TO_4PLUS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION SIGMA_ION_1_2_3(n)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_ION_1_2_3(NUM_V + 1)

    INTEGER, INTENT(IN) :: n

    REAL(KIND=DEFAULT_REAL) :: E(NUM_V + 1), &
              s(NUM_V + 1), &
              A0, &
              A(5), &
              I

    INTEGER :: J, K

    I = 13.6D00 / REAL(n ** 2,KIND=DEFAULT_REAL)

    s = 0

    DO J = 1, NUM_V

      E(J) = TEMP_0_EV * V_GRID(J) ** 2

    END DO

    E(NUM_V + 1) = TEMP_0_EV * (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1)) ** 2

    IF (n .EQ. 1) THEN

      A0 = 0.18450D00
      A(1:5) = (/ -0.032226D00, -0.034539D00, 1.4003D00, -2.8115D00, 2.2986D00 /)

    ELSE IF (n .EQ. 2) THEN

      A0 = 0.14784D00
      A(1:5) = (/ 0.0080871D00, -0.062270D00, 1.9414D00, - 2.1980D00, 0.95894D00 /)

    ELSE IF (n .EQ. 3) THEN

      A0 = 0.058463D00
      A(1:5) = (/ -0.051272D00, 0.85310D00, -0.57014D00, 0.76684D00, 0.00D00 /)

    END IF

    DO J = 1, NUM_V

      IF (E(J) .GT. I) THEN

        s(J) = 1.00D-17 * A0 * LOG(E(J)/I)/(I * E(J))

        DO K = 1, 5

          s(J) = S(J) + 1.00D-17 * A(K) * (1.00D00 - I/E(J)) ** K / (E(J) * I)

        END DO

      END IF

    END DO

    SIGMA_ION_1_2_3 = s / SIGMA_0

  END FUNCTION SIGMA_ION_1_2_3
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION SIGMA_ION_4PLUS(n)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: SIGMA_ION_4PLUS(NUM_V + 1)

    INTEGER, INTENT(IN) :: n

    REAL(KIND=DEFAULT_REAL) :: E(NUM_V + 1), &
              s(NUM_V + 1), &
              x(NUM_V + 1), &
              g0, g1, g2, &
              A, &
              r, &
              bb, &
              B, &
              I

    INTEGER :: J

    I = 13.6D00/REAL(n ** 2,KIND=DEFAULT_REAL)


    s = 0

    DO J = 1, NUM_V

      E(J) = TEMP_0_EV * V_GRID(J) ** 2

    END DO

    E(NUM_V + 1) = TEMP_0_EV * (2.00D00 * V_GRID(NUM_V) - V_GRID(NUM_V - 1)) ** 2

    x = E/I

    IF (n .eq. 1) THEN

      g0 = 1.1330D00
      g1 = -0.4059D00
      g2 = 0.0714D00

    ELSE IF (n .eq. 2) THEN

      g0 = 1.0785D00
      g1 = -0.2319D00
      g2 = 0.0295D00

    ELSE

      g0 = 0.9935D00 + 0.2328D00  / REAL(n,KIND=DEFAULT_REAL) - 0.1296D00  / REAL(n ** 2,KIND=DEFAULT_REAL)
      g1 = -(0.6282D00 - 0.5598D00/REAL(n,KIND=DEFAULT_REAL) + 0.5299D00/REAL(n ** 2,KIND=DEFAULT_REAL)) / REAL(n,KIND=DEFAULT_REAL)
      g2 = (0.3887D00 -1.181D00/REAL(n,KIND=DEFAULT_REAL) + 1.1470D00/REAL(n ** 2,KIND=DEFAULT_REAL))/REAL(n ** 2,KIND=DEFAULT_REAL)

    END IF

    A = 32.00D00 * n * (g0/3.00D00 + g1/4.00D00 + g2/5.00D00) / (3.00D00 * PI * SQRT(3.00D00))
    r = 1.94D00/ REAL(n,KIND=DEFAULT_REAL) ** 1.57D00
    bb = (4.00D00 - 18.63D00/REAL(n,KIND=DEFAULT_REAL) + 36.24D00/REAL(n ** 2 ,KIND=DEFAULT_REAL)&
    - 28.09D00/REAL(n ** 3,KIND=DEFAULT_REAL)) &
    / REAL(n,KIND=DEFAULT_REAL)
    B = 2.00D00 *REAL(n ** 2,KIND=DEFAULT_REAL) * (5.00D00 + bb)/3.00D00

    DO J = 1, NUM_V + 1

      IF (E(J) .GT. I) THEN

        s(J) = 1.76D00 * REAL(n ** 2,KIND=DEFAULT_REAL) * (1.00D00 - EXP(- r * x(J))) * &
               (A * LOG(x(J)) + (B - A * LOG(REAL(2*n ** 2,KIND=DEFAULT_REAL))) * (1.00D00 - 1.00D00/x(J)) ** 2) / x(J)

      END IF

    END DO

    SIGMA_ION_4PLUS = s * 1.00D-20 / SIGMA_0

  END FUNCTION SIGMA_ION_4PLUS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Load dipole transition probability data
  SUBROUTINE INIT_DIPOLE_TRANS_PROB

    IMPLICIT NONE

    INTEGER :: RANK, I, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    ALLOCATE(DIPOLE_TRANSITION(MAX(NUM_NEUTRALS,20),MAX(NUM_NEUTRALS,20)))

    IF (RANK .EQ. 0) CALL INPUT_DIPOL_TRANS(DIPOLE_TRANSITION,MAX(NUM_NEUTRALS,20))                   !Call dipole transition input

    DO I = 1, MAX(NUM_NEUTRALS,20)

      CALL MPI_Bcast(DIPOLE_TRANSITION(:,I),MAX(NUM_NEUTRALS,20),MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    END DO

  END SUBROUTINE INIT_DIPOLE_TRANS_PROB
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize Legendre moments of cross-section data
  SUBROUTINE INIT_SIGMA_L

    IMPLICIT NONE

    INTEGER :: I,J,K

    ALLOCATE(SIGMA_L(0:L_MAX,NUM_V+1,NUM_NEUTRALS,0:NUM_NEUTRALS))
    ALLOCATE(SIGMA_L_EL(0:L_MAX,NUM_V+1,NUM_NEUTRALS))
    ALLOCATE(SIGMA_L_DEEX(0:L_MAX,NUM_NEUTRALS,NUM_NEUTRALS,1:NUM_X,NUM_V+1))
    ALLOCATE(SIGMA_L_RECOMB(0:L_MAX,NUM_NEUTRALS,1:NUM_X,NUM_V+1))

    SIGMA_L = 0
    SIGMA_L_EL = 0
    SIGMA_L_DEEX = 0
    SIGMA_L_RECOMB = 0

    DO I = 0, L_MAX

      DO J = 1, NUM_NEUTRALS

        DO K = J + 1, NUM_NEUTRALS

          SIGMA_L(I,:,J,K) = SIGMA_GET(I, 10000 + J * 100 + K)

        END DO

        SIGMA_L(I,:,J,0) = SIGMA_GET(I, 10000 + J * 100)

        SIGMA_L_EL(I,:,J) = SIGMA_GET(I,1) * J ** 4

      END DO

    END DO

  END SUBROUTINE INIT_SIGMA_L
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update deexcitation cross section
  SUBROUTINE UPDATE_SIGMA_L_DEEX(T,ALL_X)

    IMPLICIT NONE

    INTEGER :: L,I,J,K,X_MIN, X_MAX
    REAL(KIND=DEFAULT_REAL) :: T(NUM_X)
    LOGICAL,OPTIONAL :: ALL_X

    X_MIN = MIN_X
    X_MAX = MAX_X

    IF (PRESENT(ALL_X)) THEN

      IF (ALL_X) THEN

        X_MIN = 1
        X_MAX = NUM_X

      END IF

    END IF

    DO K = X_MIN, X_MAX

      DO L = 0, L_MAX

        DO I = 1,NUM_NEUTRALS

          DO J = 1, I-1

            SIGMA_L_DEEX(L,I,J,K,:) = SIGMA_DEEX(L,I,J,T(K))

          END DO

        END DO

      END DO

    END DO

  END SUBROUTINE UPDATE_SIGMA_L_DEEX
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update recombination cross section
  SUBROUTINE UPDATE_SIGMA_L_RECOMB(T,ALL_X)

    IMPLICIT NONE

    INTEGER :: L,I,K,X_MIN, X_MAX
    REAL(KIND=DEFAULT_REAL) :: T(NUM_X)
    LOGICAL,OPTIONAL :: ALL_X

    X_MIN = MIN_X
    X_MAX = MAX_X

    IF (PRESENT(ALL_X)) THEN

      IF (ALL_X) THEN

        X_MIN = 1
        X_MAX = NUM_X

      END IF

    END IF

    DO K = X_MIN, X_MAX
      DO L = 0, L_MAX

        DO I = 1,NUM_NEUTRALS

          SIGMA_L_RECOMB(L,I,K,:) = SIGMA_RECOMB(L,I,T(K))

        END DO

      END DO

    END DO

  END SUBROUTINE UPDATE_SIGMA_L_RECOMB
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update ionization rates
  SUBROUTINE UPDATE_ION_RATES(f_lagged)

    IMPLICIT NONE

    INTEGER :: I,K
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(DIM_F)

    REAL(KIND=DEFAULT_REAL) :: F(NUM_V)

    IF (.NOT. ALLOCATED(ION_RATES)) ALLOCATE(ION_RATES(NUM_NEUTRALS,MIN_X:MAX_X))

    DO K = MIN_X, MAX_X

      F = f_lagged(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        ION_RATES(I,K) = COLL_RATE(F, 10000 + 100 * I)

      END DO

    END DO

  END SUBROUTINE UPDATE_ION_RATES
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update ionization rates
  SUBROUTINE UPDATE_TB_RECOMB_RATES(f_lagged)

    IMPLICIT NONE

    INTEGER :: I,K
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(DIM_F)

    REAL(KIND=DEFAULT_REAL) :: F(NUM_V)

    IF (.NOT. ALLOCATED(TB_RECOMB_RATES)) ALLOCATE(TB_RECOMB_RATES(NUM_NEUTRALS,MIN_X:MAX_X))

    DO K = MIN_X, MAX_X

      F = f_lagged(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        TB_RECOMB_RATES(I,K) = COLL_RATE(F, 10000 + I,K)

      END DO

    END DO

  END SUBROUTINE UPDATE_TB_RECOMB_RATES
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Update ionization rates
  SUBROUTINE UPDATE_EX_DEEX_RATES(f_lagged)

    IMPLICIT NONE

    INTEGER :: I,J,K
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(DIM_F)

    REAL(KIND=DEFAULT_REAL) :: F(NUM_V)

    IF (.NOT. ALLOCATED(EX_RATES)) ALLOCATE(EX_RATES(NUM_NEUTRALS,NUM_NEUTRALS,MIN_X:MAX_X))
    IF (.NOT. ALLOCATED(DEEX_RATES)) ALLOCATE(DEEX_RATES(NUM_NEUTRALS,NUM_NEUTRALS,MIN_X:MAX_X))

    DO K = MIN_X, MAX_X

      F = f_lagged(X_POS(K) + (NUM_H - 1) * NUM_V + 1: X_POS(K) + NUM_H  * NUM_V)

      DO I = 1,NUM_NEUTRALS

        DO J = I + 1, NUM_NEUTRALS

          EX_RATES(I,J,K) = COLL_RATE(F, 10000 + I * 100 + J)

        END DO

        DO J = 1, I - 1

          DEEX_RATES(I,J,K) = COLL_RATE(F, 10000 + I * 100 + J, K)

        END DO

      END DO

    END DO

  END SUBROUTINE UPDATE_EX_DEEX_RATES
!*************************************************************************************************************************************
!Special function routines
!*************************************************************************************************************************************
  subroutine p_polynomial_value ( m, n, x, v )

  !*****************************************************************************80
  !
  ! P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
  !
  !  Discussion:
  !
  !    P(n,1) = 1.
  !    P(n,-1) = (-1)^N.
  !    | P(n,x) | <= 1 in [-1,1].
  !
  !    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
  !    quadrature of the integral of a function F(X) with weight function 1
  !    over the interval [-1,1].
  !
  !    The Legendre polynomials are orthogonal under the inner product defined
  !    as integration from -1 to 1:
  !
  !      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
  !        = 0 if I =/= J
  !        = 2 / ( 2*I+1 ) if I = J.
  !
  !    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
  !
  !    A function F(X) defined on [-1,1] may be approximated by the series
  !      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
  !    where
  !      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
  !
  !    The formula is:
  !
  !      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
  !
  !  Differential equation:
  !
  !    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
  !
  !  First terms:
  !
  !    P( 0,x) =      1
  !    P( 1,x) =      1 X
  !    P( 2,x) = (    3 X^2 -       1)/2
  !    P( 3,x) = (    5 X^3 -     3 X)/2
  !    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
  !    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
  !    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
  !    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
  !    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
  !    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
  !    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
  !
  !  Recursion:
  !
  !    P(0,x) = 1
  !    P(1,x) = x
  !    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
  !
  !    P'(0,x) = 0
  !    P'(1,x) = 1
  !    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 March 2012
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Milton Abramowitz, Irene Stegun,
  !    Handbook of Mathematical Functions,
  !    National Bureau of Standards, 1964,
  !    ISBN: 0-486-61272-4,
  !    LC: QA47.A34.
  !
  !    Daniel Zwillinger, editor,
  !    CRC Standard Mathematical Tables and Formulae,
  !    30th Edition,
  !    CRC Press, 1996.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the number of evaluation points.
  !
  !    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
  !    Note that polynomials 0 through N will be evaluated.
  !
  !    Input, real ( kind = 8 ) X(M), the evaluation points.
  !
  !    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials
  !    of order 0 through N at the points X.
  !
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( KIND=DEFAULT_REAL ) v(m,0:n)
  real ( KIND=DEFAULT_REAL ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)

  do i = 2, n

    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )

  end do

  return
  end
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
  subroutine e1xa ( x, e1 )

  !*****************************************************************************80
  !
  ! E1XA computes the exponential integral E1(x).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
  !    they give permission to incorporate this routine into a user program
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    06 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) X, the argument.
  !
  !    Output, real ( kind = 8 ) E1, the function value.
  !
    implicit none

    real ( KIND=DEFAULT_REAL ) e1
    real ( KIND=DEFAULT_REAL ) es1
    real ( KIND=DEFAULT_REAL ) es2
    real ( KIND=DEFAULT_REAL ) x

    if ( x == 0.0D+00 ) then

      e1 = 1.0D+300

    else if ( x <= 1.0D+00 ) then

      e1 = - log ( x ) + (((( &
          1.07857D-03 * x &
        - 9.76004D-03 ) * x &
        + 5.519968D-02 ) * x &
        - 0.24991055D+00 ) * x &
        + 0.99999193D+00 ) * x &
        - 0.57721566D+00

    else

      es1 = ((( x &
        + 8.5733287401D+00 ) * x &
        +18.059016973D+00  ) * x &
        + 8.6347608925D+00 ) * x &
        + 0.2677737343D+00

      es2 = ((( x &
        +  9.5733223454D+00 ) * x &
        + 25.6329561486D+00 ) * x &
        + 21.0996530827D+00 ) * x &
        +  3.9584969228D+00

      e1 = exp ( - x ) / x * es1 / es2

    end if

    return
  end
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE COLL_CS_RATES
