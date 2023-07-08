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
!> Contains all discretization properties used in the treatment of inelastic collisions
MODULE INEL_GRID

  USE GRID
  USE NORMALIZATION
  USE NEUT_AND_HEAT
  USE MPI
  USE VAR_KIND_DATA

!> Contains inelastic collision mapping data for a point on the velocity grid and for given transition
  TYPE INEL_ABSORB

    INTEGER :: NA !< Number of absorbers for given velocity grid point

    INTEGER :: EMIT_SW !< 0 if less than 2 absorbers, 1 otherwise

    INTEGER, ALLOCATABLE :: A(:) !< List of absorbers for given grid point

    INTEGER :: P !< Number of pairs of absorber points

    INTEGER, ALLOCATABLE:: m1(:), & !< First pair point
                           m2(:) !< Second pair point

    REAL(KIND=HIGH_PREC_REAL), ALLOCATABLE :: W1(:), & !< Weights for first pair point
                                 W2(:) !< Weights for second pair point

  END TYPE INEL_ABSORB

!> Contains labels for connecting absorber and emmiter data
  TYPE LABEL

    INTEGER, ALLOCATABLE :: L(:) !< List of labels

  END TYPE LABEL

!> Contains data connecting points to their paired emmiters
  TYPE INEL_EMIT

    INTEGER :: NE !< Number of emitters for given velocity grid point

    INTEGER, ALLOCATABLE :: E(:), & !< List of emmiters for given grid point
                            P(:) !< Number of pairs for given emmiter in which grid point is present

    TYPE(LABEL), ALLOCATABLE :: PAIR_LABEL(:), & !< List of pairs in which grid point is present
                                ABS_LABEL(:) !< 1 if given grid point is the first pair point, 2 if second

  END TYPE INEL_EMIT

  TYPE (INEL_ABSORB), ALLOCATABLE :: INEL_ABS(:,:,:) !< Data structure containing emitter-absorber mapping

  TYPE (INEL_EMIT), ALLOCATABLE :: INEL_EM(:,:,:) !< Data structure containing absorber-emitter mapping

  INTEGER, ALLOCATABLE :: INEL_NNZ(:,:) !< Number of non-zero elements for given inelastic process

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes all discretization properties for inelastic collisions
  SUBROUTINE INIT_INEL_GRID

    IMPLICIT NONE

    INTEGER :: I, J, K,&
               n_1(NUM_V), m_1(NUM_V), m_2(NUM_V),& ! Closest emmiter/absorbers to ideal emmiter/absorber
               N_A(NUM_V)

    REAL(KIND=HIGH_PREC_REAL) :: E


    ALLOCATE(INEL_ABS(NUM_V,0:NUM_NEUTRALS,0:NUM_NEUTRALS))
    ALLOCATE(INEL_EM(NUM_V,0:NUM_NEUTRALS,0:NUM_NEUTRALS))

    ALLOCATE(INEL_NNZ(0:NUM_NEUTRALS,0:NUM_NEUTRALS))

    INEL_NNZ = 0;

!Loop through all inelastic processes
    DO I = 1, NUM_NEUTRALS

!Excitation loop
      DO J = I + 1, NUM_NEUTRALS

        n_1 = 0; m_1 = 0; m_2 = 0; N_A = 0

        E = ION_POT_H * (1.00D00 / REAL(I ** 2,KIND=HIGH_PREC_REAL) &
        - 1.00D00 / REAL(J ** 2,KIND=HIGH_PREC_REAL))                   !Excitation energy

!Direct sweep

        CALL INEL_SWEEP(E,n_1,ABS_COUNT = N_A)

!Inverse sweep

        CALL INEL_SWEEP(-E,m_1,SECOND_ABS=m_2)

        DO K = 1, NUM_V
!Finish count of number of absorbers

          CALL INEL_ABS_COUNT(K,N_A,m_1,m_2, IDEAL_EM = n_1)

          INEL_ABS(K,I,J)%EMIT_SW = 0

!Fill out data if 2 or more absorbers
          IF (N_A(K) .GT. 1) THEN

            CALL INEL_FILL_ABS(I,J,K,N_A,m_1,m_2,IDEAL_EM = n_1)

!Fill pairs and calculate weights
            CALL INEL_PARTITION_AND_WEIGHTS(I,J,K,N_A,m_1,m_2,E)

          END IF

        END DO

!Fill emitter data
        CALL INEL_FILL_EM(I,J)

!Calculate number of nonzero elements for particular transition
        CALL INEL_NNZ_COUNT(I,J)

      END DO
!Ionization processes

      n_1 = 0; m_1 = 0; m_2 = 0; N_A = 0
      E = ION_POT_H / REAL(I ** 2,KIND=HIGH_PREC_REAL)  + V_GRID_P(1) ** 2                                                  !Ionization potential corrected for finite first cell

!Direct sweep

        CALL INEL_SWEEP(E,n_1,ABS_COUNT = N_A)

!Inverse sweep

        CALL INEL_SWEEP(-E,m_1,SECOND_ABS=m_2)

      DO K = 1, NUM_V
!Finish count of number of absorbers
        CALL INEL_ABS_COUNT(K,N_A,m_1,m_2, IDEAL_EM = n_1)

          INEL_ABS(K,I,0)%EMIT_SW = 0

!Fill out data if 2 or more absorbers
        IF (N_A(K) .GT. 1) THEN

          CALL INEL_FILL_ABS(I,0,K,N_A,m_1,m_2,IDEAL_EM = n_1)

!Fill pairs and calculate weights
          CALL INEL_PARTITION_AND_WEIGHTS(I,0,K,N_A,m_1,m_2,E)

        END IF

      END DO

!Fill emitter data
        CALL INEL_FILL_EM(I,0)

!Calculate number of nonzero elements for particular transition
        CALL INEL_NNZ_COUNT(I,0)

!De-excitation loop
        DO J = 1, I - 1

          m_1 = 0; m_2 = 0; N_A = 0

          E = - ION_POT_H * (1.00D00 / REAL(J ** 2,KIND=HIGH_PREC_REAL) &
          - 1.00D00 / REAL(I ** 2,KIND=HIGH_PREC_REAL))                   !Excitation energy

!Inverse sweep

        CALL INEL_SWEEP(-E,m_1,SECOND_ABS=m_2)

          DO K = 1, NUM_V
!Finish count of number of absorbers
            CALL INEL_ABS_COUNT(K,N_A,m_1,m_2)

            INEL_ABS(K,I,J)%EMIT_SW = 0

!Fill out data if 2 or more absorbers
            IF (N_A(K) .GT. 1) THEN

              CALL INEL_FILL_ABS(I,J,K,N_A,m_1,m_2)

!Fill pairs and calculate weights
              CALL INEL_PARTITION_AND_WEIGHTS(I,J,K,N_A,m_1,m_2,E)

            END IF

          END DO

!Fill emitter data
          CALL INEL_FILL_EM(I,J)

!Calculate number of nonzero elements for particular transition
          CALL INEL_NNZ_COUNT(I,J)

        END DO

!Recombination

        m_1 = 0; m_2 = 0; N_A = 0

        E = - ION_POT_H / REAL(I ** 2,KIND=HIGH_PREC_REAL)  - V_GRID_P(1) ** 2                 !Ionization potential corrected for finite first cell

!Inverse sweep
        CALL INEL_SWEEP(-E,m_1,SECOND_ABS=m_2)

        DO K = 1, NUM_V

!Finish count of number of absorbers
          CALL INEL_ABS_COUNT(K,N_A,m_1,m_2)


          INEL_ABS(K,0,I)%EMIT_SW = 0

!Fill out data if 2 or more absorbers
          IF (N_A(K) .GT. 1) THEN

            CALL INEL_FILL_ABS(0,I,K,N_A,m_1,m_2)

!Fill pairs and calculate weights
            CALL INEL_PARTITION_AND_WEIGHTS(0,I,K,N_A,m_1,m_2,E)

          END IF

        END DO

!Fill emitter data
          CALL INEL_FILL_EM(0,I)

!Calculate number of nonzero elements for particular transition
          CALL INEL_NNZ_COUNT(0,I)

      END DO

  END SUBROUTINE INIT_INEL_GRID
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Performs a sweep with given transition energy and fills SWEEP_RES with obtained mapping, and SECOND_ABS with the second absorber mapping, if present
  SUBROUTINE INEL_SWEEP(TRANS_EN, SWEEP_RES, ABS_COUNT, SECOND_ABS)

    IMPLICIT NONE

    REAL(KIND=HIGH_PREC_REAL), INTENT(IN) :: TRANS_EN !< Transition energy (negative for superelastic collisions)
    INTEGER, INTENT(OUT) :: SWEEP_RES(NUM_V) !< The vector that holds the primary sweep result
    INTEGER, INTENT(INOUT), OPTIONAL :: ABS_COUNT(NUM_V)   !< Vector containing number of absorbers (to be incremented if positive transition energy)
    INTEGER, INTENT(INOUT), OPTIONAL :: SECOND_ABS(NUM_V) !< Vector that holds second absorber, if present

    INTEGER :: K, P

    REAL(KIND=HIGH_PREC_REAL) :: alpha(NUM_V)

    INTEGER :: first_res(NUM_V), &
               sec_res(NUM_V)

    alpha = 0.00D00 ; first_res = 0; sec_res = 0

    DO K = 1, NUM_V

      IF ((1.00D00 + TRANS_EN/V_GRID_P(K) ** 2) .GT. 0.00D00) &
      alpha(K) = SQRT(1.00D00 + TRANS_EN/V_GRID_P(K) ** 2)           !Calculate alpha factor of sweep

    END DO

!Perform sweep
    DO K = 1, NUM_V

      DO P = 1, NUM_V

        IF ((alpha(K) * V_GRID_P(K) .LT. V_CELL_BOUNDARY_P(P)) .AND. &
        (alpha(K) * V_GRID_P(K) .GE. V_CELL_BOUNDARY_P(P-1))) THEN

          first_res(K) = P

          EXIT

        END IF

      END DO

!Remove any mapping outside grid
      IF (alpha(K) * V_GRID_P(K) - V_GRID_P(NUM_V) .GT. 0) THEN

        first_res(K) = 0

      END IF

!Perform special case handling and second absorber mapping if applicable

      IF ((TRANS_EN .GT. 0.00D00) .AND.(PRESENT(ABS_COUNT))) THEN               !Direct sweep

          IF (first_res(K) .GT. 0) THEN

            IF (V_GRID_P(first_res(K)) ** 2 .LT. TRANS_EN)  first_res(K) = first_res(K) + 1 !Shift ideal emitter one cell in order for each cell to have a legal emitter

          END IF

          IF (first_res(K) .GT. 0) ABS_COUNT(first_res(K)) = ABS_COUNT(first_res(K)) + 1 !Increase absorber count for found ideal emitter of cell K

      END IF

      IF ((TRANS_EN .LT. 0.00D00) .AND. (PRESENT(SECOND_ABS)) .AND. (first_res(K) .GT. 0)) THEN              !Inelastic inverse sweep

        IF (alpha(K) * V_GRID_P(K) - V_GRID_P(first_res(K)) .GT. 0) THEN

          sec_res(K) = first_res(K) + 1

        ELSE

          sec_res(K) = first_res(K) - 1

        END IF

        !Forbid absorbers within one transition energy of grid edge
        IF (V_GRID_P(first_res(K)) ** 2 .GT. V_GRID_P(NUM_V) ** 2 + TRANS_EN) first_res(K) = 0

        IF (sec_res(K) .GT. 0) THEN

          IF (V_GRID_P(sec_res(K)) ** 2 .GT. V_GRID_P(NUM_V) ** 2 + TRANS_EN) sec_res(K) = 0

        END IF

      END IF

      IF ((TRANS_EN .GT. 0.00D00) .AND. (PRESENT(SECOND_ABS))) THEN              !Superelastic inverse sweep

        IF (first_res(K) .GT. 0) THEN

          IF (alpha(K) * V_GRID_P(K) - V_GRID_P(first_res(K)) .GT. 0) THEN

            sec_res(K) = first_res(K) + 1

          ELSE

            sec_res(K) = first_res(K) - 1

          END IF
          
        !(1.10.2019. Removed in order to conserve energy and maintain detailed balance)
        !  IF (V_GRID_P(first_res(K)) ** 2 .LT.  TRANS_EN) first_res(K) = 0
        !
        !  IF (sec_res(K) .GT. 0) THEN
        !
        !    IF (V_GRID_P(sec_res(K)) ** 2 .LT.  TRANS_EN) sec_res(K) = 0
        !
        !  END IF

        END IF

      END IF

    END DO

    SWEEP_RES = first_res
    IF (PRESENT(SECOND_ABS)) SECOND_ABS = sec_res

  END SUBROUTINE INEL_SWEEP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Perform final count of absorbers
  SUBROUTINE INEL_ABS_COUNT(K,ABS_COUNT, FIRST_ABS,SECOND_ABS,IDEAL_EM)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: K, FIRST_ABS(NUM_V), SECOND_ABS (NUM_V)

    INTEGER, INTENT(IN), OPTIONAL :: IDEAL_EM(NUM_V)

    INTEGER, INTENT(INOUT) :: ABS_COUNT(NUM_V)


    IF (PRESENT(IDEAL_EM)) THEN

      IF (FIRST_ABS(K) .GT. 0) THEN

        IF (IDEAL_EM(FIRST_ABS(K)) .NE. K) ABS_COUNT(K) = ABS_COUNT(K) + 1

      END IF

      IF (SECOND_ABS(K) .GT. 0) THEN

        IF (IDEAL_EM(SECOND_ABS(K)) .NE. K) ABS_COUNT(K) = ABS_COUNT(K) + 1

      END IF

      IF ((FIRST_ABS(K) .EQ. 1) .AND. (SECOND_ABS(K) .EQ. 0)) ABS_COUNT(K) = 0                 !If a cell emits below least resolved energy, ignore it

    ELSE

      IF (FIRST_ABS(K) .GT. 0) ABS_COUNT(K) = ABS_COUNT(K) + 1

      IF (SECOND_ABS(K) .GT. 0) ABS_COUNT(K) = ABS_COUNT(K) + 1

    END IF

    IF ((FIRST_ABS(K) .EQ. 0) .OR. (SECOND_ABS(K) .EQ.0)) ABS_COUNT(K) = 0

  END SUBROUTINE INEL_ABS_COUNT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fill absorber list for given transition and cell
  SUBROUTINE INEL_FILL_ABS(I,J,K,ABS_COUNT,FIRST_ABS,SECOND_ABS,IDEAL_EM)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: I, J, K, FIRST_ABS(NUM_V), SECOND_ABS (NUM_V), ABS_COUNT(NUM_V)

    INTEGER, INTENT(IN), OPTIONAL :: IDEAL_EM(NUM_V)

    INTEGER :: P, L

    INEL_ABS(K,I,J)%NA = ABS_COUNT(K)

    INEL_ABS(K,I,J)%EMIT_SW = 1

    ALLOCATE(INEL_ABS(K,I,J)%A(ABS_COUNT(K)))

    L = 1

!Fill absorber list if inelastic collision
    IF (PRESENT(IDEAL_EM)) THEN

      DO P = 1, NUM_V

        IF ((IDEAL_EM(P) .EQ. K) .OR. (P .EQ. FIRST_ABS(K)) .OR. (P .EQ. SECOND_ABS(K))) THEN

          INEL_ABS(K,I,J)%A(L) = P

          L = L + 1

        END IF

      END DO

    ELSE
!Fill absorber list if superelastic collision
      DO P = 1, NUM_V

        IF ((P .EQ. FIRST_ABS(K)) .OR. (P .EQ. SECOND_ABS(K))) THEN

          INEL_ABS(K,I,J)%A(L) = P

          L = L + 1

        END IF

      END DO

    END IF

  END SUBROUTINE INEL_FILL_ABS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fills pairs and calculates weights for given transition
  SUBROUTINE INEL_PARTITION_AND_WEIGHTS(I,J,K,ABS_COUNT,FIRST_ABS,SECOND_ABS,TRANS_EN)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: I, J, K, FIRST_ABS(NUM_V), SECOND_ABS (NUM_V), ABS_COUNT(NUM_V)

    REAL(KIND=HIGH_PREC_REAL), INTENT(IN) :: TRANS_EN

    INTEGER :: L, P, dl , dg

    REAL(KIND=HIGH_PREC_REAL), ALLOCATABLE :: b1(:), b2(:),gp(:)

    REAL(KIND=HIGH_PREC_REAL) :: b_norm

    IF (TRANS_EN .LT. 0.00D00) THEN

      INEL_ABS(K,I,J)%P = ABS_COUNT(K) - 1

      ALLOCATE(INEL_ABS(K,I,J)%m1(INEL_ABS(K,I,J)%P))
      ALLOCATE(INEL_ABS(K,I,J)%m2(INEL_ABS(K,I,J)%P))
      ALLOCATE(gp(INEL_ABS(K,I,J)%P))

      DO P = 1, INEL_ABS(K,I,J)%P

        INEL_ABS(K,I,J)%m1(P) = INEL_ABS(K,I,J)%A(P)
        INEL_ABS(K,I,J)%m2(P) = INEL_ABS(K,I,J)%A(P + 1)

      END DO

      gp = 1.00

    ELSE

      INEL_ABS(K,I,J)%P = ABS_COUNT(K) - 1

      ALLOCATE(INEL_ABS(K,I,J)%m1(INEL_ABS(K,I,J)%P))
      ALLOCATE(INEL_ABS(K,I,J)%m2(INEL_ABS(K,I,J)%P))

      DO P = 1, INEL_ABS(K,I,J)%P

        INEL_ABS(K,I,J)%m1(P) = INEL_ABS(K,I,J)%A(P)
        INEL_ABS(K,I,J)%m2(P) = INEL_ABS(K,I,J)%A(P + 1)

      END DO

!Allocate auxillary data

      ALLOCATE(gp(INEL_ABS(K,I,J)%P))
      ALLOCATE(b1(INEL_ABS(K,I,J)%P))
      ALLOCATE(b2(INEL_ABS(K,I,J)%P))

      dl = 0

      dg = 0

      L = 1

      DO P = 1, INEL_ABS(K,I,J)%NA

        IF (INEL_ABS(K,I,J)%A(P) .LT. min(FIRST_ABS(K),SECOND_ABS(K))) THEN

          INEL_ABS(K,I,J)%m1(L) = INEL_ABS(K,I,J)%A(P)
          INEL_ABS(K,I,J)%m2(L) = max(FIRST_ABS(K),SECOND_ABS(K))

          b1(L) = 2.00 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m1(L)) * V_GRID_P(INEL_ABS(K,I,J)%m1(L))
          b2(L) = 2.00 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m2(L)) * V_GRID_P(INEL_ABS(K,I,J)%m2(L))

          dl = dl + 1

          L = L + 1

        END IF

        IF (INEL_ABS(K,I,J)%A(P) .EQ. min(FIRST_ABS(K),SECOND_ABS(K))) THEN

          INEL_ABS(K,I,J)%m1(L) = min(FIRST_ABS(K),SECOND_ABS(K))
          INEL_ABS(K,I,J)%m2(L) = max(FIRST_ABS(K),SECOND_ABS(K))

          b1(L) = 2.00 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m1(L)) * V_GRID_P(INEL_ABS(K,I,J)%m1(L))
          b2(L) = 2.00 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m2(L)) * V_GRID_P(INEL_ABS(K,I,J)%m2(L))

          L = L + 1

        END IF

        IF (INEL_ABS(K,I,J)%A(P) .GT. max(FIRST_ABS(K),SECOND_ABS(K))) THEN

          INEL_ABS(K,I,J)%m1(L) = min(FIRST_ABS(K),SECOND_ABS(K))
          INEL_ABS(K,I,J)%m2(L) = INEL_ABS(K,I,J)%A(P)

          b1(L) = 2.00 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m1(L)) * V_GRID_P(INEL_ABS(K,I,J)%m1(L))
          b2(L) = 2.00 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m2(L)) * V_GRID_P(INEL_ABS(K,I,J)%m2(L))

          dg = dg + 1

          L = L + 1

        END IF

      END DO

!Calculate auxillary data

      b_norm = 0

      DO P = 1, ABS_COUNT(K)

        b_norm = b_norm + 2.00 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%A(P)) * V_GRID_P(INEL_ABS(K,I,J)%A(P))

      END DO

      b1 = b1 / b_norm
      b2 = b2 / b_norm

      DO P = 1, INEL_ABS(K,I,J)%P

        IF (INEL_ABS(K,I,J)%m1(P) .LT. min(FIRST_ABS(K),SECOND_ABS(K))) THEN

           gp(P) = b1(P) + b2(P) /REAL(dl + 1,KIND=HIGH_PREC_REAL)

         ELSE IF ((INEL_ABS(K,I,J)%m1(P) .EQ. min(FIRST_ABS(K),SECOND_ABS(K))) .AND. &
           (INEL_ABS(K,I,J)%m2(P) .EQ. max(FIRST_ABS(K),SECOND_ABS(K)))) THEN

           gp(P) = b1(P)/REAL(dg + 1,KIND=HIGH_PREC_REAL) + b2(P) /REAL(dl + 1,KIND=HIGH_PREC_REAL)

         ELSE IF (INEL_ABS(K,I,J)%m2(P) .GT. max(FIRST_ABS(K),SECOND_ABS(K))) THEN

           gp(P) = b1(P)/REAL(dg + 1,KIND=HIGH_PREC_REAL) + b2(P)

         END IF

      END DO

    END IF
!Fill weights
    ALLOCATE(INEL_ABS(K,I,J)%W1(INEL_ABS(K,I,J)%P))
    ALLOCATE(INEL_ABS(K,I,J)%W2(INEL_ABS(K,I,J)%P))

    DO P = 1, INEL_ABS(K,I,J)%P


      INEL_ABS(K,I,J)%W1(P) = gp(P) * ((V_GRID_P(K) ** 2 * V_GRID_WIDTH_P(K))&
                              / (V_GRID_P(INEL_ABS(K,I,J)%m1(P)) ** 2 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m1(P)))) * &
                              (1.00 - (V_GRID_P(K) ** 2 - V_GRID_P(INEL_ABS(K,I,J)%m1(P)) ** 2 - TRANS_EN) &
                              / (V_GRID_P(INEL_ABS(K,I,J)%m2(P)) ** 2 - V_GRID_P(INEL_ABS(K,I,J)%m1(P)) ** 2))

      INEL_ABS(K,I,J)%W2(P) = gp(P) * ((V_GRID_P(K) ** 2 * V_GRID_WIDTH_P(K))&
                              /(V_GRID_P(INEL_ABS(K,I,J)%m2(P)) ** 2 * V_GRID_WIDTH_P(INEL_ABS(K,I,J)%m2(P)))) * &
                              (V_GRID_P(K) ** 2 - V_GRID_P(INEL_ABS(K,I,J)%m1(P)) ** 2 -  TRANS_EN ) &
                              / (V_GRID_P(INEL_ABS(K,I,J)%m2(P)) ** 2 - V_GRID_P(INEL_ABS(K,I,J)%m1(P)) ** 2)

    END DO

    DEALLOCATE(gp)
    IF (ALLOCATED(b1)) DEALLOCATE(b1)
    IF (ALLOCATED(b2)) DEALLOCATE(b2)

  END SUBROUTINE INEL_PARTITION_AND_WEIGHTS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Fill data on emitters for I -> J transition
  SUBROUTINE INEL_FILL_EM(I,J)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: I, J

    INTEGER :: NE, K, P, M, N, L

    DO K = 1, NUM_V

      NE = 0

    DO P = 1, NUM_V

      IF (ALLOCATED(INEL_ABS(P,I,J)%A)) THEN

          DO L = 1, INEL_ABS(P,I,J)%NA

            IF (INEL_ABS(P,I,J)%A(L) .EQ. K) THEN

              NE = NE + 1

              EXIT

            END IF

          END DO

        END IF

      END DO

      INEL_EM(K,I,J)%NE = NE

      IF (NE .GT. 0) THEN

        ALLOCATE(INEL_EM(K,I,J)%E(NE))
        ALLOCATE(INEL_EM(K,I,J)%P(NE))
        ALLOCATE(INEL_EM(K,I,J)%PAIR_LABEL(NE))
        ALLOCATE(INEL_EM(K,I,J)%ABS_LABEL(NE))

        L = 1

        INEL_EM(K,I,J)%P = 0

        DO N = 1, NUM_V

          IF (ALLOCATED(INEL_ABS(N,I,J)%A)) THEN

            DO P = 1, INEL_ABS(N,I,J)%NA

              IF (K .EQ. INEL_ABS(N,I,J)%A(P)) THEN

                INEL_EM(K,I,J)%E(L) = N

                DO M = 1, INEL_ABS(N,I,J)%P

                  IF ((INEL_ABS(N,I,J)%m1(M) .EQ. K) .OR. (INEL_ABS(N,I,J)%m2(M) .EQ. K)) &
                    INEL_EM(K,I,J)%P(L) = INEL_EM(K,I,J)%P(L) + 1

                END DO

                L = L + 1

                EXIT

              END IF

            END DO

          END IF

        END DO

        DO P = 1, NE

          ALLOCATE(INEL_EM(K,I,J)%PAIR_LABEL(P)%L(INEL_EM(K,I,J)%P(P)))
          ALLOCATE(INEL_EM(K,I,J)%ABS_LABEL(P)%L(INEL_EM(K,I,J)%P(P)))

          M = 1

          DO N = 1, INEL_ABS(INEL_EM(K,I,J)%E(P),I,J)%P

            IF (INEL_ABS(INEL_EM(K,I,J)%E(P),I,J)%m1(N) .EQ. K) THEN

              INEL_EM(K,I,J)%PAIR_LABEL(P)%L(M) = N
              INEL_EM(K,I,J)%ABS_LABEL(P)%L(M) = 1

              M = M + 1

            ELSE IF (INEL_ABS(INEL_EM(K,I,J)%E(P),I,J)%m2(N) .EQ. K) THEN

              INEL_EM(K,I,J)%PAIR_LABEL(P)%L(M) = N
              INEL_EM(K,I,J)%ABS_LABEL(P)%L(M) = 2

              M = M + 1

            END IF

          END DO

        END DO

      END IF

    END DO

  END SUBROUTINE INEL_FILL_EM
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculate number of nonzeroes for given transition
  SUBROUTINE INEL_NNZ_COUNT(I,J)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: I, J

    INTEGER :: K, P, N, NNZ

    NNZ = 0

    DO K = 1, NUM_V

      IF (ALLOCATED(INEL_ABS(K,I,J)%A))  NNZ = NNZ + 1

      DO P = 1, NUM_V

        IF (ALLOCATED(INEL_ABS(P,I,J)%A)) THEN

          DO N = 1, INEL_ABS(P,I,J)%NA

            IF ((K .EQ. INEL_ABS(P,I,J)%A(N)) .AND. (P .NE. K)) NNZ = NNZ + 1

          END DO

        END IF

      END DO

    END DO

    INEL_NNZ(I,J) = NNZ

  END SUBROUTINE INEL_NNZ_COUNT
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE INEL_GRID
!-------------------------------------------------------------------------------------------------------------------------------------
