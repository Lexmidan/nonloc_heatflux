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
!> Contains calls to build the collision integrals and the collisional-radiative model submatrices
MODULE BUILD_0D

  USE GRID
  USE SWITCHES
  USE BUILD_COULOMB_COLL
  USE BUILD_NEUTRAL_COLL
  USE BUILD_CRM
  USE MPI
  USE VAR_KIND_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calls all collision integral and CRM submatrix builder and adds the resulting matrices
  SUBROUTINE FILL_0D(f_lagged, n_old, T_old, n_lagged, T_lagged, n_i_lagged, u_i_lagged,N_NONLIN)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(DIM_F) !< Lagged vector
    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_X), INTENT(IN) :: n_old, & !< Old density
                                             T_old, & !< Old temperature
                                             n_lagged, & !< Lagged density
                                             T_lagged !< Lagged temperature
    INTEGER, INTENT(IN) :: N_NONLIN !< Current nonlinear iteration

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: n_i_lagged(NUM_X), &  !< Lagged ion density
                                  u_i_lagged(NUM_X)  !< Lagged ion velocity

    INTEGER :: I, J, K ,&
               H_OFFSET !< Offset when calling CRM builders

    REAL(KIND=DEFAULT_REAL) :: F(NUM_V) !< f_0^0 chunk for given spatial cell

    IF (COLL_EE_0_SWITCH) CALL ALLOCATE_WEIGHTS                                 !Allocate Chang-Cooper weight arrays if not allocated

    DO I = MIN_X, MAX_X                                                         !Main spatial loop

      IF (.NOT. FULL_FLUID_MODE) THEN

        DO J = 1, NUM_H                                                           !Add all matrices that evolve the distribution function (save for recombination)

          IF (((MOD(J - 1,2) .EQ. 0) .AND. (MOD(I,2) .EQ. 1)) .OR. ((MOD(J - 1,2) .EQ. 1) .AND. (MOD(I,2) .EQ. 0))) THEN        !Evolve only even harmonics on odd spatial cells (centres) and vice-versa

            F = f_lagged(X_POS(I) + NUM_V * (NUM_H - 1) + 1: X_POS(I) + NUM_V * (NUM_H - 1) + NUM_V)  !Extract f_0^0 from total vector

            IF (J - 1 .EQ. 0) THEN                                             !Operators acting on f_0^0

              IF (COLL_EE_0_SWITCH) THEN

                CALL FILL_EE_00 (F, n_old(I), T_old(I), N_NONLIN, I)     !Call e-e collision operator matrix builder

              END IF

              IF (COLL_EN_EL_0_SWITCH) THEN

                DO K = 1, NUM_NEUTRALS

                  CALL FILL_EN_EL_00(f_lagged(X_POS(I) + NUM_H * NUM_V + K),I, K)

                END DO

              END IF

            ELSE

              IF (COLL_EI_L_SWITCH) THEN

                IF (COLD_ION_FLUID_SWITCH) THEN

                  CALL FILL_CI_EI_L(n_i_lagged(I), n_old(I), T_old(I), J - 1,I,J,u_i_lagged(I),F)   !Call e-i collision operator matrix builder

                ELSE

                  CALL FILL_EI_L(n_lagged(I)/Z_PROF(I), n_old(I), T_old(I), J - 1,I,J)   !Call e-i collision operator matrix builder

                END IF

              END IF

              IF (COLL_EE_L_SWITCH) THEN

                CALL FILL_EE_L(n_old(I), T_old(I), F, J - 1,I,J)             !Call e-e collision operator matrix builder for higher harmonics

              END IF

              IF (COLL_EN_EL_L_SWITCH) THEN

                DO K = 1, NUM_NEUTRALS

                  CALL FILL_EN_EL_L(f_lagged(X_POS(I) + NUM_H * NUM_V + K), J - 1,I,J,K)

                END DO

              END IF

            END IF

            IF (NEUTRAL_TRACK_SWITCH) THEN                                        !Inelastic collision operators

              IF (COLL_EN_EX) THEN

                CALL FILL_EN_EX(f_lagged(X_POS(I) + NUM_H * NUM_V + 1 :&     !Call excitation collision operator matrix builder
                                X_POS(I) + NUM_H * NUM_V + NUM_NEUTRALS), J - 1,I,J)


              END IF

              IF (COLL_EN_ION) THEN

              CALL FILL_EN_ION(f_lagged(X_POS(I) + NUM_H * NUM_V + 1 :&      !Call ionization collision operator matrix builder
                              X_POS(I) + NUM_H * NUM_V + NUM_NEUTRALS), J - 1,I,J)

              END IF

              IF (COLL_RECOMB) THEN

                IF (COLD_ION_FLUID_SWITCH) THEN

                  CALL FILL_EN_RECOMB(n_lagged(I),n_i_lagged(I), J - 1,I,J)

                ELSE

                  CALL FILL_EN_RECOMB(n_lagged(I),n_lagged(I), J - 1,I,J)

                END IF

              END IF

            END IF

          END IF

        END DO

      END IF

      IF (NEUTRAL_TRACK_SWITCH) THEN                                            !Add all submatrices that evolve the neutral densities due to CR model (plus the effect of recombination on the distribution function)

          H_OFFSET = X_POS(I) + NUM_V * NUM_H                                   !Neutral subvector position

          IF (MOD(I,2) .EQ. 1) THEN                                             !Evolve neutral densities only on odd spatial cells (centres)

            IF (COLL_EN_EX) THEN

                CALL FILL_CRM_EX(I)                           !Call excitation and de-excitation CR model matrix builder

            END IF

            IF (COLL_EN_ION) THEN

                CALL FILL_CRM_ION(I)                                      !Call ionization CR model matrix builder

            END IF

            IF (COLL_RECOMB .AND. (.NOT. FULL_FLUID_MODE)) THEN

              IF (COLD_ION_FLUID_SWITCH) THEN

                CALL FILL_CRM_RECOMB(T_lagged(I), n_i_lagged(I),I) !Call recombination CR model matrix builder

              ELSE

                CALL FILL_CRM_RECOMB(T_lagged(I), n_lagged(I),I) !Call recombination CR model matrix builder

              END IF

            END IF

          END IF

        END IF

    END DO


  END SUBROUTINE FILL_0D
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_0D
!-------------------------------------------------------------------------------------------------------------------------------------
