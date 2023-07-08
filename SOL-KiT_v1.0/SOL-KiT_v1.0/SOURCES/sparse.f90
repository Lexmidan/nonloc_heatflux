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
!> Contains definition of SPARSE_MAT type, along with sparse matrix operations and the identity matrix
MODULE SPARSE

USE MPI
USE VAR_KIND_DATA
!> Sparse matrix type in coordinate list representation
  TYPE SPARSE_MAT

    INTEGER :: N_ROWS, &    !< Number of rows in matrix
               N_COLUMNS, & !< Number of columns in matrix
               N_NONZ       !< Number of nonzero entries in matrix
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ROW, & !< Vector containing row numbers of nonzero elements
                                          COLUMN !< Vector containing column numbers of nonzero elements
    REAL(KIND=PETSC_DEF_REAL), ALLOCATABLE, DIMENSION(:) :: VALUE !< Vector containing values of nonzero elements

  END TYPE SPARSE_MAT

!> Sparse matrix pattern for precalculation
  TYPE SPARSITY_PATTERN

    INTEGER :: N_NZ !< Number of nonzeros in pattern
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ROW, & !< Vector containing row numbers of nonzero elements in pattern
                                          COL !< Vector containing column numbers of nonzero elements in pattern

  END TYPE SPARSITY_PATTERN

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Allocates SPARSE_MAT type vector members and sets integer properties
  SUBROUTINE ALLOCATE_SPARSE(M,N_R,N_C,N_NZ)

    IMPLICIT NONE

    TYPE (SPARSE_MAT), INTENT(OUT) :: M !< Matrix to be allocated
    INTEGER, INTENT(IN) :: N_R, & !< Number of rows to be set to M%N_ROWS
                           N_C, & !< Number of columns to be set to M%N_COLUMNS
                           N_NZ   !< Number of nonzero elements to be set to M%N_N_NONZ and as length to M%ROW, M%COLUMN, and M%VALUE

      M%N_ROWS = N_R
      M%N_COLUMNS = N_C
      M%N_NONZ = N_NZ

      ALLOCATE(M%ROW(N_NZ))
      ALLOCATE(M%COLUMN(N_NZ))
      ALLOCATE(M%VALUE(N_NZ))

  END SUBROUTINE ALLOCATE_SPARSE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> If sparse matrix member vectors are allocated, deallocate them
  SUBROUTINE DEALLOCATE_SPARSE(M)

    IMPLICIT NONE

    TYPE(SPARSE_MAT), INTENT(OUT) :: M

    IF (ALLOCATED(M%ROW)) THEN

      DEALLOCATE(M%ROW)
      DEALLOCATE(M%COLUMN)
      DEALLOCATE(M%VALUE)

    END IF

  END SUBROUTINE DEALLOCATE_SPARSE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Sets M1 = M2.
!! Deallocates all M1 vectors, and if M2 is allocated sets all M1 members to those of M2.
  SUBROUTINE SPARSE_EQ(M1,M2)

    IMPLICIT NONE

    TYPE (SPARSE_MAT), INTENT(INOUT) :: M1 !< LHS sparse matrix
    TYPE (SPARSE_MAT), INTENT(IN) :: M2 !< RHS sparse matrix

    INTEGER :: I


    CALL DEALLOCATE_SPARSE(M1)

    IF (ALLOCATED(M2%ROW)) THEN

      CALL ALLOCATE_SPARSE(M1, M2%N_ROWS,M2%N_COLUMNS,M2%N_NONZ)

      DO I = 1, M2%N_NONZ

        M1%ROW(I) = M2%ROW(I)
        M1%COLUMN(I) = M2%COLUMN(I)
        M1%VALUE(I) = M2%VALUE(I)

      END DO

    END IF

  END SUBROUTINE SPARSE_EQ
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates number of overlaping entries in sparsity patterns M1 and M2, as well as the overlap indices
  SUBROUTINE OVERLAP_SP(M1,M2,OVERLAP_NUM,OVERLAP_INDEX1,OVERLAP_INDEX2,O_START)

    IMPLICIT NONE

    TYPE (SPARSITY_PATTERN), INTENT(IN) :: M1, & !< First matrix
                                     M2    !< Second matrix
    INTEGER, INTENT(OUT) :: OVERLAP_NUM    !< Number of overlaping elements
    INTEGER, INTENT(IN) :: O_START !< First potential overlap element in M1
    INTEGER, ALLOCATABLE, INTENT(OUT) :: OVERLAP_INDEX1(:), &  !< Overlap index for M1
                                         OVERLAP_INDEX2(:)
                                                            !!! If the I-th nonzero element of M1 overlaps the J-th nonzero element of
                                                            !!! M2 then OVERLAP_INDEX(I) = J, otherwise it equals 0
    INTEGER I, J

    OVERLAP_NUM = 0

    ALLOCATE(OVERLAP_INDEX1(M1%N_NZ))
    ALLOCATE(OVERLAP_INDEX2(M2%N_NZ))

    OVERLAP_INDEX1 = 0
    OVERLAP_INDEX2 = 0

    DO I = O_START + 1, M1%N_NZ
      DO J = 1, M2%N_NZ

        IF ((M1%ROW(I) .EQ. M2%ROW(J)) .AND. (M1%COL(I) .EQ. M2%COL(J))) THEN

          OVERLAP_NUM = OVERLAP_NUM + 1
          OVERLAP_INDEX1(I) = J
          OVERLAP_INDEX2(J) = I

          EXIT
        END IF

      END DO
    END DO

  END SUBROUTINE OVERLAP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Adds sparsity pattern M_S to sparsity pattern M
  SUBROUTINE ADD_SP(M,M_S,OFFSET_R, OFFSET_C,OVERLAP_STATUS,OVERLAP_START)

    TYPE (SPARSITY_PATTERN), INTENT(INOUT) :: M!< Sparsity pattern to be added to
    TYPE (SPARSITY_PATTERN), INTENT(IN) :: M_S !< Sparsity pattern to be added
    LOGICAL, OPTIONAL, INTENT(IN) :: OVERLAP_STATUS !< .TRUE. if overlap possible, .FALSE. if not
    INTEGER, OPTIONAL, INTENT(IN) :: OVERLAP_START !< Start of overlap in M
    TYPE (SPARSITY_PATTERN) :: M_TEMP1,&  !< Temporary sparsity pattern to store result of SPARSE_ADD
                               M_TEMP2
    INTEGER, INTENT(IN) :: OFFSET_R, & !< Row offset
                           OFFSET_C !< Column offset
    LOGICAL :: O_STATUS

    INTEGER :: O_START

    IF (PRESENT(OVERLAP_STATUS)) THEN

      O_STATUS = OVERLAP_STATUS

    ELSE

      O_STATUS = .TRUE.

    END IF

    O_START = 0
    IF (PRESENT(OVERLAP_START)) O_START = OVERLAP_START

    M_TEMP1%N_NZ = M_S%N_NZ

    ALLOCATE(M_TEMP1%ROW(M_S%N_NZ))
    ALLOCATE(M_TEMP1%COL(M_S%N_NZ))

    M_TEMP1%ROW = M_S%ROW + OFFSET_R
    M_TEMP1%COL = M_S%COL + OFFSET_C

    CALL SPARSE_ADD_SP(M,M_TEMP1,M_TEMP2,O_STATUS,O_START)
    CALL SPARSE_EQ_SP(M,M_TEMP2)

    IF (ALLOCATED(M_TEMP1%ROW)) THEN

      DEALLOCATE(M_TEMP1%ROW)
      DEALLOCATE(M_TEMP1%COL)

    END IF

    IF (ALLOCATED(M_TEMP2%ROW)) THEN

      DEALLOCATE(M_TEMP2%ROW)
      DEALLOCATE(M_TEMP2%COL)

    END IF

  END SUBROUTINE ADD_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Adds two sparsity patterns and stores result in M
  SUBROUTINE SPARSE_ADD_SP(M1,M2,M, OVERLAP_STATUS, O_START)

    IMPLICIT NONE

    TYPE (SPARSITY_PATTERN), INTENT(IN) :: M1, & !< First sparsity pattern to be added
                                           M2 !< Second sparsity pattern to be added
    TYPE (SPARSITY_PATTERN), INTENT(OUT) :: M !< Output sparsity matrix

    INTEGER, INTENT(IN) :: O_START !< First potential overlap element in M1
    LOGICAL, INTENT(IN) :: OVERLAP_STATUS !< .TRUE. if overlap possible, .FALSE. if not
    INTEGER :: I, J, &
               OVERLAP_N !< Number of overlaping nonzero element positions between M1 and M2
    INTEGER, ALLOCATABLE, DIMENSION(:) :: OVERLAP_I1, & !< Overlap indices for M1
                                          OVERLAP_I2    !< Overlap indices for M2

    IF (ALLOCATED(M1%ROW) .AND. ALLOCATED(M2%ROW)) THEN                         !If both allocated

      OVERLAP_N = 0
      IF (OVERLAP_STATUS) THEN
        CALL OVERLAP_SP(M1,M2,OVERLAP_N,OVERLAP_I1,OVERLAP_I2,O_START)                                  !Calculate overlap indices if overlap present
      END IF

      IF (ALLOCATED(M%ROW)) THEN

        DEALLOCATE(M%ROW)
        DEALLOCATE(M%COL)

      END IF

      M%N_NZ = M1%N_NZ + M2%N_NZ - OVERLAP_N

      ALLOCATE(M%ROW(M%N_NZ))
      ALLOCATE(M%COL(M%N_NZ))

      IF (OVERLAP_STATUS) THEN

        DO I = 1, M1%N_NZ                                                       !Add elements of M1, including overlaping elements from M2

          M%ROW(I) = M1%ROW(I)
          M%COL(I) = M1%COL(I)

        END DO

        J = M1%N_NZ

        DO I = 1, M2%N_NZ                                                       !Add non-overlaping elements from M2

          IF (OVERLAP_I2(I) .EQ. 0) THEN

            J = J + 1

            M%ROW(J) = M2%ROW(I)
            M%COL(J) = M2%COL(I)

          END IF

        END DO

      ELSE

        DO I = 1, M1%N_NZ                                                       !Add elements of M1, including overlaping elements from M2

          M%ROW(I) = M1%ROW(I)
          M%COL(I) = M1%COL(I)

        END DO

        DO I = 1, M2%N_NZ                                                       !Add non-overlaping elements from M2

            M%ROW(M1%N_NZ + I) = M2%ROW(I)
            M%COL(M1%N_NZ + I) = M2%COL(I)

        END DO

      END IF

    ELSE IF (ALLOCATED(M1%ROW)) THEN                                            !If one allocated

      CALL SPARSE_EQ_SP(M,M1)

    ELSE IF (ALLOCATED(M2%ROW)) THEN

      CALL SPARSE_EQ_SP(M,M2)

    END IF

  END SUBROUTINE SPARSE_ADD_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Sets M1 = M2.
!! Deallocates all M1 vectors, and if M2 is allocated sets all M1 members to those of M2.
  SUBROUTINE SPARSE_EQ_SP(M1,M2)

    IMPLICIT NONE

    TYPE (SPARSITY_PATTERN), INTENT(INOUT) :: M1 !< LHS sparsity pattern
    TYPE (SPARSITY_PATTERN), INTENT(IN) :: M2 !< RHS sparsity pattern

    INTEGER :: I

    IF (ALLOCATED(M1%ROW)) THEN

      DEALLOCATE(M1%ROW)
      DEALLOCATE(M1%COL)

    END IF

    IF (ALLOCATED(M2%ROW)) THEN

      M1%N_NZ = M2%N_NZ

      ALLOCATE(M1%ROW(M2%N_NZ))
      ALLOCATE(M1%COL(M2%N_NZ))

      DO I = 1, M2%N_NZ

        M1%ROW(I) = M2%ROW(I)
        M1%COL(I) = M2%COL(I)

      END DO

    END IF

  END SUBROUTINE SPARSE_EQ_SP
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE SPARSE
!-------------------------------------------------------------------------------------------------------------------------------------
