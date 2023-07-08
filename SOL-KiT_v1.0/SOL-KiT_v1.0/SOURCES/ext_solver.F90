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
!> Manages calling PETSc
MODULE EXT_SOLVER

#include "petsc/finclude/petscksp.h"

  USE petscksp
  USE SPARSE
  USE SOLVER_PARAMS
  USE GRID
  USE MATRIX_DATA
  USE MPI
  USE VAR_KIND_DATA

!>Preallocated value buffer based on dimensions of markers for each row
  TYPE VAL_BUFFER

    REAL(KIND=PETSC_DEF_REAL), ALLOCATABLE, DIMENSION(:) :: val

  END TYPE VAL_BUFFER

  TYPE (VAL_BUFFER), ALLOCATABLE, DIMENSION(:) :: vb !< Matrix value buffer

  LOGICAL :: FIRST_PASS !<True if first time through solve (in which case objects are built and stored as module variables)

  INTEGER, ALLOCATABLE :: loc_lengths(:), displs(:) !Gatherv arguments

  Vec rhs, sol
  Mat PETSc_mat
  PC pc
  KSP solver

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Solver the matrix equation M*f_new = f_old
  SUBROUTINE SOLVE_M(M, f_old, f_new,TIMESTEP)

    IMPLICIT NONE

    TYPE(SPARSE_MAT), INTENT(INOUT) :: M

    REAL(KIND=DEFAULT_REAL), INTENT(IN), DIMENSION(DIM_F) :: f_old !< Known vector
    INTEGER, INTENT(IN) :: TIMESTEP
    REAL(KIND=DEFAULT_REAL), INTENT(INOUT), DIMENSION(DIM_F) :: f_new !< Unknown vector

    REAL(KIND=DEFAULT_REAL), DIMENSION(LOC_ROWS) :: loc_x

    INTEGER :: RANK,SIZE, IERR, I, J, IR

    KSPConvergedReason conv_r
    PetscReal norm
    PetscInt its
    PetscScalar, pointer :: xx_v(:)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

!If first time through perform all allocation and PETSc object creation
    IF (FIRST_PASS) THEN

      ALLOCATE(vb(LOC_ROWS))

      DO I = 1, LOC_ROWS

        ALLOCATE(vb(I)%val(RD(I + RANK * nd * NUM_0D)%NNZ))

      END DO

      CALL MatCreateAIJ(MPI_COMM_WORLD,LOC_ROWS,LOC_ROWS,DIM_F,DIM_F,1,D_NNZ,1,OD_NNZ,PETSc_mat,IERR)
      CHKERRQ(IERR)
      CALL VecCreateMPI(MPI_COMM_WORLD,LOC_ROWS,DIM_F,sol,IERR)
      CHKERRQ(IERR)
      CALL VecDuplicate(sol,rhs,IERR)
      CHKERRQ(IERR)

      CALL KSPCreate(MPI_COMM_WORLD,solver,IERR)
      CHKERRQ(IERR)
      CALL KSPSetType(solver,"bcgs",IERR)
      CHKERRQ(IERR)

      CALL KSPGetPC(solver,pc,IERR)
      CHKERRQ(IERR)

      !CALL PCSetType(pc,PCCHOLESKY,IERR)
      CALL PCSetType(pc,PCBJACOBI,IERR)

      CALL PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO,IERR)


      ALLOCATE(loc_lengths(0:SIZE-1))

      DO I = 0, SIZE - 1

        IF (I .EQ. SIZE - 1) THEN

          loc_lengths(I) = (NUM_X - nd * (SIZE - 1))*NUM_0D

        ELSE

          loc_lengths(I) = nd * NUM_0D

        END IF

      END DO

      ALLOCATE(displs(0:SIZE-1))

      DO I = 0, SIZE - 1

        displs(I) = nd * NUM_0D * I

      END DO

      FIRST_PASS = .FALSE.

    END IF

!Fill and assemble matrix
    DO I = 1, LOC_ROWS

      DO J = 1, RD(I + RANK * nd * NUM_0D)%NNZ

        vb(I)%val(J) = M%VALUE(RD(I + RANK * nd * NUM_0D)%MARKER(J))

      END DO

      CALL MatSetValues(PETSc_mat,1, I + RANK * nd * NUM_0D - 1, RD(I + RANK * nd * NUM_0D)%NNZ,&
                        RD(I + RANK * nd * NUM_0D)%COL - 1,&
                        vb(I)%val,INSERT_VALUES,IERR); CHKERRQ(IERR)

    END DO

    CALL MatAssemblyBegin(PETSc_mat,MAT_FINAL_ASSEMBLY,IERR)
    CHKERRQ(IERR)
    CALL MatAssemblyEnd(PETSc_mat,MAT_FINAL_ASSEMBLY,IERR)
    CHKERRQ(IERR)

!Fill and assemble vectors
    DO I = 1, LOC_ROWS

      CALL VecSetValues(rhs,1,I + RANK * nd * NUM_0D - 1,&
      REAL(f_old(I + RANK * nd * NUM_0D),KIND=PETSC_DEF_REAL),INSERT_VALUES,IERR)
      CHKERRQ(IERR)

    END DO

    CALL VecAssemblyBegin(rhs,IERR)
    CHKERRQ(IERR)
    CALL VecAssemblyEnd(rhs,IERR)
    CHKERRQ(IERR)

    DO I = 1, LOC_ROWS

      CALL VecSetValues(sol,1,I + RANK * nd * NUM_0D - 1,&
      REAL(f_new(I + RANK * nd * NUM_0D),KIND=PETSC_DEF_REAL),INSERT_VALUES,IERR)
      CHKERRQ(IERR)

    END DO

    CALL VecAssemblyBegin(sol,IERR)
    CHKERRQ(IERR)
    CALL VecAssemblyEnd(sol,IERR)
    CHKERRQ(IERR)

!Determine rhs norm and use in setting solver tolerance
    CALL VecNorm(rhs, NORM_2,norm,IERR)
    CHKERRQ(IERR)

    CALL KSPSetTolerances(solver,REAL(1.00D-20,KIND=PETSC_DEF_REAL),&
    SOLVER_TOL*norm,REAL(1.00D07,KIND=PETSC_DEF_REAL),SOLVER_MAX_ITER,IERR)
    CHKERRQ(IERR)

!Set solver operators and solve
    call KSPSetOperators(solver,PETSc_mat,PETSc_mat,IERR)
    CHKERRQ(IERR)

    call KSPSolve(solver,rhs,sol,IERR)
    CHKERRQ(IERR)

!Get converged reason and number of iterations

    CALL KSPGetConvergedReason(solver,conv_r,IERR)
    CHKERRQ(IERR)
    IF (RANK .EQ. 0) PRINT*, 'KSP converged reason: ',conv_r

    CALL KSPGetIterationNumber(solver,its,IERR)
    CHKERRQ(IERR)

    IF (conv_r .LT. 0) THEN

        IF (RANK .EQ. 0) THEN
          PRINT*, 'KSP diverged - aborting'

          OPEN(11, FILE = 'OUTPUT/STATUS.txt', ACTION = 'WRITE', IOSTAT = IR, STATUS = 'REPLACE')

          WRITE(11,'(A,I5)') 'KSP diverged at timestep ',TIMESTEP
          WRITE(11,'(A,I5)') 'KSP diverged with',conv_r

          CLOSE(11,IOSTAT=IR)
        END IF

        STOP

    END IF

    PRINT*, 'Process rank ',RANK,' converged after ',its,' iterations.'


!Fetch solution values into local vector and gather into new f_new

    CALL VecGetArrayReadF90(sol,xx_v,IERR)
    CHKERRQ(IERR)

    loc_x = REAL(xx_v,KIND=DEFAULT_REAL)

    CALL VecRestoreArrayReadF90(sol,xx_v,IERR)
    CHKERRQ(IERR)

    CALL MPI_Gatherv(loc_x,LOC_ROWS,MPI_REAL8,f_new,loc_lengths,displs,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

  END SUBROUTINE SOLVE_M
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Sets first pass to true
  SUBROUTINE SET_FIRST_PASS

    IMPLICIT NONE

    FIRST_PASS = .TRUE.

  END SUBROUTINE SET_FIRST_PASS
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE EXT_SOLVER
!-------------------------------------------------------------------------------------------------------------------------------------
