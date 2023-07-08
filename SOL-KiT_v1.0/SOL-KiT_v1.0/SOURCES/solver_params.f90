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
!>Contains all parameters used in matrix solver and nonlinear iterations
MODULE SOLVER_PARAMS

  USE INPUT
  USE MPI
  USE VAR_KIND_DATA

  IMPLICIT NONE

!KSP Parameters

  REAL(KIND=DEFAULT_REAL) :: SOLVER_TOL !< Solver tolerance
  INTEGER :: SOLVER_MAX_ITER !< Solver max iteration

!Nonlinear iteration parameters

  INTEGER :: MAX_NONLIN !< Maximum number of iterations
  REAL(KIND=DEFAULT_REAL) :: NONLIN_TOL !< Nonlinear tolerance

!Cut-off parameters

  REAL(KIND=DEFAULT_REAL) :: BIS_TOL !< Bisection tolerence for calculating divertor logical boundary condition cut-off velocity

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!>Calls input routine to
  SUBROUTINE INIT_SOLVER_PARAMS

    IMPLICIT NONE

    INTEGER :: RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    IF (RANK .EQ. 0) CALL INPUT_SOLVER_PARAMS(SOLVER_TOL, &
                                              SOLVER_MAX_ITER, &
                                              MAX_NONLIN, &
                                              NONLIN_TOL, &
                                              BIS_TOL)

    CALL MPI_Bcast(SOLVER_TOL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(SOLVER_MAX_ITER,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NONLIN_TOL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(MAX_NONLIN,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(BIS_TOL,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

  END SUBROUTINE INIT_SOLVER_PARAMS
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE SOLVER_PARAMS
!-------------------------------------------------------------------------------------------------------------------------------------
