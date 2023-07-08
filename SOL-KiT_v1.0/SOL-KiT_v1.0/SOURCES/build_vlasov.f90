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
!> Contains Vlasov submatrix builder
MODULE BUILD_VLASOV

  USE GRID
  USE SWITCHES
  USE SPARSE
  USE BUILD_E_ADV
  USE MPI
  USE VAR_KIND_DATA
  USE BUILD_X_ADV
  USE BUILD_MAXWELL
  USE MATRIX_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calls Vlasov term submatrix builders
  SUBROUTINE FILL_VLASOV(f_lagged)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_lagged(DIM_F) !< Lagged variable vectors

    IF (E_ADV_ON) THEN

      CALL FILL_E_ADV(f_lagged)                                                 !Call E-field advection submatrix builders

    END IF

  END SUBROUTINE FILL_VLASOV
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE BUILD_VLASOV
!-------------------------------------------------------------------------------------------------------------------------------------
