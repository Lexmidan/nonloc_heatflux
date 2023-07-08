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
!>Contains variable kind data
MODULE VAR_KIND_DATA

  USE MPI

  INTEGER, PARAMETER :: REAL_DOUBLE = 8, & !<Double precision kind
                        REAL_QUAD = 16     !< Quad precision kind

  INTEGER, PARAMETER :: DEFAULT_REAL = REAL_DOUBLE, &  !< Default real kind - used for most variables
                        HIGH_PREC_REAL = REAL_DOUBLE, &  !< Higher precision real kind, to be used for operators prone to round-off error accumulation
                        PETSC_DEF_REAL = REAL_DOUBLE   !< Default kind used for variables passed to PETSc

!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE VAR_KIND_DATA
!-------------------------------------------------------------------------------------------------------------------------------------
