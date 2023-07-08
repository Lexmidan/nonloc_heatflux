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
!> Contains all printing subroutines
MODULE PRINTING

  USE MPI
  USE VAR_KIND_DATA
  
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Prints starting logo
  SUBROUTINE PRINT_START

    IMPLICIT NONE

    PRINT*,  '-----------------------------------------------------------------------------------------------------------------'
    PRINT*,  '_____   _____   _                _   __  _   _____ ' // NEW_LINE('A') // &
             '/  ___| |  _  | | |              | | / / (_) |_   _|' // NEW_LINE('A') // &
             '\ `--.  | | | | | |      ______  | |/ /   _    | |  ' // NEW_LINE('A') // &
             ' `--. \ | | | | | |     |______| |    \  | |   | |  ' // NEW_LINE('A') // &
             '/\__/ / \ \_/ / | |____          | |\  \ | |   | |  ' // NEW_LINE('A') // &
             '\____/   \___/  \_____/          \_| \_/ |_|   \_/  '
    PRINT*
    PRINT*, '-----------------------------------------------------------------------------------------------------------------'
    PRINT*,'Copyright 2020 Stefan Mijin'
    PRINT*
    PRINT*,'SOL-KiT is free software: you can redistribute it and/or modify'
    PRINT*,'     it under the terms of the GNU General Public License as published by'
    PRINT*,'     the Free Software Foundation, either version 3 of the License, or'
    PRINT*,'     (at your option) any later version.'
    PRINT*
    PRINT*,'     SOL-KiT is distributed in the hope that it will be useful,'
    PRINT*,'     but WITHOUT ANY WARRANTY; without even the implied warranty of'
    PRINT*,'     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    PRINT*,'     GNU General Public License for more details.'
    PRINT*
    PRINT*,'     You should have received a copy of the GNU General Public License'
    PRINT*,'     along with SOL-KiT.  If not, see <https://www.gnu.org/licenses/>.'
    PRINT*, '-----------------------------------------------------------------------------------------------------------------'
    PRINT*
  END SUBROUTINE PRINT_START
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Prints initialization declaration
  SUBROUTINE PRINT_INITIALIZING(char)

    IMPLICIT NONE

    CHARACTER (LEN = *), INTENT(IN) :: char

    PRINT*, 'Initializing ' // char // ' ...'

  END SUBROUTINE PRINT_INITIALIZING
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Prints timestep status
  SUBROUTINE PRINT_TIMESTEP(TIMESTEP, T, step)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: TIMESTEP
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: T
    CHARACTER (LEN = *), INTENT(IN) :: step

    PRINT*, 'Finished ' // step // ' ', TIMESTEP
    WRITE(*,'(1X,A, ES14.7, A)') 'Elapsed real time: ', T ,' seconds'

  END SUBROUTINE PRINT_TIMESTEP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Generic print
  SUBROUTINE PRINT_ECHO(char)

    IMPLICIT NONE

    CHARACTER (LEN = *), INTENT(IN) :: char

    PRINT*, char

  END SUBROUTINE PRINT_ECHO
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Prints evolve start declaration
  SUBROUTINE PRINT_START_EVOLVE(TIMESTEP, char)

    IMPLICIT NONE

    CHARACTER (LEN = *), INTENT(IN) :: char
    INTEGER, INTENT(IN) :: TIMESTEP

    PRINT*, 'Starting ' // char, ' ', TIMESTEP

  END SUBROUTINE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Prints status after nonlinear step
  SUBROUTINE PRINT_END_NONLIN(N_NONLIN, DELTA_F)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N_NONLIN
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: DELTA_F

    WRITE(*,'(1X,A, I4, 10X, A, ES14.7)') 'Finished nonlinear iteration ', N_NONLIN, 'Nonlinear residual = ', DELTA_F

  END SUBROUTINE PRINT_END_NONLIN
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Prints evolve end declaration and status
  SUBROUTINE PRINT_END_EVOLVE(N_NONLIN, DELTA_F)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N_NONLIN
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: DELTA_F

    WRITE(*,'(1X,A, I4, A, ES14.7)') 'Converged at nonlinear iteration', N_NONLIN, ' and with nonlinear residual = ', DELTA_F

  END SUBROUTINE PRINT_END_EVOLVE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE PRINTING
