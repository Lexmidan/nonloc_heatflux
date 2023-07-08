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
!> Contains all grid quantities and functions
MODULE GRID

  USE INPUT
  USE SWITCHES
  USE NEUT_AND_HEAT
  USE NORMALIZATION
  USE MPI
  USE VAR_KIND_DATA

PRIVATE :: INIT_V_GRID, INIT_X_GRID, INIT_DIMS

! Time

  INTEGER :: TIMESTEP_NUM, & !< Number of full length timesteps
             PRETIMESTEP_NUM, & !< Number of small pre-timesteps
             T_SAVE !< Save interval (save every T_SAVE timesteps)

  REAL(KIND=DEFAULT_REAL) :: dt, & !< Timestep size
            pre_dt !< Pre-timestep size

! Harmonics

  INTEGER :: L_MAX, & !< L-number of maximum resolved harmonic
             NUM_H !< Total number of tracked harmonics

! Velocity

  INTEGER :: NUM_V !< Number of cells in velocity space

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: V_GRID(:), & !< Vector containing velocity cell centres
                               V_GRID_WIDTH(:), & !< Vector containing velocity grid widths
                               dvp(:), & !< Vector containing forward velocity grid widths (used in differentiaion)
                               dvm(:), & !< Vector containing backward velocity grid widths (used in differentiaion)
                               V_CELL_BOUNDARY(:), & !< Vector containing values of velocity on velocity cell boundaries
                               v_interp(:) !< Interpolation factor vector on velocity grid

! High precision velocity grid for inelastic collisions
  REAL(KIND=HIGH_PREC_REAL), ALLOCATABLE :: V_GRID_P(:), & !< Vector containing velocity cell centres (high-precision)
                               V_GRID_WIDTH_P(:), & !< Vector containing velocity grid widths (high-precision)
                               V_CELL_BOUNDARY_P(:) !< Vector containing values of velocity on velocity cell boundaries (high-precision)

  REAL(KIND=DEFAULT_REAL) :: dv, & !< Size of first velocity cell
                  V_GRID_MULT !< Velocity cell width common ratio

! x-grid

  INTEGER :: NUM_C, & !< Number of spatial cells
             NUM_X !< Length of X_GRID vector, including all cell boundaries

  REAL(KIND=DEFAULT_REAL) :: dx, & !< Size of spatial cells (or largest cell if logarithmic grid)
            SMALL_dx, & !< Size of divertor cell for logarithmic grid
            TOTAL_L !< Total simulation length in meters

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: X_GRID(:), &!< Vector containing cell centres and boundaries
                               dxc(:), & !< Vector containing grid widths (used in differentiaion)
                               dxp(:), & !< Vector containing forward grid widths (used in differentiaion)
                               dxm(:) !< Vector containing backward grid widths (used in differentiaion)

! Total and fields

  INTEGER :: DIM_F, & !< Total vector length
             NUM_FIELDS, & !< Number of fields tracked
             NUM_0D !< Number of quantities tracked per spatial cell

! Offsets

  INTEGER :: OFFSET_UP, & !< Upstream offset if fixed boundary
             OFFSET_DIV !< Divertor offset if fixed boundary

! Local processor dimensions

  INTEGER :: nd, & !< Number of spatial points per processor (rounded up)
             nd_loc, & !< Local number of spatial cells
             LOC_ROWS, & !< Local number of rows
             MIN_X, & !< First evolved spatial cell - local
             MAX_X, & !< Last evolved spatial cell - local
             LOC_NUM_C !< Local number of spatial cell centres

! Ionization profile

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: Z_PROF(:)
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes all grid quantities
  SUBROUTINE INIT_GRID

    IMPLICIT NONE

    INTEGER :: RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    IF (RANK .EQ. 0) CALL INPUT_GRID(TIMESTEP_NUM, &
                                     PRETIMESTEP_NUM, &
                                     dt, &
                                     pre_dt, &
                                     T_SAVE, &
                                     L_MAX, &
                                     NUM_V, &
                                     dv, &
                                     V_GRID_MULT, &
                                     NUM_C, &
                                     dx, &
                                     SMALL_dx &
                                     )

    CALL MPI_Bcast(TIMESTEP_NUM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(PRETIMESTEP_NUM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(pre_dt,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(T_SAVE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(L_MAX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NUM_V,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(dv,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(V_GRID_MULT,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NUM_C,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(dx,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(SMALL_dx,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    CALL INIT_V_GRID
    CALL INIT_X_GRID
    CALL INIT_DIMS

!Initialize ionization profile
    ALLOCATE(Z_PROF(NUM_X))

    IF (Z_PROFILE_FROM_FILE) THEN

      IF (RANK .EQ. 0) CALL INPUT_INIT("Z_PROFILE",Z_PROF,NUM_X)
      CALL MPI_Bcast(Z_PROF,NUM_X,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    ELSE

      Z_PROF = ION_Z

    END IF

  END SUBROUTINE INIT_GRID
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Returns position of L harmonic (from bottom up) in single spatial cell vector chunk
  INTEGER FUNCTION H_POS(L)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: L !< L-number of traced harmonic

    INTEGER :: I, &
               H

    H = 0

    DO I = L + 1, L_MAX

      H = H + 1

    END DO

    H_POS = H

  END FUNCTION H_POS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes velocity grid
  SUBROUTINE INIT_V_GRID

    IMPLICIT NONE

    INTEGER :: I

    ALLOCATE(V_GRID(NUM_V))
    ALLOCATE(V_GRID_WIDTH(NUM_V))
    ALLOCATE(dvp(NUM_V))
    ALLOCATE(dvm(NUM_V))
    ALLOCATE(V_CELL_BOUNDARY(0:NUM_V))
    ALLOCATE(v_interp(NUM_V))

    V_GRID(1) = dv/2.00D00                                                            !First cell centre
    V_GRID_WIDTH(1) = dv                                                         !First velocity cell width
    V_CELL_BOUNDARY(0) = 0.0D00
    V_CELL_BOUNDARY(1) = dv

    DO I = 2, NUM_V

      V_GRID_WIDTH(I) = V_GRID_MULT * V_GRID_WIDTH(I-1)

      V_CELL_BOUNDARY(I) = V_CELL_BOUNDARY(I-1) + V_GRID_WIDTH(I)

      V_GRID(I) = V_GRID(I - 1) + 0.50D00*(V_GRID_WIDTH(I-1)+V_GRID_WIDTH(I))

    END DO

!Initialize grid spacing data
    DO I = 2, NUM_V

      dvm(I) = V_GRID(I) - V_GRID(I-1)

    END DO

    dvm(1) = dv/2.00D00

    DO I = 1, NUM_V - 1

      dvp(I) = V_GRID(I+1) - V_GRID(I)

    END DO

    dvp(NUM_V) = 0.50D00 * V_GRID_WIDTH(NUM_V) *(1.00D00 + V_GRID_MULT)

    DO I = 1, NUM_V

      v_interp(I) = 1.00D00 - 0.50D00*V_GRID_WIDTH(I)/dvp(I)

    END DO

!High precision grid
    ALLOCATE(V_GRID_P(NUM_V))
    ALLOCATE(V_GRID_WIDTH_P(NUM_V))
    ALLOCATE(V_CELL_BOUNDARY_P(0:NUM_V))

    V_GRID_P(1) = dv/2.00D00                                                            !First cell centre
    V_GRID_WIDTH_P(1) = dv                                                         !First velocity cell width
    V_CELL_BOUNDARY_P(0) = 0.0D00
    V_CELL_BOUNDARY_P(1) = dv

    DO I = 2, NUM_V

      V_GRID_WIDTH_P(I) = V_GRID_MULT * V_GRID_WIDTH_P(I-1)

      V_CELL_BOUNDARY_P(I) = V_CELL_BOUNDARY_P(I-1) + V_GRID_WIDTH_P(I)

      V_GRID_P(I) = V_GRID_P(I - 1) + 0.50D00*(V_GRID_WIDTH_P(I-1)+V_GRID_WIDTH_P(I))

    END DO

  END SUBROUTINE INIT_V_GRID
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes spatial grid
  SUBROUTINE INIT_X_GRID

    IMPLICIT NONE

    INTEGER :: I, &
               NUM_GHOST_CELLS

    REAL(KIND=DEFAULT_REAL) :: GRID_B

    INTEGER :: RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    NUM_GHOST_CELLS = 0

    IF (PERIODIC_BOUNDARY_SWITCH) THEN                                          !Equal numbers of cell centres and boundaries

      NUM_X = 2 * NUM_C

    ELSE

      IF (FIXED_BOUNDARY_UP_SWITCH) NUM_GHOST_CELLS = NUM_GHOST_CELLS + 1       !Extra cell centres for fixed quantities
      IF (FIXED_BOUNDARY_DIV_SWITCH) NUM_GHOST_CELLS = NUM_GHOST_CELLS + 1

      NUM_X = 2 * (NUM_C + NUM_GHOST_CELLS) - 1

    END IF

    ALLOCATE(X_GRID(NUM_X))

    IF (X_GRID_FROM_FILE) THEN

      IF (RANK .EQ. 0) CALL INPUT_INIT("X_GRID",X_GRID,NUM_X)
      CALL MPI_Bcast(X_GRID,NUM_X,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

    ELSE
!Handle 0D case
      IF (NUM_X .GT. 1) THEN

        X_GRID(1) = 0

      ELSE

        X_GRID(1) = 1

      END IF

      IF (LOGARITHMIC_GRID_SWITCH) THEN

        GRID_B = LOG(dx/SMALL_dx)/(NUM_X-1)

        DO I = 2, NUM_X

          X_GRID(I) = X_GRID(I - 1) + dx * EXP(-(I - 1) * GRID_B) / 2

        END DO

      ELSE

        DO I = 2, NUM_X

          X_GRID(I) = X_GRID(I - 1) + dx/2.0

        END DO

      END IF

    END IF

!Initialize grid spacing data
    ALLOCATE(dxc(NUM_X))
    ALLOCATE(dxp(NUM_X))
    ALLOCATE(dxm(NUM_X))

    DO I = 2, NUM_X - 1

      dxc(I) = X_GRID(I + 1) - X_GRID(I - 1)

    END DO

    IF (NUM_X .GT. 1) THEN

      dxc(1) = 2 * (X_GRID(2) - X_GRID(1))
      dxc(NUM_X) = 2 * (X_GRID(NUM_X) - X_GRID(NUM_X - 1))

      DO I = 2, NUM_X

        dxm(I) = dxc(I - 1)

      END DO

      dxm(1) = dx

      DO I = 1, NUM_X - 1

        dxp(I) = dxc(I + 1)

      END DO

      dxp(NUM_X) = dx

    END IF

!Make sure periodic boundaries are consistent with input grid
    IF (PERIODIC_BOUNDARY_SWITCH .AND. X_GRID_FROM_FILE) THEN

      dxc(1) = MAX(dxc(1),dxc(NUM_X))
      dxc(NUM_X) = dxc(1)

    END IF

  END SUBROUTINE INIT_X_GRID
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes vector dimensions, offsets, and local dimensions
  SUBROUTINE INIT_DIMS

    IMPLICIT NONE

    INTEGER :: I, RANK, SIZE, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    NUM_FIELDS = 0

    IF (MAXWELL_SWITCH .OR. E_ADV_SWITCH) NUM_FIELDS = 1                        !Counts E_x field entry

    NUM_H = L_MAX + 1

    IF (FULL_FLUID_MODE .AND. NUM_H .EQ. 1) NUM_H = 2                              !Make sure f1 is present if electrons are fluid

!Initialize 0D and total dimensions
    NUM_0D = NUM_H*NUM_V + NUM_NEUTRALS + NUM_FIELDS

    IF (COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) NUM_0D = NUM_0D + 2                              !Count ion density and flow velocity

    IF (FULL_FLUID_MODE) NUM_0D = NUM_0D + 3

    DIM_F = NUM_X * NUM_0D

    OFFSET_UP = 0
    OFFSET_DIV = 0

    IF (FIXED_BOUNDARY_UP_SWITCH) OFFSET_UP = 1
    IF (FIXED_BOUNDARY_DIV_SWITCH) OFFSET_DIV = 1

!Distribute spatial cells to processors and calculate local data
    IF (MOD(NUM_X,SIZE) .GT. 0) THEN

      nd = NUM_X / SIZE + 1

      IF (nd * (SIZE - 1) .GE. NUM_X) THEN

        nd = NUM_X / SIZE

      END IF

    ELSE

      nd = NUM_X / SIZE

    END IF

    IF (RANK .EQ. SIZE - 1) THEN

      nd_loc = NUM_X - nd * (SIZE - 1)

    ELSE

      nd_loc = nd

    END IF

    LOC_ROWS = NUM_0D * nd_loc

    MIN_X = MAX(RANK * nd + 1, 1 + OFFSET_UP)
    MAX_X = MIN (RANK * nd + nd_loc, NUM_X - OFFSET_DIV)

    LOC_NUM_C = 0

    DO I = MIN_X,MAX_X

      IF (MOD(I,2) .EQ. 1) LOC_NUM_C = LOC_NUM_C + 1

    END DO
!Calculate total length of domain
    TOTAL_L = X_GRID(NUM_X) * TIME_NORM * THERM_VEL

  END SUBROUTINE INIT_DIMS
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Returns starting position of k-th spatial cell in total vector
  INTEGER FUNCTION X_POS(k)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: k

    X_POS = (k - 1) * NUM_0D

  END FUNCTION X_POS
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE GRID
!-------------------------------------------------------------------------------------------------------------------------------------
