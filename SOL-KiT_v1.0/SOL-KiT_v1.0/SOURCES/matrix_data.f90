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
!> Contains the global matrix and sparsity patterns
MODULE MATRIX_DATA

  USE GRID
  USE SPARSE
  USE NEUT_AND_HEAT
  USE BUILD_X_ADV
  USE BUILD_MAXWELL
  USE INEL_GRID
  USE SWITCHES
  USE MPI
  USE VAR_KIND_DATA

  PRIVATE :: INIT_TRI_DIAG_SP, INIT_DENSE_SP, INIT_DIAG_SP, INIT_ID_MAT_SP, INIT_DENSE_NEUT_SP, INIT_DIAG_NEUT_SP, &
             INIT_EN_ION_SP, INIT_INEL_SP, INIT_INEL_MARKER, INIT_CRM_RECOMB_MARKER, &
             INIT_E_ADV_SP, INIT_SEC_EL_MARKER, INIT_CRM_RECOMB_SP

!> Marker for elements affected by inelastic collision operators
  TYPE INEL_MARKER

    INTEGER, ALLOCATABLE :: LIST(:) !< List of positions of nonzero elements in global matrix

  END TYPE INEL_MARKER

!> Contains nonzero data for each row in local matrix band used to create PETSc matrices
  TYPE ROW_DATA

    INTEGER :: NNZ, & !< Number of nonzeros
               D_NNZ, & !< Number of nonzeros in diagonal submatrix
               OD_NNZ !< Number of nonzeros in off-diagonal submatrices

    INTEGER, ALLOCATABLE, DIMENSION(:) :: COL, & !< Vector containing column values for nonzeros in given row
                                          MARKER !< Vector containing nonzero index of nonzeros in given row (connecting to LOCAL_M)

  END TYPE ROW_DATA

!> Contains data about particle source matrix elements
  TYPE PART_SOURCE_DATA

    INTEGER :: POS, & !< Spatial position of element
               J, & !< Velocity cell index - row
               K !< Velocity cell index - column

  END TYPE PART_SOURCE_DATA


!> Contains data about inelastic collision matrix elements
  TYPE INEL_DAT

    INTEGER, ALLOCATABLE :: E(:) !< 1 if de-excitation matrix element corresponds to emmision

    REAL(KIND=HIGH_PREC_REAL), ALLOCATABLE :: W(:,:) !< Pair-summed inelastic weights


  END TYPE INEL_DAT

!> Marker for properties of nonzero elements of E field advection
  TYPE E_ADV_MARKER

    INTEGER, ALLOCATABLE, DIMENSION(:) :: H, & !< Harmonics
                                          V, & !< Velocity cell
                                          P !< Spatial cell
    REAL(KIND=DEFAULT_REAL), ALLOCATABLE, DIMENSION(:) ::  M !< Multipliers for interpolation

  END TYPE E_ADV_MARKER

!> Marker for properties of nonzero elements of ion fluid convection term
  TYPE CONV_PROP

    INTEGER, ALLOCATABLE, DIMENSION(:) :: SIGN, & !< Term sign
                                          POS !< Term position

  END TYPE CONV_PROP


  TYPE(SPARSE_MAT) :: LOCAL_M !< Global matrix

  TYPE(SPARSE_MAT) :: FIXED_MAT !< Stores fixed elements

  TYPE (SPARSITY_PATTERN) :: LOCAL_SP, & !< Global matrix sparsity pattern
                             TRI_DIAG_SP, & !< Tri-diagonal sparsity pattern in velocity space
                             DENSE_SP, &  !< Dense sparsity pattern in velocity space
                             DIAG_SP, &  !< Diagonal sparsity pattern in velocity space
                             DENSE_NEUT_SP, & !< Dense sparsity pattern for neutrals
                             DIAG_NEUT_SP, & !< Diagonal sparsity pattern for neutrals
                             EN_ION_SP, & !< Cold seconary electron formation sparsity pattern
                             CRM_RECOMB_SP, & !< Recombination sparsity pattern
                             E_ADV_SP,& !< E-field velocity space advection sparsity pattern
                             NEUTRAL_DIFF_SP, & !< Neutral diffusion sparsity pattern
                             SINK_SP, & !< Divertor sink sparsity pattern
                             REC_SP, & !< Recycling sparsity pattern
                             ID_MAT_SP, & !< Diagonal sparsity pattern for total matrix
                             MAXWELL_ION_SP, & !< Sparsity pattern for ion velocity contribution to Maxwell-Ampere law
                             ION_LORENTZ_SP, & !< Sparsity pattern for Lorentz term in ion fluid momentum equation
                             ION_CONV_SP, & !< Sparsity pattern for ion convection term
                             ION_CONT_SP, & !< Sparsity pattern for ion continuity term
                             ION_GAIN_SP, & !< Sparsity pattern for ion gain due to ionization
                             ION_LOSS_SP, & !< Sparsity pattern for ion loss due to recombination
                             ION_DRAG_SP, & !< Sparsity pattern for e-i drag in ion momentum equation
                             CI_REC_SP, & !< Sparsity pattern for recycling with cold fluid ions
                             ION_CONV_DIV_SP, & !< Sparsity pattern for the ion convection term when divertor boundary condition on
                             ION_SINK_SP, & !< Sparsity pattern for the ion sink at divertor
                             CI_EI_OD_SP, & !< Sparsity pattern for off-diagonal contribution to cold ion e-i collision operator
                             ION_MOM_SOURCE_SP,& !< Sparsity pattern for the effect of continuity source term on the ion momentum equation
                             PART_SOURCE_SP, & !< Sparsity pattern for electron particle source
                             ION_PART_SOURCE_SP, & !< Sparsity pattern for ion particle source
                             ION_PRESSURE_SP, & !< Sparsity pattern for ion pressure gradient (isothermal)
                             SIMPLE_CX_SP, & !< Sparsity pattern for simple CX model i-n friction
                             MAXWELL_EL_SP, & !< Sparsity pattern for electron velocity contribution to Maxwell-Ampere law
                             EL_LORENTZ_SP, & !< Sparsity pattern for Lorentz term in electron fluid momentum equation
                             EL_CONV_SP, & !< Sparsity pattern for electron convection term
                             EL_CONT_SP, & !< Sparsity pattern for electron continuity term
                             EL_GAIN_SP, & !< Sparsity pattern for electron gain due to ionization
                             ION_LOSS_FF_SP, & !< Sparsity pattern for ion loss due to recombination with full fluid mode
                             EL_LOSS_SP, & !< Sparsity pattern for electron loss due to recombination
                             EL_CONV_DIV_SP, & !< Sparsity pattern for the electron convection term when divertor boundary condition on
                             EL_SINK_SP, & !< Sparsity pattern for the ion sink at divertor
                             EL_MOM_SOURCE_SP,& !< Sparsity pattern for the effect of continuity source term on the electron momentum equation
                             EL_PRESSURE_SP, & !< Sparsity pattern for electron pressure gradient
                             EL_T_DRAG_SP, & !< Sparsity pattern for electron thermal drag force
                             ION_T_DRAG_SP, & !< Sparsity pattern for ion thermal drag force
                             EL_U_DRAG_SP, & !< Sparsity pattern for electron e-i drag force
                             ION_U_DRAG_SP, & !< Sparsity pattern for ion e-i drag force
                             EL_CONV_TEMP_SP, & !< Sparsity pattern for electron temperature convection term
                             EL_TEMP_PDV_SP, & !< Sparsity pattern for electron temperature pdv term
                             EL_TEMP_PDV_DIV_SP, & !< Sparsity pattern for electron divertor temperature pdv term
                             EL_TEMP_Q_U_SP, & !< Sparsity pattern for electron temperature div q_u term
                             EL_GAIN_TEMP_SP, & !< Sparsity pattern for temperature equation source term due to ionization
                             EL_LOSS_TEMP_SP, & !< Sparsity pattern for temperature equation source term due to recombination
                             EL_T_DRAG_TEMP_SP, & !< Sparsity pattern for electron thermal drag force temp eq term
                             EL_U_DRAG_TEMP_SP, & !< Sparsity pattern for e-i drag force temp eq term
                             EL_EN_DRAG_N_SP, & !< Sparsity pattern for e-n drag force (excluding recombination)
                             EL_EN_DRAG_RECOMB_SP, & !< Sparsity pattern for recombination drag
                             EL_EN_DRAG_TEMP_SP, & !< Sparsity pattern for e-n drag temperature terms
                             EL_EN_INCOLL_TEMP_SP, & !< Sparsity pattern for electron energy loss/gain term due to e-n inelastic collisions
                             EL_EN_ELCOLL_TEMP_SP, & !< Sparsity pattern for electron energy loss/gain term due to e-n inelastic collisions
                             CRM_RECOMB_FF_SP, & !< Recombination sparsity pattern for full fluid model
                             EL_HEATING_SP !< Heating pattern for full fluid model

  TYPE (SPARSITY_PATTERN), ALLOCATABLE :: INEL_SP(:,:) !< Sparsity pattern array for inelastic collisions

  TYPE (PART_SOURCE_DATA), ALLOCATABLE :: SOURCE_DATA(:), & !< Vectory holding electron source matrix element data
                                          ION_SOURCE_DATA(:) !< Vectory holding ion source matrix element data

  TYPE (INEL_MARKER), ALLOCATABLE :: INEL_M(:,:,:,:) !< Marker array for inelastic collisions
  TYPE (E_ADV_MARKER) :: PROP_E_ADV !< Properties for nonzero elements of E-advection operator

  INTEGER, ALLOCATABLE :: MARKER_EE_0(:), &  !< Markers for start of nonzeroes for e-e collision operator for l = 0
                          MARKER_EN_EL_0(:), & !< Markers for start of nonzeroes for e-n elastic collision operator for l = 0
                          MARKER_EE_L(:,:), & !< Markers for start of nonzeroes for e-e collision operators for l > 0
                          MARKER_EI_L(:,:), & !< Markers for start of nonzeroes for e-i collision operators for l > 0
                          MARKER_EN_EL_L(:,:), & !< Markers for start of nonzeroes for e-n elastic collision operators for l > 0
                          MARKER_CRM_EX(:), & !< Markers for start of nonzeroes for CR excitation model
                          MARKER_CRM_ION(:), & !< Markers for start of nonzeroes for CR excitation model
                          MARKER_RECOMB_CRM(:,:), & !< List of positions of all nonzero elements of recombination CR model in global matrix
                          MARKER_ID(:), & !< List of positions of all diagonal elements of global matrix
                          MARKER_SEC_EL(:,:), & !< List of positions of ionization cold electron production matrix elements
                          MARKER_HEATING(:), & !< Markers for the start of nonzeroes for the heating operator
                          MARKER_NEUT_DIFF(:), & !< Markers for neutral diffusion nonzeroes
                          MARKER_SINK(:), & !< Markers for divertor sink nonzeroes
                          MARKER_REC(:), & !< Markers for recycling nonzeroes
                          MARKER_MAXWELL_ION(:), & !< Marker for ion velocity contribution to Maxwell-Ampere law
                          MARKER_ION_LORENTZ(:), & !< Marker for Lorentz term in ion fluid momentum equation
                          MARKER_ION_CONV(:), & !< Marker for ion convection term
                          MARKER_ION_CONT(:), & !< Marker for ion continuity term
                          MARKER_ION_GAIN(:), & !< Marker for ion gain term due to ionization
                          MARKER_ION_LOSS(:), & !< Marker for ion loss term due to recombination
                          MARKER_CI_EI_L(:,:,:), & !< Marker for left submatrix in cold ion e-i collisions
                          MARKER_CI_EI_R(:,:,:), & !< Marker for right submatrix in cold ion e-i collisions
                          MARKER_ION_DRAG(:), & !< Marker for e-i drag in ion momentum equation
                          MARKER_CI_REC(:), & !< Marker for recycling with cold fluid ions
                          MARKER_ION_CONV_DIV(:), & !< Marker for the ion convection term when divertor boundary condition on
                          MARKER_ION_SINK(:), & !< Marker for ion sink at divertor
                          MARKER_ION_MOM_SOURCE(:), & !< Marker for source effect on ion momentum
                          MARKER_PART_SOURCE(:), & !< Marker for particle source for electrons
                          MARKER_ION_PART_SOURCE(:), & !< Marker for particle source for ions
                          MARKER_ION_PRESSURE(:), & !< Marker for ion pressure gradient (isothermal)
                          MARKER_SIMPLE_CX(:), & !< Marker for simple charge-exchange model
                          MARKER_MAXWELL_EL(:), & !< Marker for electron velocity contribution to Maxwell-Ampere law
                          MARKER_EL_LORENTZ(:), & !< Marker for Lorentz term in electron fluid momentum equation
                          MARKER_EL_CONV(:), & !< Marker for electron convection term
                          MARKER_EL_CONT(:), & !< Marker for electron continuity term
                          MARKER_EL_GAIN(:), & !< Marker for electron gain term due to ionization
                          MARKER_ION_LOSS_FF(:), & !< Marker for ion loss term due to recombination with full fluid mode
                          MARKER_EL_LOSS(:), & !< Marker for electron loss term due to recombination
                          MARKER_EL_CONV_DIV(:), & !< Marker for the ion convection term when divertor boundary condition on
                          MARKER_EL_SINK(:), & !< Marker for electron sink at divertor
                          MARKER_EL_MOM_SOURCE(:), & !< Marker for source effect on electron momentum
                          MARKER_EL_PRESSURE(:), & !< Marker for electron pressure gradient
                          MARKER_EL_T_DRAG(:), & !< Marker for electron thermal drag force
                          MARKER_ION_T_DRAG(:), & !< Marker for electron thermal drag force
                          MARKER_EL_U_DRAG(:), & !< Marker for electron e-i drag force
                          MARKER_ION_U_DRAG(:), & !< Marker for ion e-i drag force
                          MARKER_EL_CONV_TEMP(:), & !< Marker for electron convection term
                          MARKER_EL_TEMP_PDV(:), & !< Marker for electron pdv term
                          MARKER_EL_TEMP_PDV_DIV(:), & !< Marker for electron divertor pdv term
                          MARKER_EL_TEMP_Q_U(:), & !< Marker for electron temperature div q_u term
                          MARKER_EL_GAIN_TEMP(:), & !< Marker for temperature equation source term due to ionization
                          MARKER_EL_LOSS_TEMP(:), & !< Marker for temperature equation source term due to recombination
                          MARKER_EL_T_DRAG_TEMP(:), & !< Marker for electron thermal drag force temp eq term
                          MARKER_EL_U_DRAG_TEMP(:), & !< Marker for e-i drag force temp eq term
                          MARKER_EL_EN_DRAG(:), & !< Marker for e-n drag force (excluding recombination)
                          MARKER_EL_EN_DRAG_RECOMB(:), & !< Marker for recombination drag force
                          MARKER_EL_EN_DRAG_TEMP(:), & !< Marker for e-n drag temperature terms
                          MARKER_EL_EN_INCOLL_TEMP(:), & !< Marker for electron energy loss/gain term due to e-n inelastic collisions
                          MARKER_EL_EN_ELCOLL_TEMP(:), & !< Marker for electron energy loss/gain term due to e-n elastic collisions
                          MARKER_RECOMB_CRM_FF(:), & !< Marker for CRM recombination submatrix if full fluid model
                          MARKER_EL_HEATING(:) !< Marker for full fluid mode heating

  INTEGER :: MARKER_E_ADV !< Marker for start of E-advection terms

  INTEGER :: FL_INDEX, & !< Flush index
             FLUID_ION_INDEX !< Index marking start of fluid ion elements in matrix

!Local switches
  LOGICAL :: HEATING_ON, & !< True if heating pattern allocated
             E_ADV_ON, &  !< True if E-adv pattern allocated
             PLASMA_SINK_ON, & !< True if plasma sink pattern allocated
             REC_ON, & !< True if recycling pattern allocated
             NEUTRAL_DIFF_ON, & !< True if neutral diffusion pattern allocated
             ION_FLOW_ON, & !< True if ion flow operators allocated
             ION_CONT_ON, & !< True if ion continuity operators allocated
             EL_FLOW_ON, & !< True if electron flow operators allocated
             EL_CONT_ON, & !< True if electron continuity operators allocated
             PART_SOURCE_ON !< True if particle source pattern allocated

  TYPE (INEL_DAT), ALLOCATABLE :: IN_DAT(:,:)                                   !Inelastic process data

  TYPE (ROW_DATA), ALLOCATABLE :: RD(:)                                         !Local row data

  TYPE (CONV_PROP) :: ION_CONV_PROP, ION_CONT_PROP, EL_CONV_PROP, &
                      EL_CONT_PROP, EL_CONV_TEMP_PROP, EL_TEMP_PDV_PROP, &
                      EL_TEMP_Q_U_PROP   !Ion/electron convection/continuity properties

  INTEGER, ALLOCATABLE, DIMENSION(:) :: D_NNZ, OD_NNZ          !Number of nonzeros in (off-)diagonal section of given row - for use in PETSc

  INTEGER, ALLOCATABLE :: ION_GAIN_POS(:), &                      !Spatial position for each local ion gain term due to ionization
                          ION_GAIN_N(:), &                           !Neutral state being ionized for each ion gain term
                          ION_RECOMB_V(:), &                         !Velocity grid point for each ion loss term due to recombination
                          ION_RECOMB_POS(:), &                         !Spatial position for each local ion loss term due to recombination
                          SINK_SHIFT(:), &                               !Shift of sink boundary condition pre-reflected distribution sampling position
                          ION_PRESSURE_POS(:), &          !Position of column element for ion pressure gradient contribution
                          ION_PRESSURE_SIGN(:), &        !Sign of column element for ion pressure gradient contribution
                          EL_GAIN_POS(:), &                      !Spatial position for each local electron gain term due to ionization
                          EL_GAIN_N(:), &                           !Neutral state being ionized for each electron gain term
                          ION_RECOMB_FF_POS(:), &                         !Spatial position for each local ion loss term due to recombination with full fluid mode
                          EL_RECOMB_POS(:), &                         !Spatial position for each local electron loss term due to recombination
                          EL_PRESSURE_POS(:), &          !Position of column element for electron pressure gradient contribution
                          EL_PRESSURE_SIGN(:), &        !Sign of column element for electron pressure gradient contribution
                          EL_T_DRAG_SIGN(:), &        !Sign of column element for electron thermal drag force
                          ION_T_DRAG_SIGN(:), &        !Sign of column element for ion thermal drag force
                          EL_U_DRAG_SIGN(:), &        !Sign of column element for electron e-i drag force
                          ION_U_DRAG_SIGN(:), &        !Sign of column element for ion e-i drag force#
                          EL_GAIN_TEMP_POS(:), &                      !Spatial position for each local electron temperature source due to ionization
                          EL_GAIN_TEMP_N(:), &                           !Neutral state being ionized for each electron electron temperature source term
                          EL_RECOMB_TEMP_POS(:), &                         !Spatial position for each local electron temperature equation loss term due to recombination
                          EL_T_DRAG_TEMP_SIGN(:), &        !Sign of column element for electron thermal drag force temp eq term
                          EL_EN_DRAG_N(:), &                 ! Neutral state for each local e-n drag term
                          EL_EN_DRAG_N_POS(:), &             ! Position of each local e-n drag term (excluding recombination)
                          EL_EN_DRAG_REC_POS(:), &              ! Position of each local recombination drag term
                          NEUTRAL_DIFF_N(:), & ! State reference by each diffusion term
                          ION_LORENTZ_POS(:)  ! Position of each ion lorentz term

  REAL(KIND=DEFAULT_REAL), ALLOCATABLE :: EL_EN_DRAG_N_MULT(:), & ! Interpolation multiplier for each local e-n drag term (excluding recombination)
                                          EL_EN_DRAG_REC_MULT(:) !Interpolation multiplier for each local recombination drag term

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initializes all global matrix data and determines the flush index based on switches
  SUBROUTINE INIT_MATRIX_DATA

    IMPLICIT NONE

    INTEGER :: I, J, K, L, H_OFFSET, OVERLAP_START, OVERLAP_START_P

    LOGICAL :: BUILD_1D_MAT

    INTEGER :: RANK, SIZE,IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    FL_INDEX = 0
    FLUID_ION_INDEX = 0

!Initialize local switches
    HEATING_ON = .FALSE.
    E_ADV_ON = .FALSE.
    PLASMA_SINK_ON = .FALSE.
    REC_ON = .FALSE.
    ION_FLOW_ON = .FALSE.
    ION_CONT_ON = .FALSE.
    EL_FLOW_ON = .FALSE.
    EL_CONT_ON = .FALSE.
    PART_SOURCE_ON = .FALSE.
    NEUTRAL_DIFF_ON = .FALSE.

!Initialize auxillary matrices

!Make sure we're not building empty 1D local bands
    BUILD_1D_MAT = .TRUE.

    IF (nd .EQ. 1) THEN

      IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) BUILD_1D_MAT = .FALSE.
      IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) BUILD_1D_MAT = .FALSE.

    END IF

    IF (RANK .EQ. 0) PRINT*, ' Allocating matrices, sparsity patterns, and markers'

!Initialize x-advection and Ampere-Maxwell law matrices
    IF ((NUM_X .GT. 1) .AND. BUILD_1D_MAT) THEN

      IF (X_ADV_SWITCH .AND. (.NOT. FULL_FLUID_MODE)) THEN

        CALL FILL_X_ADV                                                         !Call x-advection matrix builder
        FL_INDEX = M_X_ADV%N_NONZ                                               !Increment flush index

      END IF

      IF (MAXWELL_SWITCH .AND. (.NOT. FULL_FLUID_MODE)) THEN

        IF (nd .EQ. 1 .AND. MOD(MIN_X,2) .EQ. 0) THEN

          CALL FILL_MAXWELL                                                     !Call Maxwell law matrix builder
          FL_INDEX = FL_INDEX + M_MAXWELL%N_NONZ

        ELSE IF (nd .GT. 1) THEN

          CALL FILL_MAXWELL                                                     !Call Maxwell law matrix builder
          FL_INDEX = FL_INDEX + M_MAXWELL%N_NONZ

        END IF

      END IF

      IF (E_ADV_SWITCH .AND. (.NOT. FULL_FLUID_MODE)) THEN

        CALL INIT_E_ADV_SP                                                      !Initialize E-advection sparsity pattern
        E_ADV_ON = .TRUE.

      END IF

    END IF

!Initialize generic sparisty patterns
    CALL INIT_TRI_DIAG_SP
    CALL INIT_DENSE_SP
    CALL INIT_DIAG_SP
    CALL INIT_ID_MAT_SP

!Allocate local markers
    IF (COLL_EE_0_SWITCH) THEN

      ALLOCATE(MARKER_EE_0(MIN_X:MAX_X))
      MARKER_EE_0 = 0

    END IF

    IF (COLL_EN_EL_0_SWITCH) THEN

      ALLOCATE(MARKER_EN_EL_0(MIN_X:MAX_X))
      MARKER_EN_EL_0 = 0

    END IF

    IF (COLL_EE_L_SWITCH) THEN

      ALLOCATE(MARKER_EE_L(MIN_X:MAX_X,NUM_H))
      MARKER_EE_L = 0

    END IF

    IF (COLL_EI_L_SWITCH) THEN

      ALLOCATE(MARKER_EI_L(MIN_X:MAX_X,NUM_H))
      MARKER_EI_L = 0

      IF (COLD_ION_FLUID_SWITCH) THEN

        ALLOCATE(MARKER_CI_EI_L(MIN_X:MAX_X,NUM_H,NUM_V))
        ALLOCATE(MARKER_CI_EI_R(MIN_X:MAX_X,NUM_H,NUM_V))
        MARKER_CI_EI_L = 0
        MARKER_CI_EI_R = 0

        CALL INIT_CI_EI_OD_SP

      END IF

    END IF

    IF (COLL_EN_EL_L_SWITCH) THEN

      ALLOCATE(MARKER_EN_EL_L(MIN_X:MAX_X,NUM_H))
      MARKER_EN_EL_L = 0

    END IF

    IF (COLL_EN_EX) THEN

      ALLOCATE(MARKER_CRM_EX(MIN_X:MAX_X))
      MARKER_CRM_EX = 0

    END IF

    IF (COLL_EN_ION) THEN

      ALLOCATE(MARKER_CRM_ION(MIN_X:MAX_X))
      MARKER_CRM_ION = 0

    END IF

    IF (HEATING_SWITCH) THEN

      ALLOCATE(MARKER_HEATING(MIN_X:MAX_X))
      MARKER_HEATING = 0

    END IF

!Initialize upstrem particle source sparsity pattern
    IF (PART_SOURCE_SWITCH) THEN

      CALL INIT_PART_SOURCE_SP

      IF (PART_SOURCE_ON .AND. COLD_ION_FLUID_SWITCH) CALL INIT_ION_PART_SOURCE_SP

    END IF

    IF (NEUTRAL_TRACK_SWITCH) THEN

!Initialize neutral process sparsity patterns
      CALL INIT_DENSE_NEUT_SP
      CALL INIT_DIAG_NEUT_SP
      CALL INIT_EN_ION_SP
      CALL INIT_CRM_RECOMB_SP
      CALL INIT_INEL_SP

      IF (nd .EQ. 1 .AND. MOD(MIN_X,2) .EQ. 1) THEN

        IF (NEUTRAL_DIFFUSION_SWITCH .AND. BUILD_1D_MAT) THEN

          CALL INIT_NEUT_DIFF_SP

          NEUTRAL_DIFF_ON = .TRUE.

        END IF

      ELSE IF (nd .GT. 1) THEN

        IF (NEUTRAL_DIFFUSION_SWITCH) THEN

          CALL INIT_NEUT_DIFF_SP

          NEUTRAL_DIFF_ON = .TRUE.

        END IF


      END IF


      IF (RECYCLING_SWITCH .AND. (RANK .EQ. SIZE - 1)) THEN

        CALL INIT_REC_SP
        REC_ON = .TRUE.

      END IF

    END IF

!Initialize sink sparsity pattern
    IF (PLASMA_SINK_SWITCH  .AND. (RANK .EQ. SIZE - 1)) THEN

      CALL INIT_SINK_SP
      PLASMA_SINK_ON = .TRUE.

    END IF

      IF (RANK .EQ. 0) PRINT*, ' Adding auxillary matrices'
  !Add auxillary matrix sparsity patterns before flush index
      IF (FL_INDEX .GT. 0) THEN

        LOCAL_SP%N_NZ = FL_INDEX

        ALLOCATE(LOCAL_SP%ROW(FL_INDEX))
        ALLOCATE(LOCAL_SP%COL(FL_INDEX))

        IF (X_ADV_SWITCH) THEN

          DO I = 1, M_X_ADV%N_NONZ

            LOCAL_SP%ROW(I) = M_X_ADV%ROW(I)
            LOCAL_SP%COL(I) = M_X_ADV%COLUMN(I)

          END DO

        END IF

        IF (MAXWELL_SWITCH) THEN

          DO I = FL_INDEX - M_MAXWELL%N_NONZ + 1, FL_INDEX

            LOCAL_SP%ROW(I) = M_MAXWELL%ROW(I - FL_INDEX + M_MAXWELL%N_NONZ)
            LOCAL_SP%COL(I) = M_MAXWELL%COLUMN(I - FL_INDEX + M_MAXWELL%N_NONZ)

          END DO

        END IF

      END IF

    IF (RANK .EQ. 0) PRINT*, ' Calculating local sparsity pattern'
!Main spatial loop
    DO I = MIN_X, MAX_X

      OVERLAP_START_P = LOCAL_SP%N_NZ

      IF (.NOT. FULL_FLUID_MODE) THEN
  !Loop over collision integrals for various f's
        DO J = 1, NUM_H

          IF (((MOD(J-1,2) .EQ. 0) .AND. (MOD(I,2) .EQ. 1)) .OR. ((MOD(J-1,2) .EQ. 1) .AND. (MOD(I,2) .EQ. 0))) THEN !Evolve only even harmonics on odd spatial cells (centres) and vice-versa

            H_OFFSET = X_POS(I) + NUM_V * (NUM_H - J)                             !Submatrix offset

            OVERLAP_START = LOCAL_SP%N_NZ

            IF (J - 1 .EQ. 0) THEN

              IF (COLL_EE_0_SWITCH .OR. COLL_EN_EL_0_SWITCH .OR. (HEATING_SWITCH .AND. ((I - OFFSET_UP +1) /2 .LE. N_HEATING))) THEN

                IF (COLL_EE_0_SWITCH) MARKER_EE_0(I) = LOCAL_SP%N_NZ        !Set e-e collision integral for f_0^0 marker for given spatial cell
                IF (COLL_EN_EL_0_SWITCH) MARKER_EN_EL_0(I) = LOCAL_SP%N_NZ  !Set e-n elastic collision integral for f_0^0 marker for given spatial cell
                IF (HEATING_SWITCH .AND. ((I - OFFSET_UP +1) /2 .LE. N_HEATING)) THEN

                  MARKER_HEATING(I) = LOCAL_SP%N_NZ !Set heating operator markers for given spatial cell
                  HEATING_ON = .TRUE.

                END IF

                CALL ADD_SP(LOCAL_SP,TRI_DIAG_SP,H_OFFSET,H_OFFSET,.FALSE.)       !Add tri-diagonal sparsity pattern to total sparsity pattern

              END IF

            ELSE

              IF (COLL_EE_L_SWITCH) THEN

                MARKER_EE_L(I,J) = LOCAL_SP%N_NZ                                !Set e-e collision integral for l>0 marker for given spatial cell

                IF (COLL_EI_L_SWITCH) MARKER_EI_L(I,J) = LOCAL_SP%N_NZ          !Set e-i collision integral marker for given spatial cell
                IF (COLL_EN_EL_L_SWITCH) MARKER_EN_EL_L(I,J) = LOCAL_SP%N_NZ    !Set e-n elastic collision integral marker for given spatial cell

                IF (.NOT. DIAG_EE_L_SWITCH) THEN

                  CALL ADD_SP(LOCAL_SP,DENSE_SP,H_OFFSET,H_OFFSET,.FALSE.)        !Add dense sparsity pattern to total sparsity pattern if calculating full e-e collisions

                ELSE

                  CALL ADD_SP(LOCAL_SP,TRI_DIAG_SP,H_OFFSET,H_OFFSET,.FALSE.)     !Add tri-diagonal sparsity pattern to total sparsity pattern if not calculating full e-e collisions

                END IF

              ELSE IF (COLL_EI_L_SWITCH .OR. COLL_EN_EL_L_SWITCH) THEN

                IF (COLL_EI_L_SWITCH) MARKER_EI_L(I,J) = LOCAL_SP%N_NZ          !Set e-i collision integral marker for given spatial cell
                IF (COLL_EN_EL_L_SWITCH) MARKER_EN_EL_L(I,J) = LOCAL_SP%N_NZ    !Set e-n elastic collision integral marker for given spatial cell

                IF ((COLL_EI_L_SWITCH) .AND. (COLD_ION_FLUID_SWITCH)) THEN

                  CALL ADD_SP(LOCAL_SP,TRI_DIAG_SP,H_OFFSET,H_OFFSET,.FALSE.)     !Add tri-diagonal sparsity pattern if cold ions on

                ELSE

                  CALL ADD_SP(LOCAL_SP, DIAG_SP,H_OFFSET,H_OFFSET,.FALSE.)          !Add diagonal sparisty pattern to total sparsity pattern if not dealing with e-e collisions

                END IF

              END IF

  !Cold ion fluid e-i collisions

              IF (COLD_ION_FLUID_SWITCH) THEN

                IF (COLL_EI_L_SWITCH) THEN

                  IF (MOD(J-1,2) .EQ. 0) THEN

                    CALL ADD_SP(LOCAL_SP, CI_EI_OD_SP,H_OFFSET,X_POS(I),.FALSE.)

                    CALL INIT_CI_EI_OD_MARKER(MARKER_CI_EI_R(I,J,:),H_OFFSET,X_POS(I))

                  ELSE

                    IF (I .EQ. 1) THEN

                      CALL ADD_SP(LOCAL_SP, CI_EI_OD_SP,H_OFFSET,X_POS(I+1),.FALSE.)

                      CALL INIT_CI_EI_OD_MARKER(MARKER_CI_EI_R(I,J,:),H_OFFSET,X_POS(I+1))

                      IF (PERIODIC_BOUNDARY_SWITCH) THEN

                        CALL ADD_SP(LOCAL_SP, CI_EI_OD_SP,H_OFFSET,X_POS(NUM_X),.FALSE.)

                        CALL INIT_CI_EI_OD_MARKER(MARKER_CI_EI_L(I,J,:),H_OFFSET,X_POS(NUM_X))

                      END IF

                    ELSE IF (I .EQ. NUM_X) THEN

                      CALL ADD_SP(LOCAL_SP, CI_EI_OD_SP,H_OFFSET,X_POS(I - 1),.FALSE.)

                      CALL INIT_CI_EI_OD_MARKER(MARKER_CI_EI_L(I,J,:),H_OFFSET,X_POS(I-1))

                      IF (PERIODIC_BOUNDARY_SWITCH) THEN

                        CALL ADD_SP(LOCAL_SP, CI_EI_OD_SP,H_OFFSET,0,.FALSE.)

                        CALL INIT_CI_EI_OD_MARKER(MARKER_CI_EI_R(I,J,:),H_OFFSET,0)

                      END IF

                    ELSE

                      CALL ADD_SP(LOCAL_SP, CI_EI_OD_SP,H_OFFSET,X_POS(I+1),.FALSE.)

                      CALL INIT_CI_EI_OD_MARKER(MARKER_CI_EI_R(I,J,:),H_OFFSET,X_POS(I+1))

                      CALL ADD_SP(LOCAL_SP, CI_EI_OD_SP,H_OFFSET,X_POS(I - 1),.FALSE.)

                      CALL INIT_CI_EI_OD_MARKER(MARKER_CI_EI_L(I,J,:),H_OFFSET,X_POS(I-1))

                    END IF

                  END IF

                END IF

              END IF

            END IF

  !Add inelastic collision sparisty patterns
            IF (NEUTRAL_TRACK_SWITCH) THEN

              IF (COLL_EN_EX) THEN

  !Loop over transitions
                DO K = 1, NUM_NEUTRALS

                  DO L = K + 1, NUM_NEUTRALS

                    CALL ADD_SP(LOCAL_SP, INEL_SP(K,L),H_OFFSET,H_OFFSET,.TRUE.,OVERLAP_START)         !Add inelastic collision sparisty pattern for given transition

                  END DO

                  DO L = 1, K - 1

                    CALL ADD_SP(LOCAL_SP, INEL_SP(K,L),H_OFFSET,H_OFFSET,.TRUE.,OVERLAP_START)         !Add inelastic collision sparisty pattern for given transition

                  END DO

                END DO

              END IF

              IF (COLL_EN_ION) THEN

                DO K = 1, NUM_NEUTRALS

                  CALL ADD_SP(LOCAL_SP,INEL_SP(K,0),H_OFFSET,H_OFFSET,.TRUE.,OVERLAP_START)            !Add ionization sparsity pattern

                  IF (J - 1 .EQ. 0) CALL ADD_SP(LOCAL_SP,EN_ION_SP,H_OFFSET,H_OFFSET,.TRUE.,OVERLAP_START) !If l=0 add cold secondary electron generation from ionizing collisions

                END DO

              END IF

              IF (COLL_RECOMB) THEN

                DO K = 1, NUM_NEUTRALS

                  CALL ADD_SP(LOCAL_SP,INEL_SP(0,K),H_OFFSET,H_OFFSET,.TRUE.,OVERLAP_START)            !Add recombination sparsity pattern

                END DO

              END IF

            END IF

          END IF

        END DO

      END IF

!Add CRM sparsity patterns
      IF (NEUTRAL_TRACK_SWITCH) THEN

        H_OFFSET = X_POS(I) + NUM_V * NUM_H                                     !Offset for neutral subvector position

        IF (MOD(I,2) .EQ. 1) THEN                                               !Evolve neutral density only on odd spatial cells (centres)

          IF (COLL_EN_EX) THEN

            MARKER_CRM_EX(I) = LOCAL_SP%N_NZ                                    !Set excitation CR model marker
            IF (COLL_EN_ION) MARKER_CRM_ION((I)) = LOCAL_SP%N_NZ                !Set ionization CR model marker

            CALL ADD_SP(LOCAL_SP, DENSE_NEUT_SP,H_OFFSET,H_OFFSET,.FALSE.)      !Add dense neutral sparsity pattern for excitation and de-excitation

          ELSE IF (COLL_EN_ION) THEN

            MARKER_CRM_ION(I) = LOCAL_SP%N_NZ                                   !Set ionization CR model marker if excitation not treated

            CALL ADD_SP(LOCAL_SP, DIAG_NEUT_SP, H_OFFSET, H_OFFSET,.FALSE.)     !Add ionization sparsity pattern (diagonal) if excitation not treated

          END IF

            IF (COLL_RECOMB .AND. (.NOT. FULL_FLUID_MODE)) THEN

                CALL ADD_SP(LOCAL_SP, CRM_RECOMB_SP, H_OFFSET - NUM_V, H_OFFSET - NUM_V,.TRUE.,OVERLAP_START_P) !Add recombination sparsity pattern and secondary electron loss

            END IF

        END IF

      END IF

    END DO

!Add neutral diffusion pattern
    IF (nd .EQ. 1 .AND. MOD(MIN_X,2) .EQ. 1) THEN

      IF (NEUTRAL_DIFFUSION_SWITCH .AND. BUILD_1D_MAT) CALL ADD_SP(LOCAL_SP, NEUTRAL_DIFF_SP,0,0)

    ELSE IF (nd .GT. 1) THEN

      IF (NEUTRAL_DIFFUSION_SWITCH) CALL ADD_SP(LOCAL_SP, NEUTRAL_DIFF_SP,0,0)

    END IF

    IF (.NOT. FULL_FLUID_MODE) THEN
!Add sink pattern
      IF (PLASMA_SINK_ON) CALL ADD_SP(LOCAL_SP, SINK_SP,0,0)


      IF (E_ADV_SWITCH .AND. ((NUM_X .GT. 1) .AND. BUILD_1D_MAT)) THEN

        MARKER_E_ADV = LOCAL_SP%N_NZ                                              !Set E-advection marker
        CALL ADD_SP(LOCAL_SP, E_ADV_SP,0,0,.FALSE.)                               !Add E-advection sparsity pattern

      END IF

    END IF

!Add cold fluid ion patterns
    IF ((COLD_ION_FLUID_SWITCH .OR. FULL_FLUID_MODE) .AND. ((NUM_X .GT. 1) .AND. BUILD_1D_MAT)) THEN

      FLUID_ION_INDEX = LOCAL_SP%N_NZ

      IF (nd .EQ. 1 .AND. MOD(MIN_X,2) .EQ. 0) THEN

        CALL ADD_FLUID_IONS_BOUNDARY

        ION_FLOW_ON = .TRUE.

      ELSE IF (nd .GT. 1) THEN

        CALL ADD_FLUID_IONS_BOUNDARY

        ION_FLOW_ON = .TRUE.

      END IF

      IF (nd .EQ. 1 .AND. MOD(MIN_X,2) .EQ. 1) THEN

        IF (.NOT. ION_CONT_OFF_SWITCH) THEN

          CALL ADD_FLUID_IONS_CENTRE

          ION_CONT_ON = .TRUE.

        END IF

      ELSE IF (nd .GT. 1) THEN

        IF (.NOT. ION_CONT_OFF_SWITCH) THEN

          CALL ADD_FLUID_IONS_CENTRE

          ION_CONT_ON = .TRUE.

        END IF

      END IF

    END IF

!Add fluid electron patterns
    IF ((FULL_FLUID_MODE) .AND. ((NUM_X .GT. 1) .AND. BUILD_1D_MAT)) THEN

      IF (nd .EQ. 1 .AND. MOD(MIN_X,2) .EQ. 0) THEN

        CALL ADD_FLUID_EL_BOUNDARY

        EL_FLOW_ON = .TRUE.

      ELSE IF (nd .GT. 1) THEN

        CALL ADD_FLUID_EL_BOUNDARY

        EL_FLOW_ON = .TRUE.

      END IF

      IF (nd .EQ. 1 .AND. MOD(MIN_X,2) .EQ. 1) THEN

        CALL ADD_FLUID_EL_CENTRE

        EL_CONT_ON = .TRUE.

      ELSE IF (nd .GT. 1) THEN

        CALL ADD_FLUID_EL_CENTRE

        EL_CONT_ON = .TRUE.

      END IF

    END IF

!Add recycling pattern
    IF (RECYCLING_SWITCH .AND. ((NUM_X .GT. 1) .AND. BUILD_1D_MAT) .AND. (RANK .EQ. SIZE - 1)) THEN

      IF (FULL_FLUID_MODE) THEN

        CALL ADD_SP(LOCAL_SP,REC_SP,0,0)

      ELSE

        IF (COLD_ION_FLUID_SWITCH .AND. (.NOT.SONIC_OUTFLOW_DIV_SWITCH)) THEN

          CALL INIT_CI_REC_SP
          CALL ADD_SP(LOCAL_SP,CI_REC_SP,0,0,.FALSE.)

        ELSE

          IF (COLD_ION_FLUID_SWITCH) THEN

              CALL INIT_CI_REC_SP
              CALL ADD_SP(LOCAL_SP,CI_REC_SP,0,0,.FALSE.)

          END IF
          CALL ADD_SP(LOCAL_SP,REC_SP,0,0)

        END IF

      END IF

    END IF

!Add particle source patterns
    IF (PART_SOURCE_ON) THEN

      IF (FULL_FLUID_MODE) THEN

        !Full fluid part source to be added in later version

      ELSE

        CALL ADD_SP(LOCAL_SP,PART_SOURCE_SP,0,0)
        IF (COLD_ION_FLUID_SWITCH) CALL ADD_SP(LOCAL_SP,ION_PART_SOURCE_SP,0,0)

      END IF

    END IF

    IF (FULL_FLUID_MODE .AND. HEATING_SWITCH) THEN

      IF (MIN_X .LE. 2*N_HEATING - 1) THEN

        CALL INIT_EL_HEATING_SP
        CALL ADD_SP(LOCAL_SP,EL_HEATING_SP,0,0)
        CALL INIT_SIMPLE_MARKER(EL_HEATING_SP,MARKER_EL_HEATING,FLUID_ION_INDEX)

        HEATING_ON = .TRUE.

      END IF

    END IF

    CALL ADD_SP(LOCAL_SP,ID_MAT_SP,0,0)                                         !Add main diagonal sparsity pattern

!Initialize all remaining non-trivial markers
    IF (RANK .EQ. 0) PRINT*, ' Initializing markers'
    IF (NEUTRAL_TRACK_SWITCH) THEN

      IF (.NOT. FULL_FLUID_MODE) THEN

        CALL INIT_INEL_MARKER
        CALL INIT_CRM_RECOMB_MARKER
        CALL INIT_SEC_EL_MARKER

      END IF

!Neutral diffusion marker
      IF (NEUTRAL_DIFFUSION_SWITCH .AND. ((NUM_X .GT. 1) .AND. BUILD_1D_MAT)) &
          CALL INIT_SIMPLE_MARKER(NEUTRAL_DIFF_SP,MARKER_NEUT_DIFF,FL_INDEX)

!Recycling marker
      IF (REC_ON) THEN

        IF (FULL_FLUID_MODE) THEN

          CALL INIT_SIMPLE_MARKER(REC_SP,MARKER_REC,FLUID_ION_INDEX)

        ELSE

          IF (COLD_ION_FLUID_SWITCH .AND. (.NOT.SONIC_OUTFLOW_DIV_SWITCH)) THEN

            CALL INIT_SIMPLE_MARKER(CI_REC_SP,MARKER_CI_REC,FLUID_ION_INDEX)

          ELSE

            CALL INIT_SIMPLE_MARKER(CI_REC_SP,MARKER_CI_REC,FLUID_ION_INDEX)

            CALL INIT_SIMPLE_MARKER(REC_SP,MARKER_REC,FL_INDEX)

          END IF

        END IF

      END IF

    END IF

!Plasma sink marker
    IF (PLASMA_SINK_ON .AND. (.NOT. FULL_FLUID_MODE)) &
    CALL INIT_SIMPLE_MARKER(SINK_SP,MARKER_SINK,0)

!Particle source markers
    IF (PART_SOURCE_ON) THEN

      IF (FULL_FLUID_MODE) THEN

        !Particle source full fluid markers

      ELSE

        CALL INIT_SIMPLE_MARKER(PART_SOURCE_SP,MARKER_PART_SOURCE,FL_INDEX)
        IF (COLD_ION_FLUID_SWITCH .AND. (.NOT. ION_CONT_OFF_SWITCH)) &
            CALL INIT_SIMPLE_MARKER(ION_PART_SOURCE_SP,MARKER_ION_PART_SOURCE,FLUID_ION_INDEX)

      END IF

    END IF

    CALL INIT_SIMPLE_MARKER(ID_MAT_SP,MARKER_ID,FL_INDEX)                       !Initializes identity matrix marker

!Allocate and fill main sparse matrix with auxillary sparse matrices
    CALL ALLOCATE_SPARSE(LOCAL_M,LOC_ROWS,DIM_F,LOCAL_SP%N_NZ)

    LOCAL_M%ROW = LOCAL_SP%ROW
    LOCAL_M%COLUMN = LOCAL_SP%COL

    LOCAL_M%VALUE = 0

    IF (.NOT. FULL_FLUID_MODE) THEN

      IF ((NUM_X .GT. 1) .AND. BUILD_1D_MAT) THEN

        IF (X_ADV_SWITCH) THEN

          DO I = 1, M_X_ADV%N_NONZ

            LOCAL_M%VALUE(I) = M_X_ADV%VALUE(I)

          END DO

        END IF

        IF (MAXWELL_SWITCH) THEN

          DO I = FL_INDEX - M_MAXWELL%N_NONZ + 1, FL_INDEX

            LOCAL_M%VALUE(I) = M_MAXWELL%VALUE(I - (FL_INDEX - M_MAXWELL%N_NONZ))

          END DO

        END IF

      END IF

    END IF

!Allocate the auxillary fixed sparse matrices to a reusable matrix strip
    IF (FL_INDEX .GT. 0) THEN

      CALL ALLOCATE_SPARSE(FIXED_MAT,LOC_ROWS,DIM_F,FL_INDEX)

      FIXED_MAT%ROW = LOCAL_SP%ROW(1:FL_INDEX)
      FIXED_MAT%COLUMN = LOCAL_SP%COL(1:FL_INDEX)

      FIXED_MAT%VALUE = LOCAL_M%VALUE(1:FL_INDEX)

    END IF

    IF (RANK .EQ. 0) PRINT*, ' Counting number of nonzeros in each row'
!Count numbers of nonzeros in each row
    CALL COUNT_NNZ

  END SUBROUTINE INIT_MATRIX_DATA
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize tri-diagonal sparsity pattern for velocity subspace
  SUBROUTINE INIT_TRI_DIAG_SP

    IMPLICIT NONE

    INTEGER :: K

    TRI_DIAG_SP%N_NZ = 3 * NUM_V - 2

    ALLOCATE(TRI_DIAG_SP%ROW(TRI_DIAG_SP%N_NZ))
    ALLOCATE(TRI_DIAG_SP%COL(TRI_DIAG_SP%N_NZ))

    DO K = 1, NUM_V

      TRI_DIAG_SP%ROW(K) = K
      TRI_DIAG_SP%COL(K) = K

    END DO

    DO K = NUM_V + 1, 2 * NUM_V - 1

      TRI_DIAG_SP%ROW(K) = K - NUM_V
      TRI_DIAG_SP%COL(K) = TRI_DIAG_SP%ROW(K) + 1

    END DO

    DO K = 2 * NUM_V , 3 * NUM_V - 2

      TRI_DIAG_SP%COL(K) = K - 2 * NUM_V + 1
      TRI_DIAG_SP%ROW(K) = TRI_DIAG_SP%COL(K) + 1

    END DO


  END SUBROUTINE INIT_TRI_DIAG_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initalize dense sparsity pattern for velocity subspace
  SUBROUTINE INIT_DENSE_SP

    IMPLICIT NONE

    INTEGER :: K,I,J

    DENSE_SP%N_NZ = NUM_V ** 2

    ALLOCATE(DENSE_SP%ROW(DENSE_SP%N_NZ))
    ALLOCATE(DENSE_SP%COL(DENSE_SP%N_NZ))

!Initialize tri-diagonal elements first
    DO K = 1, NUM_V

      DENSE_SP%ROW(K) = K
      DENSE_SP%COL(K) = K

    END DO

    DO K = NUM_V + 1, 2 * NUM_V - 1

      DENSE_SP%ROW(K) = K - NUM_V
      DENSE_SP%COL(K) = DENSE_SP%ROW(K) + 1

    END DO

    DO K = 2 * NUM_V , 3 * NUM_V - 2

      DENSE_SP%COL(K) = K - 2 * NUM_V + 1
      DENSE_SP%ROW(K) = DENSE_SP%COL(K) + 1

    END DO

    K = 3 * NUM_V - 1
!Initialize rest
    DO I = 1, NUM_V

      DO J = I + 2, NUM_V

        DENSE_SP%ROW(K) = I
        DENSE_SP%COL(K) = J

        K = K + 1

      END DO

    END DO

    DO I = 3, NUM_V

      DO J = 1, I - 2

        DENSE_SP%ROW(K) = I
        DENSE_SP%COL(K) = J

        K = K + 1

      END DO

    END DO

  END SUBROUTINE INIT_DENSE_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize diagonal sparsity pattern for velocity subspace
  SUBROUTINE INIT_DIAG_SP

    IMPLICIT NONE

    INTEGER :: K

    DIAG_SP%N_NZ = NUM_V

    ALLOCATE(DIAG_SP%ROW(DIAG_SP%N_NZ))
    ALLOCATE(DIAG_SP%COL(DIAG_SP%N_NZ))

    DO K = 1, NUM_V

      DIAG_SP%ROW(K) = K
      DIAG_SP%COL(K) = K

    END DO

  END SUBROUTINE INIT_DIAG_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize dense sparsity pattern for neutral subspace
  SUBROUTINE INIT_DENSE_NEUT_SP

    IMPLICIT NONE

    INTEGER :: K,I,J

    DENSE_NEUT_SP%N_NZ = NUM_NEUTRALS ** 2

    ALLOCATE(DENSE_NEUT_SP%ROW(DENSE_NEUT_SP%N_NZ))
    ALLOCATE(DENSE_NEUT_SP%COL(DENSE_NEUT_SP%N_NZ))

!Initialize diagonal elements first
    DO K = 1, NUM_NEUTRALS

      DENSE_NEUT_SP%ROW(K) = K
      DENSE_NEUT_SP%COL(K) = K

    END DO

    K = NUM_NEUTRALS + 1

!Initialize rest
    DO I = 1, NUM_NEUTRALS

      DO J = 1, I - 1

        DENSE_NEUT_SP%ROW(K) = I
        DENSE_NEUT_SP%COL(K) = J

        K = K + 1

      END DO

    END DO

    DO I = 1, NUM_NEUTRALS

      DO J = I + 1, NUM_NEUTRALS

        DENSE_NEUT_SP%ROW(K) = I
        DENSE_NEUT_SP%COL(K) = J

        K = K + 1

      END DO

    END DO

  END SUBROUTINE INIT_DENSE_NEUT_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize diagonal sparisty pattern for neutral subspace
  SUBROUTINE INIT_DIAG_NEUT_SP

    IMPLICIT NONE

    INTEGER :: K

    DIAG_NEUT_SP%N_NZ = NUM_NEUTRALS

    ALLOCATE(DIAG_NEUT_SP%ROW(DIAG_NEUT_SP%N_NZ))
    ALLOCATE(DIAG_NEUT_SP%COL(DIAG_NEUT_SP%N_NZ))

    DO K = 1, NUM_NEUTRALS

      DIAG_NEUT_SP%ROW(K) = K
      DIAG_NEUT_SP%COL(K) = K

    END DO

  END SUBROUTINE INIT_DIAG_NEUT_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize sparsity pattern for secondary electron generation from ionization
  SUBROUTINE INIT_EN_ION_SP

    IMPLICIT NONE

    INTEGER :: K

    EN_ION_SP%N_NZ = NUM_NEUTRALS
    ALLOCATE(EN_ION_SP%ROW(NUM_NEUTRALS))
    ALLOCATE(EN_ION_SP%COL(NUM_NEUTRALS))

    DO K = 1, NUM_NEUTRALS

      EN_ION_SP%ROW(K) = 1
      EN_ION_SP%COL(K) = NUM_V + K

    END DO

  END SUBROUTINE INIT_EN_ION_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize recombination sparsity pattern
  SUBROUTINE INIT_CRM_RECOMB_SP

    IMPLICIT NONE

    INTEGER :: K, I, J

    CRM_RECOMB_SP%N_NZ = NUM_V + NUM_V * NUM_NEUTRALS

    ALLOCATE(CRM_RECOMB_SP%ROW(CRM_RECOMB_SP%N_NZ))
    ALLOCATE(CRM_RECOMB_SP%COL(CRM_RECOMB_SP%N_NZ))

    DO K = 1, NUM_V

      CRM_RECOMB_SP%ROW(K) = 1
      CRM_RECOMB_SP%COL(K) = K

    END DO

    K = NUM_V + 1

    DO I = 1, NUM_NEUTRALS

      DO J = 1, NUM_V

        CRM_RECOMB_SP%ROW(K) = I + NUM_V
        CRM_RECOMB_SP%COL(K) = J

        K = K + 1

      END DO

    END DO

  END SUBROUTINE INIT_CRM_RECOMB_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize all inelastic sparsity patterns
  SUBROUTINE INIT_INEL_SP

    IMPLICIT NONE

    INTEGER :: I, J
    INTEGER ::  K, P, L

    LOGICAL, ALLOCATABLE :: W_LOG(:,:)                 !True if pair-summed weights non-zero (i.e. transition possible)


        ALLOCATE(INEL_SP(0:NUM_NEUTRALS,0:NUM_NEUTRALS))
        ALLOCATE(IN_DAT(0:NUM_NEUTRALS,0:NUM_NEUTRALS))
        ALLOCATE(W_LOG(NUM_V,NUM_V))

        DO I = 1, NUM_NEUTRALS

!Excitation
          DO J = I + 1, NUM_NEUTRALS

            !Set weight data (pair-summed)
            CALL INIT_INEL_SP_W(I,J,W_LOG)

            !Set sparisity pattern
            K = 1

            DO P = 1, NUM_V

              IF ((ALLOCATED(INEL_ABS(P,I,J)%A)) .OR. (W_LOG(P,P)) ) THEN

                INEL_SP(I,J)%ROW(K) = P
                INEL_SP(I,J)%COL(K) = P
                IF (ALLOCATED(INEL_ABS(P,I,J)%A))  IN_DAT(I,J)%E(K) = 1

                K = K + 1

              END IF

              DO L = P + 1, NUM_V

                IF (W_LOG(P,L)) THEN

                  INEL_SP(I,J)%ROW(K) = P
                  INEL_SP(I,J)%COL(K) = L
                  K = K + 1
                END IF

              END DO

            END DO

          END DO

!De-excitation
          DO J = 1, I - 1

            !Set weight data (pair-summed)
            CALL INIT_INEL_SP_W(I,J,W_LOG)

            !Set sparsity pattern
            K = 1

            DO P = 1, NUM_V

              IF ((ALLOCATED(INEL_ABS(P,I,J)%A)) .OR. W_LOG(P,P)) THEN

                INEL_SP(I,J)%ROW(K) = P
                INEL_SP(I,J)%COL(K) = P
                IF (ALLOCATED(INEL_ABS(P,I,J)%A))  IN_DAT(I,J)%E(K) = 1
                K = K + 1

              END IF

              DO L = 1, P - 1

                IF ((IN_DAT(I,J)%W(P,L) .GT. 0) .OR. W_LOG(P,L)) THEN

                  INEL_SP(I,J)%ROW(K) = P
                  INEL_SP(I,J)%COL(K) = L
                  K = K + 1

                END IF

              END DO

            END DO

          END DO
!Ionization

          !Set weight data (pair-summed)
          CALL INIT_INEL_SP_W(I,0,W_LOG)

          !Set sparisity pattern
          K = 1

          DO P = 1, NUM_V

            IF ((ALLOCATED(INEL_ABS(P,I,0)%A)) .OR. W_LOG(P,P)) THEN

              INEL_SP(I,0)%ROW(K) = P
              INEL_SP(I,0)%COL(K) = P
              IF (ALLOCATED(INEL_ABS(P,I,0)%A))  IN_DAT(I,0)%E(K) = 1
              K = K + 1

            END IF

            DO L = P + 1, NUM_V

              IF (W_LOG(P,L)) THEN

                INEL_SP(I,0)%ROW(K) = P
                INEL_SP(I,0)%COL(K) = L
                K = K + 1

              END IF

            END DO

          END DO
!Recombination

          !Set weight data (pair-summed)
          CALL INIT_INEL_SP_W(0,I,W_LOG)

          !Set sparisity pattern
          K = 1

          DO P = 1, NUM_V

            IF ((ALLOCATED(INEL_ABS(P,0,I)%A)) .OR. &
            (W_LOG(P,P))) THEN

              INEL_SP(0,I)%ROW(K) = P
              INEL_SP(0,I)%COL(K) = P
              IF (ALLOCATED(INEL_ABS(P,0,I)%A))  IN_DAT(0,I)%E(K) = 1
              K = K + 1

            END IF

            DO L = 1, P - 1

              IF (W_LOG(P,L)) THEN

                INEL_SP(0,I)%ROW(K) = P
                INEL_SP(0,I)%COL(K) = L
                K = K + 1

              END IF

            END DO

          END DO

        END DO

      END SUBROUTINE INIT_INEL_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize summed weights for given transition and set W_LOG
     SUBROUTINE INIT_INEL_SP_W(I,J,W_LOG)

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: I, J !< Transition levels

       LOGICAL, INTENT(OUT) :: W_LOG(NUM_V,NUM_V) !< Marker whether weight is nonzero

       INTEGER :: K, P, L, M

       W_LOG = .FALSE.
       INEL_SP(I,J)%N_NZ = INEL_NNZ(I,J)

       ALLOCATE (INEL_SP(I,J)%ROW(INEL_NNZ(I,J)))
       ALLOCATE (INEL_SP(I,J)%COL(INEL_NNZ(I,J)))

       ALLOCATE (IN_DAT(I,J)%E(INEL_NNZ(I,J)))
       ALLOCATE (IN_DAT(I,J)%W(NUM_V,NUM_V))

       IN_DAT(I,J)%E = 0
       IN_DAT(I,J)%W = 0

       DO K = 1, NUM_V

         DO P = 1, NUM_V

           DO L = 1, INEL_EM(K,I,J)%NE

             IF (P .EQ. INEL_EM(K,I,J)%E(L)) THEN

               DO M = 1, INEL_EM(K,I,J)%P(L)

                 IF (INEL_EM(K,I,J)%ABS_LABEL(L)%L(M) .EQ. 1) THEN

                   IN_DAT(I,J)%W(K,P) = IN_DAT(I,J)%W(K,P) + INEL_ABS(P,I,J)%W1(INEL_EM(K,I,J)%PAIR_LABEL(L)%L(M))
                   W_LOG(K,P) = .TRUE.

                 ELSE IF (INEL_EM(K,I,J)%ABS_LABEL(L)%L(M) .EQ. 2) THEN

                   IN_DAT(I,J)%W(K,P) = IN_DAT(I,J)%W(K,P) + INEL_ABS(P,I,J)%W2(INEL_EM(K,I,J)%PAIR_LABEL(L)%L(M))
                   W_LOG(K,P) = .TRUE.

                 END IF

               END DO

             END IF

           END DO

         END DO

       END DO

     END SUBROUTINE INIT_INEL_SP_W
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize E-field advection sparsity pattern and E-field advection nonzero element properties
  SUBROUTINE INIT_E_ADV_SP

    IMPLICIT NONE

    INTEGER :: I, J, K, HL, P, &
               OFFSET,& !< Offset for implicit E-field interpolation
               N_NZ, RANK, SIZE, IERR, OFFSET_0, OFFSET_SIZE

    REAL(KIND=DEFAULT_REAL) :: c, d, &
                    a_minus, a_plus  !Interpolation quantities

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

    ! Calculate number of nonzeroes
    DO I = MIN_X, MAX_X

      c = 1.00D00

      IF (((I .EQ. 1) .AND. NO_FLOW_BOUNDARY_UP_SWITCH) .OR. ((I .EQ. NUM_X) .AND. NO_FLOW_BOUNDARY_DIV_SWITCH)) THEN

        c = 0.5D00

        IF ((I .EQ. NUM_X) .AND. (PLASMA_SINK_SWITCH .OR. RECYCLING_SWITCH))  c = 1.00D00

      END IF

      DO HL = 0, L_MAX

        IF (((MOD(HL,2) .EQ. 0) .AND. (MOD(I,2) .EQ. 1)) .OR. ((MOD(HL,2) .EQ. 1) .AND. (MOD(I,2) .EQ. 0))) THEN

          d = 1.00D00

          IF (MOD(HL,2) .EQ. 0) d = 2 * c

          IF ((HL - 1 .GE. 0) .OR. (HL + 1 .LE. L_MAX)) THEN

            IF ((HL + 1 .LE. L_MAX) .OR. (HL - 1 .EQ. 0)) THEN

              N_NZ = N_NZ + INT(d) * NUM_V

            ELSE

              N_NZ = N_NZ + INT(d) * (NUM_V - 1)

            END IF

          END IF

        END IF

      END DO

    END DO

    E_ADV_SP%N_NZ = N_NZ

    ALLOCATE(E_ADV_SP%ROW(N_NZ))
    ALLOCATE(E_ADV_SP%COL(N_NZ))

    ALLOCATE(PROP_E_ADV%H(N_NZ))
    ALLOCATE(PROP_E_ADV%V(N_NZ))
    ALLOCATE(PROP_E_ADV%P(N_NZ))
    ALLOCATE(PROP_E_ADV%M(N_NZ))

    !Offsets for first and last processor (they need to fill in custom elements - need to exclude these from next loop)
    OFFSET_0 = 0
    OFFSET_SIZE = 0

    K = 1
!First first cell elements
    IF (RANK .EQ. 0) THEN

      IF (.NOT. FIXED_BOUNDARY_UP_SWITCH) THEN

        OFFSET_0 = 1

        DO I = 1, NUM_H

          IF (MOD(I - 1,2) .EQ. 0) THEN

            IF ((I - 2 .GE. 0) .OR. (I .LE. L_MAX)) THEN

              IF ((I .LE. L_MAX) .OR. (I - 2 .EQ. 0)) THEN

                OFFSET = 0

              ELSE

                OFFSET = 1

              END IF

              DO J = 1 + OFFSET, NUM_V

                E_ADV_SP%ROW(K) = NUM_V * (NUM_H - I) + J
                E_ADV_SP%COL(K) = X_POS(2) + NUM_V * NUM_H + NUM_NEUTRALS + 1

                PROP_E_ADV%H(K) = I
                PROP_E_ADV%V(K) = J
                PROP_E_ADV%P(K) = 1
                PROP_E_ADV%M(K) = 0.500D00

                K = K + 1

                IF (PERIODIC_BOUNDARY_SWITCH) THEN

                  E_ADV_SP%ROW(K) = NUM_V * (NUM_H - I) + J
                  E_ADV_SP%COL(K) = X_POS(NUM_X) + NUM_V * NUM_H + NUM_NEUTRALS + 1

                  PROP_E_ADV%H(K) = I
                  PROP_E_ADV%V(K) = J
                  PROP_E_ADV%P(K) = 1
                  PROP_E_ADV%M(K) = 0.500D00

                  K = K + 1

                END IF

              END DO

            END IF

          END IF

        END DO

      END IF

    END IF

    IF ((RANK .EQ. SIZE - 1) .AND. (.NOT.FIXED_BOUNDARY_DIV_SWITCH)) OFFSET_SIZE = 1
!Bulk elements (excluding firs/last cell)
    DO P = MIN_X + OFFSET_0, MAX_X - OFFSET_SIZE

      !Calculate interpolation quantities
      a_minus = (X_GRID(P + 1) - X_GRID(P)) / dxc(P)
      a_plus = (X_GRID(P) - X_GRID(P - 1)) / dxc(P)

      DO I = 1, NUM_H

        !Cell centres
        IF ((MOD(I - 1,2) .EQ. 0) .AND. (MOD(P,2) .EQ. 1)) THEN

          IF ((I - 2 .GE. 0) .OR. (I .LE. L_MAX)) THEN

            IF ((I .LE. L_MAX) .OR. (I - 2 .EQ. 0)) THEN

              OFFSET = 0

            ELSE

              OFFSET = 1

            END IF

            DO J = 1 + OFFSET, NUM_V

              E_ADV_SP%ROW(K) = X_POS(P) +  NUM_V * (NUM_H - I) + J
              E_ADV_SP%COL(K) = X_POS(P + 1) + NUM_V * NUM_H + NUM_NEUTRALS + 1

              PROP_E_ADV%H(K) = I
              PROP_E_ADV%V(K) = J
              PROP_E_ADV%P(K) = P
              PROP_E_ADV%M(K) = a_plus

              K = K + 1

              E_ADV_SP%ROW(K) = X_POS(P) + NUM_V * (NUM_H - I) + J
              E_ADV_SP%COL(K) = X_POS(P - 1) + NUM_V * NUM_H + NUM_NEUTRALS + 1

              PROP_E_ADV%H(K) = I
              PROP_E_ADV%V(K) = J
              PROP_E_ADV%P(K) = P
              PROP_E_ADV%M(K) = a_minus

              K = K + 1

            END DO

          END IF

        END IF

        !Cell boundaries
        IF ((MOD(I - 1,2) .EQ. 1) .AND. (MOD(P,2) .EQ. 0)) THEN

          IF ((I - 2 .GE. 0) .OR. (I .LE. L_MAX)) THEN

            IF ((I .LE. L_MAX) .OR. (I - 2 .EQ. 0)) THEN

              OFFSET = 0

            ELSE

              OFFSET = 1

            END IF

            DO J = 1 + OFFSET, NUM_V

              E_ADV_SP%ROW(K) = X_POS(P) +  NUM_V * (NUM_H - I) + J
              E_ADV_SP%COL(K) = X_POS(P) + NUM_V * NUM_H + NUM_NEUTRALS + 1

              PROP_E_ADV%H(K) = I
              PROP_E_ADV%V(K) = J
              PROP_E_ADV%P(K) = P
              PROP_E_ADV%M(K) = 1.00D00

              K = K + 1

            END DO

          END IF

        END IF

      END DO

    END DO

!Last cell
    IF (RANK .EQ. SIZE - 1) THEN

      IF (.NOT. FIXED_BOUNDARY_DIV_SWITCH) THEN

        DO I = 1, NUM_H
          !If last cell is centre
          IF ((MOD(I - 1,2) .EQ. 0) .AND. (MOD(NUM_X,2) .EQ. 1)) THEN

            IF ((I - 2 .GE. 0) .OR. (I .LE. L_MAX)) THEN

              IF ((I .LE. L_MAX) .OR. (I - 2 .EQ. 0)) THEN

                OFFSET = 0

              ELSE

                OFFSET = 1

              END IF

              IF (RECYCLING_SWITCH .OR. PLASMA_SINK_SWITCH) THEN

                DO J = 1 + OFFSET, NUM_V

                  E_ADV_SP%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                  E_ADV_SP%COL(K) = X_POS(NUM_X - 1) + NUM_V * NUM_H + NUM_NEUTRALS + 1

                  PROP_E_ADV%H(K) = I
                  PROP_E_ADV%V(K) = J
                  PROP_E_ADV%P(K) = NUM_X
                  PROP_E_ADV%M(K) = 1.00 + (X_GRID(NUM_X) - X_GRID(NUM_X - 1)) / (X_GRID(NUM_X -1) - X_GRID(NUM_X - 3))

                  K = K + 1

                  E_ADV_SP%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                  E_ADV_SP%COL(K) = X_POS(NUM_X - 3) + NUM_V * NUM_H + NUM_NEUTRALS + 1

                  PROP_E_ADV%H(K) = I
                  PROP_E_ADV%V(K) = J
                  PROP_E_ADV%P(K) = NUM_X
                  PROP_E_ADV%M(K) = - (X_GRID(NUM_X) - X_GRID(NUM_X - 1)) / (X_GRID(NUM_X -1) - X_GRID(NUM_X - 3))

                  K = K + 1


                END DO

              ELSE

                DO J = 1 + OFFSET, NUM_V

                  E_ADV_SP%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                  E_ADV_SP%COL(K) = X_POS(NUM_X - 1) + NUM_V * NUM_H + NUM_NEUTRALS + 1

                  PROP_E_ADV%H(K) = I
                  PROP_E_ADV%V(K) = J
                  PROP_E_ADV%P(K) = NUM_X
                  PROP_E_ADV%M(K) = 0.50D00

                  K = K + 1

                  IF (PERIODIC_BOUNDARY_SWITCH) THEN

                    E_ADV_SP%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                    E_ADV_SP%COL(K) = + NUM_V * NUM_H + NUM_NEUTRALS + 1

                    PROP_E_ADV%H(K) = I
                    PROP_E_ADV%V(K) = J
                    PROP_E_ADV%P(K) = NUM_X
                    PROP_E_ADV%M(K) = 0.50D00

                    K = K + 1

                  END IF

                END DO

              END IF

            END IF

          END IF
          !If last cell is boundary
          IF ((MOD(I - 1,2) .EQ. 1) .AND. (MOD(NUM_X,2) .EQ. 0)) THEN

            IF ((I - 2 .GE. 0) .OR. (I .LE. L_MAX)) THEN

              IF ((I .LE. L_MAX) .OR. (I - 2 .EQ. 0)) THEN

                OFFSET = 0

              ELSE

                OFFSET = 1

              END IF

              DO J = 1 + OFFSET, NUM_V

                E_ADV_SP%ROW(K) = X_POS(NUM_X) + NUM_V * (NUM_H - I) + J
                E_ADV_SP%COL(K) = X_POS(NUM_X) +  NUM_V * NUM_H + NUM_NEUTRALS + 1

                PROP_E_ADV%H(K) = I
                PROP_E_ADV%V(K) = J
                PROP_E_ADV%P(K) = NUM_X
                PROP_E_ADV%M(K) = 1.00D00

                K = K + 1

              END DO

            END IF

          END IF

        END DO

      END IF

    END IF

  END SUBROUTINE INIT_E_ADV_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize neutral diffusion sparsity pattern
  SUBROUTINE INIT_NEUT_DIFF_SP

    IMPLICIT NONE

    INTEGER :: I, K,J, NUM_BOUNDARIES, RANK, SIZE, IERR

    INTEGER :: MIN_C, MAX_C !First and last bulk cell centres in local band

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

!Calculate number of nonzeroes
    NUM_BOUNDARIES = 0
    IF (RANK .EQ. 0) NUM_BOUNDARIES = NUM_BOUNDARIES + 1
    IF (RANK .EQ. SIZE - 1) NUM_BOUNDARIES = NUM_BOUNDARIES + 1

    NEUTRAL_DIFF_SP%N_NZ = (3 * LOC_NUM_C - NUM_BOUNDARIES)*NUM_NEUTRALS

    ALLOCATE(NEUTRAL_DIFF_SP%ROW(NEUTRAL_DIFF_SP%N_NZ))
    ALLOCATE(NEUTRAL_DIFF_SP%COL(NEUTRAL_DIFF_SP%N_NZ))
    ALLOCATE(NEUTRAL_DIFF_N(NEUTRAL_DIFF_SP%N_NZ))

!Determine local cell centres
    IF (RANK .EQ. 0) THEN

      MIN_C = 3

    ELSE

      MIN_C = MIN_X + MOD(MIN_X + 1, 2)

    END IF

    IF (RANK .EQ. SIZE - 1) THEN

      MAX_C = NUM_X - 2

    ELSE

      MAX_C = MAX_X - MOD(MAX_X + 1, 2)

    END IF

    K = 1

    DO J = 1, NUM_NEUTRALS
  !First cell centre
      IF (RANK .EQ. 0) THEN

        NEUTRAL_DIFF_SP%ROW(K) = NUM_H * NUM_V + J
        NEUTRAL_DIFF_SP%COL(K) = NUM_H * NUM_V + J

        NEUTRAL_DIFF_N(K) = J

        K = K + 1

        NEUTRAL_DIFF_SP%ROW(K) = NUM_H * NUM_V + J
        NEUTRAL_DIFF_SP%COL(K) = X_POS(3) +  NUM_H * NUM_V + J

        NEUTRAL_DIFF_N(K) = J

        K = K + 1

      END IF

      DO I = MIN_C, MAX_C, 2

        NEUTRAL_DIFF_SP%ROW(K) = X_POS(I) + NUM_H * NUM_V + J
        NEUTRAL_DIFF_SP%COL(K) = X_POS(I) + NUM_H * NUM_V + J

        NEUTRAL_DIFF_N(K) = J

        K = K + 1

        NEUTRAL_DIFF_SP%ROW(K) = X_POS(I) + NUM_H * NUM_V + J
        NEUTRAL_DIFF_SP%COL(K) = X_POS(I - 2) + NUM_H * NUM_V + J

        NEUTRAL_DIFF_N(K) = J

        K = K + 1

        NEUTRAL_DIFF_SP%ROW(K) = X_POS(I) + NUM_H * NUM_V + J
        NEUTRAL_DIFF_SP%COL(K) = X_POS(I + 2) + NUM_H * NUM_V + J

        NEUTRAL_DIFF_N(K) = J

        K = K + 1

      END DO
  !Last cell centre
      IF (RANK .EQ. SIZE - 1) THEN

        NEUTRAL_DIFF_SP%ROW(K) = X_POS(NUM_X) + NUM_H * NUM_V + J
        NEUTRAL_DIFF_SP%COL(K) = X_POS(NUM_X) + NUM_H * NUM_V + J

        NEUTRAL_DIFF_N(K) = J

        K = K + 1

        NEUTRAL_DIFF_SP%ROW(K) = X_POS(NUM_X) + NUM_H * NUM_V + J
        NEUTRAL_DIFF_SP%COL(K) = X_POS(NUM_X - 2) + NUM_H * NUM_V + J

        NEUTRAL_DIFF_N(K) = J

        K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_NEUT_DIFF_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize divertor sink sparsity pattern
  SUBROUTINE INIT_SINK_SP

    IMPLICIT NONE

    INTEGER :: I, J, K, L

    SINK_SP%N_NZ = 0

    DO J = 1, NUM_H

      IF (MOD(J,2) .EQ. 1) THEN

          SINK_SP%N_NZ = SINK_SP%N_NZ + NUM_V * NUM_H

      END IF

    END DO

    ALLOCATE(SINK_SP%ROW(SINK_SP%N_NZ))
    ALLOCATE(SINK_SP%COL(SINK_SP%N_NZ))
    ALLOCATE(SINK_SHIFT(SINK_SP%N_NZ))

    K = 1

    DO I = 1, NUM_V

      DO J = 1, NUM_H

        IF (MOD(J,2) .EQ. 1) THEN

          DO L = 1, NUM_H

            !Set shifts (odd/even harmonic - centre or boundary sampling)
            IF (MOD(L,2) .EQ. 1) THEN

              SINK_SHIFT(K) = 0

            ELSE

              SINK_SHIFT(K) = 1

            END IF

            SINK_SP%ROW(K) = X_POS(NUM_X) + NUM_V*(NUM_H - J) + I
            SINK_SP%COL(K) = X_POS(NUM_X - SINK_SHIFT(K)) + NUM_V*(NUM_H - L) + I

            K = K + 1

          END DO

        END IF

      END DO

    END DO


  END SUBROUTINE INIT_SINK_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize recycling sparsity pattern
  SUBROUTINE INIT_REC_SP

    IMPLICIT NONE

    INTEGER :: I

    IF (FULL_FLUID_MODE) THEN

      REC_SP%N_NZ = 1

      ALLOCATE(REC_SP%ROW(1))
      ALLOCATE(REC_SP%COL(1))

      REC_SP%ROW(1) = X_POS(NUM_X) + NUM_V * NUM_H + 1
      REC_SP%COL(1) = X_POS(NUM_X) + NUM_0D - 4

    ELSE

      IF (SONIC_OUTFLOW_DIV_SWITCH .AND. COLD_ION_FLUID_SWITCH .AND. (.NOT. ION_CONT_OFF_SWITCH)) THEN

        REC_SP%N_NZ = 1

        ALLOCATE(REC_SP%ROW(1))
        ALLOCATE(REC_SP%COL(1))

        REC_SP%ROW(1) = X_POS(NUM_X) + NUM_V * NUM_H + 1
        REC_SP%COL(1) = X_POS(NUM_X) + NUM_0D - 1

      ELSE

        REC_SP%N_NZ = NUM_V

        ALLOCATE(REC_SP%ROW(NUM_V))
        ALLOCATE(REC_SP%COL(NUM_V))

        DO I = 1, NUM_V

          REC_SP%ROW(I) = X_POS(NUM_X) + NUM_V * NUM_H + 1
          REC_SP%COL(I) = X_POS(NUM_X) + NUM_V * (NUM_H - 1) + I

        END DO

      END IF

    END IF

  END SUBROUTINE INIT_REC_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize electron source sparsity pattern
  SUBROUTINE INIT_PART_SOURCE_SP

    IMPLICIT NONE

    INTEGER :: I, J, K, N

    PART_SOURCE_SP%N_NZ = 0
!Count number of nonzeroes
    DO I = MIN_X, MAX_X

      IF ((MOD(I,2) .EQ. 1) .AND. (I - OFFSET_UP + 1)/2 .LE. N_PART_SOURCE) PART_SOURCE_SP%N_NZ = PART_SOURCE_SP%N_NZ + NUM_V ** 2

    END DO

    IF (PART_SOURCE_SP%N_NZ .GT. 0) THEN

      PART_SOURCE_ON = .TRUE.

      N = 1

      ALLOCATE(PART_SOURCE_SP%ROW(PART_SOURCE_SP%N_NZ))
      ALLOCATE(PART_SOURCE_SP%COL(PART_SOURCE_SP%N_NZ))
      ALLOCATE(SOURCE_DATA(PART_SOURCE_SP%N_NZ))

      DO I = MIN_X, MAX_X

        IF ((MOD(I,2) .EQ. 1) .AND. (I - OFFSET_UP + 1)/2 .LE. N_PART_SOURCE) THEN

          DO J = 1, NUM_V

            DO K = 1, NUM_V

              PART_SOURCE_SP%ROW(N) = X_POS(I) + (NUM_H - 1) * NUM_V + J
              PART_SOURCE_SP%COL(N) = X_POS(I) + (NUM_H - 1) * NUM_V + K

              SOURCE_DATA(N)%POS = I
              SOURCE_DATA(N)%J = J
              SOURCE_DATA(N)%K = K

              N = N + 1

            END DO

          END DO

        END IF

      END DO


    END IF

  END SUBROUTINE INIT_PART_SOURCE_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize diagonal element sparsity pattern for main diagonal
  SUBROUTINE INIT_ID_MAT_SP

    IMPLICIT NONE

    INTEGER :: K, RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    ID_MAT_SP%N_NZ = LOC_ROWS

    ALLOCATE(ID_MAT_SP%ROW(ID_MAT_SP%N_NZ))
    ALLOCATE(ID_MAT_SP%COL(ID_MAT_SP%N_NZ))

    DO K = 1, LOC_ROWS

      ID_MAT_SP%ROW(K) = K + NUM_0D * nd * RANK
      ID_MAT_SP%COL(K) = K + NUM_0D * nd * RANK

    END DO

  END SUBROUTINE INIT_ID_MAT_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid ion Maxwell_Ampere contribution
  SUBROUTINE INIT_MAXWELL_ION_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, IERR, OFFSET, SIZE

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    MAXWELL_ION_SP%N_NZ = nd_loc - (LOC_NUM_C + OFFSET)

    ALLOCATE(MAXWELL_ION_SP%ROW(MAXWELL_ION_SP%N_NZ))
    ALLOCATE(MAXWELL_ION_SP%COL(MAXWELL_ION_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        MAXWELL_ION_SP%ROW(K) = X_POS(I) + NUM_H*NUM_V + NUM_NEUTRALS + 1
        MAXWELL_ION_SP%COL(K) = X_POS(I) + NUM_0D

        K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_MAXWELL_ION_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid ion Lorentz force contribution to convection
  SUBROUTINE INIT_ION_LORENTZ_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, IERR, OFFSET,SIZE

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    ION_LORENTZ_SP%N_NZ = nd_loc - (LOC_NUM_C + OFFSET)

    ALLOCATE(ION_LORENTZ_SP%ROW(ION_LORENTZ_SP%N_NZ))
    ALLOCATE(ION_LORENTZ_SP%COL(ION_LORENTZ_SP%N_NZ))
    ALLOCATE(ION_LORENTZ_POS(ION_LORENTZ_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        ION_LORENTZ_SP%ROW(K) = X_POS(I) + NUM_0D
        ION_LORENTZ_SP%COL(K) = X_POS(I) + NUM_H*NUM_V + NUM_NEUTRALS + 1
        ION_LORENTZ_POS(K) = I

        K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_ION_LORENTZ_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid ion convection contribution
  SUBROUTINE INIT_ION_CONV_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, OFFSET, CONV_MIN_X, CONV_MAX_X, NUM_BOUNDARIES_CONV,NNZ,OFFSET_CONV

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET_CONV = 0

    OFFSET = 0

    NUM_BOUNDARIES_CONV = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      OFFSET_CONV = 1 !Take into account that for periodic boundary last point is a cell boundary
      NNZ = 3 * (nd_loc - LOC_NUM_C)

    ELSE

      IF (MIN_X .LT. 3) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      IF (MAX_X .GT. NUM_X - 2) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      NNZ = 3 * (nd_loc - LOC_NUM_C - OFFSET) - NUM_BOUNDARIES_CONV

    END IF

    ION_CONV_SP%N_NZ = NNZ

    CONV_MIN_X = MAX(MIN_X,4)
    CONV_MAX_X = MIN(MAX_X, NUM_X - 3 + OFFSET_CONV)

    ALLOCATE(ION_CONV_SP%ROW(ION_CONV_SP%N_NZ))
    ALLOCATE(ION_CONV_SP%COL(ION_CONV_SP%N_NZ))
    ALLOCATE(ION_CONV_PROP%SIGN(ION_CONV_SP%N_NZ))
    ALLOCATE(ION_CONV_PROP%POS(ION_CONV_SP%N_NZ))

    K = 1

    IF (MIN_X .LT. 3) THEN

      ION_CONV_SP%ROW(K) = X_POS(2) + NUM_0D
      ION_CONV_SP%COL(K) = X_POS(4) + NUM_0D
      ION_CONV_PROP%SIGN(K) = - 1
      ION_CONV_PROP%POS(K) = 2

      K = K + 1

      ION_CONV_SP%ROW(K) = X_POS(2) + NUM_0D
      ION_CONV_SP%COL(K) = X_POS(2) + NUM_0D
      ION_CONV_PROP%SIGN(K) = 0
      ION_CONV_PROP%POS(K) = 2

      K = K + 1

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        ION_CONV_SP%ROW(K) = X_POS(2) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(NUM_X) + NUM_0D
        ION_CONV_PROP%SIGN(K) = 1
        ION_CONV_PROP%POS(K) = 2

        K = K + 1

      END IF

    END IF

    DO I = CONV_MIN_X, CONV_MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        ION_CONV_SP%ROW(K) = X_POS(I) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(I+2) + NUM_0D
        ION_CONV_PROP%SIGN(K) = - 1
        ION_CONV_PROP%POS(K) = I

        K = K + 1

        ION_CONV_SP%ROW(K) = X_POS(I) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(I) + NUM_0D
        ION_CONV_PROP%SIGN(K) = 0
        ION_CONV_PROP%POS(K) = I

        K = K + 1

        ION_CONV_SP%ROW(K) = X_POS(I) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(I-2) + NUM_0D
        ION_CONV_PROP%SIGN(K) = 1
        ION_CONV_PROP%POS(K) = I

        K = K + 1

      END IF

    END DO

    IF (MAX_X .GT. NUM_X - 2) THEN

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        ION_CONV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(2) + NUM_0D
        ION_CONV_PROP%SIGN(K) = - 1
        ION_CONV_PROP%POS(K) = NUM_X

        K = K + 1

        ION_CONV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(NUM_X) + NUM_0D
        ION_CONV_PROP%SIGN(K) = 0
        ION_CONV_PROP%POS(K) = NUM_X

        K = K + 1

        ION_CONV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(NUM_X-2) + NUM_0D
        ION_CONV_PROP%SIGN(K) = 1
        ION_CONV_PROP%POS(K) = NUM_X

        K = K + 1

      ELSE

        ION_CONV_SP%ROW(K) = X_POS(NUM_X-1) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D
        ION_CONV_PROP%SIGN(K) = 0
        ION_CONV_PROP%POS(K) = NUM_X - 1

        K = K + 1

        ION_CONV_SP%ROW(K) = X_POS(NUM_X - 1) + NUM_0D
        ION_CONV_SP%COL(K) = X_POS(NUM_X - 3) + NUM_0D
        ION_CONV_PROP%SIGN(K) = 1
        ION_CONV_PROP%POS(K) = NUM_X - 1

        K = K + 1

      END IF

    END IF

  END SUBROUTINE INIT_ION_CONV_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid ion continuity contribution
  SUBROUTINE INIT_ION_CONT_SP

    IMPLICIT NONE

    INTEGER :: I, K,CONV_MIN_X, CONV_MAX_X, NUM_BOUNDARIES_CONV,NNZ

    NUM_BOUNDARIES_CONV = 0

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      NNZ = 2 * LOC_NUM_C

    ELSE

      IF ((MIN_X .EQ. 1) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH))  NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      NNZ = 2 * LOC_NUM_C - NUM_BOUNDARIES_CONV

    END IF

    ION_CONT_SP%N_NZ = NNZ

    CONV_MIN_X = MAX(MIN_X,2)
    CONV_MAX_X = MAX_X
    IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) CONV_MAX_X = NUM_X - 2

    ALLOCATE(ION_CONT_SP%ROW(ION_CONT_SP%N_NZ))
    ALLOCATE(ION_CONT_SP%COL(ION_CONT_SP%N_NZ))
    ALLOCATE(ION_CONT_PROP%SIGN(ION_CONT_SP%N_NZ))
    ALLOCATE(ION_CONT_PROP%POS(ION_CONT_SP%N_NZ))

    K = 1

    IF (MIN_X .EQ. 1) THEN

      ION_CONT_SP%ROW(K) = X_POS(1) + NUM_0D - 1
      ION_CONT_SP%COL(K) = X_POS(2) + NUM_0D
      ION_CONT_PROP%SIGN(K) = - 1
      ION_CONT_PROP%POS(K) = 1

      K = K + 1

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        ION_CONT_SP%ROW(K) = X_POS(1) + NUM_0D - 1
        ION_CONT_SP%COL(K) = X_POS(NUM_X) + NUM_0D
        ION_CONT_PROP%SIGN(K) = 1
        ION_CONT_PROP%POS(K) = 1

        K = K + 1

      END IF

    END IF

    DO I = CONV_MIN_X, CONV_MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

        ION_CONT_SP%ROW(K) = X_POS(I) + NUM_0D - 1
        ION_CONT_SP%COL(K) = X_POS(I+1) + NUM_0D
        ION_CONT_PROP%SIGN(K) = - 1
        ION_CONT_PROP%POS(K) = I

        K = K + 1

        ION_CONT_SP%ROW(K) = X_POS(I) + NUM_0D - 1
        ION_CONT_SP%COL(K) = X_POS(I-1) + NUM_0D
        ION_CONT_PROP%SIGN(K) = 1
        ION_CONT_PROP%POS(K) = I

        K = K + 1

      END IF

    END DO

    IF (MAX_X .EQ. NUM_X) THEN

      IF (MOD(NUM_X,2) .EQ. 1) THEN

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

          ION_CONT_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 1
          ION_CONT_SP%COL(K) = X_POS(1) + NUM_0D
          ION_CONT_PROP%SIGN(K) = - 1
          ION_CONT_PROP%POS(K) = NUM_X

          K = K + 1

          ION_CONT_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 1
          ION_CONT_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D
          ION_CONT_PROP%SIGN(K) = 1
          ION_CONT_PROP%POS(K) = NUM_X

          K = K + 1

        ELSE

          ION_CONT_SP%ROW(K) = X_POS(NUM_X) + NUM_0D -1
          ION_CONT_SP%COL(K) = X_POS(NUM_X - 1) + NUM_0D
          ION_CONT_PROP%SIGN(K) = 1
          ION_CONT_PROP%POS(K) = NUM_X

          K = K + 1

        END IF

      END IF

    END IF

  END SUBROUTINE INIT_ION_CONT_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize ion gain term (due to ionization) sparsity pattern
  SUBROUTINE INIT_ION_GAIN_SP

    IMPLICIT NONE

    INTEGER :: I, K, J


    ION_GAIN_SP%N_NZ = NUM_NEUTRALS * LOC_NUM_C

    ALLOCATE(ION_GAIN_SP%ROW(ION_GAIN_SP%N_NZ))
    ALLOCATE(ION_GAIN_SP%COL(ION_GAIN_SP%N_NZ))
    ALLOCATE(ION_GAIN_POS(ION_GAIN_SP%N_NZ))
    ALLOCATE(ION_GAIN_N(ION_GAIN_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

          DO J = 1, NUM_NEUTRALS

            ION_GAIN_SP%ROW(K) = X_POS(I) + NUM_0D - 1
            ION_GAIN_SP%COL(K) = X_POS(I) + NUM_H*NUM_V + J
            ION_GAIN_POS(K) = I
            ION_GAIN_N(K) = J

            K = K + 1

          END DO

      END IF

    END DO

  END SUBROUTINE INIT_ION_GAIN_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize ion loss term (due to recombination) sparsity pattern
  SUBROUTINE INIT_ION_LOSS_SP

    IMPLICIT NONE

    INTEGER :: I, K, J


    ION_LOSS_SP%N_NZ = NUM_V * LOC_NUM_C

    ALLOCATE(ION_LOSS_SP%ROW(ION_LOSS_SP%N_NZ))
    ALLOCATE(ION_LOSS_SP%COL(ION_LOSS_SP%N_NZ))
    ALLOCATE(ION_RECOMB_POS(ION_LOSS_SP%N_NZ))
    ALLOCATE(ION_RECOMB_V(ION_LOSS_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

          DO J = 1, NUM_V

            ION_LOSS_SP%ROW(K) = X_POS(I) + NUM_0D - 1
            ION_LOSS_SP%COL(K) = X_POS(I) + (NUM_H-1)*NUM_V + J
            ION_RECOMB_POS(K) = I
            ION_RECOMB_V(K) = J

            K = K + 1

          END DO

      END IF

    END DO

  END SUBROUTINE INIT_ION_LOSS_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize ion loss term (due to recombination) sparsity pattern with full fluid mode
  SUBROUTINE INIT_ION_LOSS_FF_SP

    IMPLICIT NONE

    INTEGER :: I, K


    ION_LOSS_FF_SP%N_NZ = LOC_NUM_C

    ALLOCATE(ION_LOSS_FF_SP%ROW(ION_LOSS_FF_SP%N_NZ))
    ALLOCATE(ION_LOSS_FF_SP%COL(ION_LOSS_FF_SP%N_NZ))
    ALLOCATE(ION_RECOMB_FF_POS(ION_LOSS_FF_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

          ION_LOSS_FF_SP%ROW(K) = X_POS(I) + NUM_0D - 1
          ION_LOSS_FF_SP%COL(K) = X_POS(I) + NUM_0D - 4
          ION_RECOMB_FF_POS(K) = I

          K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_ION_LOSS_FF_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize e-i drag sparsity pattern
  SUBROUTINE INIT_ION_DRAG_SP

    IMPLICIT NONE

    INTEGER :: I, K, J, OFFSET, RANK,SIZE,IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    ION_DRAG_SP%N_NZ =  (NUM_V + 2) * (nd_loc - (LOC_NUM_C + OFFSET))

    ALLOCATE(ION_DRAG_SP%ROW(ION_DRAG_SP%N_NZ))
    ALLOCATE(ION_DRAG_SP%COL(ION_DRAG_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        IF ((I .EQ. NUM_X) .AND. (PERIODIC_BOUNDARY_SWITCH)) THEN

          DO J = 1, NUM_V

            ION_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_DRAG_SP%COL(K) = X_POS(I) + (NUM_H-2)*NUM_V + J

            K = K + 1

          END DO

          ION_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
          ION_DRAG_SP%COL(K) = NUM_0D - 1

          K = K + 1

          ION_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
          ION_DRAG_SP%COL(K) = X_POS(I-1) + NUM_0D - 1

          K = K + 1

        ELSE

          DO J = 1, NUM_V

            ION_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_DRAG_SP%COL(K) = X_POS(I) + (NUM_H-2)*NUM_V + J

            K = K + 1

          END DO

          ION_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
          ION_DRAG_SP%COL(K) = X_POS(I+1) + NUM_0D - 1

          K = K + 1

          ION_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
          ION_DRAG_SP%COL(K) = X_POS(I-1) + NUM_0D - 1

          K = K + 1

        END IF

      END IF

    END DO

  END SUBROUTINE INIT_ION_DRAG_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid ion momentum equation pressure gradient contribution
  SUBROUTINE INIT_ION_PRESSURE_SP

    IMPLICIT NONE

    INTEGER :: I, J, K, RANK, SIZE, IERR, N_NZ, FF_OFFSET

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0
    FF_OFFSET = 0
    IF (FULL_FLUID_MODE) FF_OFFSET = 3

      IF ((.NOT. ION_CONT_OFF_SWITCH) .OR. FULL_FLUID_MODE) THEN

        DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2

        END DO

        ION_PRESSURE_SP%N_NZ = N_NZ

        ALLOCATE(ION_PRESSURE_SP%ROW(N_NZ))
        ALLOCATE(ION_PRESSURE_SP%COL(N_NZ))
        ALLOCATE(ION_PRESSURE_POS(N_NZ))
        ALLOCATE(ION_PRESSURE_SIGN(N_NZ))

        K = 1

        DO I = MIN_X, MAX_X

          IF ((I .EQ. NUM_X) .AND. (PERIODIC_BOUNDARY_SWITCH)) THEN

            ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_PRESSURE_SP%COL(K) = X_POS(I-1) + NUM_0D - 1 - FF_OFFSET
            ION_PRESSURE_POS(K) = I - 1
            ION_PRESSURE_SIGN(K) = 1

            K = K + 1

            ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_PRESSURE_SP%COL(K) = NUM_0D - 1 - FF_OFFSET
            ION_PRESSURE_POS(K) = 1
            ION_PRESSURE_SIGN(K) = - 1
            K = K + 1

          ELSE

            IF (MOD(I,2) .EQ. 0) THEN

              ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
              ION_PRESSURE_SP%COL(K) = X_POS(I-1) + NUM_0D - 1 - FF_OFFSET
              ION_PRESSURE_POS(K) = I - 1
              ION_PRESSURE_SIGN(K) = 1

              K = K + 1

              ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
              ION_PRESSURE_SP%COL(K) = X_POS(I+1) + NUM_0D - 1 - FF_OFFSET
              ION_PRESSURE_POS(K) = I + 1
              ION_PRESSURE_SIGN(K) = - 1

              K = K + 1

            END IF

          END IF

        END DO

      ELSE

        DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2 * NUM_V

        END DO

        ION_PRESSURE_SP%N_NZ = N_NZ

        ALLOCATE(ION_PRESSURE_SP%ROW(N_NZ))
        ALLOCATE(ION_PRESSURE_SP%COL(N_NZ))
        ALLOCATE(ION_PRESSURE_POS(N_NZ))
        ALLOCATE(ION_PRESSURE_SIGN(N_NZ))

        K = 1

        DO I = MIN_X, MAX_X

          IF ((I .EQ. NUM_X) .AND. (PERIODIC_BOUNDARY_SWITCH)) THEN

            DO J = 1, NUM_V

              ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
              ION_PRESSURE_SP%COL(K) = X_POS(I-1) + NUM_V*(NUM_H-1) + J
              ION_PRESSURE_POS(K) = I - 1
              ION_PRESSURE_SIGN(K) = 1

              K = K + 1

              ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
              ION_PRESSURE_SP%COL(K) = NUM_V*(NUM_H-1) + J
              ION_PRESSURE_POS(K) = 1
              ION_PRESSURE_SIGN(K) = - 1

              K = K + 1

            END DO

          ELSE

            IF (MOD(I,2) .EQ. 0) THEN

              DO J = 1, NUM_V

                ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
                ION_PRESSURE_SP%COL(K) = X_POS(I-1) + NUM_V*(NUM_H-1) + J
                ION_PRESSURE_POS(K) = I - 1
                ION_PRESSURE_SIGN(K) = 1

                K = K + 1

                ION_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D
                ION_PRESSURE_SP%COL(K) = X_POS(I+1) + NUM_V*(NUM_H-1) + J
                ION_PRESSURE_POS(K) = I + 1
                ION_PRESSURE_SIGN(K) = - 1

                K = K + 1

              END DO

            END IF

          END IF

        END DO

      END IF

  END SUBROUTINE INIT_ION_PRESSURE_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize recycling sparsity pattern (with cold fluid ions)
  SUBROUTINE INIT_CI_REC_SP

    IMPLICIT NONE

    CI_REC_SP%N_NZ = 2

    ALLOCATE(CI_REC_SP%ROW(CI_REC_SP%N_NZ))
    ALLOCATE(CI_REC_SP%COL(CI_REC_SP%N_NZ))

    CI_REC_SP%ROW(1) = X_POS(NUM_X) + NUM_H * NUM_V + 1
    CI_REC_SP%COL(1) = X_POS(NUM_X - 1) + NUM_0D

    CI_REC_SP%ROW(2) = X_POS(NUM_X) + NUM_H * NUM_V + 1
    CI_REC_SP%COL(2) = X_POS(NUM_X - 3) + NUM_0D

  END SUBROUTINE INIT_CI_REC_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize divertor ion convection pattern
  SUBROUTINE INIT_ION_CONV_DIV_SP

    IMPLICIT NONE

      ION_CONV_DIV_SP%N_NZ = 2

      ALLOCATE(ION_CONV_DIV_SP%ROW(ION_CONV_DIV_SP%N_NZ))
      ALLOCATE(ION_CONV_DIV_SP%COL(ION_CONV_DIV_SP%N_NZ))

      ION_CONV_DIV_SP%ROW(1) = X_POS(NUM_X - 1) + NUM_0D
      ION_CONV_DIV_SP%COL(1) = X_POS(NUM_X - 1) + NUM_0D

      ION_CONV_DIV_SP%ROW(2) = X_POS(NUM_X - 1) + NUM_0D
      ION_CONV_DIV_SP%COL(2) = X_POS(NUM_X - 3) + NUM_0D


  END SUBROUTINE INIT_ION_CONV_DIV_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize divertor ion sink pattern
  SUBROUTINE INIT_ION_SINK_SP

    IMPLICIT NONE

      ION_SINK_SP%N_NZ = 3

      ALLOCATE(ION_SINK_SP%ROW(ION_SINK_SP%N_NZ))
      ALLOCATE(ION_SINK_SP%COL(ION_SINK_SP%N_NZ))

      ION_SINK_SP%ROW(1) = X_POS(NUM_X) + NUM_0D - 1
      ION_SINK_SP%COL(1) = X_POS(NUM_X) + NUM_0D - 1

      ION_SINK_SP%ROW(2) = X_POS(NUM_X) + NUM_0D - 1
      ION_SINK_SP%COL(2) = X_POS(NUM_X - 1) + NUM_0D

      ION_SINK_SP%ROW(3) = X_POS(NUM_X) + NUM_0D - 1
      ION_SINK_SP%COL(3) = X_POS(NUM_X - 3) + NUM_0D


  END SUBROUTINE INIT_ION_SINK_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize e-i collision off-diagonal components when using cold ions
  SUBROUTINE INIT_CI_EI_OD_SP

    IMPLICIT NONE

    INTEGER :: I

    CI_EI_OD_SP%N_NZ = NUM_V

    ALLOCATE(CI_EI_OD_SP%ROW(CI_EI_OD_SP%N_NZ))
    ALLOCATE(CI_EI_OD_SP%COL(CI_EI_OD_SP%N_NZ))

    DO I = 1, NUM_V

      CI_EI_OD_SP%ROW(I) = I
      CI_EI_OD_SP%COL(I) = NUM_0D - 1

    END DO

  END SUBROUTINE INIT_CI_EI_OD_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize source term in ion momentum equation
  SUBROUTINE INIT_ION_MOM_SOURCE_SP

    IMPLICIT NONE

    INTEGER :: I,K, RANK, SIZE, IERR, OFFSET


    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1

    ION_MOM_SOURCE_SP%N_NZ = nd_loc - (LOC_NUM_C + OFFSET)

    ALLOCATE(ION_MOM_SOURCE_SP%ROW(ION_MOM_SOURCE_SP%N_NZ))
    ALLOCATE(ION_MOM_SOURCE_SP%COL(ION_MOM_SOURCE_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

          ION_MOM_SOURCE_SP%ROW(K) = X_POS(I) + NUM_0D
          ION_MOM_SOURCE_SP%COL(K) = X_POS(I) + NUM_0D

          K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_ION_MOM_SOURCE_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize ion particle source sparsity pattern
  SUBROUTINE INIT_ION_PART_SOURCE_SP

    IMPLICIT NONE

    INTEGER :: I, N, J

    ION_PART_SOURCE_SP%N_NZ = 0

!Count number of nonzeroes
    DO I = MIN_X, MAX_X

      IF ((MOD(I,2) .EQ. 1) .AND. (I - OFFSET_UP + 1)/2 .LE. N_PART_SOURCE) &
      ION_PART_SOURCE_SP%N_NZ = ION_PART_SOURCE_SP%N_NZ + NUM_V

    END DO

    N = 1

    ALLOCATE(ION_PART_SOURCE_SP%ROW(ION_PART_SOURCE_SP%N_NZ))
    ALLOCATE(ION_PART_SOURCE_SP%COL(ION_PART_SOURCE_SP%N_NZ))
    ALLOCATE(ION_SOURCE_DATA(PART_SOURCE_SP%N_NZ))

    DO I = MIN_X, MAX_X

      IF ((MOD(I,2) .EQ. 1) .AND. (I - OFFSET_UP + 1)/2 .LE. N_PART_SOURCE) THEN

        DO J = 1, NUM_V

          ION_PART_SOURCE_SP%ROW(N) = X_POS(I) + NUM_0D - 1
          ION_PART_SOURCE_SP%COL(N) = X_POS(I) + NUM_V*(NUM_H-1) + J

          ION_SOURCE_DATA(N)%K = J
          ION_SOURCE_DATA(N)%POS = I
          ION_SOURCE_DATA(N)%J = J

          N = N + 1

        END DO

      END IF

    END DO

  END SUBROUTINE INIT_ION_PART_SOURCE_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize CX friction sparsity pattern
  SUBROUTINE INIT_SIMPLE_CX_SP

    IMPLICIT NONE

    INTEGER :: I, N

    SIMPLE_CX_SP%N_NZ = 0

!Count number of nonzeroes
    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) &
      SIMPLE_CX_SP%N_NZ = SIMPLE_CX_SP%N_NZ + 1

    END DO

    N = 1

    ALLOCATE(SIMPLE_CX_SP%ROW(SIMPLE_CX_SP%N_NZ))
    ALLOCATE(SIMPLE_CX_SP%COL(SIMPLE_CX_SP%N_NZ))

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

          SIMPLE_CX_SP%ROW(N) = X_POS(I) + NUM_0D
          SIMPLE_CX_SP%COL(N) = X_POS(I) + NUM_0D

          N = N + 1

      END IF

    END DO

  END SUBROUTINE INIT_SIMPLE_CX_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid elextron Maxwell_Ampere contribution
  SUBROUTINE INIT_MAXWELL_EL_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, IERR, OFFSET, SIZE

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    MAXWELL_EL_SP%N_NZ = nd_loc - (LOC_NUM_C + OFFSET)

    ALLOCATE(MAXWELL_EL_SP%ROW(MAXWELL_EL_SP%N_NZ))
    ALLOCATE(MAXWELL_EL_SP%COL(MAXWELL_EL_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        MAXWELL_EL_SP%ROW(K) = X_POS(I) + NUM_H*NUM_V + NUM_NEUTRALS + 1
        MAXWELL_EL_SP%COL(K) = X_POS(I) + NUM_0D - 3

        K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_MAXWELL_EL_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron Lorentz force contribution to convection
  SUBROUTINE INIT_EL_LORENTZ_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, IERR, OFFSET,SIZE

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    EL_LORENTZ_SP%N_NZ = nd_loc - (LOC_NUM_C + OFFSET)

    ALLOCATE(EL_LORENTZ_SP%ROW(EL_LORENTZ_SP%N_NZ))
    ALLOCATE(EL_LORENTZ_SP%COL(EL_LORENTZ_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        EL_LORENTZ_SP%ROW(K) = X_POS(I) + NUM_0D - 3
        EL_LORENTZ_SP%COL(K) = X_POS(I) + NUM_H*NUM_V + NUM_NEUTRALS + 1

        K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_EL_LORENTZ_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron convection contribution
  SUBROUTINE INIT_EL_CONV_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, OFFSET, CONV_MIN_X, CONV_MAX_X, NUM_BOUNDARIES_CONV,NNZ,OFFSET_CONV

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET_CONV = 0

    OFFSET = 0

    NUM_BOUNDARIES_CONV = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_DIV_SWITCH) OFFSET = OFFSET + 1

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      OFFSET_CONV = 1 !Take into account that for periodic boundary last point is a cell boundary
      NNZ = 3 * (nd_loc - LOC_NUM_C)

    ELSE

      IF (MIN_X .LT. 3) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      IF (MAX_X .GT. NUM_X - 2) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      NNZ = 3 * (nd_loc - LOC_NUM_C - OFFSET) - NUM_BOUNDARIES_CONV

    END IF

    EL_CONV_SP%N_NZ = NNZ

    CONV_MIN_X = MAX(MIN_X,4)
    CONV_MAX_X = MIN(MAX_X, NUM_X - 3 + OFFSET_CONV)

    ALLOCATE(EL_CONV_SP%ROW(EL_CONV_SP%N_NZ))
    ALLOCATE(EL_CONV_SP%COL(EL_CONV_SP%N_NZ))
    ALLOCATE(EL_CONV_PROP%SIGN(EL_CONV_SP%N_NZ))
    ALLOCATE(EL_CONV_PROP%POS(EL_CONV_SP%N_NZ))

    K = 1

    IF (MIN_X .LT. 3) THEN

      EL_CONV_SP%ROW(K) = X_POS(2) + NUM_0D - 3
      EL_CONV_SP%COL(K) = X_POS(4) + NUM_0D - 3
      EL_CONV_PROP%SIGN(K) = - 1
      EL_CONV_PROP%POS(K) = 2

      K = K + 1

      EL_CONV_SP%ROW(K) = X_POS(2) + NUM_0D - 3
      EL_CONV_SP%COL(K) = X_POS(2) + NUM_0D - 3
      EL_CONV_PROP%SIGN(K) = 0
      EL_CONV_PROP%POS(K) = 2

      K = K + 1

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        EL_CONV_SP%ROW(K) = X_POS(2) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = 1
        EL_CONV_PROP%POS(K) = 2

        K = K + 1

      END IF

    END IF

    DO I = CONV_MIN_X, CONV_MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

        EL_CONV_SP%ROW(K) = X_POS(I) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(I+2) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = - 1
        EL_CONV_PROP%POS(K) = I

        K = K + 1

        EL_CONV_SP%ROW(K) = X_POS(I) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(I) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = 0
        EL_CONV_PROP%POS(K) = I

        K = K + 1

        EL_CONV_SP%ROW(K) = X_POS(I) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(I-2) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = 1
        EL_CONV_PROP%POS(K) = I

        K = K + 1

      END IF

    END DO

    IF (MAX_X .GT. NUM_X - 2) THEN

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        EL_CONV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(2) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = - 1
        EL_CONV_PROP%POS(K) = NUM_X

        K = K + 1

        EL_CONV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = 0
        EL_CONV_PROP%POS(K) = NUM_X

        K = K + 1

        EL_CONV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(NUM_X-2) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = 1
        EL_CONV_PROP%POS(K) = NUM_X

        K = K + 1

      ELSE

        EL_CONV_SP%ROW(K) = X_POS(NUM_X-1) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = 0
        EL_CONV_PROP%POS(K) = NUM_X - 1

        K = K + 1

        EL_CONV_SP%ROW(K) = X_POS(NUM_X - 1) + NUM_0D - 3
        EL_CONV_SP%COL(K) = X_POS(NUM_X - 3) + NUM_0D - 3
        EL_CONV_PROP%SIGN(K) = 1
        EL_CONV_PROP%POS(K) = NUM_X - 1

        K = K + 1

      END IF

    END IF

  END SUBROUTINE INIT_EL_CONV_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron continuity contribution
  SUBROUTINE INIT_EL_CONT_SP

    IMPLICIT NONE

    INTEGER :: I, K,CONV_MIN_X, CONV_MAX_X, NUM_BOUNDARIES_CONV,NNZ

    NUM_BOUNDARIES_CONV = 0

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      NNZ = 2 * LOC_NUM_C

    ELSE

      IF ((MIN_X .EQ. 1) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH))  NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      NNZ = 2 * LOC_NUM_C - NUM_BOUNDARIES_CONV

    END IF

    EL_CONT_SP%N_NZ = NNZ

    CONV_MIN_X = MAX(MIN_X,2)
    CONV_MAX_X = MAX_X
    IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) CONV_MAX_X = NUM_X - 2

    ALLOCATE(EL_CONT_SP%ROW(EL_CONT_SP%N_NZ))
    ALLOCATE(EL_CONT_SP%COL(EL_CONT_SP%N_NZ))
    ALLOCATE(EL_CONT_PROP%SIGN(EL_CONT_SP%N_NZ))
    ALLOCATE(EL_CONT_PROP%POS(EL_CONT_SP%N_NZ))

    K = 1

    IF (MIN_X .EQ. 1) THEN

      EL_CONT_SP%ROW(K) = X_POS(1) + NUM_0D - 4
      EL_CONT_SP%COL(K) = X_POS(2) + NUM_0D - 3
      EL_CONT_PROP%SIGN(K) = - 1
      EL_CONT_PROP%POS(K) = 1

      K = K + 1

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        EL_CONT_SP%ROW(K) = X_POS(1) + NUM_0D - 4
        EL_CONT_SP%COL(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_CONT_PROP%SIGN(K) = 1
        EL_CONT_PROP%POS(K) = 1

        K = K + 1

      END IF

    END IF

    DO I = CONV_MIN_X, CONV_MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

        EL_CONT_SP%ROW(K) = X_POS(I) + NUM_0D - 4
        EL_CONT_SP%COL(K) = X_POS(I+1) + NUM_0D - 3
        EL_CONT_PROP%SIGN(K) = - 1
        EL_CONT_PROP%POS(K) = I

        K = K + 1

        EL_CONT_SP%ROW(K) = X_POS(I) + NUM_0D - 4
        EL_CONT_SP%COL(K) = X_POS(I-1) + NUM_0D - 3
        EL_CONT_PROP%SIGN(K) = 1
        EL_CONT_PROP%POS(K) = I

        K = K + 1

      END IF

    END DO

    IF (MAX_X .EQ. NUM_X) THEN

      IF (MOD(NUM_X,2) .EQ. 1) THEN

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

          EL_CONT_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 4
          EL_CONT_SP%COL(K) = X_POS(1) + NUM_0D - 3
          EL_CONT_PROP%SIGN(K) = - 1
          EL_CONT_PROP%POS(K) = NUM_X

          K = K + 1

          EL_CONT_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 4
          EL_CONT_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D - 3
          EL_CONT_PROP%SIGN(K) = 1
          EL_CONT_PROP%POS(K) = NUM_X

          K = K + 1

        ELSE

          EL_CONT_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 4
          EL_CONT_SP%COL(K) = X_POS(NUM_X - 1) + NUM_0D - 3
          EL_CONT_PROP%SIGN(K) = 1
          EL_CONT_PROP%POS(K) = NUM_X

          K = K + 1

        END IF

      END IF

    END IF

  END SUBROUTINE INIT_EL_CONT_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize electron gain term (due to ionization) sparsity pattern
  SUBROUTINE INIT_EL_GAIN_SP

    IMPLICIT NONE

    INTEGER :: I, K, J


    EL_GAIN_SP%N_NZ = NUM_NEUTRALS * LOC_NUM_C

    ALLOCATE(EL_GAIN_SP%ROW(EL_GAIN_SP%N_NZ))
    ALLOCATE(EL_GAIN_SP%COL(EL_GAIN_SP%N_NZ))
    ALLOCATE(EL_GAIN_POS(EL_GAIN_SP%N_NZ))
    ALLOCATE(EL_GAIN_N(EL_GAIN_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

          DO J = 1, NUM_NEUTRALS

            EL_GAIN_SP%ROW(K) = X_POS(I) + NUM_0D - 4
            EL_GAIN_SP%COL(K) = X_POS(I) + NUM_H*NUM_V + J
            EL_GAIN_POS(K) = I
            EL_GAIN_N(K) = J

            K = K + 1

          END DO

      END IF

    END DO

  END SUBROUTINE INIT_EL_GAIN_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize electron loss term (due to recombination) sparsity pattern
  SUBROUTINE INIT_EL_LOSS_SP

    IMPLICIT NONE

    INTEGER :: I, K


    EL_LOSS_SP%N_NZ = LOC_NUM_C

    ALLOCATE(EL_LOSS_SP%ROW(EL_LOSS_SP%N_NZ))
    ALLOCATE(EL_LOSS_SP%COL(EL_LOSS_SP%N_NZ))
    ALLOCATE(EL_RECOMB_POS(EL_LOSS_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

          EL_LOSS_SP%ROW(K) = X_POS(I) + NUM_0D - 4
          EL_LOSS_SP%COL(K) = X_POS(I) + NUM_0D - 4
          EL_RECOMB_POS(K) = I

          K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_EL_LOSS_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize divertor electron convection pattern
  SUBROUTINE INIT_EL_CONV_DIV_SP

    IMPLICIT NONE

      EL_CONV_DIV_SP%N_NZ = 2

      ALLOCATE(EL_CONV_DIV_SP%ROW(EL_CONV_DIV_SP%N_NZ))
      ALLOCATE(EL_CONV_DIV_SP%COL(EL_CONV_DIV_SP%N_NZ))

      EL_CONV_DIV_SP%ROW(1) = X_POS(NUM_X - 1) + NUM_0D - 3
      EL_CONV_DIV_SP%COL(1) = X_POS(NUM_X - 1) + NUM_0D - 3

      EL_CONV_DIV_SP%ROW(2) = X_POS(NUM_X - 1) + NUM_0D - 3
      EL_CONV_DIV_SP%COL(2) = X_POS(NUM_X - 3) + NUM_0D - 3


  END SUBROUTINE INIT_EL_CONV_DIV_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize divertor electron sink pattern
  SUBROUTINE INIT_EL_SINK_SP

    IMPLICIT NONE

      EL_SINK_SP%N_NZ = 3

      ALLOCATE(EL_SINK_SP%ROW(EL_SINK_SP%N_NZ))
      ALLOCATE(EL_SINK_SP%COL(EL_SINK_SP%N_NZ))

      EL_SINK_SP%ROW(1) = X_POS(NUM_X) + NUM_0D - 4
      EL_SINK_SP%COL(1) = X_POS(NUM_X) + NUM_0D - 4

      EL_SINK_SP%ROW(2) = X_POS(NUM_X) + NUM_0D - 4
      EL_SINK_SP%COL(2) = X_POS(NUM_X - 1) + NUM_0D - 3

      EL_SINK_SP%ROW(3) = X_POS(NUM_X) + NUM_0D - 4
      EL_SINK_SP%COL(3) = X_POS(NUM_X - 3) + NUM_0D - 3


  END SUBROUTINE INIT_EL_SINK_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize source term in electron momentum equation
  SUBROUTINE INIT_EL_MOM_SOURCE_SP

    IMPLICIT NONE

    INTEGER :: I,K, RANK, SIZE, IERR, OFFSET


    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    OFFSET = 0

    IF ((RANK .EQ. 0) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1
    IF ((RANK .EQ. SIZE - 1) .AND. FIXED_BOUNDARY_UP_SWITCH) OFFSET = OFFSET + 1

    EL_MOM_SOURCE_SP%N_NZ = nd_loc - (LOC_NUM_C + OFFSET)

    ALLOCATE(EL_MOM_SOURCE_SP%ROW(EL_MOM_SOURCE_SP%N_NZ))
    ALLOCATE(EL_MOM_SOURCE_SP%COL(EL_MOM_SOURCE_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 0) THEN

          EL_MOM_SOURCE_SP%ROW(K) = X_POS(I) + NUM_0D - 3
          EL_MOM_SOURCE_SP%COL(K) = X_POS(I) + NUM_0D - 3

          K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_EL_MOM_SOURCE_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron momentum equation pressure gradient contribution
  SUBROUTINE INIT_EL_PRESSURE_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, N_NZ

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2

      END DO

      EL_PRESSURE_SP%N_NZ = N_NZ

      ALLOCATE(EL_PRESSURE_SP%ROW(N_NZ))
      ALLOCATE(EL_PRESSURE_SP%COL(N_NZ))
      ALLOCATE(EL_PRESSURE_POS(N_NZ))
      ALLOCATE(EL_PRESSURE_SIGN(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

        IF ((I .EQ. NUM_X) .AND. (PERIODIC_BOUNDARY_SWITCH)) THEN

          EL_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D - 3
          EL_PRESSURE_SP%COL(K) = X_POS(I-1) + NUM_0D - 4
          EL_PRESSURE_POS(K) = I - 1
          EL_PRESSURE_SIGN(K) = 1

          K = K + 1

          EL_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D - 3
          EL_PRESSURE_SP%COL(K) = NUM_0D - 4
          EL_PRESSURE_POS(K) = 1
          EL_PRESSURE_SIGN(K) = - 1
          K = K + 1

        ELSE

          IF (MOD(I,2) .EQ. 0) THEN

            EL_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_PRESSURE_SP%COL(K) = X_POS(I-1) + NUM_0D - 4
            EL_PRESSURE_POS(K) = I - 1
            EL_PRESSURE_SIGN(K) = 1

            K = K + 1

            EL_PRESSURE_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_PRESSURE_SP%COL(K) = X_POS(I+1) + NUM_0D - 4
            EL_PRESSURE_POS(K) = I + 1
            EL_PRESSURE_SIGN(K) = - 1

            K = K + 1

          END IF

        END IF

      END DO


  END SUBROUTINE INIT_EL_PRESSURE_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron momentum equation thermal force contribution
  SUBROUTINE INIT_EL_T_DRAG_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, N_NZ

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2

      END DO

      EL_T_DRAG_SP%N_NZ = N_NZ

      ALLOCATE(EL_T_DRAG_SP%ROW(N_NZ))
      ALLOCATE(EL_T_DRAG_SP%COL(N_NZ))
      ALLOCATE(EL_T_DRAG_SIGN(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

        IF ((I .EQ. NUM_X) .AND. (PERIODIC_BOUNDARY_SWITCH)) THEN

          EL_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D - 3
          EL_T_DRAG_SP%COL(K) = X_POS(I-1) + NUM_0D - 2
          EL_T_DRAG_SIGN(K) = 1

          K = K + 1

          EL_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D - 3
          EL_T_DRAG_SP%COL(K) = NUM_0D - 2
          EL_T_DRAG_SIGN(K) = - 1
          K = K + 1

        ELSE

          IF (MOD(I,2) .EQ. 0) THEN

            EL_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_T_DRAG_SP%COL(K) = X_POS(I-1) + NUM_0D - 2
            EL_T_DRAG_SIGN(K) = 1

            K = K + 1

            EL_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_T_DRAG_SP%COL(K) = X_POS(I+1) + NUM_0D - 2
            EL_T_DRAG_SIGN(K) = - 1

            K = K + 1

          END IF

        END IF

      END DO


  END SUBROUTINE INIT_EL_T_DRAG_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid ion momentum equation thermal force contribution
  SUBROUTINE INIT_ION_T_DRAG_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, N_NZ

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2

      END DO

      ION_T_DRAG_SP%N_NZ = N_NZ

      ALLOCATE(ION_T_DRAG_SP%ROW(N_NZ))
      ALLOCATE(ION_T_DRAG_SP%COL(N_NZ))
      ALLOCATE(ION_T_DRAG_SIGN(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

        IF ((I .EQ. NUM_X) .AND. (PERIODIC_BOUNDARY_SWITCH)) THEN

          ION_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
          ION_T_DRAG_SP%COL(K) = X_POS(I-1) + NUM_0D - 2
          ION_T_DRAG_SIGN(K) = 1

          K = K + 1

          ION_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
          ION_T_DRAG_SP%COL(K) = NUM_0D - 2
          ION_T_DRAG_SIGN(K) = - 1
          K = K + 1

        ELSE

          IF (MOD(I,2) .EQ. 0) THEN

            ION_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_T_DRAG_SP%COL(K) = X_POS(I-1) + NUM_0D - 2
            ION_T_DRAG_SIGN(K) = 1

            K = K + 1

            ION_T_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_T_DRAG_SP%COL(K) = X_POS(I+1) + NUM_0D - 2
            ION_T_DRAG_SIGN(K) = - 1

            K = K + 1

          END IF

        END IF

      END DO


  END SUBROUTINE INIT_ION_T_DRAG_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron momentum equation friction contribution
  SUBROUTINE INIT_EL_U_DRAG_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, N_NZ

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2

      END DO

      EL_U_DRAG_SP%N_NZ = N_NZ

      ALLOCATE(EL_U_DRAG_SP%ROW(N_NZ))
      ALLOCATE(EL_U_DRAG_SP%COL(N_NZ))
      ALLOCATE(EL_U_DRAG_SIGN(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 0) THEN

            EL_U_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_U_DRAG_SP%COL(K) = X_POS(I) + NUM_0D
            EL_U_DRAG_SIGN(K) = 1

            K = K + 1

            EL_U_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_U_DRAG_SP%COL(K) = X_POS(I) + NUM_0D - 3
            EL_U_DRAG_SIGN(K) = - 1

            K = K + 1

          END IF

      END DO


  END SUBROUTINE INIT_EL_U_DRAG_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid ion momentum equation friction contribution
  SUBROUTINE INIT_ION_U_DRAG_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, N_NZ

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2

      END DO

      ION_U_DRAG_SP%N_NZ = N_NZ

      ALLOCATE(ION_U_DRAG_SP%ROW(N_NZ))
      ALLOCATE(ION_U_DRAG_SP%COL(N_NZ))
      ALLOCATE(ION_U_DRAG_SIGN(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 0) THEN

            ION_U_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_U_DRAG_SP%COL(K) = X_POS(I) + NUM_0D
            ION_U_DRAG_SIGN(K) = 1

            K = K + 1

            ION_U_DRAG_SP%ROW(K) = X_POS(I) + NUM_0D
            ION_U_DRAG_SP%COL(K) = X_POS(I) + NUM_0D - 3
            ION_U_DRAG_SIGN(K) = - 1

            K = K + 1

          END IF

      END DO


  END SUBROUTINE INIT_ION_U_DRAG_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron temperature convection contribution
  SUBROUTINE INIT_EL_CONV_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, CONV_MIN_X, CONV_MAX_X, NUM_BOUNDARIES_CONV,NNZ

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    NUM_BOUNDARIES_CONV = 0

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      NNZ = 3 * LOC_NUM_C

    ELSE

      IF ((MIN_X .LT. 2) .AND. (.NOT. FIXED_BOUNDARY_UP_SWITCH)) &
      NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      IF ((MAX_X .GT. NUM_X - 1) .AND. (.NOT. FIXED_BOUNDARY_DIV_SWITCH)) &
      NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      NNZ = 3 * LOC_NUM_C - NUM_BOUNDARIES_CONV

    END IF

    EL_CONV_TEMP_SP%N_NZ = NNZ

    CONV_MIN_X = MAX(MIN_X,3)
    CONV_MAX_X = MIN(MAX_X, NUM_X - 2)

    ALLOCATE(EL_CONV_TEMP_SP%ROW(EL_CONV_TEMP_SP%N_NZ))
    ALLOCATE(EL_CONV_TEMP_SP%COL(EL_CONV_TEMP_SP%N_NZ))
    ALLOCATE(EL_CONV_TEMP_PROP%SIGN(EL_CONV_TEMP_SP%N_NZ))
    ALLOCATE(EL_CONV_TEMP_PROP%POS(EL_CONV_TEMP_SP%N_NZ))

    K = 1

    IF (MIN_X .LT. 2) THEN

      IF (.NOT. FIXED_BOUNDARY_UP_SWITCH) THEN

        EL_CONV_TEMP_SP%ROW(K) = NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(3) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = - 1
        EL_CONV_TEMP_PROP%POS(K) = 1

        K = K + 1

        EL_CONV_TEMP_SP%ROW(K) = NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = 0
        EL_CONV_TEMP_PROP%POS(K) = 1

        K = K + 1

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

          EL_CONV_TEMP_SP%ROW(K) = NUM_0D - 2
          EL_CONV_TEMP_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D - 2
          EL_CONV_TEMP_PROP%SIGN(K) = 1
          EL_CONV_TEMP_PROP%POS(K) = 1

          K = K + 1

        END IF

      END IF

    END IF

    DO I = CONV_MIN_X, CONV_MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

        EL_CONV_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(I+2) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = - 1
        EL_CONV_TEMP_PROP%POS(K) = I

        K = K + 1

        EL_CONV_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(I) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = 0
        EL_CONV_TEMP_PROP%POS(K) = I

        K = K + 1

        EL_CONV_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(I-2) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = 1
        EL_CONV_TEMP_PROP%POS(K) = I

        K = K + 1

      END IF

    END DO

    IF (MAX_X .GT. NUM_X - 2) THEN

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        EL_CONV_TEMP_SP%ROW(K) = X_POS(NUM_X - 1) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) =  NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = - 1
        EL_CONV_TEMP_PROP%POS(K) = NUM_X - 1

        K = K + 1

        EL_CONV_TEMP_SP%ROW(K) = X_POS(NUM_X - 1) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(NUM_X - 1) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = 0
        EL_CONV_TEMP_PROP%POS(K) = NUM_X - 1

        K = K + 1

        EL_CONV_TEMP_SP%ROW(K) = X_POS(NUM_X - 1) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(NUM_X - 3) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = 1
        EL_CONV_TEMP_PROP%POS(K) = NUM_X - 1

        K = K + 1

      ELSE IF (.NOT. FIXED_BOUNDARY_DIV_SWITCH) THEN

        EL_CONV_TEMP_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(NUM_X) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = 0
        EL_CONV_TEMP_PROP%POS(K) = NUM_X

        K = K + 1

        EL_CONV_TEMP_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
        EL_CONV_TEMP_SP%COL(K) = X_POS(NUM_X - 2) + NUM_0D - 2
        EL_CONV_TEMP_PROP%SIGN(K) = 1
        EL_CONV_TEMP_PROP%POS(K) = NUM_X

        K = K + 1

      END IF

    END IF

  END SUBROUTINE INIT_EL_CONV_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron pdv temperature equation contribution
  SUBROUTINE INIT_EL_TEMP_PDV_SP

    IMPLICIT NONE

    INTEGER :: I, K,CONV_MIN_X, CONV_MAX_X, NUM_BOUNDARIES_CONV,NNZ

    NUM_BOUNDARIES_CONV = 0

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      NNZ = 2 * LOC_NUM_C

    ELSE

      IF ((MIN_X .EQ. 1) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH))  NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      NNZ = 2 * LOC_NUM_C - NUM_BOUNDARIES_CONV

    END IF

    EL_TEMP_PDV_SP%N_NZ = NNZ

    CONV_MIN_X = MAX(MIN_X,2)
    CONV_MAX_X = MAX_X
    IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) CONV_MAX_X = NUM_X - 2

    ALLOCATE(EL_TEMP_PDV_SP%ROW(EL_TEMP_PDV_SP%N_NZ))
    ALLOCATE(EL_TEMP_PDV_SP%COL(EL_TEMP_PDV_SP%N_NZ))
    ALLOCATE(EL_TEMP_PDV_PROP%SIGN(EL_TEMP_PDV_SP%N_NZ))
    ALLOCATE(EL_TEMP_PDV_PROP%POS(EL_TEMP_PDV_SP%N_NZ))

    K = 1

    IF (MIN_X .EQ. 1) THEN

      EL_TEMP_PDV_SP%ROW(K) = X_POS(1) + NUM_0D - 2
      EL_TEMP_PDV_SP%COL(K) = X_POS(2) + NUM_0D - 3
      EL_TEMP_PDV_PROP%SIGN(K) = - 1
      EL_TEMP_PDV_PROP%POS(K) = 1

      K = K + 1

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        EL_TEMP_PDV_SP%ROW(K) = X_POS(1) + NUM_0D - 2
        EL_TEMP_PDV_SP%COL(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_TEMP_PDV_PROP%SIGN(K) = 1
        EL_TEMP_PDV_PROP%POS(K) = 1

        K = K + 1

      END IF

    END IF

    DO I = CONV_MIN_X, CONV_MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

        EL_TEMP_PDV_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_TEMP_PDV_SP%COL(K) = X_POS(I+1) + NUM_0D - 3
        EL_TEMP_PDV_PROP%SIGN(K) = - 1
        EL_TEMP_PDV_PROP%POS(K) = I

        K = K + 1

        EL_TEMP_PDV_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_TEMP_PDV_SP%COL(K) = X_POS(I-1) + NUM_0D - 3
        EL_TEMP_PDV_PROP%SIGN(K) = 1
        EL_TEMP_PDV_PROP%POS(K) = I

        K = K + 1

      END IF

    END DO

    IF (MAX_X .EQ. NUM_X) THEN

      IF (MOD(NUM_X,2) .EQ. 1) THEN

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

          EL_TEMP_PDV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_PDV_SP%COL(K) = X_POS(1) + NUM_0D - 3
          EL_TEMP_PDV_PROP%SIGN(K) = - 1
          EL_TEMP_PDV_PROP%POS(K) = NUM_X

          K = K + 1

          EL_TEMP_PDV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_PDV_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D - 3
          EL_TEMP_PDV_PROP%SIGN(K) = 1
          EL_TEMP_PDV_PROP%POS(K) = NUM_X

          K = K + 1

        ELSE

          EL_TEMP_PDV_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_PDV_SP%COL(K) = X_POS(NUM_X - 1) + NUM_0D - 3
          EL_TEMP_PDV_PROP%SIGN(K) = 1
          EL_TEMP_PDV_PROP%POS(K) = NUM_X

          K = K + 1

        END IF

      END IF

    END IF

  END SUBROUTINE INIT_EL_TEMP_PDV_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize divertor electron convection pattern
  SUBROUTINE INIT_EL_TEMP_PDV_DIV_SP

    IMPLICIT NONE

      EL_TEMP_PDV_DIV_SP%N_NZ = 1

      ALLOCATE(EL_TEMP_PDV_DIV_SP%ROW(EL_TEMP_PDV_DIV_SP%N_NZ))
      ALLOCATE(EL_TEMP_PDV_DIV_SP%COL(EL_TEMP_PDV_DIV_SP%N_NZ))

      EL_TEMP_PDV_DIV_SP%ROW(1) = X_POS(NUM_X) + NUM_0D - 2
      EL_TEMP_PDV_DIV_SP%COL(1) = X_POS(NUM_X) + NUM_0D - 2


  END SUBROUTINE INIT_EL_TEMP_PDV_DIV_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize div q_u fluid electron temperature equation contribution
  SUBROUTINE INIT_EL_TEMP_Q_U_SP

    IMPLICIT NONE

    INTEGER :: I, K,CONV_MIN_X, CONV_MAX_X, NUM_BOUNDARIES_CONV,NNZ

    NUM_BOUNDARIES_CONV = 0

    IF (PERIODIC_BOUNDARY_SWITCH) THEN

      NNZ = 4 * LOC_NUM_C

    ELSE

      IF ((MIN_X .EQ. 1) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH))  NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) NUM_BOUNDARIES_CONV = NUM_BOUNDARIES_CONV + 1

      NNZ = 4 * LOC_NUM_C - 2 * NUM_BOUNDARIES_CONV

    END IF

    EL_TEMP_Q_U_SP%N_NZ = NNZ

    CONV_MIN_X = MAX(MIN_X,2)
    CONV_MAX_X = MAX_X
    IF ((MAX_X .EQ. NUM_X) .AND. (NO_FLOW_BOUNDARY_DIV_SWITCH)) CONV_MAX_X = NUM_X - 2

    ALLOCATE(EL_TEMP_Q_U_SP%ROW(EL_TEMP_Q_U_SP%N_NZ))
    ALLOCATE(EL_TEMP_Q_U_SP%COL(EL_TEMP_Q_U_SP%N_NZ))
    ALLOCATE(EL_TEMP_Q_U_PROP%SIGN(EL_TEMP_Q_U_SP%N_NZ))
    ALLOCATE(EL_TEMP_Q_U_PROP%POS(EL_TEMP_Q_U_SP%N_NZ))

    K = 1

    IF (MIN_X .EQ. 1) THEN

      EL_TEMP_Q_U_SP%ROW(K) = X_POS(1) + NUM_0D - 2
      EL_TEMP_Q_U_SP%COL(K) = X_POS(2) + NUM_0D - 3
      EL_TEMP_Q_U_PROP%SIGN(K) = - 1
      EL_TEMP_Q_U_PROP%POS(K) = 1

      K = K + 1

      EL_TEMP_Q_U_SP%ROW(K) = X_POS(1) + NUM_0D - 2
      EL_TEMP_Q_U_SP%COL(K) = X_POS(2) + NUM_0D
      EL_TEMP_Q_U_PROP%SIGN(K) = 1
      EL_TEMP_Q_U_PROP%POS(K) = 1

      K = K + 1

      IF (PERIODIC_BOUNDARY_SWITCH) THEN

        EL_TEMP_Q_U_SP%ROW(K) = X_POS(1) + NUM_0D - 2
        EL_TEMP_Q_U_SP%COL(K) = X_POS(NUM_X) + NUM_0D - 3
        EL_TEMP_Q_U_PROP%SIGN(K) = 1
        EL_TEMP_Q_U_PROP%POS(K) = 1

        K = K + 1

        EL_TEMP_Q_U_SP%ROW(K) = X_POS(1) + NUM_0D - 2
        EL_TEMP_Q_U_SP%COL(K) = X_POS(NUM_X) + NUM_0D
        EL_TEMP_Q_U_PROP%SIGN(K) = - 1
        EL_TEMP_Q_U_PROP%POS(K) = 1

        K = K + 1

      END IF

    END IF

    DO I = CONV_MIN_X, CONV_MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

        EL_TEMP_Q_U_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_TEMP_Q_U_SP%COL(K) = X_POS(I+1) + NUM_0D - 3
        EL_TEMP_Q_U_PROP%SIGN(K) = - 1
        EL_TEMP_Q_U_PROP%POS(K) = I

        K = K + 1

        EL_TEMP_Q_U_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_TEMP_Q_U_SP%COL(K) = X_POS(I-1) + NUM_0D - 3
        EL_TEMP_Q_U_PROP%SIGN(K) = 1
        EL_TEMP_Q_U_PROP%POS(K) = I

        K = K + 1

        EL_TEMP_Q_U_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_TEMP_Q_U_SP%COL(K) = X_POS(I+1) + NUM_0D
        EL_TEMP_Q_U_PROP%SIGN(K) = 1
        EL_TEMP_Q_U_PROP%POS(K) = I

        K = K + 1

        EL_TEMP_Q_U_SP%ROW(K) = X_POS(I) + NUM_0D - 2
        EL_TEMP_Q_U_SP%COL(K) = X_POS(I-1) + NUM_0D
        EL_TEMP_Q_U_PROP%SIGN(K) = - 1
        EL_TEMP_Q_U_PROP%POS(K) = I

        K = K + 1

      END IF

    END DO

    IF (MAX_X .EQ. NUM_X) THEN

      IF (MOD(NUM_X,2) .EQ. 1) THEN

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

          EL_TEMP_Q_U_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_Q_U_SP%COL(K) = X_POS(1) + NUM_0D - 3
          EL_TEMP_Q_U_PROP%SIGN(K) = - 1
          EL_TEMP_Q_U_PROP%POS(K) = NUM_X

          K = K + 1

          EL_TEMP_Q_U_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_Q_U_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D - 3
          EL_TEMP_Q_U_PROP%SIGN(K) = 1
          EL_TEMP_Q_U_PROP%POS(K) = NUM_X

          K = K + 1

          EL_TEMP_Q_U_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_Q_U_SP%COL(K) = X_POS(1) + NUM_0D
          EL_TEMP_Q_U_PROP%SIGN(K) = 1
          EL_TEMP_Q_U_PROP%POS(K) = NUM_X

          K = K + 1

          EL_TEMP_Q_U_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_Q_U_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D
          EL_TEMP_Q_U_PROP%SIGN(K) = - 1
          EL_TEMP_Q_U_PROP%POS(K) = NUM_X

          K = K + 1

        ELSE

          EL_TEMP_Q_U_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_Q_U_SP%COL(K) = X_POS(NUM_X - 1) + NUM_0D - 3
          EL_TEMP_Q_U_PROP%SIGN(K) = 1
          EL_TEMP_Q_U_PROP%POS(K) = NUM_X

          K = K + 1

          EL_TEMP_Q_U_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_TEMP_Q_U_SP%COL(K) = X_POS(NUM_X - 1) + NUM_0D
          EL_TEMP_Q_U_PROP%SIGN(K) = - 1
          EL_TEMP_Q_U_PROP%POS(K) = NUM_X

          K = K + 1

        END IF

      END IF

    END IF

  END SUBROUTINE INIT_EL_TEMP_Q_U_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize electron temperature equation gain term (due to ionization) sparsity pattern
  SUBROUTINE INIT_EL_GAIN_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K, J


    EL_GAIN_TEMP_SP%N_NZ = NUM_NEUTRALS * LOC_NUM_C

    ALLOCATE(EL_GAIN_TEMP_SP%ROW(EL_GAIN_TEMP_SP%N_NZ))
    ALLOCATE(EL_GAIN_TEMP_SP%COL(EL_GAIN_TEMP_SP%N_NZ))
    ALLOCATE(EL_GAIN_TEMP_POS(EL_GAIN_TEMP_SP%N_NZ))
    ALLOCATE(EL_GAIN_TEMP_N(EL_GAIN_TEMP_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

          DO J = 1, NUM_NEUTRALS

            EL_GAIN_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
            EL_GAIN_TEMP_SP%COL(K) = X_POS(I) + NUM_H*NUM_V + J
            EL_GAIN_TEMP_POS(K) = I
            EL_GAIN_TEMP_N(K) = J

            K = K + 1

          END DO

      END IF

    END DO

  END SUBROUTINE INIT_EL_GAIN_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize electron temp loss term (due to recombination) sparsity pattern
  SUBROUTINE INIT_EL_LOSS_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K


    EL_LOSS_TEMP_SP%N_NZ = LOC_NUM_C

    ALLOCATE(EL_LOSS_TEMP_SP%ROW(EL_LOSS_TEMP_SP%N_NZ))
    ALLOCATE(EL_LOSS_TEMP_SP%COL(EL_LOSS_TEMP_SP%N_NZ))
    ALLOCATE(EL_RECOMB_TEMP_POS(EL_LOSS_TEMP_SP%N_NZ))

    K = 1

    DO I = MIN_X, MAX_X

      IF (MOD(I,2) .EQ. 1) THEN

          EL_LOSS_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
          EL_LOSS_TEMP_SP%COL(K) = X_POS(I) + NUM_0D - 4
          EL_RECOMB_TEMP_POS(K) = I

          K = K + 1

      END IF

    END DO

  END SUBROUTINE INIT_EL_LOSS_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron temperature equation thermal force contribution
  SUBROUTINE INIT_EL_T_DRAG_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, N_NZ, OFFSET_DRAG_UP, OFFSET_DRAG_DIV

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 1) N_NZ = N_NZ + 2

      END DO

      EL_T_DRAG_TEMP_SP%N_NZ = N_NZ

      ALLOCATE(EL_T_DRAG_TEMP_SP%ROW(N_NZ))
      ALLOCATE(EL_T_DRAG_TEMP_SP%COL(N_NZ))
      ALLOCATE(EL_T_DRAG_TEMP_SIGN(N_NZ))

      OFFSET_DRAG_UP = 0
      OFFSET_DRAG_DIV = 0

      IF ((MIN_X .EQ. 1) .AND. (.NOT. FIXED_BOUNDARY_UP_SWITCH)) OFFSET_DRAG_UP = 1

      IF ((MAX_X .GT. NUM_X - 2).AND. (.NOT. FIXED_BOUNDARY_UP_SWITCH)) OFFSET_DRAG_DIV = 2

      K = 1

      IF (MIN_X .EQ. 1) THEN

        EL_T_DRAG_TEMP_SP%ROW(K) = NUM_0D - 2
        EL_T_DRAG_TEMP_SP%COL(K) = X_POS(3) + NUM_0D - 2
        EL_T_DRAG_TEMP_SIGN(K) = - 1

        K = K + 1

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

          EL_T_DRAG_TEMP_SP%ROW(K) = NUM_0D - 2
          EL_T_DRAG_TEMP_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D - 2
          EL_T_DRAG_TEMP_SIGN(K) = 1

          K = K + 1

        ELSE

          EL_T_DRAG_TEMP_SP%ROW(K) = NUM_0D - 2
          EL_T_DRAG_TEMP_SP%COL(K) = NUM_0D - 2
          EL_T_DRAG_TEMP_SIGN(K) = 1

          K = K + 1

        END IF

      END IF

      DO I = MIN_X+OFFSET_DRAG_UP, MAX_X-OFFSET_DRAG_DIV

          IF (MOD(I,2) .EQ. 1) THEN

            EL_T_DRAG_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
            EL_T_DRAG_TEMP_SP%COL(K) = X_POS(I-2) + NUM_0D - 2
            EL_T_DRAG_TEMP_SIGN(K) = 1

            K = K + 1

            EL_T_DRAG_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
            EL_T_DRAG_TEMP_SP%COL(K) = X_POS(I+2) + NUM_0D - 2
            EL_T_DRAG_TEMP_SIGN(K) = - 1

            K = K + 1

          END IF

      END DO

      IF (MAX_X .GT. NUM_X - 2) THEN

        IF (PERIODIC_BOUNDARY_SWITCH) THEN

          EL_T_DRAG_TEMP_SP%ROW(K) = X_POS(NUM_X-1) + NUM_0D - 2
          EL_T_DRAG_TEMP_SP%COL(K) = NUM_0D - 2
          EL_T_DRAG_TEMP_SIGN(K) = - 1

          K = K + 1

          EL_T_DRAG_TEMP_SP%ROW(K) = X_POS(NUM_X-1) + NUM_0D - 2
          EL_T_DRAG_TEMP_SP%COL(K) = X_POS(NUM_X-3) + NUM_0D - 2
          EL_T_DRAG_TEMP_SIGN(K) = 1

          K = K + 1

        ELSE IF (NO_FLOW_BOUNDARY_DIV_SWITCH) THEN

          EL_T_DRAG_TEMP_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_T_DRAG_TEMP_SP%COL(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_T_DRAG_TEMP_SIGN(K) = - 1

          K = K + 1

          EL_T_DRAG_TEMP_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 2
          EL_T_DRAG_TEMP_SP%COL(K) = X_POS(NUM_X-2) + NUM_0D - 2
          EL_T_DRAG_TEMP_SIGN(K) = 1

          K = K + 1

        END IF

      END IF

  END SUBROUTINE INIT_EL_T_DRAG_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron temperature equation e-i friction contribution
  SUBROUTINE INIT_EL_U_DRAG_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K, RANK, SIZE, IERR, N_NZ

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD,SIZE,IERR)

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 1) N_NZ = N_NZ + 1

      END DO

      EL_U_DRAG_TEMP_SP%N_NZ = N_NZ

      ALLOCATE(EL_U_DRAG_TEMP_SP%ROW(N_NZ))
      ALLOCATE(EL_U_DRAG_TEMP_SP%COL(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 1) THEN

            EL_U_DRAG_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
            EL_U_DRAG_TEMP_SP%COL(K) = X_POS(I) + NUM_0D - 4

            K = K + 1

          END IF

      END DO


  END SUBROUTINE INIT_EL_U_DRAG_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron e-n drag (excluding recombination) sparsity pattern
  SUBROUTINE INIT_EL_EN_DRAG_N_SP

    IMPLICIT NONE

    INTEGER :: I, K, N_NZ , OFFSET, J

    REAL(KIND=DEFAULT_REAL) :: a_minus, a_plus

    N_NZ = 0

    OFFSET = 0

    IF ((PERIODIC_BOUNDARY_SWITCH) .AND. (MAX_X .EQ. NUM_X)) OFFSET = 1

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2*NUM_NEUTRALS

      END DO

      EL_EN_DRAG_N_SP%N_NZ = N_NZ

      ALLOCATE(EL_EN_DRAG_N_SP%ROW(N_NZ))
      ALLOCATE(EL_EN_DRAG_N_SP%COL(N_NZ))
      ALLOCATE(EL_EN_DRAG_N(N_NZ))
      ALLOCATE(EL_EN_DRAG_N_POS(N_NZ))
      ALLOCATE(EL_EN_DRAG_N_MULT(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X - OFFSET

          IF (MOD(I,2) .EQ. 0) THEN

            a_minus = (X_GRID(I + 1) - X_GRID(I)) / dxc(I)
            a_plus = (X_GRID(I) - X_GRID(I - 1)) / dxc(I)

            DO J = 1, NUM_NEUTRALS

              EL_EN_DRAG_N_SP%ROW(K) = X_POS(I) + NUM_0D - 3
              EL_EN_DRAG_N_SP%COL(K) = X_POS(I-1) + NUM_V*NUM_H + J
              EL_EN_DRAG_N(K) = J
              EL_EN_DRAG_N_POS(K) = I
              EL_EN_DRAG_N_MULT(K) = a_minus

              K = K + 1

              EL_EN_DRAG_N_SP%ROW(K) = X_POS(I) + NUM_0D - 3
              EL_EN_DRAG_N_SP%COL(K) = X_POS(I+1) + NUM_V*NUM_H + J
              EL_EN_DRAG_N(K) = J
              EL_EN_DRAG_N_POS(K) = I
              EL_EN_DRAG_N_MULT(K) = a_plus

              K = K + 1

            END DO

          END IF

      END DO

      IF (OFFSET .EQ. 1) THEN

        DO J = 1, NUM_NEUTRALS

          EL_EN_DRAG_N_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 3
          EL_EN_DRAG_N_SP%COL(K) = X_POS(NUM_X-1) + NUM_V*NUM_H + J
          EL_EN_DRAG_N(K) = J
          EL_EN_DRAG_N_POS(K) = NUM_X
          EL_EN_DRAG_N_MULT(K) = 0.5D00

          K = K + 1

          EL_EN_DRAG_N_SP%ROW(K) = X_POS(I) + NUM_0D - 3
          EL_EN_DRAG_N_SP%COL(K) = NUM_V*NUM_H + J
          EL_EN_DRAG_N(K) = J
          EL_EN_DRAG_N_POS(K) = NUM_X
          EL_EN_DRAG_N_MULT(K) = 0.5D00

          K = K + 1

        END DO

      END IF

  END SUBROUTINE INIT_EL_EN_DRAG_N_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron recombination drag sparsity pattern
  SUBROUTINE INIT_EL_EN_DRAG_RECOMB_SP

    IMPLICIT NONE

    INTEGER :: I, K, N_NZ , OFFSET

    REAL(KIND=DEFAULT_REAL) :: a_minus, a_plus

    N_NZ = 0

    OFFSET = 0

    IF ((PERIODIC_BOUNDARY_SWITCH) .AND. (MAX_X .EQ. NUM_X)) OFFSET = 1

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 0) N_NZ = N_NZ + 2

      END DO

      EL_EN_DRAG_RECOMB_SP%N_NZ = N_NZ

      ALLOCATE(EL_EN_DRAG_RECOMB_SP%ROW(N_NZ))
      ALLOCATE(EL_EN_DRAG_RECOMB_SP%COL(N_NZ))
      ALLOCATE(EL_EN_DRAG_REC_POS(N_NZ))
      ALLOCATE(EL_EN_DRAG_REC_MULT(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X - OFFSET

          IF (MOD(I,2) .EQ. 0) THEN

            a_minus = (X_GRID(I + 1) - X_GRID(I)) / dxc(I)
            a_plus = (X_GRID(I) - X_GRID(I - 1)) / dxc(I)

            EL_EN_DRAG_RECOMB_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_EN_DRAG_RECOMB_SP%COL(K) = X_POS(I-1) + NUM_0D - 4
            EL_EN_DRAG_REC_POS(K) = I
            EL_EN_DRAG_REC_MULT(K) = a_minus

            K = K + 1

            EL_EN_DRAG_RECOMB_SP%ROW(K) = X_POS(I) + NUM_0D - 3
            EL_EN_DRAG_RECOMB_SP%COL(K) = X_POS(I+1) + NUM_0D - 4
            EL_EN_DRAG_REC_POS(K) = I
            EL_EN_DRAG_REC_MULT(K) = a_plus

            K = K + 1

          END IF

      END DO

      IF (OFFSET .EQ. 1) THEN

          EL_EN_DRAG_RECOMB_SP%ROW(K) = X_POS(NUM_X) + NUM_0D - 3
          EL_EN_DRAG_RECOMB_SP%COL(K) = X_POS(NUM_X-1) + NUM_0D - 4
          EL_EN_DRAG_REC_POS(K) = NUM_X
          EL_EN_DRAG_REC_MULT(K) = 0.5D00

          K = K + 1

          EL_EN_DRAG_RECOMB_SP%ROW(K) = X_POS(I) + NUM_0D - 3
          EL_EN_DRAG_RECOMB_SP%COL(K) = NUM_V*NUM_H + NUM_0D - 4
          EL_EN_DRAG_REC_POS(K) = NUM_X
          EL_EN_DRAG_REC_MULT(K) = 0.5D00

          K = K + 1

      END IF

  END SUBROUTINE INIT_EL_EN_DRAG_RECOMB_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron temperature equation drag term sparsity pattern
  SUBROUTINE INIT_EL_EN_DRAG_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K, N_NZ , J

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 1) THEN

          IF (COLL_RECOMB) N_NZ = N_NZ + 1
          IF (COLL_EN_EL_L_SWITCH .OR. COLL_EN_EX .OR. COLL_EN_ION) N_NZ = N_NZ + NUM_NEUTRALS

        END IF

      END DO

      EL_EN_DRAG_TEMP_SP%N_NZ = N_NZ

      ALLOCATE(EL_EN_DRAG_TEMP_SP%ROW(N_NZ))
      ALLOCATE(EL_EN_DRAG_TEMP_SP%COL(N_NZ))

      K = 1

      IF (COLL_EN_EL_L_SWITCH .OR. COLL_EN_EX .OR. COLL_EN_ION) THEN

        DO I = MIN_X, MAX_X

            IF (MOD(I,2) .EQ. 1) THEN

              DO J = 1, NUM_NEUTRALS

                EL_EN_DRAG_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
                EL_EN_DRAG_TEMP_SP%COL(K) = X_POS(I) + NUM_V*NUM_H + J

                K = K + 1

              END DO

            END IF

        END DO

      END IF

      IF (COLL_RECOMB) THEN

        DO I = MIN_X, MAX_X

            IF (MOD(I,2) .EQ. 1) THEN

              EL_EN_DRAG_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
              EL_EN_DRAG_TEMP_SP%COL(K) = X_POS(I) + NUM_0D - 4

              K = K + 1

            END IF

        END DO

      END IF

  END SUBROUTINE INIT_EL_EN_DRAG_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron electron energy loss/gain term due to e-n inelastic collisions SP
  SUBROUTINE INIT_EL_EN_INCOLL_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K, N_NZ , J

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 1) THEN

           N_NZ = N_NZ + NUM_NEUTRALS + 1

        END IF

      END DO

      EL_EN_INCOLL_TEMP_SP%N_NZ = N_NZ

      ALLOCATE(EL_EN_INCOLL_TEMP_SP%ROW(N_NZ))
      ALLOCATE(EL_EN_INCOLL_TEMP_SP%COL(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 1) THEN

            DO J = 1, NUM_NEUTRALS

              EL_EN_INCOLL_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
              EL_EN_INCOLL_TEMP_SP%COL(K) = X_POS(I) + NUM_V*NUM_H + J

              K = K + 1

            END DO

            EL_EN_INCOLL_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
            EL_EN_INCOLL_TEMP_SP%COL(K) = X_POS(I) + NUM_0D - 4

            K = K + 1

          END IF

      END DO


  END SUBROUTINE INIT_EL_EN_INCOLL_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron electron energy loss/gain term due to e-n elastic collisions SP
  SUBROUTINE INIT_EL_EN_ELCOLL_TEMP_SP

    IMPLICIT NONE

    INTEGER :: I, K, N_NZ , J

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 1) THEN

           N_NZ = N_NZ + NUM_NEUTRALS

        END IF

      END DO

      EL_EN_ELCOLL_TEMP_SP%N_NZ = N_NZ

      ALLOCATE(EL_EN_ELCOLL_TEMP_SP%ROW(N_NZ))
      ALLOCATE(EL_EN_ELCOLL_TEMP_SP%COL(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 1) THEN

            DO J = 1, NUM_NEUTRALS

              EL_EN_ELCOLL_TEMP_SP%ROW(K) = X_POS(I) + NUM_0D - 2
              EL_EN_ELCOLL_TEMP_SP%COL(K) = X_POS(I) + NUM_V*NUM_H + J

              K = K + 1

            END DO

          END IF

      END DO


  END SUBROUTINE INIT_EL_EN_ELCOLL_TEMP_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize full fluid mode CRM recombination sparsity pattern
  SUBROUTINE INIT_CRM_RECOMB_FF_SP

    IMPLICIT NONE

    INTEGER :: I, K, N_NZ , J

    N_NZ = 0

      DO I = MIN_X, MAX_X

        IF (MOD(I,2) .EQ. 1) THEN

           N_NZ = N_NZ + NUM_NEUTRALS

        END IF

      END DO

      CRM_RECOMB_FF_SP%N_NZ = N_NZ

      ALLOCATE(CRM_RECOMB_FF_SP%ROW(N_NZ))
      ALLOCATE(CRM_RECOMB_FF_SP%COL(N_NZ))

      K = 1

      DO I = MIN_X, MAX_X

          IF (MOD(I,2) .EQ. 1) THEN

            DO J = 1, NUM_NEUTRALS

              CRM_RECOMB_FF_SP%ROW(K) = X_POS(I) + NUM_V*NUM_H + J
              CRM_RECOMB_FF_SP%COL(K) = X_POS(I) + NUM_0D - 4

              K = K + 1

            END DO

          END IF

      END DO


  END SUBROUTINE INIT_CRM_RECOMB_FF_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize fluid electron electron energy loss/gain term due to e-n elastic collisions SP
  SUBROUTINE INIT_EL_HEATING_SP

    IMPLICIT NONE

    INTEGER :: I, K, N_NZ

    N_NZ = 0

      DO I = MIN_X, MIN(MAX_X,2*N_HEATING - 1)

        IF (MOD(I,2) .EQ. 1) THEN

           N_NZ = N_NZ + 1

        END IF

      END DO

      EL_HEATING_SP%N_NZ = N_NZ

      ALLOCATE(EL_HEATING_SP%ROW(N_NZ))
      ALLOCATE(EL_HEATING_SP%COL(N_NZ))

      K = 1

      DO I = MIN_X, MIN(MAX_X,2*N_HEATING - 1)

          IF (MOD(I,2) .EQ. 1) THEN

              EL_HEATING_SP%ROW(K) = X_POS(I) + NUM_0D - 2
              EL_HEATING_SP%COL(K) = X_POS(I) + NUM_0D - 2

              K = K + 1

          END IF

      END DO


  END SUBROUTINE INIT_EL_HEATING_SP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize all inelastic process markers
  SUBROUTINE INIT_INEL_MARKER

    IMPLICIT NONE

    INTEGER :: P, H, L, K, I, J, OFFSET

    ALLOCATE(INEL_M(MIN_X:MAX_X,NUM_H,0:NUM_NEUTRALS,0:NUM_NEUTRALS))
!Allocate nonzero lists based on the data from inelastic mapping
    DO P = MIN_X, MAX_X

      DO H = 1, NUM_H

        DO K = 1, NUM_NEUTRALS

          DO L = K + 1, NUM_NEUTRALS

            ALLOCATE(INEL_M(P,H,K,L)%LIST(INEL_NNZ(K,L)))

            INEL_M(P,H,K,L)%LIST = 0

          END DO

          DO L = 1, K - 1

            ALLOCATE(INEL_M(P,H,K,L)%LIST(INEL_NNZ(K,L)))

            INEL_M(P,H,K,L)%LIST = 0

          END DO

          ALLOCATE(INEL_M(P,H,K,0)%LIST(INEL_NNZ(K,0)))

          INEL_M(P,H,K,0)%LIST = 0

          ALLOCATE(INEL_M(P,H,0,K)%LIST(INEL_NNZ(0,K)))

          INEL_M(P,H,0,K)%LIST = 0

        END DO

      END DO

    END DO
!Set markers
    DO P = MIN_X, MAX_X

      DO H = 1, NUM_H
        !Markers for odd harmonics on boundaries, even in centres
        IF (((MOD(H - 1,2) .EQ. 0) .AND. (MOD(P,2) .EQ. 1)) .OR. ((MOD(H - 1,2) .EQ. 1) .AND. (MOD(P,2) .EQ. 0))) THEN

          OFFSET = X_POS(P) + NUM_V * (NUM_H - H)

          DO K = 1, NUM_NEUTRALS
            !Excitation
            DO L = K + 1, NUM_NEUTRALS

              DO I = 1, INEL_NNZ(K,L)

                DO J = FL_INDEX + 1, LOCAL_SP%N_NZ

                   IF ((INEL_SP(K,L)%ROW(I) + OFFSET .EQ. LOCAL_SP%ROW(J)) .AND. &
                    (INEL_SP(K,L)%COL(I) + OFFSET .EQ. LOCAL_SP%COL(J))) THEN

                     INEL_M(P,H,K,L)%LIST(I) = J
                     EXIT

                   END IF

                 END DO

               END DO

             END DO
             !De-excitation
             DO L = 1, K - 1

               DO I = 1, INEL_NNZ(K,L)

                 DO J = FL_INDEX + 1, LOCAL_SP%N_NZ

                    IF ((INEL_SP(K,L)%ROW(I) + OFFSET .EQ. LOCAL_SP%ROW(J)) .AND. &
                     (INEL_SP(K,L)%COL(I) + OFFSET .EQ. LOCAL_SP%COL(J))) THEN

                      INEL_M(P,H,K,L)%LIST(I) = J
                      EXIT

                    END IF

                  END DO

                END DO

              END DO
              !Ionization
             DO I = 1, INEL_NNZ(K,0)

               DO J = FL_INDEX + 1, LOCAL_SP%N_NZ

                  IF ((INEL_SP(K,0)%ROW(I) + OFFSET .EQ. LOCAL_SP%ROW(J)) .AND. &
                   (INEL_SP(K,0)%COL(I) + OFFSET .EQ. LOCAL_SP%COL(J))) THEN

                    INEL_M(P,H,K,0)%LIST(I) = J

                    EXIT

                  END IF

                END DO

              END DO
              !Recombination
              DO I = 1, INEL_NNZ(0,K)

                DO J = FL_INDEX + 1, LOCAL_SP%N_NZ

                   IF ((INEL_SP(0,K)%ROW(I) + OFFSET .EQ. LOCAL_SP%ROW(J)) .AND. &
                    (INEL_SP(0,K)%COL(I) + OFFSET .EQ. LOCAL_SP%COL(J))) THEN

                     INEL_M(P,H,0,K)%LIST(I) = J

                     EXIT

                   END IF

                 END DO

               END DO

           END DO

         END IF

       END DO

     END DO

   END SUBROUTINE INIT_INEL_MARKER
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize CRM recombination markers
  SUBROUTINE INIT_CRM_RECOMB_MARKER

    IMPLICIT NONE

    INTEGER :: P,I,J,OFFSET

    ALLOCATE(MARKER_RECOMB_CRM(MIN_X:MAX_X,CRM_RECOMB_SP%N_NZ))
    MARKER_RECOMB_CRM = 0

    DO P = MIN_X, MAX_X

      OFFSET = X_POS(P) + NUM_V * (NUM_H - 1)

      IF (MOD(P,2) .EQ. 1) THEN

        DO I = 1, CRM_RECOMB_SP%N_NZ

         DO J = FL_INDEX + 1, LOCAL_SP%N_NZ

           IF ((CRM_RECOMB_SP%ROW(I) + OFFSET .EQ. LOCAL_SP%ROW(J)) .AND. &
              (CRM_RECOMB_SP%COL(I) + OFFSET .EQ. LOCAL_SP%COL(J))) THEN

            MARKER_RECOMB_CRM(P, I) = J

            EXIT

           END IF

         END DO

       END DO

      END IF

    END DO

  END SUBROUTINE INIT_CRM_RECOMB_MARKER
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize secondary electron generation markers
  SUBROUTINE INIT_SEC_EL_MARKER

    IMPLICIT NONE

    INTEGER :: P, I, J,OFFSET

    ALLOCATE(MARKER_SEC_EL(MIN_X:MAX_X,NUM_NEUTRALS))
    MARKER_SEC_EL = 0

    DO P = MIN_X, MAX_X

      IF (MOD(P,2) .EQ. 1) THEN

        OFFSET = X_POS(P) + NUM_V * (NUM_H - 1)

        DO I = 1, NUM_NEUTRALS

          DO J = FL_INDEX + 1, LOCAL_SP%N_NZ

            IF ((EN_ION_SP%ROW(I) + OFFSET .EQ. LOCAL_SP%ROW(J)) .AND. &
             (EN_ION_SP%COL(I) + OFFSET .EQ. LOCAL_SP%COL(J))) THEN

             MARKER_SEC_EL(P,I) = J

             EXIT

            END IF

          END DO

        END DO

      END IF

    END DO

 END SUBROUTINE INIT_SEC_EL_MARKER
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize e-i collision terms in cold ion case
  SUBROUTINE INIT_CI_EI_OD_MARKER(MARKER,OFFSET_ROW,OFFSET_COL)

    IMPLICIT NONE

    INTEGER :: I, J

    INTEGER, DIMENSION(:), INTENT(INOUT) :: MARKER

    INTEGER, INTENT(IN) :: OFFSET_ROW, OFFSET_COL

    DO I = 1, CI_EI_OD_SP%N_NZ

      DO J = 1, LOCAL_SP%N_NZ

        IF ((CI_EI_OD_SP%ROW(I) + OFFSET_ROW .EQ. LOCAL_SP%ROW(J)) .AND. &
         (CI_EI_OD_SP%COL(I) + OFFSET_COL .EQ. LOCAL_SP%COL(J))) THEN

         MARKER(I) = J

         EXIT

        END IF

      END DO

    END DO

  END SUBROUTINE INIT_CI_EI_OD_MARKER
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize simple marker (just using sparsity pattern)
  SUBROUTINE INIT_SIMPLE_MARKER(PATTERN,MARKER,START_INDEX)

    IMPLICIT NONE

    TYPE(SPARSITY_PATTERN), INTENT(IN) :: PATTERN !< Sparsity patter to be used for marker

    INTEGER, DIMENSION(:), INTENT(INOUT), ALLOCATABLE :: MARKER !< Marker to fill

    INTEGER, INTENT(IN) :: START_INDEX !< Number of local matrix strip nonzeroes to skip

    INTEGER :: I, J

    ALLOCATE(MARKER(PATTERN%N_NZ))

    DO I = 1, PATTERN%N_NZ

      DO J = START_INDEX + 1, LOCAL_SP%N_NZ

        IF ((PATTERN%ROW(I) .EQ. LOCAL_SP%ROW(J)) .AND. &
         (PATTERN%COL(I) .EQ. LOCAL_SP%COL(J))) THEN

         MARKER(I) = J

         EXIT

        END IF

      END DO

    END DO

  END SUBROUTINE INIT_SIMPLE_MARKER
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Add all fluid ion terms on cell boundaries and initialize markers
  SUBROUTINE ADD_FLUID_IONS_BOUNDARY

    IMPLICIT NONE

    IF (MAXWELL_SWITCH) THEN

      CALL INIT_MAXWELL_ION_SP
      CALL ADD_SP(LOCAL_SP, MAXWELL_ION_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(MAXWELL_ION_SP,MARKER_MAXWELL_ION,FLUID_ION_INDEX)

    END IF

    IF (E_ADV_SWITCH) THEN

      CALL INIT_ION_LORENTZ_SP
      CALL ADD_SP(LOCAL_SP, ION_LORENTZ_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(ION_LORENTZ_SP,MARKER_ION_LORENTZ,FLUID_ION_INDEX)

    END IF

    IF (X_ADV_SWITCH) THEN

      CALL INIT_ION_CONV_SP
      CALL ADD_SP(LOCAL_SP, ION_CONV_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(ION_CONV_SP,MARKER_ION_CONV,FLUID_ION_INDEX)

      IF (ION_EL_TEMP_SWITCH) THEN

        CALL INIT_ION_PRESSURE_SP
        CALL ADD_SP(LOCAL_SP, ION_PRESSURE_SP,0,0)
        CALL INIT_SIMPLE_MARKER(ION_PRESSURE_SP,MARKER_ION_PRESSURE,FLUID_ION_INDEX)

      END IF

    END IF

    IF (COLL_EI_L_SWITCH) THEN

      IF (FULL_FLUID_MODE) THEN

        CALL INIT_ION_T_DRAG_SP
        CALL ADD_SP(LOCAL_SP, ION_T_DRAG_SP,0,0)
        CALL INIT_SIMPLE_MARKER(ION_T_DRAG_SP,MARKER_ION_T_DRAG,FLUID_ION_INDEX)

        CALL INIT_ION_U_DRAG_SP
        CALL ADD_SP(LOCAL_SP, ION_U_DRAG_SP,0,0)
        CALL INIT_SIMPLE_MARKER(ION_U_DRAG_SP,MARKER_ION_U_DRAG,FLUID_ION_INDEX)

      ELSE

        CALL INIT_ION_DRAG_SP
        CALL ADD_SP(LOCAL_SP, ION_DRAG_SP,0,0)
        CALL INIT_SIMPLE_MARKER(ION_DRAG_SP,MARKER_ION_DRAG,FLUID_ION_INDEX)

      END IF

    END IF

    IF ((MIN_X .LE. NUM_X - 1) .AND. (MAX_X .GE. NUM_X - 1) .AND. (PLASMA_SINK_SWITCH)) THEN

      CALL INIT_ION_CONV_DIV_SP
      CALL ADD_SP(LOCAL_SP, ION_CONV_DIV_SP,0,0)
      CALL INIT_SIMPLE_MARKER(ION_CONV_DIV_SP,MARKER_ION_CONV_DIV,FLUID_ION_INDEX)

    END IF

    IF (COLL_RECOMB .OR. COLL_EN_ION) THEN

      CALL INIT_ION_MOM_SOURCE_SP
      CALL ADD_SP(LOCAL_SP, ION_MOM_SOURCE_SP,0,0)
      CALL INIT_SIMPLE_MARKER(ION_MOM_SOURCE_SP,MARKER_ION_MOM_SOURCE,FLUID_ION_INDEX)

    END IF

    IF (SIMPLE_CX_SWITCH) THEN

      CALL INIT_SIMPLE_CX_SP
      CALL ADD_SP(LOCAL_SP, SIMPLE_CX_SP, 0, 0)
      CALL INIT_SIMPLE_MARKER(SIMPLE_CX_SP,MARKER_SIMPLE_CX,FLUID_ION_INDEX)

    END IF

  END SUBROUTINE ADD_FLUID_IONS_BOUNDARY
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Add all fluid ion terms in cell centres and initialize markers
  SUBROUTINE ADD_FLUID_IONS_CENTRE

    IMPLICIT NONE

    IF (X_ADV_SWITCH) THEN

      CALL INIT_ION_CONT_SP
      CALL ADD_SP(LOCAL_SP, ION_CONT_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(ION_CONT_SP,MARKER_ION_CONT,FLUID_ION_INDEX)

    END IF

    IF (COLL_EN_ION) THEN

      CALL INIT_ION_GAIN_SP
      CALL ADD_SP(LOCAL_SP, ION_GAIN_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(ION_GAIN_SP,MARKER_ION_GAIN,FLUID_ION_INDEX)

    END IF

    IF (COLL_RECOMB) THEN

      IF (FULL_FLUID_MODE) THEN

        CALL INIT_ION_LOSS_FF_SP
        CALL ADD_SP(LOCAL_SP, ION_LOSS_FF_SP,0,0)
        CALL INIT_SIMPLE_MARKER(ION_LOSS_FF_SP,MARKER_ION_LOSS_FF,FLUID_ION_INDEX)

      ELSE

        CALL INIT_ION_LOSS_SP
        CALL ADD_SP(LOCAL_SP, ION_LOSS_SP,0,0,.FALSE.)
        CALL INIT_SIMPLE_MARKER(ION_LOSS_SP,MARKER_ION_LOSS,FLUID_ION_INDEX)

      END IF

    END IF

    IF (PLASMA_SINK_ON) THEN

      CALL INIT_ION_SINK_SP
      CALL ADD_SP(LOCAL_SP, ION_SINK_SP,0,0)
      CALL INIT_SIMPLE_MARKER(ION_SINK_SP,MARKER_ION_SINK,FLUID_ION_INDEX)

    END IF

  END SUBROUTINE ADD_FLUID_IONS_CENTRE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Add all fluid electron terms on cell boundaries and initialize markers
  SUBROUTINE ADD_FLUID_EL_BOUNDARY

    IMPLICIT NONE

    IF (MAXWELL_SWITCH) THEN

      CALL INIT_MAXWELL_EL_SP
      CALL ADD_SP(LOCAL_SP, MAXWELL_EL_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(MAXWELL_EL_SP,MARKER_MAXWELL_EL,FLUID_ION_INDEX)

    END IF

    IF (E_ADV_SWITCH) THEN

      CALL INIT_EL_LORENTZ_SP
      CALL ADD_SP(LOCAL_SP, EL_LORENTZ_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(EL_LORENTZ_SP,MARKER_EL_LORENTZ,FLUID_ION_INDEX)

    END IF

    IF (X_ADV_SWITCH) THEN

      CALL INIT_EL_CONV_SP
      CALL ADD_SP(LOCAL_SP, EL_CONV_SP,0,0,.FALSE.)
      CALL INIT_SIMPLE_MARKER(EL_CONV_SP,MARKER_EL_CONV,FLUID_ION_INDEX)


      CALL INIT_EL_PRESSURE_SP
      CALL ADD_SP(LOCAL_SP, EL_PRESSURE_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_PRESSURE_SP,MARKER_EL_PRESSURE,FLUID_ION_INDEX)

    END IF

    IF (COLL_EI_L_SWITCH) THEN

        CALL INIT_EL_T_DRAG_SP
        CALL ADD_SP(LOCAL_SP, EL_T_DRAG_SP,0,0)
        CALL INIT_SIMPLE_MARKER(EL_T_DRAG_SP,MARKER_EL_T_DRAG,FLUID_ION_INDEX)

        CALL INIT_EL_U_DRAG_SP
        CALL ADD_SP(LOCAL_SP, EL_U_DRAG_SP,0,0)
        CALL INIT_SIMPLE_MARKER(EL_U_DRAG_SP,MARKER_EL_U_DRAG,FLUID_ION_INDEX)


    END IF

    IF ((MIN_X .LE. NUM_X - 1) .AND. (MAX_X .GE. NUM_X - 1) .AND. (PLASMA_SINK_SWITCH)) THEN

      CALL INIT_EL_CONV_DIV_SP
      CALL ADD_SP(LOCAL_SP, EL_CONV_DIV_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_CONV_DIV_SP,MARKER_EL_CONV_DIV,FLUID_ION_INDEX)

    END IF

    IF (COLL_RECOMB .OR. COLL_EN_ION) THEN

      CALL INIT_EL_MOM_SOURCE_SP
      CALL ADD_SP(LOCAL_SP, EL_MOM_SOURCE_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_MOM_SOURCE_SP,MARKER_EL_MOM_SOURCE,FLUID_ION_INDEX)

    END IF

    IF (COLL_EN_EX .OR. COLL_EN_ION .OR. COLL_EN_EL_L_SWITCH) THEN

      CALL INIT_EL_EN_DRAG_N_SP
      CALL ADD_SP(LOCAL_SP, EL_EN_DRAG_N_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_EN_DRAG_N_SP,MARKER_EL_EN_DRAG,FLUID_ION_INDEX)

    END IF

    IF (COLL_RECOMB) THEN

      CALL INIT_EL_EN_DRAG_RECOMB_SP
      CALL ADD_SP(LOCAL_SP, EL_EN_DRAG_RECOMB_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_EN_DRAG_RECOMB_SP,MARKER_EL_EN_DRAG_RECOMB,FLUID_ION_INDEX)

    END IF

  END SUBROUTINE ADD_FLUID_EL_BOUNDARY
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Add all fluid electron terms in cell centres and initialize markers
  SUBROUTINE ADD_FLUID_EL_CENTRE

    IMPLICIT NONE

    IF (X_ADV_SWITCH) THEN

      CALL INIT_EL_CONT_SP
      CALL ADD_SP(LOCAL_SP, EL_CONT_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_CONT_SP,MARKER_EL_CONT,FLUID_ION_INDEX)

      CALL INIT_EL_CONV_TEMP_SP
      CALL ADD_SP(LOCAL_SP, EL_CONV_TEMP_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_CONV_TEMP_SP,MARKER_EL_CONV_TEMP,FLUID_ION_INDEX)

      CALL INIT_EL_TEMP_PDV_SP
      CALL ADD_SP(LOCAL_SP, EL_TEMP_PDV_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_TEMP_PDV_SP,MARKER_EL_TEMP_PDV,FLUID_ION_INDEX)

      IF (COLL_EI_L_SWITCH) THEN

        CALL INIT_EL_TEMP_Q_U_SP
        CALL ADD_SP(LOCAL_SP, EL_TEMP_Q_U_SP,0,0)
        CALL INIT_SIMPLE_MARKER(EL_TEMP_Q_U_SP,MARKER_EL_TEMP_Q_U,FLUID_ION_INDEX)

        CALL INIT_EL_T_DRAG_TEMP_SP
        CALL ADD_SP(LOCAL_SP, EL_T_DRAG_TEMP_SP,0,0)
        CALL INIT_SIMPLE_MARKER(EL_T_DRAG_TEMP_SP,MARKER_EL_T_DRAG_TEMP,FLUID_ION_INDEX)

        CALL INIT_EL_U_DRAG_TEMP_SP
        CALL ADD_SP(LOCAL_SP, EL_U_DRAG_SP,0,0)
        CALL INIT_SIMPLE_MARKER(EL_U_DRAG_SP,MARKER_EL_U_DRAG_TEMP,FLUID_ION_INDEX)

      END IF

    END IF

    IF (COLL_EN_ION) THEN

      CALL INIT_EL_GAIN_SP
      CALL ADD_SP(LOCAL_SP, EL_GAIN_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_GAIN_SP,MARKER_EL_GAIN,FLUID_ION_INDEX)

      CALL INIT_EL_GAIN_TEMP_SP
      CALL ADD_SP(LOCAL_SP, EL_GAIN_TEMP_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_GAIN_TEMP_SP,MARKER_EL_GAIN_TEMP,FLUID_ION_INDEX)

    END IF

    IF (COLL_RECOMB) THEN

      CALL INIT_EL_LOSS_SP
      CALL ADD_SP(LOCAL_SP, EL_LOSS_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_LOSS_SP,MARKER_EL_LOSS,FLUID_ION_INDEX)

      CALL INIT_EL_LOSS_TEMP_SP
      CALL ADD_SP(LOCAL_SP, EL_LOSS_TEMP_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_LOSS_TEMP_SP,MARKER_EL_LOSS_TEMP,FLUID_ION_INDEX)

      CALL INIT_CRM_RECOMB_FF_SP
      CALL ADD_SP(LOCAL_SP, CRM_RECOMB_FF_SP,0,0)
      CALL INIT_SIMPLE_MARKER(CRM_RECOMB_FF_SP,MARKER_RECOMB_CRM_FF,FL_INDEX)

    END IF

    IF (PLASMA_SINK_ON) THEN

      CALL INIT_EL_SINK_SP
      CALL ADD_SP(LOCAL_SP, EL_SINK_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_SINK_SP,MARKER_EL_SINK,FLUID_ION_INDEX)

      IF (X_ADV_SWITCH) THEN

        CALL INIT_EL_TEMP_PDV_DIV_SP
        CALL ADD_SP(LOCAL_SP, EL_TEMP_PDV_DIV_SP,0,0)
        CALL INIT_SIMPLE_MARKER(EL_TEMP_PDV_DIV_SP,MARKER_EL_TEMP_PDV_DIV,FLUID_ION_INDEX)

      END IF

    END IF

    IF (COLL_EN_EX .OR. COLL_EN_ION .OR. COLL_EN_EL_L_SWITCH .OR. COLL_RECOMB) THEN

      CALL INIT_EL_EN_DRAG_TEMP_SP
      CALL ADD_SP(LOCAL_SP, EL_EN_DRAG_TEMP_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_EN_DRAG_TEMP_SP,MARKER_EL_EN_DRAG_TEMP,FLUID_ION_INDEX)

    END IF

    IF (COLL_EN_EX .OR. COLL_EN_ION .OR. COLL_RECOMB) THEN

      CALL INIT_EL_EN_INCOLL_TEMP_SP
      CALL ADD_SP(LOCAL_SP, EL_EN_INCOLL_TEMP_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_EN_INCOLL_TEMP_SP,MARKER_EL_EN_INCOLL_TEMP,FLUID_ION_INDEX)

    END IF

    IF (COLL_EN_EL_L_SWITCH) THEN

      CALL INIT_EL_EN_ELCOLL_TEMP_SP
      CALL ADD_SP(LOCAL_SP, EL_EN_ELCOLL_TEMP_SP,0,0)
      CALL INIT_SIMPLE_MARKER(EL_EN_ELCOLL_TEMP_SP,MARKER_EL_EN_ELCOLL_TEMP,FLUID_ION_INDEX)

    END IF


  END SUBROUTINE ADD_FLUID_EL_CENTRE
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Prepare matrix for solver M_SOLVER = I - dt * M
  SUBROUTINE MAT_PREP(M,M_SOLVER,dt)

    IMPLICIT NONE

    TYPE(SPARSE_MAT), INTENT(INOUT) :: M, &
                                       M_SOLVER

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: dt

    INTEGER :: I

    CALL SPARSE_EQ(M_SOLVER,M)

    M_SOLVER%VALUE = - dt * M_SOLVER%VALUE                                      !Multiply all elements with negative timestep

    DO I = 1, LOC_ROWS

      M_SOLVER%VALUE(MARKER_ID(I)) = M_SOLVER%VALUE(MARKER_ID(I)) + 1.00D00     !Add identity matrix

    END DO

  END SUBROUTINE MAT_PREP
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Sets all changing matrix elements to 0, and if sink on reset pre-flush elements
  SUBROUTINE MAT_FLUSH

    IMPLICIT NONE

    INTEGER :: I

    DO I = FL_INDEX + 1, LOCAL_SP%N_NZ

      LOCAL_M%VALUE(I) = 0.00D00

    END DO

!Reset fixed matrix elements (from x-advection/Maxwell-Ampere law)
    IF (PLASMA_SINK_ON .AND. (FL_INDEX .GT. 0)) THEN

      DO I = 1, FL_INDEX

        LOCAL_M%VALUE(I) = FIXED_MAT%VALUE(I)

      END DO

    END IF

  END SUBROUTINE MAT_FLUSH
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Counts numbers of nonzeros in each row for later PETSc use
  SUBROUTINE COUNT_NNZ

    IMPLICIT NONE

    INTEGER :: I, RANK, IERR
    INTEGER, ALLOCATABLE :: COUNTER(:)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    ALLOCATE(RD(RANK * nd * NUM_0D + 1:RANK * nd * NUM_0D + LOC_ROWS))

!Initialize row data to 0
    DO I = 1, LOC_ROWS

      RD(I + RANK * nd * NUM_0D)%NNZ = 0
      RD(I + RANK * nd * NUM_0D)%D_NNZ = 0
      RD(I + RANK * nd * NUM_0D)%OD_NNZ = 0

    END DO
!Count number of nonzeroes in each row
    DO I = 1, LOCAL_SP%N_NZ

      RD(LOCAL_SP%ROW(I))%NNZ = RD(LOCAL_SP%ROW(I))%NNZ + 1

!Count diagonal block nonzeroes in each row
      IF ((LOCAL_SP%COL(I) .GT. RANK * nd * NUM_0D) .AND. LOCAL_SP%COL(I) .LE. RANK * nd * NUM_0D + LOC_ROWS) &
          RD(LOCAL_SP%ROW(I))%D_NNZ = RD(LOCAL_SP%ROW(I))%D_NNZ + 1

    END DO
!Initializes markers and calculate off-diagonal number of nonzeroes
    DO I = 1 + RANK * nd * NUM_0D,LOC_ROWS + RANK * nd * NUM_0D

      ALLOCATE(RD(I)%COL(RD(I)%NNZ))
      ALLOCATE(RD(I)%MARKER(RD(I)%NNZ))
      RD(I)%OD_NNZ = RD(I)%NNZ - RD(I)%D_NNZ

    END DO
!Allocate counter for each row
    ALLOCATE(COUNTER(1 + RANK * nd * NUM_0D:LOC_ROWS + RANK * nd * NUM_0D))

    COUNTER = 1

!Set column value for each nonzero in each row, and link it to corresponding sparse matrix element via marker
    DO I = 1, LOCAL_SP%N_NZ

      RD(LOCAL_SP%ROW(I))%COL(COUNTER(LOCAL_SP%ROW(I))) = LOCAL_SP%COL(I)
      RD(LOCAL_SP%ROW(I))%MARKER(COUNTER(LOCAL_SP%ROW(I))) = I

      COUNTER(LOCAL_SP%ROW(I)) = COUNTER(LOCAL_SP%ROW(I)) + 1

    END DO

    ALLOCATE(D_NNZ(LOC_ROWS))
    ALLOCATE(OD_NNZ(LOC_ROWS))

!Set local (off-)diagonal number of nonzeroes for PETSc
    DO I = 1, LOC_ROWS

      D_NNZ(I) = RD(I + RANK * nd * NUM_0D)%D_NNZ
      OD_NNZ(I) = RD(I + RANK * nd * NUM_0D)%OD_NNZ

    END DO

  END SUBROUTINE COUNT_NNZ
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE MATRIX_DATA
