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
!> Contains all switches
MODULE SWITCHES

  USE INPUT
  USE MPI
  USE VAR_KIND_DATA

! Vlasov and field switches
  LOGICAL :: MAXWELL_SWITCH, & !< Evolve E-field using Maxwell's equations (Ampere-Maxwell's law)
             X_ADV_SWITCH, &   !< Include spatial advection
             E_ADV_SWITCH      !< Include velocity space advection due to E-field

! Collisional and radiative switches
  LOGICAL :: COLL_EE_0_SWITCH, & !< Include e-e collisions for \f_0^0
             COLL_EE_L_SWITCH, & !< Include e-e collisions for  f_{l>0}^m
             DIAG_EE_L_SWITCH, & !< Use only tridiagonal matrix elements in e-e collisions for  f_{l>0}^m
             COLL_EI_L_SWITCH, & !< Include pitch-angle scattering with ions
             COLL_EN_EL_0_SWITCH, & !< Include e-n elastic scattering for f_0^0
             COLL_EN_EL_L_SWITCH, & !< Include e-n elastic scattering for  f_{l>0}^m
             COLL_EN_EX, & !< Include e-n excitation and radiative de-excitation
             COLL_EN_ION, & !< Include e-n ionization
             COLL_RECOMB !< Include radiative and 3-body recombination

! Grid and boundary switches
  LOGICAL :: LOGARITHMIC_GRID_SWITCH, & !< Logarithmic grid with largest cell equal to dx and smallest to SMALL_dx
             PERIODIC_BOUNDARY_SWITCH, & !< Turn on periodic boundaries
             FIXED_BOUNDARY_UP_SWITCH, & !< Fix upstream boundary cell centre
             FIXED_BOUNDARY_DIV_SWITCH, & !< Fix divertor cell centre
             NO_FLOW_BOUNDARY_UP_SWITCH, & !< Assume no electric field or flows on upstream boundary
             NO_FLOW_BOUNDARY_DIV_SWITCH !< Assume no electric field or flows on downstream boundary

! Neutral switches
  LOGICAL :: NEUTRAL_TRACK_SWITCH, & !< Track neutral state densities
             NEUTRAL_DIFFUSION_SWITCH, & !< Turn on effective neutral diffusions
             RECYCLING_SWITCH, & !< Turn on neutral recycling at divertor
             FAST_DETAILED_BALANCE_SWITCH !< Update detailed balance cross-sections only once every timestep using previous timestep data

! Heating and sink switches
  LOGICAL :: HEATING_SWITCH, & !< Turn on effective heating term
             PLASMA_SINK_SWITCH !< Turn on plasma sink at divertor

! Output switches
  LOGICAL :: OUTPUT_DENSITY, & !< Output electron/ion density
             OUTPUT_TEMP, &    !< Output electron temperature
             OUTPUT_FLOW_VEL, & !< Output electron/ion flow velocity
             OUTPUT_HEAT_FLOW, & !< Output electron conduction heat flux
             OUTPUT_E_FIELD, & !< Output E-field
             OUTPUT_SH_TEST, & !< Output Spitzer-Harm test results for heat flow
             OUTPUT_RATE_DATA, & !< Output ionization rates
             OUTPUT_NEUTRAL_DATA, & !< Output neutral state densities
             OUTPUT_ATOMIC_EN_TEST_SWITCH, & !< Output total energy contained in atomic states and thermal motion (for hydrogen)
             OUTPUT_QN_TEST, & !< Output quasineutrality tests
             OUTPUT_CURRENT_TEST !< Output current (ambipolarity) test

! Initialization switches (if not periodic)
  LOGICAL :: LINEAR_INIT, & !< Initialize temperature and density as linear
             DROP_INIT, & !< Initialize temperature and density with an optional exponential drop/jump
             TWO_POINT_M_INIT, & !< Initialize temperature and density form Two-Point model
             UNIFORM_NEUTRAL_INIT, & !< Initialize neutrals in ground state uniformly
             NEUTRAL_CLOUD_INIT, &  !< Initialize neutrals as cloud at divertor (with optional exponential ramp maintining constant density if one used for DROP_INIT)
             DENSITY_FROM_FILE_INIT, & !< Load electron density from file
             TEMPERATURE_FROM_FILE_INIT, & !< Load electron temperature from file
             NEUTRAL_GROUND_DENS_FROM_FILE_INIT,& !< Load neutral ground state density from file
             ION_VEL_FROM_FILE_INIT, & !< Load ion velocity from file if cold ions on
             LOCAL_INIT_SWITCH, & !< Initialize f_1^0 and E-field as local values (works for periodic boundaries as well)
             LOCAL_SAHA_BOLTZMANN_INIT_SWITCH !< Initialize ionization degree and excited state populations assuming local initial electron density is total density

! Restart switches

  LOGICAL :: RESTART_SWITCH, & !< Initialize function from saved restart point
             SAVE_RESTART_SWITCH, & !<Save restart point after every timestep
             CONTINUE_RUN_SWITCH, & !< Continue previous run from saved restart point
             ADAPTIVE_RESTART_SWITCH !< Input restart vector updated for current (larger) number of neutral states and/or harmonics

! Fluid Ions

  LOGICAL :: COLD_ION_FLUID_SWITCH !< Include continuity and momentum equations for cold ions
  LOGICAL :: ION_CONT_OFF_SWITCH !< Turn off ion continuity equation and force n_i = n_e
  LOGICAL :: ION_CONV_UPWINDING_SWITCH !< Use first order upwinding for ion convection term
  LOGICAL :: SIMPLE_CX_SWITCH !< Include simple (cold ion, cold neutral, constant cross-section) charge exchange
  LOGICAL :: ION_EL_TEMP_SWITCH !< Include pressure gradient term in ion momentum equation with T_i = T_e
  LOGICAL :: NO_EXTRAPOLATION_SWITCH !< Turn off extrapolation for Bohm criterion
  LOGICAL :: SONIC_OUTFLOW_DIV_SWITCH !< Fix divertor outflow to ion sound speed in last cell

! Particle source switches

  LOGICAL :: PART_SOURCE_SWITCH !< Turn on upstream particle source
  LOGICAL :: PART_SOURCE_BACKGROUND_TEMP_SWITCH !< Set particle source temperature to background electron temperature

! Other switches

  LOGICAL :: ADAPTIVE_TIMESTEP_SWITCH !< Use adaptive timestep (rescale dt with min(T^(3/2)/n))
  LOGICAL :: FULL_FLUID_MODE !< Both ions and electrons as fluid, turn off kinetic operators, force everything to Maxwellian
  LOGICAL :: Z_PROFILE_FROM_FILE !< Load ionization profile from file
  LOGICAL :: X_GRID_FROM_FILE !< Load the spatial grid from file
!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Initialize all switches though input
  SUBROUTINE INIT_SWITCHES

    IMPLICIT NONE

    INTEGER :: RANK, IERR

    CALL MPI_Comm_rank(MPI_COMM_WORLD,RANK,IERR)

    IF (RANK .EQ. 0) CALL INPUT_SWITCHES(MAXWELL_SWITCH, &
                                         X_ADV_SWITCH, &
                                         E_ADV_SWITCH,&
                                         COLL_EE_0_SWITCH, &
                                         COLL_EE_L_SWITCH, &
                                         DIAG_EE_L_SWITCH, &
                                         COLL_EI_L_SWITCH, &
                                         COLL_EN_EL_0_SWITCH, &
                                         COLL_EN_EL_L_SWITCH, &
                                         COLL_EN_EX, &
                                         COLL_EN_ION, &
                                         COLL_RECOMB, &
                                         LOGARITHMIC_GRID_SWITCH, &
                                         PERIODIC_BOUNDARY_SWITCH, &
                                         FIXED_BOUNDARY_UP_SWITCH, &
                                         FIXED_BOUNDARY_DIV_SWITCH, &
                                         NO_FLOW_BOUNDARY_UP_SWITCH, &
                                         NO_FLOW_BOUNDARY_DIV_SWITCH, &
                                         NEUTRAL_TRACK_SWITCH, &
                                         NEUTRAL_DIFFUSION_SWITCH, &
                                         RECYCLING_SWITCH, &
                                         FAST_DETAILED_BALANCE_SWITCH, &
                                         HEATING_SWITCH, &
                                         PLASMA_SINK_SWITCH, &
                                         OUTPUT_DENSITY, &
                                         OUTPUT_TEMP, &
                                         OUTPUT_FLOW_VEL, &
                                         OUTPUT_HEAT_FLOW, &
                                         OUTPUT_E_FIELD, &
                                         OUTPUT_SH_TEST, &
                                         OUTPUT_RATE_DATA, &
                                         OUTPUT_NEUTRAL_DATA, &
                                         OUTPUT_ATOMIC_EN_TEST_SWITCH, &
                                         OUTPUT_QN_TEST, &
                                         OUTPUT_CURRENT_TEST, &
                                         LINEAR_INIT, &
                                         DROP_INIT, &
                                         TWO_POINT_M_INIT, &
                                         UNIFORM_NEUTRAL_INIT, &
                                         NEUTRAL_CLOUD_INIT, &
                                         DENSITY_FROM_FILE_INIT, &
                                         TEMPERATURE_FROM_FILE_INIT, &
                                         NEUTRAL_GROUND_DENS_FROM_FILE_INIT, &
                                         ION_VEL_FROM_FILE_INIT, &
                                         LOCAL_INIT_SWITCH, &
                                         LOCAL_SAHA_BOLTZMANN_INIT_SWITCH, &
                                         RESTART_SWITCH, &
                                         SAVE_RESTART_SWITCH, &
                                         CONTINUE_RUN_SWITCH, &
                                         ADAPTIVE_RESTART_SWITCH, &
                                         COLD_ION_FLUID_SWITCH, &
                                         ION_CONT_OFF_SWITCH, &
                                         ION_CONV_UPWINDING_SWITCH, &
                                         SIMPLE_CX_SWITCH, &
                                         ION_EL_TEMP_SWITCH, &
                                         NO_EXTRAPOLATION_SWITCH, &
                                         SONIC_OUTFLOW_DIV_SWITCH, &
                                         PART_SOURCE_SWITCH, &
                                         PART_SOURCE_BACKGROUND_TEMP_SWITCH, &
                                         ADAPTIVE_TIMESTEP_SWITCH, &
                                         FULL_FLUID_MODE, &
                                         Z_PROFILE_FROM_FILE, &
                                         X_GRID_FROM_FILE &
                                         )

    CALL MPI_Bcast(MAXWELL_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(X_ADV_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(E_ADV_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_EE_0_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_EE_L_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(DIAG_EE_L_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_EI_L_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_EN_EL_0_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_EN_EL_L_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_EN_EX,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_EN_ION,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLL_RECOMB,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(LOGARITHMIC_GRID_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(PERIODIC_BOUNDARY_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(FIXED_BOUNDARY_UP_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(FIXED_BOUNDARY_DIV_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NO_FLOW_BOUNDARY_UP_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NO_FLOW_BOUNDARY_DIV_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NEUTRAL_TRACK_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NEUTRAL_DIFFUSION_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(RECYCLING_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(FAST_DETAILED_BALANCE_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(HEATING_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(PLASMA_SINK_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_DENSITY,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_TEMP,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_FLOW_VEL,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_HEAT_FLOW,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_E_FIELD,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_SH_TEST,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_RATE_DATA,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_NEUTRAL_DATA,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_ATOMIC_EN_TEST_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_QN_TEST,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(OUTPUT_CURRENT_TEST,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(LINEAR_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(DROP_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(TWO_POINT_M_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(UNIFORM_NEUTRAL_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NEUTRAL_CLOUD_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(DENSITY_FROM_FILE_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(TEMPERATURE_FROM_FILE_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NEUTRAL_GROUND_DENS_FROM_FILE_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(ION_VEL_FROM_FILE_INIT,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(LOCAL_INIT_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(LOCAL_SAHA_BOLTZMANN_INIT_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(RESTART_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(SAVE_RESTART_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(CONTINUE_RUN_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(ADAPTIVE_RESTART_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(COLD_ION_FLUID_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(ION_CONT_OFF_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(ION_CONV_UPWINDING_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(SIMPLE_CX_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(ION_EL_TEMP_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(NO_EXTRAPOLATION_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(SONIC_OUTFLOW_DIV_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(PART_SOURCE_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(PART_SOURCE_BACKGROUND_TEMP_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(ADAPTIVE_TIMESTEP_SWITCH,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(FULL_FLUID_MODE,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(Z_PROFILE_FROM_FILE,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)
    CALL MPI_Bcast(X_GRID_FROM_FILE,1,MPI_LOGICAL,0, MPI_COMM_WORLD,IERR)

  END SUBROUTINE INIT_SWITCHES
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE SWITCHES
!-------------------------------------------------------------------------------------------------------------------------------------
