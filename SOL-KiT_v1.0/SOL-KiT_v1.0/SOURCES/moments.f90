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
!> Contains all moment calculation routines
MODULE MOMENTS

  USE GRID
  USE NORMALIZATION
  USE MPI
  USE VAR_KIND_DATA

!-------------------------------------------------------------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------------------------------------------------------------
!> Generic MOM_N-th moment multiplied by 4 PI
  REAL(KIND=DEFAULT_REAL) FUNCTION MOMENT_V(f,MOM_N)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(NUM_V) !< Integrand distribution
    INTEGER, INTENT (IN) :: MOM_N !< Moment number
    INTEGER :: I
    REAL(KIND=DEFAULT_REAL) :: m

    m = 0

    DO I = 1, NUM_V

      m = m + 4.0 * PI * f(I) * V_GRID(I) ** (2 + MOM_N) * V_GRID_WIDTH(I)

    END DO

    MOMENT_V = m

  END FUNCTION MOMENT_V
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates density moment
  REAL(KIND=DEFAULT_REAL) FUNCTION DENSITY_MOMENT(f)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f(NUM_V) !< Input f_0

    DENSITY_MOMENT = MOMENT_V(f,0)

  END FUNCTION DENSITY_MOMENT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates flow velocity in x-direction
  REAL(KIND=DEFAULT_REAL) FUNCTION FLOW_V_X_MOMENT(f_1,f_0)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_1 (NUM_V), & !< Input f_10
                          f_0(NUM_V) !< Input f_0

    FLOW_V_X_MOMENT = MOMENT_V(f_1,1)/(3*DENSITY_MOMENT(f_0))

  END FUNCTION FLOW_V_X_MOMENT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates flow velocity in y-direction
  REAL(KIND=DEFAULT_REAL) FUNCTION FLOW_V_Y_MOMENT(f_1,f_0)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_1 (NUM_V), & !< Input Re(f_11)
                          f_0(NUM_V) !< Input f_0

    FLOW_V_Y_MOMENT = 2.00D00 * MOMENT_V(f_1,1)/(3*DENSITY_MOMENT(f_0))

  END FUNCTION FLOW_V_Y_MOMENT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!>Calculates flow velocity in z-direction
  REAL(KIND=DEFAULT_REAL) FUNCTION FLOW_V_Z_MOMENT(f_1,f_0)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_1 (NUM_V), & !< Input Im(f_11)
                          f_0(NUM_V) !< Input f_0

    FLOW_V_Z_MOMENT = - 2.00D00 * MOMENT_V(f_1,1)/(3*DENSITY_MOMENT(f_0))

  END FUNCTION FLOW_V_Z_MOMENT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculate temperature
  REAL(KIND=DEFAULT_REAL) FUNCTION TEMPERATURE_MOMENT(f_0,f_10, Ref_11, Imf_11)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_0(NUM_V) !< Input f_0
    REAL(KIND=DEFAULT_REAL), OPTIONAL, DIMENSION(NUM_V), INTENT(IN) :: f_10, & !< Input f_10
                                                      Ref_11, & !< Input Re(f_11)
                                                      Imf_11 !< Input Re(f_11)
    REAL(KIND=DEFAULT_REAL) :: u_x, &
              u_y, &
              u_z

    u_x = 0
    u_y = 0
    u_z = 0

    IF (PRESENT(f_10)) u_x = FLOW_V_X_MOMENT(f_10,f_0)
    IF (PRESENT(Ref_11)) u_y = FLOW_V_Y_MOMENT(Ref_11,f_0)
    IF (PRESENT(Imf_11)) u_z = FLOW_V_Z_MOMENT(Imf_11,f_0)

    TEMPERATURE_MOMENT = 2.00D00 * MOMENT_V(f_0,2)/(3.00D00*DENSITY_MOMENT(f_0)) &
                       - 2.00D00 * (u_x * u_x + u_y * u_y + u_z * u_z)/3.00D00                   !Calculate total energy moment ad subtract kinetic term

  END FUNCTION TEMPERATURE_MOMENT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculate pressure tensor
  FUNCTION PRESSURE_MOMENT(f_0, f_10, Ref_11, Imf_11, f_20, Ref_21, Imf_21, Ref_22, Imf_22)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: PRESSURE_MOMENT(3,3), &
              P(3,3)
    REAL(KIND=DEFAULT_REAL), INTENT(IN) :: f_0(NUM_V) !< Input f_0
    REAL(KIND=DEFAULT_REAL), OPTIONAL, DIMENSION(NUM_V), INTENT(IN) :: f_10, &  !< Input f_10
                                                      Ref_11, & !< Input Re(f_11)
                                                      Imf_11, & !< Input Im(f_11)
                                                      f_20, & !< Input f_20
                                                      Ref_21, & !< Input Re(f_21)
                                                      Imf_21, & !< Input Im(f_21)
                                                      Ref_22, & !< Input Re(f_22)
                                                      Imf_22 !< Input Im(f_22)

    REAL(KIND=DEFAULT_REAL) :: ap_11, &
              ap_12, &
              ap_13, &
              ap_22, &
              ap_23, &
              ap_33, &
              u_x, &
              u_y, &
              u_z, &
              ip, &
              n

    ap_11 = 0 ; ap_12 = 0 ; ap_13 = 0 ; ap_22 = 0 ; ap_23 = 0 ; ap_33 = 0
    u_x = 0 ; u_y = 0; u_z = 0

    ip = MOMENT_V(f_0,2)/3.0                                                    !Isotropic part of total pressure tensor
    n = DENSITY_MOMENT(f_0)

    IF (PRESENT(f_10)) u_x = FLOW_V_X_MOMENT(f_10,f_0)
    IF (PRESENT(Ref_11)) u_y = FLOW_V_Y_MOMENT(Ref_11,f_0)
    IF (PRESENT(Imf_11)) u_z = FLOW_V_Z_MOMENT(Imf_11,f_0)

    IF (PRESENT(f_20)) THEN                                                     !Anisotropic parts of the tensor

      ap_11 = 2.00D00 * MOMENT_V(f_20,2)/15.00D00
      ap_22 = ap_22 - MOMENT_V(f_20,2)/15.00D00
      ap_33 = ap_33 - MOMENT_V(f_20,2)/15.00D00

    END IF

    IF (PRESENT(Ref_21)) THEN

      ap_12 = 6.00D00 * MOMENT_V(Ref_21,2)/15.00D00
      ap_13 = - 6.00D00 * MOMENT_V(Imf_21,2)/15.00D00

    END IF

    IF (PRESENT(Ref_22)) THEN

      ap_22 = ap_22 + 12.00D00 * MOMENT_V(Ref_22,2)/15.00D00
      ap_33 = ap_33 - 12.00D00 * MOMENT_V(Ref_22,2)/15.00D00
      ap_23 = -12.00D00 * MOMENT_V(Imf_22,2)/15.00D00

    END IF

    P(1,1) = ip + ap_11 - n * u_x * u_x/3                                         !Add isotropic and anisotropic parts and subtract dynamic pressure
    P(1,2) = ap_12 - n * u_x * u_y/3
    P(1,3) = ap_13 - n * u_x * u_z/3
    P(2,1) = P(1,2)
    P(2,2) = IP + ap_22 - n * u_y * u_y/3
    P(2,3) = ap_23 - n * u_y * u_z/3
    P(3,1) = P(1,3)
    P(3,2) = P(2,3)
    P(3,3) = ap_33 + IP - n * u_z * u_z/3

    PRESSURE_MOMENT = P

  END FUNCTION PRESSURE_MOMENT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!> Calculates heat flow vector
  FUNCTION HEAT_FLOW_MOMENT(f_0, f_10, Ref_11, Imf_11, f_20, Ref_21, Imf_21, Ref_22, Imf_22)

    IMPLICIT NONE

    REAL(KIND=DEFAULT_REAL) :: HEAT_FLOW_MOMENT(3)

    REAL(KIND=DEFAULT_REAL), DIMENSION(NUM_V), INTENT(IN) :: f_0, & !< Input f_0
                                            f_10 !< Input f_10

    REAL(KIND=DEFAULT_REAL), OPTIONAL, DIMENSION(NUM_V), INTENT(IN) :: Ref_11, & !< Input Re(f_11)
                                                      Imf_11, & !< Input Im(f_11)
                                                      f_20, & !< Input f_20
                                                      Ref_21, & !< Input Re(f_21)
                                                      Imf_21, & !< Input Im(_21)
                                                      Ref_22, & !< Input Re(f_22)
                                                      Imf_22 !< Input Im(f_22)

    REAL(KIND=DEFAULT_REAL) :: u_x, &
              u_y, &
              u_z, &
              n, &
              T, &
              p(3,3), &
              q_T(3), &
              q(3)


    n = DENSITY_MOMENT(f_0)
    q_T(1) = MOMENT_V(f_10,3)/6.00D00                                           !Calculate total heat flow
    q_T(2) = 0
    q_T(3) = 0
    p = 0
    u_x = 0; u_y =0; u_z = 0
    IF (PRESENT(Ref_11)) THEN

      q_T(2) = 2.00D00 * MOMENT_V(Ref_11,3)/6.00D00
      q_T(3) = - 2.00D00* MOMENT_V(Imf_11,3)/6.00D00

    END IF

    u_x = FLOW_V_X_MOMENT(f_10,f_0)

    IF (PRESENT(Ref_11)) THEN

      u_y = FLOW_V_Y_MOMENT(Ref_11,f_0)
      u_z = FLOW_V_Z_MOMENT(Imf_11,f_0)
      T = TEMPERATURE_MOMENT(f_0, f_10, Ref_11, Imf_11)

    ELSE

      T = TEMPERATURE_MOMENT(f_0,f_10)

    END IF

    IF (PRESENT(Ref_22)) THEN

      p = PRESSURE_MOMENT(f_0, f_10, Ref_11, Imf_11, f_20, Ref_21, Imf_21, Ref_22, Imf_22)

    ELSE IF (PRESENT(Ref_21)) THEN

      p = PRESSURE_MOMENT(f_0, f_10, Ref_11, Imf_11, f_20, Ref_21, Imf_21)

    ELSE IF (PRESENT(f_20)) THEN

      IF (PRESENT(Ref_21)) THEN

        p = PRESSURE_MOMENT(f_0, f_10, Ref_11, Imf_11, f_20)

      ELSE

        p = PRESSURE_MOMENT(f_0, f_10, f_20)

      END IF

    ELSE IF (PRESENT(Ref_11)) THEN

      p = PRESSURE_MOMENT(f_0, f_10, Ref_11, Imf_11)

    ELSE

      p = PRESSURE_MOMENT(f_0, f_10)

    END IF

    q(1) = q_T(1) - n * (u_x * u_x + u_y * u_y + u_z * u_z) * u_x/2.00D00 - 3.00D00 * n * T * u_x/4.00D00 - &
          (u_x * p(1,1) + u_y * p(2,1) + u_z * p(3,1))                                                                               ! Subtract convective flows from total flow


     q(2) = q_T(2) - n * (u_x * u_x + u_y * u_y + u_z * u_z) * u_y/2.00D00 - 3.00D00 * n * T * u_y/4.00D00 - &
           (u_x * p(1,2) + u_y * p(2,2) + u_z* p(3,2))

     q(3) = q_T(3) - n * (u_x * u_x + u_y * u_y + u_z * u_z) * u_z/2.00D00 - 3.00D00 * n * T * u_z/4.00D00 -  &
           (u_x * p(1,3) + u_y * p(2,3) + u_z* p(3,3))

    HEAT_FLOW_MOMENT = q

  END FUNCTION HEAT_FLOW_MOMENT
!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
END MODULE MOMENTS
!-------------------------------------------------------------------------------------------------------------------------------------
