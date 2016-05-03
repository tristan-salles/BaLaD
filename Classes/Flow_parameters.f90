! ============================================================================
! Name        : Flow_parameters.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Flow_parameters.f90
!
! Description :  This module encapsulates hydrological parameters used to simulate flow.
!
!<
! ============================================================================
!> Module flow_data
!<
module flow_data

    use parallel

    implicit none

   ! Gravity
   real( tkind ), parameter :: gravity = 9.81_8

    ! von Karman constant
    real( tkind ), parameter :: vonKarman = 0.41_8

    ! Kinematic viscosity
    real( tkind ), parameter :: kinematic_viscosity = 2.e-6_8

    ! Zero check accuracy
    real( tkind ), parameter :: tor = 1.0e-8

    ! Pi number
    real( tkind ), parameter :: pi = 4.0_8 * atan( 1.0_8 )

    ! Rain event duration
    real( tkind ) :: rain_time

    ! Depressionless step fill value
    real( tkind ) :: step_fill

    ! Declaration of manning parameters
    real( tkind ) :: manning_open, manning_hypo, manning_hyper

    ! Parameters to define channel width based on flow rate and slope
    real( tkind ) :: kwidth
    real( tkind ) :: awidth
    real( tkind ) :: bwidth

    ! Parameters to define erosion / deposition model
    real( tkind ) :: mdetach
    real( tkind ) :: ndetach
    real( tkind ) :: mtransport
    real( tkind ) :: ntransport

    ! Erosional threshold to trigger flow walkers (m)
    real( tkind ) :: erosion_trigger

    ! Number of rain events during the simulation
    integer :: rain_event

    ! Current rain event map ID
    integer :: rainID

    ! Starting time of the considered rainfall period
    real( tkind ),dimension(:),allocatable :: rain_tstart

    ! Ending time of the considered rainfall period
    real( tkind ),dimension(:),allocatable :: rain_tend

    ! Average duration of a rainfall period
    real( tkind ),dimension(:),allocatable :: rain_duration

    ! Precipitation rate in metres per year
    real( tkind ),dimension(:),allocatable :: precipitation

    ! Steady state flow height
    real( tkind ),dimension(:),allocatable :: hflow

    ! Steady state flow width
    real( tkind ),dimension(:),allocatable :: wflow

    ! Steady state flow rate
    real( tkind ),dimension(:),allocatable :: qflow

    ! Steady state flow velocity
    real( tkind ),dimension(:),allocatable :: vflow

    !------------------
    ! Debris flows
    !------------------

    ! Maximum critical slope
    real( tkind ) :: critical_slope_marine, critical_slope_aerial

    ! Debris flow max concentration aerial
    real( tkind ) :: debris_max_conc_aerial

    ! Debris flow max concentration marine
    real( tkind ) :: debris_max_conc_marine

    ! Debris flow max correction factor aerial
    real( tkind ) :: debris_factor_aerial

    ! Debris flow max correction factor marine
    real( tkind ) :: debris_factor_marine

    !------------------
    ! Flow walkers for overland flows
    !------------------

    logical :: rainfws

    ! Flow walker face ID
    integer :: fw_face

    ! Maximum number of FWs
    integer :: maxFWs

    ! Local number of FWs
    integer :: procFWs

    ! Position values
    real( tkind ) :: fw_xpos, fw_ypos, fw_zpos, fw_zstrat

    ! Slope values
    real( tkind ) :: fw_xslp, fw_yslp

    ! Flow height and discharge values
    real( tkind ) :: fw_h, fw_q

    ! Friction values
    real( tkind ) :: fw_Cfric

    ! Velocity values
    real( tkind ) :: fw_xv, fw_yv

    ! Velocity new values
    real( tkind ) :: fw_xv_new, fw_yv_new

    ! Acceleration values
    real( tkind ) :: fw_accx, fw_accy

    ! Flow walker initial sediment thickness
    real( tkind ), dimension(:), allocatable :: fw_hs0

    ! Flow walker current sediment thickness
    real( tkind ), dimension(:), allocatable :: fw_hs

    !------------------
    ! River / turbidity flows declaration
    !------------------

    ! Number of sources
    integer :: num_src

    ! Number of FWs per source
    integer :: nbfw_src

    ! Souce splitting number
    integer,dimension(:),allocatable :: fws_num

    ! Source flow type
    integer,dimension(:),allocatable :: fws_type

    ! Source start time
    real( tkind ),dimension(:),allocatable :: fws_tstrt

    ! Source end time
    real( tkind ),dimension(:),allocatable :: fws_tend

    ! Source x position
    real( tkind ),dimension(:),allocatable :: fws_xposition

    ! Source y position
    real( tkind ),dimension(:),allocatable :: fws_yposition

    ! Source x range
    real( tkind ),dimension(:),allocatable :: fws_xrange

    ! Source y range
    real( tkind ),dimension(:),allocatable :: fws_yrange

    ! Source x velocity
    real( tkind ),dimension(:),allocatable :: fws_xvel

    ! Source y velocity
    real( tkind ),dimension(:),allocatable :: fws_yvel

    ! Source flow discharge
    real( tkind ),dimension(:),allocatable :: fws_qfl

    ! Source water volume
    real( tkind ),dimension(:),allocatable :: fws_volume

    ! Source sediment concentration
    real( tkind ),dimension(:),allocatable :: fws_sedconc

    ! Source sediment percentage
    real( tkind ),dimension(:,:),allocatable :: fws_sedperc

    ! Source sediment volume (m3)
    real( tkind ),dimension(:,:),allocatable :: fws_sedcharge

    !------------------
    ! Flow walkers for river / turbidity flows
    !------------------

    ! Position values
    real( tkind ) :: fr_xpos, fr_ypos, fr_zpos, fr_zstrat

    ! Slope values
    real( tkind ) :: fr_xslp, fr_yslp

    ! Flow height and discharge values
    real( tkind ) :: fr_h, fr_w, fr_q

    ! Friction values
    real( tkind ) :: fr_Cfric, fr_vol, fr_vol0

    ! Velocity values
    real( tkind ) :: fr_xv, fr_yv

    ! Velocity new values
    real( tkind ) :: fr_xv_new, fr_yv_new

    ! Acceleration values
    real( tkind ) :: fr_accx, fr_accy

    ! Flow walker sediment concentration
    real( tkind ), dimension(:), allocatable :: fr_sedcharge

    ! Flow walker density
    real( tkind ) :: fr_density

    ! Flow walker face ID
    integer :: fr_face
    integer, dimension( 50000 ) :: fr_recface

    !------------------
    ! Record flow walkers for visualisation
    !------------------

    integer :: max_rec = 5000000

    integer :: rec_nb

    ! FW x position
    real( tkind ),dimension(:),allocatable :: rfw_x

    ! FW y position
    real( tkind ),dimension(:),allocatable :: rfw_y

    ! FW z position
    real( tkind ),dimension(:),allocatable :: rfw_z

    ! FW velocity
    real( tkind ),dimension(:),allocatable :: rfw_v

    ! FW source height
    real( tkind ),dimension(:),allocatable :: rfw_h

    ! FW source width
    integer,dimension(:),allocatable :: rfw_id

end module flow_data
! ============================================================================
