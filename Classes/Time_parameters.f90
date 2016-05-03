! ============================================================================
! Name        : Time_parameters.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Time_parameters.f90
!
! Description :  This module encapsules parameters used to setup time.
!
!<
! ============================================================================
!> Module time_data
!<
module time_data

    use parallel

    real( tkind ) :: time_tolerance
    real( tkind ) :: elapsed_time

    ! Start and end time of the simulation
    real( tkind ) :: time_start, time_end

    ! Layer time step
    real( tkind ) :: time_layer

    ! Output time step
    real( tkind ) :: time_output

    ! Overland flow sampling interval
    real( tkind ) :: rain_int

    ! River/turbidity flow sampling interval
    real( tkind ) :: flow_int

   ! Sampling interval
   real( tkind ) :: sampling_int

   ! Flexure interval
   real( tkind ) :: flex_int

   ! Current simulation time
   real( tkind ) :: tnow

   ! Next sampling time
   real( tkind ) :: tsampling

   ! Next output time
   real( tkind ) :: toutput

   ! Next layer time
   real( tkind ) :: tlayer

   ! Next flexure time
   real( tkind ) :: tflex

   ! Next flow time for rain
   real( tkind ) :: tflow

   ! Next flow time for river / turbidity
   real( tkind ) :: trflow

   ! Next sea time
   real( tkind ) :: tsea

   ! FW time step for overland flows
   real( tkind ) :: fw_time
   real( tkind ) :: fw_dt
   real( tkind ) :: fw_dtnext

   ! FW time step for river / turbidity flows
   real( tkind ) :: fr_time
   real( tkind ) :: fr_dt
   real( tkind ) :: fr_dtnext

   ! Number of seconds per year
   real( tkind ), parameter :: secyear = 31536000.0_8

end module time_data
! ============================================================================
