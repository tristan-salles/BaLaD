! ============================================================================
! Name        : RunProcesses.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file RunProcesses.f90
!
! Description :  This module  is called when wathever main subroutines stopped. It
! then determines which subroutine will be call next.
!
!<
! ============================================================================
!> Module RunProcesses 
!<
module run_phase

    use tin
    use lem
    use mass
    use hillslp
    use isoflex
    use steady
    use fwpath
    use parallel
    use regolith
    use diffusion
    use file_data
    use flow_swe
    use flow_algo
    use flow_data
    use init_phase
    use strata_out
    use time_data
    use UDWplugin
    use checkpoint
    use compaction
    use forces_data
    use hemipelagite

    implicit none

    logical :: output_sim

    ! Event to perform
    integer :: nextevent, oldevent

    ! Get CPU time
    real( tkind ) :: ts_tp1, ts_tp2

    ! Undefined flag
    integer, parameter :: SPM_UNDEFINED = 0

    ! Quit simulation flag
    integer, parameter :: SPM_QUIT = 1

    ! Displacement event flag
    integer, parameter :: SPM_DISPLACEMENT = 2

    ! Display event flag
    integer, parameter :: SPM_LAYER = 3

    ! Flow walker rain flag
    integer, parameter :: SPM_RAIN = 4

    ! Flow walker rain flag
    integer, parameter :: SPM_FLOW = 5

    ! Sea level fluctuation time step flag
    integer, parameter :: SPM_SEALEVEL = 6

contains

    ! ============================================================================
    !> Subroutine badlands_initialise
    !! This subroutine performs the initialisation of badlands run.
    !<
    ! ============================================================================
    subroutine badlands_initialise

        real( tkind ) :: th

        ts_tp1 = mpi_wtime( )

        ! Initialisation phase
        call initialisation_phase

        ! Build TIN grid
        if( transport_mode == 0 .or. num_src > 0 ) call generate_TIN_surface

        ! Performs steady state
        call compute_steady_state

        ! Visualisation of initial settings
        call xdmf_output( 1 )
        if( rain_event > 0 .and. num_src > 0 )then
            rain_int = min( rain_int, flow_int )
        endif

        ! Update next process time
        toutput = time_start + time_output
        tlayer = time_start + time_layer
        tsampling = time_start
        tsea = time_start
        tflex = time_start + flex_int
        tflow = time_start
        trflow = time_start
        elapsed_time = time_start
        if( rain_event == 0 ) tflow = time_end + time_layer
        if( num_src == 0 ) trflow = time_end + time_layer
        layID = layID + 1

        ! Compute soil production
        elapsed_time = time_start - min( rain_int, flow_int )
        if( transport_mode /= 1 )then
            call soil_production
            call compute_hillslope( 2 )
        endif

        ! Flexural initialisation
        if( flexureon == 1 ) call initialise_isostatic_flexure( th )

        ! Allocate record arrays
        if( fwsvis == 1 )then
            if( allocated( rfw_x ) ) deallocate( rfw_x )
            if( allocated( rfw_y ) ) deallocate( rfw_y )
            if( allocated( rfw_z ) ) deallocate( rfw_z )
            if( allocated( rfw_v ) ) deallocate( rfw_v )
            if( allocated( rfw_h ) ) deallocate( rfw_h )
            if( allocated( rfw_id ) ) deallocate( rfw_id )
            allocate( rfw_x( max_rec ) )
            allocate( rfw_y( max_rec ) )
            allocate( rfw_z( max_rec ) )
            allocate( rfw_v( max_rec ) )
            allocate( rfw_h( max_rec ) )
            allocate( rfw_id( max_rec ) )
        endif

        ! Build the kdtree for flow trajectories values interpolation
        if( transport_mode == 0 .or. num_src > 0 ) call build_DEM_kdtree

        ts_tp2 = mpi_wtime( )
        if( iam == 0 )then
            write(6,*)'Initialisation phase done ...'
            write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
            write(6,*)
        endif

100   format(A18,F12.4)

        return

    end subroutine badlands_initialise
    ! ============================================================================
    !> Function shortest_interval
    !! Function shortest_interval determines the shortest interval time
    !<
    ! ============================================================================
    function shortest_interval( ) result( eval )

        real( tkind ) :: eval

        eval = 1.e16_8
        if( gdisp%event  > 0 ) eval = min( eval, tsampling - tnow )
        eval = min( eval, tlayer - tnow )
        eval = min( eval, tsea - tnow )
        eval = min( eval, tflow - tnow )
        eval = min( eval, trflow - tnow )
!        if( eval < time_tolerance ) eval = time_tolerance

    end function shortest_interval
    ! ============================================================================
    !> Subroutine control_process_events
    !! Subroutine control_process_events controls SP Model execution by generating events for main
    !! processing loop.
    !<
    ! ============================================================================
    subroutine control_process_events

        integer    :: period

        real( tkind ) :: elapsedtime
        real( tkind ), save  :: previous_time

        elapsedtime = tnow - previous_time
        previous_time = tnow
        oldevent = nextevent
        nextevent = SPM_UNDEFINED
        period = gdisp%actual

        ! Next displacement event
        if( tnow >= tsampling .and. gdisp%event > 0 .and. tnow < time_end )then
            tnow = tsampling
            tsampling = sampling_int + tsampling

            ! In case the displacement finishes before the sampling interval
            ! adjust the next displacement time
            if(gdisp_time( period, 2 ) < tsampling ) tsampling = gdisp_time( period , 2 )

            ! In case the displacement event ends before the current time
            ! upload the next event
            if( gdisp_time( period, 2 ) <= tnow .and. tnow < time_end ) then
                new_disp = .true.
                ! In case the displacement finishes before the sampling interval
                ! adjust the next displacement time
                if( gdisp_time( period + 1, 2 ) > time_start .and. gdisp_time( period + 1, 2 ) < tsampling ) &
                    tsampling = gdisp_time( period + 1, 2 )
            endif

            if( tnow >= time_end ) tsampling = tsampling + time_end
            nextevent = SPM_DISPLACEMENT

        elseif( tnow >= tsampling .and. gdisp%event == 0 .and. tnow < time_end .and. flexureon == 1 )then
            tnow = tsampling
            tsampling = sampling_int + tsampling

            if( tnow >= time_end ) tsampling = tsampling + time_end
            nextevent = SPM_DISPLACEMENT

        ! Next sea level event
        elseif( tnow >= tsea .and. tnow < time_end )then
            tnow = tsea
            tsea = sampling_int + tsea
            nextevent = SPM_SEALEVEL

        ! Next layer event
        elseif( tnow >= tlayer )then
            tnow = tlayer
            tlayer = time_layer + tlayer
            output_sim = .false.
            if( tnow >= toutput )then
                output_sim = .true.
                toutput = toutput + time_output
            endif
            nextevent = SPM_LAYER

        ! Next river / turbidity flow event
        elseif( tnow >= trflow .and. tnow < time_end )then
            tnow = trflow
            trflow = trflow + flow_int
            nextevent = SPM_FLOW

        ! Next rain flow event
        elseif( tnow >= tflow .and. tnow < time_end )then
            tnow = tflow
            tflow = tflow + rain_int
            nextevent = SPM_RAIN

        ! Final event
        elseif( tnow >= time_end )then
            nextevent = SPM_QUIT

        endif

        return

    end subroutine control_process_events
    ! ============================================================================
    !> Subroutine badlands_run
    !! Subroutine badlands_run controls the execution phase of a considered simulation.
    !! processing loop.
    !<
    ! ============================================================================
    subroutine badlands_run

        integer :: iter

        ! Output iteration
        iter = 1

        ! Allocate displacement fields for UDW plugin
        if( tnow == time_start .and. udw_plug )then
            gdisp%disptime = 0.0_8
            call create_vtk_topsurface( 0 )
            call udw_plugin_wait_function
            new_disp = .true.
            ! Assign displacement rate
            call assign_vertical_displacements_rate( 1 )
            ! Update displacement
            call update_local_displacements
        endif

        ! Manage events procedure
        event_control: do

            call control_process_events

            if( nextevent == SPM_QUIT )then

                ! Compute flexural deformation
                if( flexureon == 1 ) call compute_isostatic_flexure

                ! Exit the event control loop
                exit event_control

            elseif( nextevent == SPM_LAYER )then

                ! Get compaction layers
                if( gporo%compaction ) call cmpt_compaction

                ! Update the creep sedimentary layer arrays
                top_sedprev = 0.0_8
                top_sedh = 0.0_8

                ! In case we need to visualise the ouput
                if( output_sim .or. tnow == time_end )then
                    iter = iter + 1
                    ts_tp1 = mpi_wtime( )
                    call xdmf_output( iter )
                    ts_tp2 = mpi_wtime( )
                    if( iam == 0 )then
                        write(6,*)'----------------------------------------'
                        write(6,*)'Visualisation step done ...'
                        write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
                        write(6,101)'Simulation time (yr):',tnow
                        write(6,*)'----------------------------------------'
                        write(6,*)
                    endif
                endif
                
                ! If check pointing on, records sedimentary layers
                if( checkpointing .and. mod( layID, checkfreq ) == 0 )then
                    if( iam == checkproc ) call write_checkpoint
                ! If last time step, records sedimentary layers
                elseif( tnow == time_end )then
                    if( iam == checkproc ) call write_checkpoint
                endif

                ! Add a new sedimentary layer on top
                layID = layID + 1
                rec_nb = 0
                if( fwsvis == 1 ) rfw_id = -1

            elseif( nextevent == SPM_DISPLACEMENT )then

                if( gdisp%event > 0 )then
                    ! Apply displacements
                    call apply_displacements

                    ! Update displacement rates if required
                    if( new_disp .and. tnow < time_end )then
                        if( udw_plug .and. tnow > time_start )then
                            call create_vtk_topsurface( 0 )
                            call udw_plugin_wait_function
                        endif
                        call load_displacements( 0 )
                    endif

                    ! Update displacement time
                    gdisp%lastdisp = tnow
                endif

                ! Performs flexural isostatic deformation
                if( flexureon == 1 .and. tflex <= tnow )then
                    ts_tp1 = mpi_wtime( )
                    call compute_isostatic_flexure
                    tflex = tnow + flex_int
                    ts_tp2 = mpi_wtime( )
                    if( iam == 0 )then
                        write(6,*)'Flexural isostatic deformation done ...'
                        write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
                        write(6,101)'Simulation time (yr):',tnow
                        write(6,*)
                    endif
                endif

            elseif( nextevent == SPM_SEALEVEL )then

                call eustatism

                ! Deposit hemipelagites
                if( hemi_flag ) call hemipelagites

                ! Mass wasting
                if( rain_event == 0 .and. num_src == 0 )then
                    if( masson == 1 )then
                        ts_tp1 = mpi_wtime( )
                        call compute_steady_state
                        if( transport_mode /= 1 )then
                            call soil_production
                            call compute_hillslope( 2 )
                        endif
                        call build_debris_flow_depth
                        ts_tp2 = mpi_wtime( )
                        if( iam == 0 )then
                            write(6,*)'Steady state hydraulics and debris flow done ...'
                            write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
                            write(6,101)'Simulation time (yr):',tnow
                            write(6,*)
                        endif
                    endif
                endif

            elseif( nextevent == SPM_RAIN )then

                ! Performs steady state
                if( tnow > time_start )then
                    ts_tp1 = mpi_wtime( )
                    call compute_steady_state
                    if( transport_mode /= 1 )then
                        call soil_production
                        call compute_hillslope( 2 )
                    endif
                    ! Perform debris flow
                    if( masson == 1 ) call build_debris_flow_depth
                    ts_tp2 = mpi_wtime( )
                    if( iam == 0 )then
                        write(6,*)'Steady state hydraulics and soil production for rain done ...'
                        write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
                        write(6,101)'Simulation time (yr):',tnow
                        write(6,*)
                    endif
                endif

                ! Compute rain flow walkers evolution
                rainfws = .true.
                ts_tp1 = mpi_wtime( )
                if( transport_mode < 2 )then
                    call compute_overland_flow_walkers
                   ! For purely erosive mode add the thickness transport
                   ! due to creeping on top layer
                    if( transport_mode == 1 ) call define_soil_creep_thickness
                else
                    call advance_landscape_evolution
                endif

                ts_tp2 = mpi_wtime( )
                if( iam == 0 )then
                    write(6,*)'Overland flow evolution done ...'
                    write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
                    write(6,101)'Simulation time (yr):',tnow
                    write(6,*)
                endif

            elseif( nextevent == SPM_FLOW )then

                ! Performs steady state
                if( tnow > time_start )then
                    ts_tp1 = mpi_wtime( )
                    call compute_steady_state
                    if( transport_mode /= 1 )then
                        call soil_production
                        call compute_hillslope( 2 )
                    endif
                    ! Perform debris flow and creep transport
                    if( rain_event == 0 .and. masson == 1 ) call build_debris_flow_depth
                    ts_tp2 = mpi_wtime( )
                    if( iam == 0 )then
                        write(6,*)'Steady state hydraulics and soil production for rivers done ...'
                        write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
                        write(6,101)'Simulation time (yr):',tnow
                        write(6,*)
                    endif
                endif

                ! Compute river / turbidity current flow walkers evolution
                rainfws = .false.
                ts_tp1 = mpi_wtime( )
                call perform_rivers_evolution

                ts_tp2 = mpi_wtime( )
                if( iam == 0 .and. num_src > 0 )then
                    write(6,*)'Stream flow evolution done ...'
                    write(6,100)'Time elapse (s):',ts_tp2 - ts_tp1
                    write(6,101)'Simulation time (yr):',tnow
                    write(6,*)
                endif

            else

                ! Next process event
                tnow = tnow + shortest_interval( )

            endif

        enddo event_control

100   format(A18,F12.4)
101   format(A21,F12.2)

        return

    end subroutine badlands_run
    ! ============================================================================
    !> Subroutine badlands_finalise
    !! Subroutine badlands_finalise deallocates arrays.
    !<
    ! ============================================================================
    subroutine badlands_finalise

        ! Destroy the kdtree for flow trajectories values interpolation
        if( transport_mode == 0 .or. num_src > 0 ) call destroy_DEM_kdtree

        ! File data
        if( allocated( frainmap ) ) deallocate( frainmap )
        if( allocated( fdisp ) ) deallocate( fdisp )
        if( allocated( fdep ) ) deallocate( fdep )

        ! Sediment data
        if( allocated( material_name ) ) deallocate( material_name )
        if( allocated( sediment ) ) deallocate( sediment )

        ! Mesh data
        if( allocated( coordX ) ) deallocate( coordX )
        if( allocated( coordY ) ) deallocate( coordY )
        if( allocated( coordY ) ) deallocate( coordZ )
        if( allocated( lcoordX ) ) deallocate( lcoordX )
        if( allocated( lcoordY ) ) deallocate( lcoordY )
        if( allocated( lcoordZ ) ) deallocate( lcoordZ )
        if( allocated( lbase ) ) deallocate( lbase )
        if( allocated( lngbID ) ) deallocate( lngbID )
        if( allocated( ngbID ) ) deallocate( ngbID )
        if( allocated( fcoordX ) ) deallocate( fcoordX )
        if( allocated( fcoordY ) ) deallocate( fcoordY )
        if( allocated( ngbID ) ) deallocate( ngbID )
        if( allocated( fptIDs ) ) deallocate( fptIDs )
        if( allocated( lfptIDs ) ) deallocate( lfptIDs )
        if( allocated( lvertexID ) ) deallocate( lvertexID )
        if( allocated( fPid ) ) deallocate( fPid )
        if( allocated( global_nid ) ) deallocate( global_nid )
        if( allocated( v_nbLays ) ) deallocate( v_nbLays )
        if( allocated( v_LaysID ) ) deallocate( v_LaysID )
        if( allocated( stratLays ) ) deallocate( stratLays )
        if( allocated( locf_points ) ) deallocate( locf_points )
        if( allocated( part_rows ) ) deallocate( part_rows )
        if( allocated( Fdata ) ) deallocate( Fdata )
        if( allocated( porosityLays ) ) deallocate( porosityLays )
        if( allocated( hemipelagic ) ) deallocate( hemipelagic )

        ! Hydraulic data
        if( allocated( lfilldem ) ) deallocate( lfilldem )
        if( allocated( filldem ) ) deallocate( filldem )
        if( allocated( lfacc ) ) deallocate( lfacc )
        if( allocated( facc ) ) deallocate( facc )
        if( allocated( lfacr ) ) deallocate( lfacr )
        if( allocated( facr ) ) deallocate( facr )
        if( allocated( slp ) ) deallocate( slp )
        if( allocated( vcurv ) ) deallocate( vcurv )
        if( allocated( hcurv ) ) deallocate( hcurv )
        if( allocated( orient ) ) deallocate( orient )
        if( allocated( rain_tstart ) ) deallocate( rain_tstart )
        if( allocated( rain_tend ) ) deallocate( rain_tend )
        if( allocated( rain_duration ) ) deallocate( rain_duration )
        if( allocated( precipitation ) ) deallocate( precipitation )
        if( allocated( hflow ) ) deallocate( hflow )
        if( allocated( wflow ) ) deallocate( wflow )
        if( allocated( qflow ) ) deallocate( qflow )
        if( allocated( vflow ) ) deallocate( vflow )

        ! Creep transport data
        if( allocated( creepz ) ) deallocate( creepz )
        if( allocated( top_sedprev ) ) deallocate( top_sedprev )
        if( allocated( top_sedh ) ) deallocate( top_sedh )
        if( allocated( depo ) ) deallocate( depo )
        if( allocated( dstart ) ) deallocate( dstart )
        if( allocated( cdif ) ) deallocate( cdif )
        if( allocated( difo ) ) deallocate( difo )
        if( allocated( difp ) ) deallocate( difp )

        ! Rain and river flow walkers
        if( allocated( fw_hs ) ) deallocate( fw_hs )
        if( allocated( fw_hs0 ) ) deallocate( fw_hs0 )
        if( allocated( fr_sedcharge ) ) deallocate( fr_sedcharge )
        if( allocated( rfw_x ) ) deallocate( rfw_x )
        if( allocated( rfw_y ) ) deallocate( rfw_y )
        if( allocated( rfw_z ) ) deallocate( rfw_z )
        if( allocated( rfw_v ) ) deallocate( rfw_v )
        if( allocated( rfw_h ) ) deallocate( rfw_h )
        if( allocated( rfw_id ) ) deallocate( rfw_id )

        ! Checkpointing and restart parameters
        if( allocated( check_NbLay ) ) deallocate( check_NbLay )
        if( allocated( check_layID ) ) deallocate( check_layID )
        if( allocated( check_sed ) ) deallocate( check_sed )
        if( allocated( check_base ) ) deallocate( check_base )
        if( allocated( soil_thickness ) ) deallocate( soil_thickness )
        if( allocated( check_porosity ) ) deallocate( check_porosity )

        ! Flow accumulation
        if( allocated( degreeMatrix ) ) deallocate( degreeMatrix )
        if( allocated( directMatrix ) ) deallocate( directMatrix )
        if( allocated( reverseMatrix ) ) deallocate( reverseMatrix )
        if( allocated( oldfacc ) ) deallocate( oldfacc )
        if( allocated( newfacc ) ) deallocate( newfacc )
        if( allocated( oldfacr ) ) deallocate( oldfacr )
        if( allocated( newfacr ) ) deallocate( newfacr )
        if( allocated( weight ) ) deallocate( weight )

        return

    end subroutine badlands_finalise
    ! ============================================================================

end module run_phase
! ============================================================================
