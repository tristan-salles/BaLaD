! ============================================================================
! Name        : SteadyState.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file SteadyState.f90
!
! Description :  This module computes morphometrics based on steady state.
!
!<
! ============================================================================
!> Module Steady_state 
!<
module steady

    use tin
    use parallel
    use morpho
    use interpol
    use tin_data
    use file_data
    use flow_algo
    use time_data
    use flow_data
    use mesh_data
    use forces_data

    implicit none

contains

    ! ============================================================================
    !> Subroutine compute_steady_state
    !! This subroutine is the general function encapsulating subroutines used to perform steady state.
    !<
    ! ============================================================================
    subroutine compute_steady_state

        ! Load rain map
        call find_current_rain_map

        ! Depression fill algorithm
        call planchon_dem_fill_algorithm

        ! Compute flow accumulation grid
        call compute_MFD

        ! Prepare for elevation Lagrangian interpolation based on filled DEM
        call lagrangian_topofield_declaration
        
        ! Update DEM parameters
        call update_DEM_parameters

        ! Compute morphometrics based on circular neighborhood statistic
        call circular_neighborhood_statistic
        call mpi_allreduce( mpi_in_place, slp, nbPts, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, hcurv, nbPts, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, vcurv, nbPts, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, orient, nbPts, dbl_type, max_type, badlands_comm_world, ierr )
        
        ! Compute velocity fields and hydrological parameters
        call relief_hydrological_features
        call mpi_allreduce( mpi_in_place, qflow, nbPts, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, wflow, nbPts, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, hflow, nbPts, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, vflow, nbPts, dbl_type, max_type, badlands_comm_world, ierr )

        ! Prepare for velocity field Lagrangian interpolation based on filled DEM
        call lagrangian_velocityfield_declaration

        ! Prepare for flow rate Lagrangian interpolation based on filled DEM
        call lagrangian_qflowfield_declaration

        ! Prepare for flow height Lagrangian interpolation based on filled DEM
        call lagrangian_flowheightfield_declaration

        ! Prepare for flow width Lagrangian interpolation based on filled DEM
        call lagrangian_flowwidthfield_declaration

        if( transport_mode >= 1 .and. num_src == 0 ) return

        ! Update TIN
        call update_TIN_surface

        return

    end subroutine compute_steady_state
    ! ============================================================================
    !> Subroutine find_current_rain_map
    !! This subroutine checks if a rain map needs to be loaded.
    !<
    ! ============================================================================
    subroutine find_current_rain_map

        integer :: n, event_nb

        ! Find the rain event if it exists and update it if required
        event_nb = 0
        rain_time = 24.0_8
        find_event: do n = 1, rain_event
            if( rain_tstart( n ) <= tnow .and. rain_tend( n ) >= tnow )then
                event_nb = n
                rain_time = rain_duration( n )
                exit find_event
            endif
        enddo find_event

        ! Is there a new map to load
        rainID = event_nb
        call loading_new_rain_map

        return

    end subroutine find_current_rain_map
    ! ============================================================================
    !> Subroutine loading_new_rain_map
    !! This subroutine loads a new rain map if required and write it as a ASCII file.
    !! Rain maps are in metres per year
    !<
    ! ============================================================================
    subroutine loading_new_rain_map

        integer :: p, iu, ios

        if( .not. allocated( precipitation ) ) allocate( precipitation( lnbPts ) )

        if( rainID > 0 )then
            iu = 19
            open( iu, file=frainmap( rainID ), status="old", action="read", iostat=ios )
            do p = 1, lnbPts
                read( iu, * ) precipitation( p )
                if( precipitation( p ) > 0.0_8 .and. lcoordZ( p ) < gsea%actual_sea ) &
                    precipitation( p ) = 0.0_8
            enddo
            close( iu )
        else
            rain_time = 0.0_8
            precipitation = 0.0_8
        endif

        return

    end subroutine loading_new_rain_map
    ! ============================================================================
    !> Subroutine update_DEM_parameters
    !! From the stratigraphic grid define DEM morphometrics.
    !<
    ! ============================================================================
    subroutine update_DEM_parameters

        integer :: k

        do k = 1, nbPts
            coordZ( k ) = lcoordZ( lvertexID( k ) )
            filldem( k ) = lfilldem( lvertexID( k ) )
            facc( k ) = lfacc( lvertexID( k ) )
            facr( k ) = lfacr( lvertexID( k ) )
        enddo

        return

    end subroutine update_DEM_parameters
    ! ============================================================================
    !> Subroutine relief_hydrological_features
    !! Based on steady state computes hydraulic parameters.
    !<
    ! ============================================================================
    subroutine relief_hydrological_features

        integer :: gid, p

        real( tkind ) :: q, width, coeff_a

        ! Check that arrays for curvature and slope have been initialised
        if( .not. allocated( vflow ) ) allocate( vflow( nbPts ) )
        if( .not. allocated( hflow ) ) allocate( hflow( nbPts ) )
        if( .not. allocated( qflow ) ) allocate( qflow( nbPts ) )
        if( .not. allocated( wflow ) ) allocate( wflow( nbPts ) )

        ! Initialise hydraulic parameters
        qflow = 0.0_8
        wflow = 0.0_8
        hflow = 0.0_8
        vflow = 0.0_8
        slpmax = 0.0_8
        slpmin = 100.0_8

        do p = 1, local_nbPts

            gid = global_nid( p )

            ! Rain accumulation
            if( facr( gid ) < 0.0_8 ) facr( gid ) = 0.0_8
            if( rain_time == 0.0_8 ) facr( gid ) = 0.0_8

            ! Compute the flow rate based on flow accumulation (cubic metres per second)
            q = facr( gid ) * strat_dx**2.0_8 / ( secyear )
            if( q < 1.e-8_8 ) q = 0.0_8

            slpmin = min( slpmin, slp( gid ) )
            slpmax = max( slpmax, slp( gid ) )

            if( q > 0.0_8 )then
                ! Compute channel width following Attal 2008
                width = kwidth * q**awidth *  (slp( gid ))**bwidth
                if( width < 0.1_8 ) width = 0.1_8

                ! Get stream flow depth
                coeff_a = manning_open * q / ( sqrt( slp( gid ) ) * width )
                hflow( gid ) = coeff_a**( 3.0_8 / 5.0_8 )

                if( width > 0.0_8 .and. hflow( gid ) > 0.0_8 )then
                    wflow( gid ) = width
                    qflow( gid ) = q
                    ! Velocity for open-channel uniform steady flow based on Manning equation
                    vflow( gid ) = q / ( hflow( gid ) * width )
                endif
            endif

        enddo

        call mpi_allreduce( mpi_in_place, slpmax, 1, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, slpmin, 1, dbl_type, min_type, badlands_comm_world, ierr )

        return

    end subroutine relief_hydrological_features
    ! ============================================================================

end module steady
! ============================================================================
