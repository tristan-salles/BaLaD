! ============================================================================
! Name        : ChannelFWs.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file ChannelFWs.f90
!
! Description :  This module computes river/turbidity flows based on shallow water equation.
!
!<
! ============================================================================
!> Module flow_swe
!<
module flow_swe

    use tin
    use stream
    use fwpath
    use parallel
    use morpho
    use interpol
    use tin_data
    use diffusion
    use file_data
    use ode_data
    use flow_data
    use time_data
    use mesh_data
    use forces_data
    use stream_erodep
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

contains

    ! ============================================================================
    !> Subroutine perform_rivers_evolution
    !! Subroutine perform_rivers_evolution introduces rivers by allocating stream flow walkers into
    !! the simulation area.
    !<
    ! ============================================================================
    subroutine perform_rivers_evolution

        integer :: k, p, ks, seed, srcID, nlay, layerid

        real( tkind ) :: slsed, pt( 2 ), randomVal, ranges, fv, slope, volsed

        if( num_src == 0 .or. transport_mode > 0 ) return

        ! Allocate stream flow walkers arrays
        if( .not. allocated( fr_sedcharge ) ) allocate( fr_sedcharge( totgrn ) )

        srcID = 0
        ! Loop over the river sources
        do p = 1, nbfw_src

            ! Get the top layer thickness
            do k = 1, lnbPts
                nlay = v_nbLays( k )
                layerid = v_LaysID( k, nlay )
                do ks = 1, totgrn
                    if( layID == layerid )then
                        top_sedh( k, ks ) = stratLays( nlay, k, ks )
                        top_sedprev( k, ks ) = stratLays( nlay, k, ks )
                    else
                        top_sedh( k, ks ) = 0.0_8
                        top_sedprev( k, ks ) = 0.0_8
                    endif
                enddo
            enddo

            do k = 1, num_src

                ! Find the ones active for the current time step
                if( tnow >= fws_tstrt( k ) .and. tnow < fws_tend( k ) )then

                    srcID = srcID + 1
                    fr_recface = -1

                    ! Move randomly the river position based on defined range
                    if( iam == 0 )then
                        call random_seed( size = seed )
                        call random_number( harvest = randomVal )
                    endif
                    call mpi_bcast( randomVal, 1, dbl_type, 0, badlands_comm_world, ierr )

                    ranges = 0.0_8
                    pt( 1 ) = fws_xposition( k ) + fws_xrange( k ) * ( 2.0_8 * randomVal - 1.0_8 )
                    pt( 2 ) = fws_yposition( k ) + fws_yrange( k ) * ( 2.0_8 * randomVal - 1.0_8 )
                    if( pt( 1 ) <= strat_xm .and. pt( 1 ) >= strat_xo .and. &
                        pt( 2 ) <= strat_ym .and. pt( 2 ) >= strat_yo )then

                        fr_xpos = pt( 1 )
                        fr_ypos = pt( 2 )
                        call search_all_TIN_face( pt, fr_face )

                        ! Calculate the elevation and local slope of the flow walker
                        call find_elevation_slope_TIN
                        if( fr_zpos < gsea%actual_sea .and. fws_type( k ) == 0 ) goto 30

                        ! Calculate the elevation in regards to stratal grid
                        call find_Mesh_elevation
                        if( fr_zstrat < gsea%actual_sea .and. fws_type( k ) == 0 ) goto 30

                        ! Get closest slope value
                    	call find_closest_slope_value( pt, slope )
                        glb_slp = slope

                    	! Compute channel width
                        fr_q = fws_qfl( k )
                    	fr_w = kwidth * fr_q**awidth * slope**bwidth
                    	if( fr_w < 0.1_8 ) fr_w = 0.1_8

                        ! Add a random component to the velocities
                        if( iam == 0 )then
                            call random_seed( size = seed )
                            call random_number( harvest = randomVal )
                        endif
                        call mpi_bcast( randomVal, 1, dbl_type, 0, badlands_comm_world, ierr )
                        fr_xv = fws_xvel( k ) * ( 1.0_8 + 0.5_8 * randomVal )
                        if( iam == 0 )then
                            call random_seed( size = seed )
                            call random_number( harvest = randomVal )
                        endif
                        call mpi_bcast( randomVal, 1, dbl_type, 0, badlands_comm_world, ierr )
                        fr_yv = fws_yvel( k ) * ( 1.0_8 + 0.5_8 * randomVal )
                        fv = sqrt( fr_xv**2 + fr_yv**2 )
                        fr_vol0 = fws_volume( k )

                        ! Define sediment properties and flow density
                        slsed = 0.0_8
                        volsed = 0.0_8
                        fr_sedcharge = 0.0_8
                        do ks = 1, silgrn
                            volsed = volsed + fr_sedcharge( ks )
                            fr_sedcharge( ks ) = fws_sedcharge( k, ks )
                            slsed = slsed + fr_sedcharge( ks ) * sediment( ks )%density
                        enddo
                        fr_vol = fr_vol0 + volsed
                        fr_density = slsed / fr_vol + fluid_density

                        ! Find flow height recursively based on Manning coefficient
                        fr_h = fr_q / ( fv * fr_w )
                        call hflow_recursive_search

                        ! Advance stream
                        call compute_river_flow_walkers( srcID )

                    endif

                endif

30 continue

            enddo

            ! Diffuse freshly deposited sediments over the topography
            call soil_diffusion_transport

            ! Update TIN surface
            if( p < nbfw_src )  call update_TIN_surface

        enddo

        return

    end subroutine perform_rivers_evolution
    ! ============================================================================
    !> Subroutine hflow_recursive_search
    !! Finds flow depth (considering a rectangular shape channel) using an iterative process.
    !<
    ! ============================================================================
    subroutine hflow_recursive_search

        real( tkind ) :: coeff_a, xnn

        ! Calculate the Manning coefficient of the FW.
        xnn = manning_fct( gsea%actual_sea - fr_zstrat, fr_density )

        coeff_a = xnn * fr_q / ( sqrt( glb_slp ) * fr_w )
        fr_h = coeff_a**( 3.0_8 / 5.0_8 )

        return

    end subroutine hflow_recursive_search
    ! ============================================================================
    !> Subroutine compute_river_flow_walkers
    !! This subroutine solves river flow walkers trajectories based on shallow water equation.
    !<
    ! ============================================================================
    subroutine compute_river_flow_walkers( srcID )

        integer :: id, pp, tt, k, ks, n, srcID, step

        real( tkind ) :: dh, fv

        step = 0
        id = 0
        tt = 0
        fr_loop: do

            if( fr_xpos < strat_xo .or. fr_xpos > strat_xm ) exit fr_loop
            if( fr_ypos < strat_yo .or. fr_ypos > strat_ym ) exit fr_loop

            if( step == 1 .and. iam == 0 .and. fwsvis == 1 ) call record_flow_walkers( srcID )
            if( step == 1 )then
                if( id > 50000 ) id = 0
                id = id + 1
                fr_recface( id ) = fr_face
                if( id > 5000 )then
                    n = 0
                    face_count: do k = 1, 50000
                        if( fr_recface( id ) == fr_face ) n = n + 1
                        if( fr_recface( id ) == -1 ) exit face_count
                    enddo face_count
                    if( n > 500 ) tt = 1
                endif
            endif
            step = step + 1
            if( step > 5 ) step = 1

            ! Update the acceleration vectors of the local flow walkers.
            call calculate_acceleration

            ! Apply Cash-Karp Runge Kutta method, and update timestep parameters.
            ! This also updates the FWs' acceleration and velocity components.
            call advance_flow_runge_kutta

            if( fr_xpos < strat_xo .or. fr_xpos > strat_xm ) exit fr_loop
            if( fr_ypos < strat_yo .or. fr_ypos > strat_ym ) exit fr_loop

            ! Apply deposition/erosion rules
            call apply_erosion_deposition_rules

            ! Update current time
            fr_time = fr_time + fr_dt

            ! Force deposition
            if( fr_zstrat < gsea%actual_sea .and. fr_density < sea_density )then
                call force_deposit_flow_walkers_sediment
                exit fr_loop
            endif

            ! Limit evolution based on velocity
            if( sqrt( fr_xv**2.0_8 + fr_yv**2.0_8 ) < 1.e-3 ) pp = pp + 1
            if( sqrt( fr_xv**2.0_8 + fr_yv**2.0_8 ) > 1.e-3 ) pp = 0

            ! In case there is no more sediment to transport
            dh = 0.0_8
            do ks = 1, totgrn
                dh = dh + fr_sedcharge( ks ) / strat_dx**2
            enddo
            if( dh < 1.e-3 ) pp = 500

            ! Force deposition
            fv = sqrt( fr_xv**2 + fr_yv**2 )
            if( pp >= 500 .or. tt == 1 .or. fv < 1.e-6 .or. fr_h < 1.e-6 .or. fr_q < 1.e-6 )then
               call force_deposit_flow_walkers_sediment
               exit fr_loop
            endif

        enddo fr_loop

        return

    end subroutine compute_river_flow_walkers
    ! ============================================================================
    !> Subroutine advance_flow_runge_kutta
    !! Subroutine advance_flow_runge_kutta repeatedly applies a Runge-Kutta solver,
    !! and adapts the timestep width according to observed error. The subroutine
    !! also handles timestep evolution for the plume case.
    !<
    ! ============================================================================
    subroutine advance_flow_runge_kutta

        real(tkind) :: large_dt
        real( tkind ) :: sc1, sc2, fv

        ! Initialise delta-t values.
        large_dt = 1.0_8 * secyear
        if( fr_time == 0.0_8 ) fr_dtnext = large_dt
        fr_dt = 5.0e2_8

        ! Run adaptive Runge-Kutta solver and get the next timestep width.
        fo_dt = fr_dt
        fo_xv = fr_xv
        fo_yv = fr_yv
        fo_Cfric = fr_Cfric
        fo_accx = fr_accx
        fo_accy = fr_accy
        fr_dtnext = runge_kutta( )
        fr_dt = fo_dt
        fr_xv_new = fo_xvnew
        fr_yv_new = fo_yvnew
        if( fr_dt > large_dt )then
            sc1 = ( fr_dt - large_dt ) / fr_dt
            sc2 = large_dt / fr_dt
            fr_dt = large_dt
            fr_xv_new = fr_xv_new * sc2 + fr_xv * sc1
            fr_yv_new = fr_yv_new * sc2 + fr_yv * sc1
        endif

        ! Using new delta-t and updated velocity, find the new position of the flow walker
        fr_xpos = fr_xpos + fr_xv_new * fr_dt
        fr_ypos = fr_ypos + fr_yv_new * fr_dt

        ! Using the new location find elevation and position on the TIN
        call calculate_local_slope_elevation
        
        ! Calculate the elevation in regards to stratal grid
        call find_Mesh_elevation

        fv = sqrt( fr_xv_new**2.0_8 + fr_yv_new**2.0_8 )

		! Compute channel width following Attal 2008
        fr_w = kwidth * fr_q**awidth *  glb_slp**bwidth
        if( fr_w < 0.1_8 ) fr_w = 0.1_8
        
        ! Find flow height recursively based on Manning coefficient
!        fr_h = fr_q / ( fv * fr_w )
        call hflow_recursive_search

        fr_q = fv * fr_h * fr_w

        call water_entrainment( fv )

        ! Finally update flow velocity
        fr_xv = fr_xv_new
        fr_yv = fr_yv_new

        return

    end subroutine advance_flow_runge_kutta
    ! ============================================================================
    !> Subroutine calculate_acceleration
    !! Subroutine calculate_acceleration considers gravity, Manning coefficient,
    !! and plume status to calculate the x/y components of each local flow walker's
    !! acceleration vector.
    !<
    ! ============================================================================
    subroutine calculate_acceleration

        real( tkind ) :: xnn, ga

        ! Calculate the Manning coefficient of the FW.
        xnn = manning_fct( gsea%actual_sea - fr_zstrat, &
            fr_density )

        ! Calculate the effect of gravity, scaled by the ratio of stream density
        ! to seawater density.
        ga = gravity
        if( fr_zstrat <=  gsea%actual_sea ) &
            ga = ga * ( fr_density - sea_density ) / fr_density

        ! Calculate the FW's coefficient of friction.
        fr_Cfric = gravity * xnn**2.0_8 / fr_h**(4.0_8/3.0_8)

        ! If the FW represents a flow above sea-level, or if the FW has
        ! density greater than or equal to sea-water density.
        if( fr_zstrat > gsea%actual_sea .or. ( fr_zstrat <= gsea%actual_sea &
            .and. fr_density >= sea_density ) ) then
            call slope_SF_direction

        ! Otherwise, assume zero slope.
        else
            fr_xslp = 0.0_8
            fr_yslp = 0.0_8
        endif

        ! Set the acceleration components according to slope and gravity.
        fr_accx = -ga * fr_xslp
        fr_accy = -ga * fr_yslp

        return

    end subroutine calculate_acceleration
    ! ============================================================================
    !> Subroutine water_entrainment
    !! Subroutine water_entrainment computes water entrainment for turbidity currents.
    !<
    ! ============================================================================
    subroutine water_entrainment( fv2 )

        real( tkind ) :: fv2, Ri, ew

        ! Water entrainment for turbidity currents
        if( sea_density < fr_density .and. fv2 > 0.0_8 )then
            Ri = ( fr_density - sea_density ) * gravity / fr_density
            Ri = Ri * fr_h / fv2**2
            ew = 1.0_8 + 718.0_8 * Ri**2.4_8
            ew = sqrt( ew )
            ew = 0.075_8 / ew
            ew = ew * fr_dt * sqrt( fv2 )

            fr_vol0 = fr_vol0 + ew * strat_dx**2.0_8
        endif

        return

    end subroutine water_entrainment
    ! ============================================================================
    !> Subroutine deposit_flow_walkers_sediment
    !! Searches a kdtree to allocate deposition.
    !<
    ! ============================================================================
    subroutine force_deposit_flow_walkers_sediment

        integer :: nid, k, ks, dFcS, mFcS
        integer, dimension( 4 ) :: idM, idD

        real( tkind ) :: sumdist, sumwgth, sumdep, th, sumsed, dis
        real( tkind ), dimension( 4 ) :: dist, wgth
        real( tkind ), dimension( 4, totgrn ) :: dep


        ! Find DEM face containing the flow walker
        call check_DEMface_containing_point( fr_xpos, fr_ypos, dFcS )
        idD( 1:4 ) = fptIDs( dFcS, 1:4 )

        ! Find stratigraphic face containing the flow walker
        call check_Meshface_containing_point( fr_xpos, fr_ypos, mFcS )
        idM( 1:4 ) = lfptIDs( mFcS, 1:4 )

        dep = 0.0_8
        sumdist = 0.0_8
        ngb_loop: do k = 1, 4
            dis = ( fr_xpos - lcoordX( idM( k ) ) )**2 + ( fr_ypos - lcoordY( idM( k ) ) )**2
            if( dis == 0.0_8 )then
                dep = 0.0_8
                do ks = 1, totgrn
                    dep( k, ks ) = fr_sedcharge( ks ) / strat_dx**2
                enddo
                goto 45
            endif
            dist( k ) = 1 / dis
            sumdist = sumdist + dist( k )
        enddo ngb_loop

        ! Distribute deposits based on inverse weighted averaged distance
        sumwgth = 0.0_8
        do k = 1, 4
            wgth( k ) = dist( k ) / sumdist
            sumwgth = wgth( k ) + sumwgth
        enddo

        sumdep = 0.0_8
        do k = 1, 4
            sumsed = 0.0_8
            wgth( k ) = wgth( k ) / sumwgth
            do ks = 1, totgrn
                dep( k, ks ) = wgth( k ) * fr_sedcharge( ks ) / strat_dx**2
                sumdep = sumdep + dep( k, ks )
                sumsed = sumsed + fr_sedcharge( ks ) / strat_dx**2
            enddo
        enddo
        if( abs( sumdep - sumsed ) > tor )then
          print*,'Problem when forcing deposition',sumdep,sumsed,iam
          stop
        endif

45 continue

        do k = 1, 4
            nid = idM( k )
            th = 0.0_8
            if( layID == v_LaysID( nid, v_nbLays( nid ) ) )then
                do ks = 1, totgrn
                    if( gporo%compaction )then
                        dep( k, ks ) = dep( k, ks ) / ( 1.0_8 - porosity( ks, 1 ) )
                        porosityLays( v_nbLays( nid ), nid, ks ) = porosity( ks, 1 )
                    endif
                    stratLays( v_nbLays( nid ), nid, ks ) = stratLays( v_nbLays( nid ), nid, ks ) + dep( k, ks )
                    th = th + dep( k, ks )
                enddo
                lcoordZ( nid ) = lcoordZ( nid ) + th
                coordZ( idD( k ) ) = lcoordZ( nid )

                soil_thick( nid ) = soil_thick( nid ) + th

               if( soil_thick( nid ) > dtb_marine .and. lcoordZ( nid ) < gsea%actual_sea ) &
                    soil_thick( nid ) = dtb_marine
               if( soil_thick( nid ) > dtb_aerial .and. lcoordZ( nid ) >= gsea%actual_sea ) &
                    soil_thick( nid ) = dtb_aerial

            else
                v_nbLays( nid ) = v_nbLays( nid ) + 1
                v_LaysID( nid, v_nbLays( nid ) ) = layID
                do ks = 1, totgrn
                    if( gporo%compaction )then
                        dep( k, ks ) = dep( k, ks ) / ( 1.0_8 - porosity( ks, 1 ) )
                        porosityLays( v_nbLays( nid ), nid, ks ) = porosity( ks, 1 )
                    endif
                    stratLays( v_nbLays( nid ), nid, ks ) = dep( k, ks )
                    th = th + dep( k, ks )
                enddo
                lcoordZ( nid ) = lcoordZ( nid ) + th
                coordZ( idD( k ) ) = lcoordZ( nid )

                soil_thick( nid ) = soil_thick( nid ) + th

               if( soil_thick( nid ) > dtb_marine .and. lcoordZ( nid ) < gsea%actual_sea ) &
                    soil_thick( nid ) = dtb_marine
               if( soil_thick( nid ) > dtb_aerial .and. lcoordZ( nid ) >= gsea%actual_sea ) &
                    soil_thick( nid ) = dtb_aerial

           endif
        enddo

        return

    end subroutine force_deposit_flow_walkers_sediment
    ! ============================================================================

end module flow_swe
! ============================================================================















