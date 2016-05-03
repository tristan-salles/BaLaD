! ============================================================================
! Name        : OverlandFWs.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file OverlandFWs.f90
!
! Description :  This module computes overland flows based on shallow water equation.
!
!<
! ============================================================================
!> Module fwpath
!<
module fwpath

    use stream
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

    ! Define the new deposits from flow
    real( tkind ), dimension(:,:), allocatable :: deposit

    ! Define the new deposits from flow
    real( tkind ), dimension(:,:), allocatable :: mergedeposit

    ! New topographic values
    real( tkind ), dimension( : ), allocatable :: newz

contains

    ! ============================================================================
    !> Subroutine compute_overland_flow_walkers
    !! This subroutine solves flow walkers trajectories based on shallow water equation.
    !<
    ! ============================================================================
    subroutine compute_overland_flow_walkers

        integer :: k, ks

        real( tkind ) :: th

        if( .not. allocated( deposit ) ) allocate( deposit( totgrn, lnbPts ) )
        if( .not. allocated( mergedeposit ) ) allocate( mergedeposit( totgrn, lnbPts ) )
        deposit = 0.0_8
     
        ! Define flow walkers for considered time step and compute trajectories
        call advance_overland_flow_walkers

        ! If purely erosive or diffusive-dependent mode no need to perform
        ! the erosion / deposition flow walkers transport
        if( transport_mode == 1 ) return

        ! Merge all deposit thicknesses
        do ks = 1, totgrn
            call mpi_allreduce( deposit( ks, 1:lnbPts ), mergedeposit( ks, 1:lnbPts ), lnbPts, dbl_type, &
                sum_type, badlands_comm_world, ierr )
        enddo

        ! Update stratigraphic layer
        do k = 1, lnbPts
            ! In case the layer already exists
            th = 0.0_8
            if( layID == v_LaysID( k, v_nbLays( k ) ) )then
                do ks = 1, totgrn
                    ! Prepare diffusion arrays
                    top_sedh( k, ks ) = stratLays( v_nbLays( k ), k, ks )
                    top_sedprev( k, ks ) = stratLays( v_nbLays( k ), k, ks )
                    ! Compute sedimentary changes
                    if( gporo%compaction )then
                        mergedeposit( ks, k ) = mergedeposit( ks, k ) / ( 1.0_8 - porosity( ks, 1 ) )
                        porosityLays( v_nbLays( k ), k, ks ) = porosity( ks, 1 )
                    endif
                    stratLays( v_nbLays( k ), k, ks ) = stratLays( v_nbLays( k ), k, ks ) + mergedeposit( ks, k )
                    th = th + mergedeposit( ks, k )
                enddo
                lcoordZ( k ) = lcoordZ( k ) + th
                soil_thick( k ) = soil_thick( k ) + th
                if( soil_thick( k ) > dtb_marine .and. lcoordZ( k ) < gsea%actual_sea ) &
                    soil_thick( k ) = dtb_marine
                if( soil_thick( k ) > dtb_aerial .and. lcoordZ( k ) >= gsea%actual_sea ) &
                    soil_thick( k ) = dtb_aerial
            else
                do ks = 1, totgrn
                    ! Prepare diffusion arrays
                    top_sedh( k, ks ) = 0.0_8
                    top_sedprev( k, ks ) = 0.0_8
                    ! Compute sedimentary changes
                    th = th + mergedeposit( ks, k )
                enddo
                if( th > 0.0_8 )then
                    v_nbLays( k ) = v_nbLays( k ) + 1
                    v_LaysID( k, v_nbLays( k ) ) = layID
                    do ks = 1, totgrn
                        if( gporo%compaction )then
                            mergedeposit( ks, k ) = mergedeposit( ks, k ) / ( 1.0_8 - porosity( ks, 1 ) )
                            porosityLays( v_nbLays( k ), k, ks ) = porosity( ks, 1 )
                        endif
                        stratLays( v_nbLays( k ), k, ks ) = mergedeposit( ks, k )
                    enddo
                    lcoordZ( k ) = lcoordZ( k ) + th
                    soil_thick( k ) = soil_thick( k ) + th
                    if( soil_thick( k ) > dtb_marine .and. lcoordZ( k ) < gsea%actual_sea ) &
                        soil_thick( k ) = dtb_marine
                    if( soil_thick( k ) > dtb_aerial .and. lcoordZ( k ) >= gsea%actual_sea ) &
                        soil_thick( k ) = dtb_aerial
                endif
            endif

        enddo

        ! Update DEM grid topography
        do k = 1, nbPts
            coordZ( k ) = lcoordZ( lvertexID( k ) )
        enddo

31 continue

        ! Diffuse freshly deposited sediments over the topography
        call soil_diffusion_transport

        return

    end subroutine compute_overland_flow_walkers
    ! ============================================================================
    !> Subroutine advance_overland_flow_walkers
    !! This subroutine defines the number and parameters of the active overland flow walkers
    !! and advances them through time.
    !<
    ! ============================================================================
    subroutine advance_overland_flow_walkers

        integer :: k, ks, id, p, nlay, pp, step, mxstp, n

        integer, dimension( nbPts ) :: gid

        real( tkind ) :: stillsed, dh, pt( 2 ), eromax
        real( tkind ) :: totero, dtb, frac, discharge, minz

        real( tkind ), dimension( nbPts, totgrn ) :: qs
        real( tkind ), dimension( totgrn ) :: sedh, remain, erode, hs, poro, prop !, prop2

        maxFWs = 0
        gid = -1

        if( .not. allocated( newz ) ) allocate( newz( lnbPts ) )
        newz = lcoordZ

        do k = 1, nbPts

            qs( k, 1:totgrn ) = 0.0_8

            if( vflow( k ) > 0.0_8 .and. coordZ( k ) > gsea%actual_sea )then
                
                ! Get soil cover composition
!                dtb = soil_thick( id )
!                call get_regolith_composition(  lvertexID( k ), dtb, dh, hs, poro )
!                if( dh == 0.0_8 )then
!                    prop( 1 ) = 1.0
!                else
!                    do ks = 1, totgrn
!                        prop( ks ) = hs( ks ) / dh
!                    enddo
!                endif
!                dh = 0.99_8 * dh
!                hs = 0.99_8 * hs

                ! Get maximum erosion thickness to prevent formation of hole
                call get_erosion_thickness( k, eromax )
                eromax = eromax * 0.95_8

                call get_regolith_composition(  lvertexID( k ), eromax, dh, hs, poro )
                if( dh == 0.0_8 )then
                    prop( 1 ) = 1.0
                else
                    do ks = 1, totgrn
                        prop( ks ) = hs( ks ) / dh
                    enddo
                endif

                if( eromax == 0.0_8 )then
                    prop = 0.0_8
                    dh = 0.0_8
                    hs = 0.0_8
                elseif( dh > eromax )then
                    do ks = 1, totgrn
                        hs( ks ) = hs( ks ) * eromax / dh
                    enddo
                    dh = eromax
                endif
                dh = 0.99_8 * dh
                hs = 0.99_8 * hs

                ! In case the depth to bedrock is greater than the maximum erosion change it to match
                if( dtb > eromax .or. dtb == 0.0_8 ) dtb = eromax

                ! Fluvial discharge
                if( newz( lvertexID( k ) ) >= gsea%actual_sea .and. vflow( k ) > 0.0_8 )then
                    discharge = qflow( k ) * secyear
                else
                    discharge = 0.0_8
                endif

                ! Compute erosion fluxes
                totero = 0.0_8
                qs( k, 1:totgrn ) = 0.0_8
                if( dh > 0 )then

                    do ks = 1, totgrn
                        
                        ! Detachment-limited mode
                        qs( k, ks ) = prop( ks ) * sediment( ks )%ero_coeff * ( discharge )**( mdetach ) * &
                            rain_int * ( slp( k ) )**( ndetach )

                        if( qs( k, ks ) > hs( ks ) ) qs( k, ks ) = hs( ks )

                        ! Update erosion value
                        totero = totero + qs( k, ks )

                    enddo

                endif

                if( ( totero > 1.e-6_8 .and. transport_mode == 1 ) .or. ( totero > erosion_trigger .and. transport_mode == 0 ) )then

                    maxFWs = maxFWs + 1

                    gid( maxFWs ) = k
                    id = lvertexID( k )
                    hs( 1:totgrn ) = qs( k, 1:totgrn )

                    ! Erode stratigraphic layer
                    nlay = v_nbLays( id )

                    erosion_loop: do p = v_nbLays( id ), 1, -1
                        stillsed = 0.0_8
                        remain = 0.0_8
                        erode = 0.0_8
                        sedh = 0.0_8
                        do ks = 1, totgrn

                            ! Get the amount to erode for the considered layer
                            if( gporo%compaction )then
                                sedh( ks ) = stratLays( p, id, ks ) * ( 1.0_8 - porosityLays( p, id, ks ) )
                                sedh( ks ) = sedh( ks ) - hs( ks )
                            else
                                sedh( ks ) = stratLays( p, id, ks ) - hs( ks )
                            endif

                            ! In case there isn't enough sediment in the layer
                            ! store the sediment that needs to be eroded in
                            ! underlying layers
                            if( sedh( ks ) <= 0.0_8 )then
                                sedh( ks ) = 0.0_8
                                if( gporo%compaction )then
                                    erode( ks ) = stratLays( p, id, ks ) * ( 1.0_8 - porosityLays( p, id, ks ) )
                                    porosityLays( p, id, ks ) = 0.0_8
                                else
                                    erode( ks ) = stratLays( p, id, ks )
                                endif
                                remain( ks ) = 0.0_8
                                hs( ks ) = hs( ks ) - erode( ks )

                            ! In case there is still sediments update the remaining
                            ! thickness accordingly
                            else
                                if( gporo%compaction )then
                                    remain( ks ) = sedh( ks ) / ( 1.0_8 - porosityLays( p, id, ks ) )
                                else
                                    remain( ks ) =  sedh( ks )
                                endif
                                hs( ks ) = 0.0_8
                            endif

                            ! Get the total thickness of sediment present in
                            ! the current layer
                            stillsed = stillsed + remain( ks )
                        enddo

                        ! If there is still some sediment on the layer
                        ! stop the erosion
                        if( stillsed > 0.0_8 )then
                            dh = 0.0_8
                            do ks = 1, totgrn
                                if( stratLays( p, id, ks ) < remain( ks ) )then
                                    print*,'Something went wrong when eroding stratigraphy for overland FWs.',iam
                                    print*,ks,id,p,remain( ks ),v_nbLays( id ),stratLays( p, id, 1:totgrn ),qs( k,ks )
                                    stop
                                endif
                                dh = dh + stratLays( p, id, ks ) - remain( ks )
                                stratLays( p, id, ks ) = remain( ks )
                                qs( k, ks ) = qs( k, ks ) - hs( ks )
                            enddo
                            lcoordZ( id ) = lcoordZ( id ) - dh
                            coordZ( k ) = lcoordZ( id )
                            if( transport_mode == 0 )then
                                soil_thick( id ) = soil_thick( id ) - dh
                                if( soil_thick( id ) < 0.0_8 ) soil_thick( id ) = 0.0_8
                            endif
                            exit erosion_loop
                        endif

                        ! If there is no more sediment
                        if( stillsed == 0.0_8 )then
                            dh = 0.0_8
                            do ks = 1, totgrn
                                dh = dh + erode( ks )
                                stratLays( p, id, ks ) = 0.0_8
                                if( gporo%compaction ) porosityLays( p, id, ks ) = 0.0_8
                            enddo
                            lcoordZ( id ) = lcoordZ( id ) - dh
                            coordZ( k ) = lcoordZ( id )
                            if( transport_mode == 0 )then
                                soil_thick( id ) = soil_thick( id ) - dh
                                if( soil_thick( id ) < 0.0_8 ) soil_thick( id ) = 0.0_8
                            endif
                            nlay = nlay - 1
                        endif
                    enddo erosion_loop
                    v_nbLays( id ) = nlay

                endif
            endif
        enddo

        if( transport_mode == 1 ) return

        ! Initialise FWs
        if( .not. allocated( fw_hs ) ) allocate( fw_hs( totgrn ) )
        if( .not. allocated( fw_hs0 ) ) allocate( fw_hs0( totgrn ) )

        p = iam + 1
        do k = 1, maxFWs
            id = gid( k )
            if( id > 0 )then

                if( p == k )then

                    ! Define initial rain flow walker parameter

                    ! Position
                    fw_xpos = coordX( id )
                    fw_ypos = coordY( id )
                    pt( 1 ) = coordX( id )
                    pt( 2 ) = coordY( id )

                    call search_all_TIN_face( pt, fw_face )

                    ! Calculate the elevation and local slope of the flow walker
                    call find_elevation_slope_TIN

                    ! Calculate the elevation in regards to stratal grid
                    call find_Rain_Mesh_elevation

                    ! Velocity
                    fw_xv = vflow( id ) * sin( orient( id ) )
                    fw_yv = vflow( id ) * cos( orient( id ) )

                    ! Hydraulics
                    fw_h = hflow( id )
                    fw_q = qflow( id )

                    ! Sediment charge (metres of sediment transported)
                    fw_hs0( 1:totgrn ) = qs( id, 1:totgrn )
                    fw_hs( 1:totgrn ) = 0.0_8

                    ! Solve flow walkers trajectories based on shallow water equation.
                    pp = 0
                    step = 0
                    mxstp = 0
                    fw_time = 0.0_8

                    fw_loop: do

                        mxstp = mxstp + 1
                        if( fw_xpos < strat_xo .or. fw_xpos > strat_xm ) exit fw_loop
                        if( fw_ypos < strat_yo .or. fw_ypos > strat_ym ) exit fw_loop

                        ! Force deposition
                        if( fw_zpos < gsea%actual_sea .or. fw_zstrat < gsea%actual_sea )then
                            call deposit_flow_walkers_sediments( fw_hs( 1:totgrn ) )
                            exit fw_loop
                        endif

                        if( step == 10  .and. fwsvis == 1 ) call record_flow_walkers( k )
                        step = step + 1
                        if( step > 100 ) step = 0

                        if( fw_q >= 1.e-3_8 .or. fw_h >= 5.e-3_8 )then
                            call advance_fw_through_time
                        else
                            pp = 101
                        endif

                        if( sqrt( fw_xv**2.0_8 + fw_yv**2.0_8 ) < 1.e-3 ) pp = pp + 1
                        if( sqrt( fw_xv**2.0_8 + fw_yv**2.0_8 ) > 1.e-3 ) pp = 0

                        ! Force deposition
                        if( pp > 100 .or. fw_q < 1.e-3_8 .or. fw_h < 5.e-3_8  .or. mxstp > 1e5 )then
                            call deposit_flow_walkers_sediments( fw_hs( 1:totgrn ) )
                            exit fw_loop
                        endif

                        ! In case there is no more sediment to transport
                        dh = 0.0_8
                        do ks = 1, totgrn
                            dh = dh + fw_hs( ks )
                        enddo
                        if( dh == 0.0_8 ) exit fw_loop

                    enddo fw_loop

                    ! Go to next flow walker
                    p = nproc + p

                endif
            endif
        enddo

        return

    end subroutine advance_overland_flow_walkers
    ! ============================================================================
    !> Subroutine get_erosion_thickness
    !! Subroutine get_erosion_thickness get the maximum erosion thickness
    !<
    ! ============================================================================
    subroutine get_erosion_thickness( id, eromax )

        integer :: idD, id, n, k

        real( tkind ) :: elev, minz, eromax

        eromax = 0.0_8

        idD = lvertexID( id )

        ! Current layer elevation
        elev = newz( idD )
        minz = elev

        ! Find neighborhood minimum elevation
        do n = 1, 8
            k = lngbID( idD, n )
            if(  k > 0 ) minz = min( minz, newz( k ) )
        enddo

        ! Maximum erosion to prevent formation of holes
        eromax = ( elev - minz )
        if( eromax < 0.0_8 ) eromax = 0.0_8

        if( elev > gsea%actual_sea .and. minz < gsea%actual_sea ) eromax = ( elev - gsea%actual_sea )

        return

    end subroutine get_erosion_thickness
    ! ============================================================================
    !> Subroutine find_Rain_Mesh_elevation
    !! Subroutine find_Rain_Mesh_elevation is used to determine the elevation of a point in
    !!  the mesh grid.
    !<
    ! ============================================================================
    subroutine find_Rain_Mesh_elevation

        integer :: Xid, Yid, fcS, id( 4 ), k

        real( tkind ) :: wght( 4 )

        ! Define the probable row and column number
        Xid = int( ( fw_xpos - lstrat_xo ) / strat_dx ) + 1
        Yid = int( ( fw_ypos - lstrat_yo ) / strat_dx ) + 1

        ! Define the face id accordingly
        if( Yid == 1 .and. Xid < lstrat_X )then
            fcS = Xid
        elseif( Yid == 1 .and. Xid == lstrat_X )then
            fcS = Xid - 1
        elseif( Xid == 1 .and. Yid < lstrat_Y )then
            fcS = ( Yid - 1 )*( lstrat_X - 1 ) + 1
        elseif( Yid == lstrat_Y .and. Xid < lstrat_X )then
            fcS = ( Yid - 2 )*( lstrat_X - 1 ) + Xid
        elseif( Yid == lstrat_Y .and. Xid == lstrat_X )then
            fcS = ( Yid - 1 )*( Xid - 1 )
        elseif( Yid < lstrat_Y .and. Xid == lstrat_X )then
            fcS = Yid * ( Xid - 1 )
        else
            fcS = ( lstrat_X - 1 )* ( Yid - 1 ) + Xid
        endif

        ! Check that picked face is containing the point
        id( 1:4 ) = lfptIDs( fcS, 1:4 )

        ! Weigth
        call inverse_distance_weighting_dem( id, wght )

        fw_zstrat = 0.0_8
        do k = 1, 4
            fw_zstrat = fw_zstrat  + wght( k ) * lcoordZ( id( k ) )
        enddo

        return

    end subroutine find_Rain_Mesh_elevation
    ! ============================================================================
    !> Subroutine inverse_distance_weighting_dem
    !! Subroutine inverse_distance_weighting_dem get weight of the fw on each points
    !<
    ! ============================================================================
    subroutine inverse_distance_weighting_dem( id, weights )

        integer :: k
        integer, dimension( 4 ) :: id

        real( tkind ) :: dis, sumdist, sumwgth
        real( tkind ), dimension( 4 ) :: weights, dist

        weights = 0.0_8
        sumdist = 0.0_8
        face_loop: do k = 1, 4
            dis = ( fw_xpos - lcoordX( id( k ) ) )**2 + ( fw_ypos - lcoordY( id( k ) ) )**2
            if( dis == 0.0_8 )then
                weights = 0.0_8
                weights( k ) = 1.0_8
                return
            endif
            dist( k ) = 1 / dis
            sumdist = sumdist + dist( k )
        enddo face_loop

        sumwgth = 0.0_8
        do k = 1, 4
            weights( k ) = dist( k ) / sumdist
            sumwgth = weights( k ) + sumwgth
        enddo

        do k = 1, 4
            weights( k ) = weights( k ) / sumwgth
        enddo

        sumwgth = 0.0_8
        do k = 1, 4
            sumwgth = weights( k ) + sumwgth
        enddo

        if( ( sumwgth - 1.0_8 ) > tor )then
            print*,'Problem when computing inverse weighting values in rain',sumwgth
            stop
        endif

        return

    end subroutine inverse_distance_weighting_dem
    ! ============================================================================
    !> Subroutine advance_fw_through_time
    !! Flow walkers algorithm.
    !<
    ! ============================================================================
    subroutine advance_fw_through_time

        integer :: ks
        real( tkind ) :: ftime, depo( totgrn ), prevh, trspt, xnn

        ! Calculate the FW's coefficient of friction using manning.
        if( fw_zpos > gsea%actual_sea .or. fw_zstrat > gsea%actual_sea )then
            xnn = manning_open
        else
            xnn = manning_hypo
        endif

        fw_Cfric = gravity * xnn**2.0_8 / fw_h**(4.0_8/3.0_8)

        ! If the FW represents a flow above sea-level, or if the FW has
        ! density greater than or equal to sea-water density.
        if( fw_zpos > gsea%actual_sea .or. fw_zstrat > gsea%actual_sea ) then
            call slope_SF_direction

        ! Otherwise, assume zero slope.
        else
            fw_xslp = 0.0_8
            fw_yslp = 0.0_8
        endif

        ! Set the acceleration components according to slope and gravity.
        fw_accx =  - gravity * fw_xslp
        fw_accy =  - gravity * fw_yslp

        ! Apply Cash-Karp Runge Kutta method, and update timestep parameters.
        ! This also updates the FWs' acceleration and velocity components.
        call advance_overland_runge_kutta

        ! Update current time
        fw_time = fw_time + fw_dt
        ftime = fw_time / 3600.0_8
        depo = 0.0_8
        do ks = 1, totgrn
            trspt = sediment(ks)%transport * rain_time
            prevh = fw_hs( ks )
            if( fw_time - fw_dt == 0.0_8 ) prevh = fw_hs0( ks )
            fw_hs( ks ) = fw_hs0( ks ) * ( 1 - ftime * exp( ftime - trspt ) / trspt )
            if( fw_hs( ks ) < 0.0_8 ) fw_hs( ks ) = 0.0_8
            depo( ks ) = prevh - fw_hs( ks )
            if( depo( ks ) < 0.0_8 )then
                print*,'problem with deposits'
                print*,depo(ks),iam,fw_hs0( ks )-fw_hs( ks ) ,ftime
            endif
        enddo

        call deposit_flow_walkers_sediments( depo )

        return

    end subroutine advance_fw_through_time
    ! ============================================================================
    !> Subroutine advance_overland_runge_kutta
    !! Subroutine advance_overland_runge_kutta repeatedly applies a Runge-Kutta solver,
    !! and adapts the timestep width according to observed error. The subroutine
    !! also handles timestep evolution for the plume case.
    !<
    ! ============================================================================
    subroutine advance_overland_runge_kutta

        real(tkind) :: large_dt, sc1, sc2, npt( 2 ), hh, qq

        ! Initialise delta-t values.
        large_dt = 1.0_8 * secyear
        if( fw_time == 0.0_8 ) fw_dtnext = large_dt
        fw_dt = 5.0e2_8

        ! Run adaptive Runge-Kutta solver and get the next timestep width.
        fo_dt = fw_dt
        fo_xv = fw_xv
        fo_yv = fw_yv
        fo_Cfric = fw_Cfric
        fo_accx = fw_accx
        fo_accy = fw_accy
        hh = fw_h
        qq = fw_q
        fw_dtnext = runge_kutta( )
        fw_dt = fo_dt
        fw_xv_new = fo_xvnew
        fw_yv_new = fo_yvnew

        if( fw_dt > large_dt )then
            sc1 = ( fw_dt - large_dt ) / fw_dt
            sc2 = large_dt / fw_dt
            fw_dt = large_dt
            fw_xv_new = dble( fw_xv_new * sc2 + fw_xv * sc1 )
            fw_yv_new = dble( fw_yv_new * sc2 + fw_yv * sc1 )
        endif
        
        ! Using new delta-t and updated velocity, find the new position of the flow walker
        fw_xpos = fw_xpos + fw_xv_new * fw_dt
        fw_ypos = fw_ypos + fw_yv_new * fw_dt
        
        ! Update flow height
        npt( 1 ) = fw_xpos
        npt( 2 ) = fw_ypos
        call get_lagrangian_interpolated_flowheightfield_values( npt, fw_h )
        if( fw_h < 5.e-3_8 ) fw_h = 0.0_8 !1.e-3

        ! Update flow discharge
        call get_lagrangian_interpolated_qflowfield_values( npt, fw_q )
        if( fw_q < 0 )then
            call search_DEM_kdtree_flowdischarge( npt, fw_q )
            if( fw_q < 1.e-3_8 ) fw_q = 0.0_8 !1.e-5
        endif

        ! Using the new location find elevation and position on the TIN
        call calculate_local_slope_elevation

        ! Calculate the elevation in regards to stratal grid
        call find_Rain_Mesh_elevation

        ! Update flow velocity
        fw_xv = fw_xv_new
        fw_yv = fw_yv_new

        return

    end subroutine advance_overland_runge_kutta
    ! ============================================================================
    !> Subroutine build_DEM_kdtree
    !! Creates a kdtree for FW search.
    !<
    ! ============================================================================
    subroutine build_DEM_kdtree

        integer :: k

        if( .not. allocated( Fdata ) ) allocate( Fdata( 2, lnbPts ) )

        ! Create the kd-tree
        do k = 1, lnbPts
            Fdata( 1, k ) =  lcoordX( k )
            Fdata( 2, k ) =  lcoordY( k )
        enddo
        Ftree => kdtree2_create( Fdata, sort = .true., rearrange = .true. )

        return

    end subroutine build_DEM_kdtree
    ! ============================================================================
    !> Subroutine destroy_DEM_kdtree
    !! Destroys a kdtree for FW search.
    !<
    ! ============================================================================
    subroutine destroy_DEM_kdtree

        call kdtree2_destroy(Ftree)

        deallocate( Fdata )

        return

    end subroutine destroy_DEM_kdtree
    ! ============================================================================
    !> Subroutine deposit_flow_walkers_sediments
    !! Searches a kdtree to allocate deposition.
    !<
    ! ============================================================================
    subroutine deposit_flow_walkers_sediments( depo )

        integer :: ngid, k, ks

        real( tkind ) :: sumdist, sumwgth, dh
        real( tkind ), dimension( 2 ) :: pt
        real( tkind ), dimension( 4 ) :: dist, wgth
        real( tkind ), dimension( totgrn ) :: depo

        type(kdtree2_result), dimension( 4 ) :: FRslt

        pt( 1 ) = fw_xpos
        pt( 2 ) = fw_ypos

        ! Find face containing the flow walker
        call kdtree2_n_nearest(Ftree, pt, nn=4, results=FRslt)

        sumdist = 0.0_8
        ngid = -1
        ngb_loop: do k = 1, 4
            if( sqrt( FRslt( k )%dis ) <= 0.01_8 )then
                do ks = 1, totgrn
                    ngid = FRslt( k )%idx
                    deposit( ks, FRslt( k )%idx  ) = deposit( ks, FRslt( k )%idx ) + depo( ks )
                enddo
                exit ngb_loop
            endif
            dist( k ) = 1 / FRslt( k )%dis
            sumdist = sumdist + dist( k )
        enddo ngb_loop

        if( ngid >= 1 ) return

        ! Distribute deposits based on inverse weighted averaged distance
        sumwgth = 0.0_8
        do k = 1, 4
            wgth( k ) = dist( k ) / sumdist
            sumwgth = wgth( k ) + sumwgth
        enddo

        dh = 0.0_8
        do ks = 1, totgrn
            dh = depo( ks ) + dh
        enddo

        do k = 1, 4
            wgth( k ) = wgth( k ) / sumwgth
            do ks = 1, totgrn
                deposit( ks, FRslt( k )%idx ) = deposit( ks, FRslt( k )%idx ) + depo( ks ) * wgth( k )
            enddo
        enddo

        return

    end subroutine deposit_flow_walkers_sediments
    ! ============================================================================
    !> Subroutine search_DEM_kdtree_flowdischarge
    !! Searches a kdtree.
    !<
    ! ============================================================================
    subroutine search_DEM_kdtree_flowdischarge( pt, q )

        integer :: k, gid, ngid

        real( tkind ) :: q, sumdist, minq, maxq, sumwgth
        real( tkind ), dimension( 2 ) :: pt
        real( tkind ), dimension( 4 ) :: qnode, dist, wgth

        type(kdtree2_result), dimension( 4 ) :: FRslt

        call kdtree2_n_nearest(Ftree, pt, nn=4, results=FRslt)

        sumdist = 0.0_8
        minq = 1.e8
        maxq = -1.e8
        ngid = -1

        ngb_loop: do k = 1, 4
            gid = FRslt( k )%idx
            qnode( k ) = lfacr( gid ) * strat_dx**2.0_8 / ( secyear )
            dist( k ) = sqrt( FRslt( k )%dis )
            minq = min( minq, qnode( k ) )
            maxq = max( maxq, qnode( k ) )
            if( sqrt( FRslt( k )%dis ) <= 0.01_8 )then
                ngid = gid
                q = qnode( k )
                exit ngb_loop
            endif
            dist( k ) = 1 / FRslt( k )%dis
            sumdist = sumdist + dist( k )
        enddo ngb_loop

        if( ngid >= 1 ) return


        ! Compute elevation based on inverse weighted averaged distance
        q = 0.0_8
        sumwgth = 0.0_8
        do k = 1, 4
            wgth( k ) = (  dist( k ) / sumdist )
            q = wgth( k ) * qnode( k ) + q
            sumwgth = wgth( k ) + sumwgth
        enddo
        q = q / sumwgth

        if( minq > q .or. maxq < q )then
            if( abs( q - minq ) > abs( q - maxq ) )then
                q = maxq
            else
                q = minq
            endif
        endif

        return

    end subroutine search_DEM_kdtree_flowdischarge
    ! ============================================================================

end module fwpath
! ============================================================================
