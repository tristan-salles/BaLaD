! ============================================================================
! Name        : LandscapeEvol.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file LandscapeEvol.f90
!
! Description :  This module performs landscape evolution based on diffusion.
!
!<
! ============================================================================
!> Module lem
!<
module lem

    use stream
    use fwpath
    use parallel
    use morpho
    use interpol
    use regolith
    use tin_data
    use diffusion
    use file_data
    use ode_data
    use flow_algo
    use flow_data
    use time_data
    use mesh_data
    use forces_data
    use stream_erodep
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    logical :: dcht_mode, trspt_mode

contains

    ! ============================================================================
    !> Subroutine advance_landscape_evolution
    !! This subroutine computes landscape evolution through time
    !<
    ! ============================================================================
    subroutine advance_landscape_evolution

        integer :: s, k, ks, p, id, nlay, sd

        integer, dimension( nbPts ) :: orderlist

        real( tkind ) :: influx, detachment, transport, discharge, maxh, eromax, mardiff
        real( tkind ) :: frac, dtb, totero, dh, dh2, stillsed, deposit_hlimit, rval

        real( tkind ), dimension( totgrn ) :: prop, prop2, poro, sedh, remain, erode, hs, hs2, erodeflux, depositflux
        real( tkind ), dimension( nbPts ) :: temparray

        if( .not. allocated( newz ) ) allocate( newz( lnbPts ) )
        newz = lcoordZ

        if( .not. allocated( sedin ) ) allocate( sedin( nbPts, totgrn ) )
        sedin = 0.0_8

        ! First we sort the flow accumulation array by increasing order
        temparray = filldem
        call quick_sort( temparray, orderlist )

        call getFluxWeight

        ! Start from the highest elevation to the lowest one.
        do s = nbPts, 1, -1

            ! Take the sorted point ID
            k = orderlist( s )

            totero = 0.0_8
            influx = 0.0_8
            erodeflux = 0.0_8
            depositflux = 0.0_8
            id = lvertexID( k )
            prop = 0.0_8

            ! Get soil cover composition
            dtb = soil_thick( id )
            call get_regolith_composition(  lvertexID( k ), dtb, dh, hs, poro )
            if( dh == 0.0_8 )then
                prop( 1 ) = 1.0
            else
                do ks = 1, totgrn
                    prop( ks ) = hs( ks ) / dh
                enddo
            endif
            dh = 0.9_8 * dh
            hs = 0.9_8 * hs

            ! Get maximum erosion thickness to prevent formation of hole
            call get_erosion_thickness( k, eromax )
            eromax = eromax * 0.8_8
            call get_regolith_composition(  lvertexID( k ), eromax, dh2, hs2, poro )
            if( dh2 == 0.0_8 )then
                prop2( 1 ) = 1.0
            else
                do ks = 1, totgrn
                    prop2( ks ) = hs2( ks ) / dh2
                enddo
            endif

            if( newz( id ) <= gsea%actual_sea ) eromax = 0.0_8

            if( eromax == 0.0_8 )then
                prop2 = 0.0_8
                dh2 = 0.0_8
                hs2 = 0.0_8
                prop = 0.0_8
                dh = 0.0_8
                hs = 0.0_8
            elseif( dh2 > eromax )then
                do ks = 1, totgrn
                    hs2( ks ) = hs2( ks ) * eromax / dh2
                enddo
                dh2 = eromax
            endif
            dh2 = 0.8_8 * dh2
            hs2 = 0.8_8 * hs2
            
            ! In case the depth to bedrock is greater than the maximum erosion change it to match
            if( dtb > eromax .or. dtb == 0.0_8 ) dtb = eromax

            ! Fluvial discharge
            if( newz( id ) >= gsea%actual_sea .and. vflow( k ) > 0.0_8 )then
                discharge = qflow( k ) * secyear
            else
                discharge = 0.0_8
            endif

            if( discharge > 0.0_8 )then
                ! Find if we are in transport or detachment limited mode
                dcht_mode = .false.
                trspt_mode = .false.
                tpt_mode: do ks = 1, totgrn

                    ! Transport capacity
                    transport = prop( ks ) * sediment( ks )%transport * discharge**( mtransport ) * ( slp( k ) )**( ntransport ) &
                        * rain_int / ( wflow( k ) * strat_dx )

                    ! Detachment law
                    detachment = prop2( ks ) * sediment( ks )%ero_coeff * ( discharge )**( mdetach) * &
                        ( slp( k ) )**( ndetach ) * rain_int

                    ! Are we in detachment or transport mode
                    if( soil_thick( id ) == 0.0_8 )then
                        dcht_mode = .true.
                    elseif( transport - sedin( k, ks ) > detachment )then
                        dcht_mode = .true.
                    else
                        trspt_mode = .true.
                        exit tpt_mode
                    endif
                enddo tpt_mode

                ! Based on type of fluvial transport get the erosion fluxes
                do ks = 1, totgrn

                    influx = sedin( k, ks ) + influx

                    ! For detachment limited mode
                    if( .not. trspt_mode )then

                        ! Detachment law
                        detachment = prop2( ks ) * sediment( ks )%ero_coeff * ( discharge )**( mdetach) * &
                            ( slp( k ) )**( ndetach ) * rain_int

                        erodeflux( ks ) = detachment
                        if( hs2( ks ) <  erodeflux( ks ) ) erodeflux( ks ) = hs2( ks )

                    ! For transport limited mode
                    else

                        ! Transport capacity
                        transport = prop( ks ) * sediment( ks )%transport * discharge**( mtransport ) * &
                            ( slp( k ) )**( ntransport )  * rain_int / ( wflow( k ) * strat_dx ) !strat_dx**2

                        if( transport > sedin( k, ks ) )then
                            erodeflux( ks ) = transport - sedin( k, ks )
                        endif

                        if( erodeflux( ks ) < 0.0_8 ) erodeflux( ks ) = 0.0_8
                        if( hs2( ks ) <  erodeflux( ks ) ) erodeflux( ks ) = hs2( ks )

                    endif

                    ! Update erosion value
                    totero = totero + erodeflux( ks )

                enddo
            endif
            
            ! Perform erosion
            if( totero > 1.e-6_8 )then

                hs = erodeflux

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
                                print*,'Something went wrong when eroding stratigraphy in LEM.',iam
                                print*,ks,id,p,remain( ks ),v_nbLays( id ),stratLays( p, id, 1:totgrn ),erodeflux( ks )
                                stop
                            endif
                            dh = dh + stratLays( p, id, ks ) - remain( ks )
                            stratLays( p, id, ks ) = remain( ks )
                            erodeflux( ks ) = erodeflux( ks ) - hs( ks )
                        enddo
                        soil_thick( id ) = soil_thick( id ) - dh
                        if( soil_thick( id ) < 0.0_8 ) soil_thick( id ) = 0.0_8
                        lcoordZ( id ) = lcoordZ( id ) - dh
                        coordZ( k ) = lcoordZ( id )
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
                        soil_thick( id ) = soil_thick( id ) - dh
                        if( soil_thick( id ) < 0.0_8 ) soil_thick( id ) = 0.0_8
                        lcoordZ( id ) = lcoordZ( id ) - dh
                        coordZ( k ) = lcoordZ( id )
                        nlay = nlay - 1
                    endif
                enddo erosion_loop

                v_nbLays( id ) = nlay
            else
                erodeflux = 0.0_8
            endif

            ! Perform deposition and update sediment flux
            depositflux = 0.0_8
            if( influx > 0.0_8 )then

                ! Get maximum deposition for stability based
                ! on neighborhood elevations
                deposit_hlimit = 1.0e6_8
                if( newz( id ) >= gsea%actual_sea )then
                    maxh = newz( id ) + 1.e6_8
                    do p = 1, 8
                        if( lcoordZ( lngbID( id, p ) ) > newz( id ) ) &
                            maxh = min( maxh, lcoordZ( lngbID( id, p ) ) )
                    enddo
                    deposit_hlimit = ( maxh - newz( id ) ) * 0.8_8
                    if( deposit_hlimit < 0.0_8 ) deposit_hlimit = 0.0_8
                else
                    deposit_hlimit = ( gsea%actual_sea - newz( id ) ) * 0.8_8
!                    deposit_hlimit = - 0.1e6_8 * rain_int / strat_dx**2 * &
!                        ( newz( id + 1 ) - 2.0_8 * newz( id ) + newz( id - 1 ) + &
!                        newz( id + lstrat_X ) - 2.0_8 * newz( id ) + newz( id - lstrat_X ) )
!                    maxh = 1.0e6_8
!                    do p = 1, 8
!                        if( newz( lngbID( id, p ) ) > gsea%actual_sea ) &
!                            maxh = min( maxh, newz( lngbID( id, p ) ) )
!                    enddo
!                    if( maxh < 1.e6_8 )then
!                        if( iam == 0 )then
!                            call random_seed( size = sd )
!                            call random_number( harvest = rval )
!                        endif
!                        call mpi_bcast( rval, 1, dbl_type, 0, badlands_comm_world, ierr )
!                        maxh = ( maxh - gsea%actual_sea ) * rval
!                        maxh = maxh + gsea%actual_sea
!                        deposit_hlimit = max( deposit_hlimit, ( maxh - newz( id ) ) )
!                    endif
                    !if( hcurv( k ) + vcurv( k ) < 0.0_8 )then
!                        rval = ( hcurv( k ) + vcurv( k ) ) * rain_int * 0.01e6
!                        if( rval < deposit_hlimit ) deposit_hlimit = rval
!                    else
!                        deposit_hlimit = 0.0_8
!                    endif
                    if( deposit_hlimit < 0.0_8 ) deposit_hlimit = 0.0_8
                endif

                ! Deposit sediment in respect to their diameter sizes
                ! from coarser to finer
                deposit_loop: do ks = 1, totgrn

                    if( erodeflux( ks ) == 0.0_8 )then

                        if( sedin( k, ks ) > deposit_hlimit )then
                            depositflux( ks ) = deposit_hlimit
                            deposit_hlimit = 0.0_8
                            exit deposit_loop
                        else
                            depositflux( ks ) = sedin( k, ks )
                            deposit_hlimit = deposit_hlimit - sedin( k, ks )
                        endif
                        if( deposit_hlimit < 1.e-6_8 ) exit deposit_loop

                    endif

                enddo deposit_loop

                ! Update stratigraphic layer and topography
                ! with incoming sediment fluxes
                dh = 0.0_8
                if( layID == v_LaysID( id, v_nbLays( id ) ) )then
                    do ks = 1, totgrn
                        ! Prepare diffusion arrays
                        top_sedh( id, ks ) = stratLays( v_nbLays( id ), id, ks )
                        top_sedprev( id, ks ) = stratLays( v_nbLays( id ), id, ks )
                        ! Compute sedimentary changes
                        if( gporo%compaction )then
                            depositflux( ks ) = depositflux( ks ) / ( 1.0_8 - porosity( ks, 1 ) )
                            porosityLays( v_nbLays( id ), id, ks ) = porosity( ks, 1 )
                        endif
                        stratLays( v_nbLays( id ), id, ks ) = stratLays( v_nbLays( id ), id, ks ) + depositflux( ks )
                        dh = dh + depositflux( ks )
                    enddo
                    lcoordZ( id ) = lcoordZ( id ) + dh
                    soil_thick( id ) = soil_thick( id ) + dh
                    if( soil_thick( id ) > dtb_marine .and. lcoordZ( id ) < gsea%actual_sea ) &
                        soil_thick( id ) = dtb_marine
                    if( soil_thick( id ) > dtb_aerial .and. lcoordZ( id ) >= gsea%actual_sea ) &
                        soil_thick( id ) = dtb_aerial
                    coordZ( k ) = lcoordZ( id )
                else
                    do ks = 1, totgrn
                        ! Prepare diffusion arrays
                        top_sedh( id, ks ) = 0.0_8
                        top_sedprev( id, ks ) = 0.0_8
                        ! Compute sedimentary changes
                        dh = dh + depositflux( ks )
                    enddo
                    if( dh > 0.0_8 )then
                        v_nbLays( id ) = v_nbLays( id ) + 1
                        v_LaysID( id, v_nbLays( id ) ) = layID
                        do ks = 1, totgrn
                            if( gporo%compaction )then
                                depositflux( ks ) = depositflux( ks ) / ( 1.0_8 - porosity( ks, 1 ) )
                                porosityLays( v_nbLays( id ), id, ks ) = porosity( ks, 1 )
                            endif
                            stratLays( v_nbLays( id ), id, ks ) = depositflux( ks )
                        enddo
                        lcoordZ( id ) = lcoordZ( id ) + dh
                        soil_thick( id ) = soil_thick( id ) + dh
                        if( soil_thick( id ) > dtb_marine .and. lcoordZ( id ) < gsea%actual_sea ) &
                            soil_thick( id ) = dtb_marine
                        if( soil_thick( id ) > dtb_aerial .and. lcoordZ( id ) >= gsea%actual_sea ) &
                            soil_thick( id ) = dtb_aerial
                        coordZ( k ) = lcoordZ( id )
                    endif
                endif

            endif

            ! Find outgoing sediment fluxes
            do p = 1, 8
                if( fluxweight( id, p ) > 0.0_8 .and. ngbID( k, p ) > 0 )then
                    do ks = 1, totgrn
                        ! Erosion
                        if( erodeflux( ks ) > 0.0_8 )then
                            sedin( ngbID( k, p ), ks ) = sedin( ngbID( k, p ), ks ) + &
                                fluxweight( id, p ) * ( sedin( k, ks ) + erodeflux( ks ) )
                        ! Deposition
                        elseif( depositflux( ks ) > 0.0_8 )then
                            sedin( ngbID( k, p ), ks ) = sedin( ngbID( k, p ), ks ) + &
                                fluxweight( id, p ) * ( sedin( k, ks ) - depositflux( ks ) )
                        ! Transport
                        else
                            sedin( ngbID( k, p ), ks ) = sedin( ngbID( k, p ), ks ) + &
                                fluxweight( id, p ) * sedin( k, ks )
                        endif
                    enddo
!                elseif( fluxweight( id, p ) > 0.0_8 .and. ngbID( k, p ) <= 0 )then
!                   outsed = outsed + fluxweight( id, p ) * sedin( k, 1 )
                endif
           enddo

        enddo

        ! Diffuse freshly deposited sediments over the topography
!        call soil_diffusion_transport

        return

    end subroutine advance_landscape_evolution
    ! ============================================================================
    !> Subroutine quadratic_polynomial_slope
    !! This subroutine computes slope value using quadratic polynomial interpolation.
    !<
    ! ============================================================================
    subroutine quadratic_polynomial_slope( id, slp_val )

        integer :: id

        real( tkind ) :: z( 9 ), p, q, slp_val

        ! Elevation
        z( 1 ) = lcoordZ( lngbID( id, 8 ) )
        z( 2 ) = lcoordZ( lngbID( id, 1 ) )
        z( 3 ) = lcoordZ( lngbID( id, 2 ) )
        z( 4 ) = lcoordZ( lngbID( id, 7 ) )
        z( 5 ) = lcoordZ( id )
        z( 6 ) = lcoordZ( lngbID( id, 3 ) )
        z( 7 ) = lcoordZ( lngbID( id, 6 ) )
        z( 8 ) = lcoordZ( lngbID( id, 5 ) )
        z( 9 ) = lcoordZ( lngbID( id, 4 ) )

        ! Quadratic polynomial coefficient
        p = ( z(3) + z(6) + z(9) - z(1) - z(4) - z(7) ) / ( 6.0_8 * strat_dx )
        q = ( z(1) + z(2) + z(3) - z(7) - z(8) - z(9) ) / ( 6.0_8 * strat_dx )

        ! Get the gradient
        slp_val = atan( sqrt( p*p + q*q ) )

    end subroutine quadratic_polynomial_slope
    ! ============================================================================

end module lem
! ============================================================================

