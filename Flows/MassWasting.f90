! ============================================================================
! Name        : MassWasting.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file MassWasting.f90
!
! Description :  This module performs mass wasting computation.
!
!<
! ============================================================================
!> Module MassWasting 
!<
module mass

    use hillslp
    use fwpath
    use parallel
    use morpho
    use interpol
    use diffusion
    use file_data
    use ode_data
    use flow_data
    use time_data
    use mesh_data
    use forces_data
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

contains

    ! ============================================================================
    !> Subroutine build_debris_flow_depth
    !! This subroutine defines the depth of sediment transported by debris flow walkers.
    !<
    ! ============================================================================
    subroutine build_debris_flow_depth

        integer :: k, id, ngid, nlay, p, ks, layerid

        real( tkind ) :: SuA, limit_slope, unstable_soil_depth, remain_thick, prop
        real( tkind ) :: minz, dmz, deltah, total_thick, volume_sed, volume_water

        real( tkind ), dimension( totgrn ) :: sedh, lay_thick

        call update_creep_borders

        if( transport_mode == 0 .or. transport_mode == 2 )then
            do id = 1, nbPts
                k = lvertexID( id )
                nlay = v_nbLays( k )
                layerid = v_LaysID(  k, nlay )
                ! Get the updated thickness
                do ks = 1, totgrn
                    if( layID == layerid )then
                        top_sedh( k ,ks ) = stratLays( nlay, k, ks )
                        top_sedprev( k ,ks ) = stratLays( nlay, k, ks )
                    else
                        top_sedh( k ,ks ) = 0.0_8
                        top_sedprev( k ,ks ) = 0.0_8
                    endif
                enddo
            enddo

            ! Update borders
            do k = 1, lnbPts
                nlay = v_nbLays( k )
                layerid = v_LaysID(  k, nlay )
                if( k <= lstrat_X )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k + 2 * lstrat_X, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k + 2 * lstrat_X, 1:totgrn )
                elseif( k <= 2 * lstrat_X  )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k + lstrat_X, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k + lstrat_X, 1:totgrn )
                elseif( k >  lnbPts - 2 * lstrat_X )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k - lstrat_X, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k - lstrat_X, 1:totgrn )
                elseif( k >  lnbPts - lstrat_X )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k - 2 * lstrat_X, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k - 2 * lstrat_X, 1:totgrn )
                elseif( lcoordX( k ) == lstrat_xo )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k + 2, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k + 2, 1:totgrn )
                elseif( lcoordX( k ) == lstrat_xo + strat_dx )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k + 1, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k + 1, 1:totgrn )
                elseif( lcoordX( k ) == lstrat_xm )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k - 2, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k - 2, 1:totgrn )
                elseif( lcoordX( k ) == lstrat_xm - strat_dx )then
                    top_sedh( k ,1:totgrn ) = top_sedh( k - 1, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( k - 1, 1:totgrn )
                endif
                if( lcoordX( k ) < strat_xo .and. lcoordY( k ) < strat_yo )then
                    top_sedh( k ,1:totgrn ) = top_sedh( 1, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( 1, 1:totgrn )
                endif
                if( lcoordX( k ) > strat_xm .and. lcoordY( k ) < strat_yo )then
                    top_sedh( k ,1:totgrn ) = top_sedh( strat_X, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( strat_X, 1:totgrn )
                endif
                if( lcoordX( k ) < strat_xo .and. lcoordY( k ) > strat_ym )then
                    top_sedh( k ,1:totgrn ) = top_sedh( nbPts - strat_X + 1, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( nbPts - strat_X + 1, 1:totgrn )
                endif
                if( lcoordX( k ) > strat_xm .and. coordY( k ) > lstrat_ym )then
                    top_sedh( k ,1:totgrn ) = top_sedh( nbPts, 1:totgrn )
                    top_sedprev( k ,1:totgrn ) = top_sedprev( nbPts, 1:totgrn )
                endif
            enddo
        endif

        ! Triggering areas detection using Horton et al. (2008)
        do k = 1, nbPts

            ! Surface of the upslope contributing area
            SuA = facc( k ) * ( strat_dx / 1000.0_8 )**2.0_8

            ! Compute limiting slope gradient
            if( SuA >= 2.5_8 )then
                limit_slope = 0.26_8
            else
                limit_slope = 0.31_8 * SuA**(-0.15_8)
            endif

            unstable_soil_depth = 0.0_8

            ! In case a triggering area is detected
            if( slp( k ) > limit_slope .and. hcurv( k ) < 0.0_8 )then

                ! Compute the depth of material requires to reach the slope limit
                unstable_soil_depth = strat_dx * ( slp( k ) - limit_slope )

                ! Get the drained water volume in the area
                volume_water = facc( k ) * strat_dx**2.0_8

                ! Compute the maximum volume of sediment available for transport
                if( coordZ( k ) >= gsea%actual_sea )then
                    volume_sed = debris_factor_aerial * debris_max_conc_aerial * &
                        volume_water / ( 1.0_8 - debris_max_conc_aerial )
                else
                    volume_sed = debris_factor_marine * debris_max_conc_marine * &
                        volume_water / ( 1.0_8 - debris_max_conc_marine )
                endif

                ! Limit sediment availability
                unstable_soil_depth = min( unstable_soil_depth, volume_sed / strat_dx**2.0_8 )

                ! Perform a second limitation based on neighbors elevation
                id = lvertexID( k )
                minz = lcoordZ( id )
                do ks = 1, 8
                    ngid =  lngbID( id, ks )
                    if( ngid > 0 )then
                        dmz = lcoordZ( ngid )
                        ! South
                        if( lcoordY( ngid ) < strat_yo .and. boundcond( 1 ) == 2 )then
                            dmz = 1.e6_8
                        endif
                        ! North
                        if( lcoordY( ngid ) > strat_ym .and. boundcond( 2 ) == 2 )then
                            dmz = 1.e6_8
                        endif
                        ! West
                        if( lcoordX( ngid ) < strat_xo .and. boundcond( 3 ) == 2 )then
                            dmz = 1.e6_8
                        endif
                        ! East
                        if( lcoordX( ngid ) > strat_xm .and. boundcond( 4 ) == 2 )then
                            dmz = 1.e6_8
                        endif
                        minz = min( dmz, minz )
                    endif
                enddo
                deltah = max( lcoordZ( id ) - minz, 0.0_8 )

                unstable_soil_depth = min( unstable_soil_depth, deltah )
                if( unstable_soil_depth == 0.0_8 ) goto 19

            elseif( ( coordZ( k ) >= gsea%actual_sea .and. slp( k ) > critical_slope_aerial ) .or. &
                ( coordZ( k ) < gsea%actual_sea .and. slp( k ) > critical_slope_marine ) )then

                ! Compute the depth of material requires to reach the critical slope
                if( coordZ( k ) >= gsea%actual_sea )then
                    unstable_soil_depth = strat_dx * ( slp( k ) - critical_slope_aerial )
                else
                    unstable_soil_depth = strat_dx * ( slp( k ) - critical_slope_marine )
                endif

                ! Perform a second limitation based on neighbors elevation
                id = lvertexID( k )
                minz = lcoordZ( id )
                do ks = 1, 8
                    ngid =  lngbID( id, ks )
                    if( ngid > 0 )then
                        minz = min( lcoordZ( ngid ), minz )
                    endif
                enddo
                deltah = max( lcoordZ( id ) - minz, 0.0_8 )

                unstable_soil_depth = min( unstable_soil_depth, deltah )
                if( unstable_soil_depth == 0.0_8 ) goto 19

            endif

            if( unstable_soil_depth > 0.0_8 )then

                ! Erode stratigraphic layer
                nlay = v_nbLays( id )
                layerid = v_LaysID( id, nlay )
                if( layerid == layID )then
                    total_thick = 0.0_8
                    do ks = 1, totgrn
                        total_thick = total_thick + stratLays( nlay, id, ks )
                    enddo
                    if( total_thick >= unstable_soil_depth ) goto 19
                endif

                lay_thick = 0.0_8
                total_thick = 0.0_8
                remain_thick = unstable_soil_depth
                erosion_loop: do p = v_nbLays( id ), 1, -1

                    if( remain_thick == 0.0_8 ) exit erosion_loop

                    ! Get the amount to erode for the considered layer
                    sedh( 1:totgrn ) = 0.0_8
                    do ks = 1, totgrn
                        sedh( ks ) = stratLays( p, id, ks )
                        total_thick = total_thick + stratLays( p, id, ks )
                    enddo

                    ! In case the layer contains enough materials
                    if( remain_thick < total_thick )then

                        prop = remain_thick / total_thick
                        do ks = 1, totgrn
                            sedh( ks ) = stratLays( p, id, ks ) * prop
                            stratLays( p, id, ks ) = stratLays( p, id, ks ) - sedh( ks )
                            lay_thick( ks )  = lay_thick( ks ) + sedh( ks )
                        enddo
                        exit erosion_loop

                    ! Otherwise grab all the layer
                    else

                        do ks = 1, totgrn
                            lay_thick( ks ) = lay_thick( ks ) + sedh( ks )
                            stratLays( p, id, ks ) = 0.0_8
                            if( gporo%compaction ) porosityLays( p, id, ks ) = 0.0_8
                        enddo
                        nlay = nlay - 1
                        remain_thick = remain_thick - total_thick
                        if( remain_thick < tor ) remain_thick = 0.0_8

                    endif


                enddo erosion_loop

                ! Update the top layer
                if( nlay < v_nbLays( id ) .or. layerid < layID ) nlay = nlay + 1
                v_nbLays( id ) = nlay
                v_LaysID( id, nlay ) = layID
                do ks = 1, totgrn
                    stratLays( nlay, id, ks ) = lay_thick( ks )
                    if( gporo%compaction ) porosityLays( nlay, id, ks ) = porosity( ks, 1 )
                enddo
            endif

 19 continue

        enddo

        call update_creep_borders

        ! Diffuse freshly deposited sediments over the topography
        call soil_diffusion_transport

        call update_creep_borders

        return

    end subroutine build_debris_flow_depth
    ! ============================================================================
    !> Subroutine define_soil_creep_thickness
    !! This subroutine compute soil creep process based on diffusion rule.
    !<
    ! ============================================================================
    subroutine define_soil_creep_thickness

        integer :: k, id, nlay, p, ks, layerid

        real( tkind ) :: remain_thick, prop, dh
        real( tkind ) :: total_thick, minz, Cdiff
        real( tkind ) :: eromax, depth_creep

        real( tkind ), dimension( lnbPts ) :: oldz 
        real( tkind ), dimension( totgrn ) :: sedh, lay_thick

        oldz = lcoordZ

        do k = 1, nbPts

            depth_creep = 0.0_8

            if( hcurv( k ) + vcurv( k ) > 0.0_8 )then
                Cdiff = 0.0_8
                do ks = 1, totgrn
                    if( oldz( k ) > gsea%actual_sea )then
                        Cdiff = Cdiff + sediment( ks )%creep_marine
                    else
                        Cdiff = Cdiff + sediment( ks )%creep_aerial
                    endif
                enddo
                Cdiff = Cdiff / totgrn

                depth_creep = Cdiff * ( hcurv( k ) + vcurv( k ) ) * rain_int

                ! Get maximum erosion thickness to prevent formation of hole
                id = lvertexID( k )
                minz = oldz( id )

                ! Find neighborhood minimum elevation
                do p = 1, 8
                    if(  lngbID( id, p ) > 0 ) minz = min( minz, oldz( lngbID( id, p ) ) )
                enddo

                ! Maximum erosion to prevent formation of holes
                eromax = 0.8_8 * ( oldz( id ) - minz )
                if( eromax < 0.0_8 ) eromax = 0.0_8

                depth_creep = min( depth_creep, eromax )

                if( depth_creep > 1.e-6_8 )then
                    nlay = v_nbLays( id )
                    depth_creep = depth_creep / totgrn
                    do ks = 1, totgrn
                        stratLays( nlay, id, ks ) = stratLays( nlay, id, ks ) - depth_creep
                        if( stratLays( nlay, id, ks ) < 0.0_8 )  stratLays( nlay, id, ks ) = 0.0_8
                    enddo
                    lcoordZ( id ) = lcoordZ( id ) - depth_creep * totgrn
                    coordZ( k ) = lcoordZ( id )
                    soil_thick( id ) = soil_thick( id ) - depth_creep * totgrn
                    if( soil_thick( id ) < 0.0_8 ) soil_thick( id ) = 0.0_8
                endif

            endif

        enddo


        return

    end subroutine define_soil_creep_thickness
    ! ============================================================================

end module mass
! ============================================================================
