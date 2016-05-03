! ============================================================================
! Name        : HillSlope.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file HillSlope.f90
!
! Description :  This module performs upland hillslope evolution.
!
!<
! ============================================================================
!> Module HillSlope 
!<
module hillslp

    use fwpath
    use parallel
    use morpho
    use interpol
    use regolith
    use file_data
    use sedtrans
    use ode_data
    use flow_data
    use time_data
    use mesh_data
    use forces_data

    implicit none

    real( tkind ), dimension( : ), allocatable :: inflx
    real( tkind ), dimension( : ), allocatable :: nsoil

contains

    ! ============================================================================
    !> Subroutine compute_hillslope
    !! This subroutine computes hillsope processes based on soil creep formulation.
    !! Erosion and deposition in the model is computed using the Forward?Time?Centered?Space
    !! (FTCS) method based on depth? and slope?dependent model.
    !<
    ! ============================================================================
    subroutine compute_hillslope( creep_mtd )

        integer :: s, k, id, p, ks, creep_mtd, nlay

        integer, dimension( nbPts ) :: orderlist

        real( tkind ) :: creep, qyhalf1, qyhalf2, qxhalf1, qxhalf2
        real( tkind ) :: total_curvature, potential_change

        real( tkind ) :: depositflx, erodeflx, outflx
        real( tkind ), dimension( nbPts ) :: temparray

        if( tnow < elapsed_time )then
            elapsed_time = tnow
            return
        endif

        call update_creep_borders

        ! Define temporary arrays
        if( .not. allocated( inflx ) ) allocate( inflx( nbPts ) )
        if( .not. allocated( nsoil ) ) allocate( nsoil( lnbPts ) )
        if( .not. allocated( newz ) ) allocate( newz( lnbPts ) )
        nsoil = soil_thick
        newz = lcoordZ

        do ks = totgrn, 1, -1

            inflx = 0.0_8

            ! Adjust for close boundary
            do k = 1, lnbPts
                ! South
                if( lcoordY( k ) < strat_yo .and. boundcond( 1 ) == 2 )then
                    newz( k ) = 1.e6_8
                endif
                ! North
                if( lcoordY( k ) > strat_ym .and. boundcond( 2 ) == 2 )then
                    newz( k ) = 1.e6_8
                endif
                ! West
                if( lcoordX( k ) < strat_xo .and. boundcond( 3 ) == 2 )then
                    newz( k ) = 1.e6_8
                endif
                ! East
                if( lcoordX( k ) > strat_xm .and. boundcond( 4 ) == 2 )then
                    newz( k ) = 1.e6_8
                endif
            enddo

            ! Sort the grid elevation array by increasing order
            temparray = coordZ
            call quick_sort( temparray, orderlist )

            ! Routing of colluvial material is calculated via MFD method
            ! Get MFD fluxes weight
            call getFluxWeightTopo

            ! Start from the highest elevation to the lowest one.
            do s = nbPts, 1, -1

                ! Take the sorted point ID
                k = orderlist( s )
                id = lvertexID( k )

                ! Define the creep diffusion coefficient
                erodeflx = 0.0_8
                depositflx = 0.0_8

                creep = sediment( ks )%creep_marine
                if( newz( id ) > gsea%actual_sea ) creep = sediment( ks )%creep_aerial

                ! Compute x and y- components for the slope- dependent model.
                if( creep_mtd >= 3 )then
                    qxhalf1 = ( newz( id + 1 ) - newz( id ) )
                    qxhalf2 = ( newz( id - 1 ) - newz( id ) )
                    qyhalf1 = ( newz( id + lstrat_X ) - newz( id ) )
                    qyhalf2 = ( newz( id - lstrat_X ) - newz( id ) )

                    if( creep_mtd == 4 )then
                        ! Use the curvature as a diffusion factor
                        ! way to take into account the biomass
                        total_curvature = hcurv( k ) + vcurv( k )
                        if( total_curvature > -0.04_8 .and. total_curvature <= 0.0_8 )then
                            qxhalf1 = 1.5_8 * qxhalf1
                            qxhalf2 = 1.5_8 * qxhalf2
                            qyhalf1 = 1.5_8 * qyhalf1
                            qyhalf2 = 1.5_8 * qyhalf2
                        elseif( total_curvature < -0.04_8 )then
                            qxhalf1 = 2.0_8 * qxhalf1
                            qxhalf2 = 2.0_8 * qxhalf2
                            qyhalf1 = 2.0_8 * qyhalf1
                            qyhalf2 = 2.0_8 * qyhalf2
                        endif
                    endif
                else
                    ! Compute x and y- components for the depth- and slope- dependent model.
                    qxhalf1 = 0.5_8 * ( prop( id + 1, ks ) * nsoil( id + 1 ) + prop( id, ks ) * nsoil( id ) ) * &
                        ( newz( id + 1 ) - newz( id ) )
                    qxhalf2 = 0.5_8 * ( prop( id - 1, ks ) * nsoil( id - 1 ) + prop( id, ks ) * nsoil( id ) ) * &
                        ( newz( id - 1 ) - newz( id ) )
                    qyhalf1 = 0.5_8 * ( prop( id + lstrat_X, ks ) * nsoil( id + lstrat_X ) + prop( id, ks ) * nsoil( id ) ) * &
                        ( newz( id + lstrat_X ) - newz( id ) )
                    qyhalf2 = 0.5_8 * ( prop( id - lstrat_X, ks ) * nsoil( id - lstrat_X ) + prop( id, ks ) * nsoil( id ) ) * &
                        ( newz( id - lstrat_X ) - newz( id ) )

                    if( creep_mtd == 2 )then
                        ! Use the curvature as a diffusion factor
                        ! way to take into account the biomass
                        total_curvature = hcurv( k ) + vcurv( k )
                        if( total_curvature > -0.04_8 .and. total_curvature <= 0.0_8 )then
                            qxhalf1 = 1.5_8 * qxhalf1
                            qxhalf2 = 1.5_8 * qxhalf2
                            qyhalf1 = 1.5_8 * qyhalf1
                            qyhalf2 = 1.5_8 * qyhalf2
                        elseif( total_curvature < -0.04_8 )then
                            qxhalf1 = 2.0_8 * qxhalf1
                            qxhalf2 = 2.0_8 * qxhalf2
                            qyhalf1 = 2.0_8 * qyhalf1
                            qyhalf2 = 2.0_8 * qyhalf2
                        endif
                    endif

                endif

                ! Potential erosion/deposition
                potential_change = -( qxhalf1 + qxhalf2 + qyhalf1 + qyhalf2 ) * &
                    ( tnow - elapsed_time ) * creep / strat_dx**2.0_8
                if( potential_change > 0.0_8 )then
                    erodeflx = potential_change
                    depositflx = 0.0_8
                elseif( potential_change < 0.0_8 )then
                    erodeflx = 0.0_8
                    depositflx = -potential_change
                endif

                call manage_sedfluxes( k, ks, erodeflx, depositflx, inflx( k ), .true., nsoil( id ), outflx )

                ! Update outgoing sediment fluxes
                ! Loop over the neighborhood
                do p = 1, 8
                    ! Flux weights computed based on MFD approach
                    if( fluxweightdem( id, p ) > 0.0_8 .and. ngbID( k, p ) > 0 )then
                        inflx( ngbID( k, p ) ) = inflx( ngbID( k, p ) ) + fluxweightdem( id, p ) * outflx

                    ! In case sediment flux exits the simulation
                    elseif( ngbID( k, p ) <= 0 )then

                        ! South border
                        if( coordY( k ) == strat_yo )then
                            ! Reflective
                            if( boundcond( 1 ) == 0 )then
                                if( p == 4 )then
                                    inflx( k + 1 ) = inflx( k + 1 ) + fluxweightdem( id, p ) * outflx
                                elseif( p == 5 .and. fluxweightdem( id, p ) > 0.0_8 )then
                                    print*,'There is an issue with the reflective south boundary'
                                elseif( p == 6 )then
                                    inflx( k - 1 ) = inflx( k - 1 ) + fluxweightdem( id, p ) * outflx
                                endif
                            endif
                        ! North border
                        elseif( coordY( k ) == strat_ym )then
                            ! Reflective
                            if( boundcond( 2 ) == 0 )then
                                if( p == 2 )then
                                    inflx( k + 1 ) = inflx( k + 1 ) + fluxweightdem( id, p ) * outflx
                                elseif( p == 1 .and. fluxweightdem( id, p ) > 0.0_8  )then
                                    print*,'There is an issue with the reflective north boundary'
                                elseif( p == 8 )then
                                    inflx( k - 1 ) = inflx( k - 1 ) + fluxweightdem( id, p ) * outflx
                                endif
                            endif
                        ! West border
                        elseif( coordX( k ) == strat_xo )then
                            ! Reflective
                            if( boundcond( 3 ) == 0 )then
                                if( p == 6 )then
                                    inflx( k - strat_X ) = inflx( k - strat_X ) + fluxweightdem( id, p ) * outflx
                                elseif( p == 7 .and. fluxweightdem( id, p ) > 0.0_8  )then
                                    print*,'There is an issue with the reflective west boundary'
                                elseif( p == 8 )then
                                    inflx( k + strat_X ) = inflx( k + strat_X ) + fluxweightdem( id, p ) * outflx
                                endif
                            endif
                        ! East border
                        elseif( coordX( k ) == strat_xm )then
                            ! Reflective
                            if( boundcond( 4 ) == 0 )then
                                if( p == 2 )then
                                    inflx( k + strat_X ) = inflx( k + strat_X ) + fluxweightdem( id, p ) * outflx
                                elseif( p == 3 .and. fluxweightdem( id, p ) > 0.0_8  )then
                                    print*,'There is an issue with the reflective west boundary'
                                elseif( p == 4 )then
                                    inflx( k - strat_X ) = inflx( k - strat_X ) + fluxweightdem( id, p ) * outflx
                                endif
                            endif
                        endif
                    endif
                enddo
            enddo

            ! Update elevation
            if( ks < totgrn ) call update_creep_borders_elevation

        enddo

        elapsed_time = tnow

        call update_creep_borders

        return

    end subroutine compute_hillslope
    ! ============================================================================
    !> Subroutine update_creep_borders
    !!
    !<
    ! ============================================================================
    subroutine update_creep_borders

        integer :: k, p, ks, id

        ! Update borders large DEM
        do k = 1, lnbPts
            p = v_nbLays( k )
            ! South
            if( k <= lstrat_X )then
                lcoordZ( k ) = lcoordZ( k + 2 * lstrat_X )
                soil_thick( k ) = soil_thick( k + 2 * lstrat_X )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k + 2 * lstrat_X, ks )
                enddo
            elseif( k <= 2 * lstrat_X  )then
                lcoordZ( k ) = lcoordZ( k + lstrat_X )
                soil_thick( k ) = soil_thick( k + lstrat_X )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k + lstrat_X, ks )
                enddo
            ! North
            elseif( k >  lnbPts - 2 * lstrat_X )then
                lcoordZ( k ) = lcoordZ( k - lstrat_X )
                soil_thick( k ) = soil_thick( k - lstrat_X )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k - lstrat_X, ks )
                enddo
            elseif( k >  lnbPts - lstrat_X )then
                lcoordZ( k ) = lcoordZ( k - 2 * lstrat_X )
                soil_thick( k ) = soil_thick( k - 2 * lstrat_X )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k - 2 * lstrat_X, ks )
                enddo
            ! West
            elseif( lcoordX( k ) == lstrat_xo )then
                lcoordZ( k ) = lcoordZ( k + 2 )
                soil_thick( k ) = soil_thick( k + 2 )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k + 2, ks )
                enddo
            elseif( lcoordX( k ) == lstrat_xo + strat_dx )then
                lcoordZ( k ) = lcoordZ( k + 1 )
                soil_thick( k ) = soil_thick( k + 1 )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k + 1, ks )
                enddo
            ! East
            elseif( lcoordX( k ) == lstrat_xm )then
                lcoordZ( k ) = lcoordZ( k - 2 )
                soil_thick( k ) = soil_thick( k - 2 )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k - 2, ks )
                enddo
            elseif( lcoordX( k ) == lstrat_xm - strat_dx )then
                lcoordZ( k ) = lcoordZ( k - 1 )
                soil_thick( k ) = soil_thick( k - 1 )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, k - 1, ks )
                enddo
            endif
            ! SW corner
            if( lcoordX( k ) < strat_xo .and. lcoordY( k ) < strat_yo )then
                lcoordZ( k ) = coordZ( 1 )
                id = lvertexID( 1 )
                soil_thick( k ) = soil_thick( id )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, id, ks )
                enddo
            endif
            ! SE corner
            if( lcoordX( k ) > strat_xm .and. lcoordY( k ) < strat_yo )then
                lcoordZ( k ) = coordZ( strat_X )
                id = lvertexID( strat_X )
                soil_thick( k ) = soil_thick( id )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, id, ks )
                enddo
            endif
            ! NW corner
            if( lcoordX( k ) < strat_xo .and. lcoordY( k ) > strat_ym )then
                lcoordZ( k ) = coordZ( nbPts - strat_X + 1 )
                id = lvertexID( nbPts - strat_X + 1 )
                soil_thick( k ) = soil_thick( id )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, id, ks )
                enddo
            endif
            ! NE corner
            if( lcoordX( k ) > strat_xm .and. coordY( k ) > lstrat_ym )then
                lcoordZ( k ) = coordZ( nbPts )
                id = lvertexID( nbPts )
                soil_thick( k ) = soil_thick( id )
                do ks  = 1, totgrn
                    stratLays( p, k, ks ) = stratLays( p, id, ks )
                enddo
            endif
        enddo

        return

    end subroutine update_creep_borders
    ! ============================================================================
    !> Subroutine update_creep_borders_elevation
    !!
    !<
    ! ============================================================================
    subroutine update_creep_borders_elevation

        integer :: k

        ! Update borders large DEM
        do k = 1, lnbPts
            ! South
            if( k <= lstrat_X )then
                lcoordZ( k ) = lcoordZ( k + 2 * lstrat_X )
            elseif( k <= 2 * lstrat_X  )then
                lcoordZ( k ) = lcoordZ( k + lstrat_X )
            ! North
            elseif( k >  lnbPts - 2 * lstrat_X )then
                lcoordZ( k ) = lcoordZ( k - lstrat_X )
            elseif( k >  lnbPts - lstrat_X )then
                lcoordZ( k ) = lcoordZ( k - 2 * lstrat_X )
            ! West
            elseif( lcoordX( k ) == lstrat_xo )then
                lcoordZ( k ) = lcoordZ( k + 2 )
            elseif( lcoordX( k ) == lstrat_xo + strat_dx )then
                lcoordZ( k ) = lcoordZ( k + 1 )
            ! East
            elseif( lcoordX( k ) == lstrat_xm )then
                lcoordZ( k ) = lcoordZ( k - 2 )
            elseif( lcoordX( k ) == lstrat_xm - strat_dx )then
                lcoordZ( k ) = lcoordZ( k - 1 )
            endif
            ! SW corner
            if( lcoordX( k ) < strat_xo .and. lcoordY( k ) < strat_yo )then
                lcoordZ( k ) = coordZ( 1 )
            endif
            ! SE corner
            if( lcoordX( k ) > strat_xm .and. lcoordY( k ) < strat_yo )then
                lcoordZ( k ) = coordZ( strat_X )
            endif
            ! NW corner
            if( lcoordX( k ) < strat_xo .and. lcoordY( k ) > strat_ym )then
                lcoordZ( k ) = coordZ( nbPts - strat_X + 1 )
            endif
            ! NE corner
            if( lcoordX( k ) > strat_xm .and. coordY( k ) > lstrat_ym )then
                lcoordZ( k ) = coordZ( nbPts )
            endif
            newz( k ) = lcoordZ( k )
        enddo

        return

    end subroutine update_creep_borders_elevation
    ! ============================================================================

end module hillslp
! ============================================================================
