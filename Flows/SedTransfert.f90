! ============================================================================
! Name        : SedTransfert.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file SedTransfert.f90
!
! Description :  This module records changes in stratigraphy from sediment fluxes.
!
!<
! ============================================================================
!> Module sedtrans
!<
module sedtrans

    use fwpath
    use parallel
    use morpho
    use interpol
    use regolith
    use file_data
    use ode_data
    use flow_data
    use time_data
    use mesh_data
    use forces_data

    implicit none

contains

    ! ============================================================================
    !> Subroutine manage_sedfluxes
    !! This subroutine performs sedimentary flux on the stratigraphic grid
    !!
    !! k: ID of the stratigraphic node
    !! pp: sediment ID number
    !! erodeflx: thickness of sediments that can be eroded
    !! depositflx: thickness of sediments that can be deposited
    !! inflx: thickness of sediments that is passing through the cell
    !! soil_limit: is regolith thickness limiting the erosion?
    !! hrego: thickness of regolith material
    !!
    !! outflx: sediment composition leaving the considered node
    !<
    ! ============================================================================
    subroutine manage_sedfluxes( k, pp, erodeflx, depositflx, inflx, soil_limit, hrego, outflx )

        logical :: soil_limit

        integer :: k, ks, id, p, nlay, pp

        real( tkind ) :: eromax, depomax, hrego, maxh
        real( tkind ) :: totero, totdep, frac, dh, stillsed
        real( tkind ) :: depositflx, erodeflx, inflx, outflx
        real( tkind ) :: hs, erode, remain, sedh

        id = lvertexID( k )

        outflx = 0.0_8

        ! First get the maximum possible erosion to :
        ! 1- prevent formation of hole
        ! 2- preclude inversing landscape
        call get_erosion_height( k, eromax )
        eromax =  eromax * 0.8_8

        if( soil_limit .and. eromax > hrego )then
            eromax = hrego
        endif

        ! Then get the maximum possible deposition to ensure
        ! landscape stability based on neighborhood elevations
        maxh = 1.0e6_8
        maxh = newz( id ) + 1.e6_8
        do p = 1, 8
            if( ( newz( lngbID( id, p ) ) - newz( id ) ) > 0.01_8 ) &
                maxh = min( maxh, newz( lngbID( id, p ) ) )
        enddo
        depomax = ( maxh - newz( id ) ) * 0.8_8
        if( depomax < 0.0_8 ) depomax = 0.0_8

        ! Sum total fluxes
        totero = erodeflx
        totdep = depositflx

        ! Limit erosion
        if( totero > eromax )then
            frac = eromax / totero
            erodeflx = frac * erodeflx
        endif
        
        ! Update node evolution from erosion fluxes
        if( totero > 0.0_8 )then

            ! Get stratigraphic layers nb
            nlay = v_nbLays( id )

            stillsed = 0.0_8
            hs = erodeflx
            erodeflx = 0.0_8

            remain = 0.0_8
            erode = 0.0_8
            sedh = 0.0_8

            ! Get the amount to erode for the considered layer
            if( gporo%compaction )then
                sedh = stratLays( nlay, id, pp ) * ( 1.0_8 - porosityLays( nlay, id, pp ) )
                sedh = sedh - hs
            else
                sedh = stratLays( nlay, id, pp ) - hs
            endif

            ! In case there isn't enough sediment in the layer adjust the
            ! sediment that needs to be eroded for the considered layers
            if( sedh <= 0.0_8 )then
                sedh = 0.0_8
                if( gporo%compaction )then
                    erode = stratLays( nlay, id, pp ) * ( 1.0_8 - porosityLays( nlay, id, pp ) )
                    porosityLays( nlay, id, pp ) = 0.0_8
                else
                    erode = stratLays( nlay, id, pp )
                endif
                remain = 0.0_8
                ! Thickness of sediment still needs to be eroded
                ! in underlying layer
                hs = hs - erode

            ! In case there is still sediments in the layer update the
            ! remaining thickness accordingly
            else
                if( gporo%compaction )then
                    remain = sedh / ( 1.0_8 - porosityLays( nlay, id, pp ) )
                else
                    remain =  sedh
                endif
                ! No more sediment to erode
                hs = 0.0_8
            endif

            ! Get the total thickness of sediment present in
            ! the current layer
            do ks = 1, totgrn
                if( ks /= pp )then
                    stillsed = stillsed + stratLays( nlay, id, ks )
                else
                    stillsed = stillsed + remain
                endif
            enddo

            ! If there is still some sediment on the layer, consider that they will
            ! have an armoring effect on the underlying sedimentary layers and
            ! stop the erosion
            dh = 0.0_8
            if( stillsed > 0.0_8 )then

                if( stratLays( nlay, id, pp ) < remain )then
                    print*,'Something went wrong when eroding stratigraphy in LEM.',iam
                    print*,pp,id,nlay,remain,v_nbLays( id ),stratLays( nlay, id, pp ),erodeflx
                    stop
                endif
                dh = dh + stratLays( nlay, id, pp ) - remain
                erodeflx = erodeflx + stratLays( nlay, id, pp ) - remain
                stratLays( nlay, id, pp ) = remain

            ! If there is no more sediment, the regolith mantle has been totally eroded
            elseif( stillsed <= 0.0_8 )then

                dh = dh + erode
                erodeflx = erodeflx + stratLays( nlay, id, pp ) - remain
                stratLays( nlay, id, 1:totgrn ) = 0.0_8
                if( gporo%compaction ) porosityLays( nlay, id, 1:totgrn ) = 0.0_8

            endif

            soil_thick( id ) = soil_thick( id ) - dh
            if( soil_thick( id ) < 0.0_8 ) soil_thick( id ) = 0.0_8
            lcoordZ( id ) = lcoordZ( id ) - dh
            coordZ( k ) = lcoordZ( id )
            if( stillsed <= 0.0_8 ) nlay = nlay - 1
            v_nbLays( id ) = nlay

        endif

        ! From influxes check maximum deposition values
        ! and get outfluxes
        dh = 0.0_8
        if( depositflx > 0.0_8 )then
            if( depositflx > inflx )then
                depositflx = inflx
                outflx = 0.0_8
            else
                outflx = inflx - depositflx
            endif
            dh = dh + depositflx
        elseif( erodeflx > 0.0_8 )then
            outflx = inflx + erodeflx
        else
            outflx = inflx
        endif

        ! Update stratigraphic layer and topography from deposit thickness
        if( dh > 0.0_8 )then
            ! The top layer corresponds to the current time layer
            if( layID == v_LaysID( id, v_nbLays( id ) ) )then
                ! Compute sedimentary changes
                if( gporo%compaction )then
                    depositflx = depositflx / ( 1.0_8 - porosity( pp, 1 ) )
                    porosityLays( v_nbLays( id ), id, pp ) = porosity( pp, 1 )
                endif
                stratLays( v_nbLays( id ), id, pp ) = stratLays( v_nbLays( id ), id, pp ) + depositflx

            ! Otherwise create a new layer
            else

                v_nbLays( id ) = v_nbLays( id ) + 1
                v_LaysID( id, v_nbLays( id ) ) = layID
                if( gporo%compaction )then
                    depositflx = depositflx / ( 1.0_8 - porosity( pp, 1 ) )
                    porosityLays( v_nbLays( id ), id, ks ) = porosity( pp, 1 )
                endif
                stratLays( v_nbLays( id ), id, pp ) = depositflx

            endif
            lcoordZ( id ) = lcoordZ( id ) + dh
            soil_thick( id ) = soil_thick( id ) + dh
            coordZ( k ) = lcoordZ( id )
        endif

       ! Update regolith thickness
        if( soil_thick( id ) > dtb_marine .and. lcoordZ( id ) < gsea%actual_sea ) &
            soil_thick( id ) = dtb_marine
        if( soil_thick( id ) > dtb_aerial .and. lcoordZ( id ) >= gsea%actual_sea ) &
            soil_thick( id ) = dtb_aerial

        return

    end subroutine manage_sedfluxes
    ! ============================================================================
    !> Subroutine get_erosion_height
    !! Subroutine get_erosion_height get the maximum erosion thickness
    !<
    ! ============================================================================
    subroutine get_erosion_height( id, eromax )

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

        return

    end subroutine get_erosion_height
    ! ============================================================================

end module sedtrans
! ============================================================================
