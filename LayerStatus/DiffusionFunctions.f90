! ============================================================================
! Name        : DiffusionFunctions.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file DiffusionFunctions.f90
!!
!! DiffusionFunctions set of functions used for diffusion transport.
!!
!<
! ============================================================================
module diffusionfct

    use parallel
    use flow_data
    use mesh_data
    use forces_data

    implicit none

    logical, parameter :: parallel_diff = .false.

    real( tkind ), dimension( : ), allocatable :: depdif

contains

    ! ============================================================================
    !> Subroutine build_current_layer_creep_grid
    !! Subroutine build_current_layer_creep_grid determines nodes with deposition.
    !<
    ! ============================================================================
    subroutine build_current_layer_diffusion_grid

        integer :: i, ks, lid, lb, ld, p

        real( tkind ) :: dh, ddc, tth

        if( parallel_diff .and. .not. allocated( depdif ) ) allocate( depdif( lnbPts ) )

        ! Update local diffusion grid top elevation
        creepz = lcoordZ
        do i = 1 , nbPts

            ! Get the mesh ID
            lid = lvertexID( i )

            if( lid > 0 )then

                lb = v_nbLays( lid )
                ld = v_LaysID( lid, lb )
                creepz( lid ) = lcoordZ( lid )
                ! South
                if( lcoordY( lid ) < strat_yo .and. boundcond( 1 ) == 2 )then
                    creepz( lid ) = 1.e6_8
                endif
                ! North
                if( lcoordY( lid ) > strat_ym .and. boundcond( 2 ) == 2 )then
                    creepz( lid ) = 1.e6_8
                endif
                ! West
                if( lcoordX( lid ) < strat_xo .and. boundcond( 3 ) == 2 )then
                    creepz( lid ) = 1.e6_8
                endif
                ! East
                if( lcoordX( lid ) > strat_xm .and. boundcond( 4 ) == 2 )then
                    creepz( lid ) = 1.e6_8
                endif

                dh = 0.0_8

                ! In case the current layer is active
                if(  layID == ld )then

                    do ks = 1, totgrn
                        dstart( lid, ks ) = 0.0_8
                        dh = dh + stratLays( lb, lid, ks )
                        ! Get the amount deposited since last creep transport
                        ddc = stratLays( lb, lid, ks ) - top_sedprev( lid, ks )
                        ! If something was deposited put it as possible material for
                        ! creep transport
                        if( ddc > 0.0_8 )then
                            dstart( lid, ks ) = ddc / dble( diff_nb )
                            top_sedh( lid, ks ) = top_sedprev( lid,ks )
                        ! Otherwise update creep sediment stack
                        else
                            top_sedh( lid, ks ) = stratLays( lb, lid, ks )
                            top_sedprev( lid ,ks ) = top_sedh( lid, ks )
                        endif
                    enddo

                    ! Now update creep grid elevation
                    tth = 0.0_8
                    do ks = 1, totgrn
                        tth = tth + top_sedh( lid, ks )
                    enddo
                    creepz( lid ) = creepz( lid ) - ( dh - tth )
                ! If the top layer is not the active one
                else
                    top_sedprev( lid, 1:totgrn ) = 0.0_8
                    top_sedh( lid, 1:totgrn ) = 0.0_8
                    creepz( lid ) = lcoordZ( lid )
                    ! South
                    if( lcoordY( lid ) < strat_yo .and. boundcond( 1 ) == 2 )then
                        creepz( lid ) = 1.e6_8
                    endif
                    ! North
                    if( lcoordY( lid ) > strat_ym .and. boundcond( 2 ) == 2 )then
                        creepz( lid ) = 1.e6_8
                    endif
                    ! West
                    if( lcoordX( lid ) < strat_xo .and. boundcond( 3 ) == 2 )then
                        creepz( lid ) = 1.e6_8
                    endif
                    ! East
                    if( lcoordX( lid ) > strat_xm .and. boundcond( 4 ) == 2 )then
                        creepz( lid ) = 1.e6_8
                    endif
                endif

            endif
        enddo

        do lid = 1, lnbPts
            ! South
            if( lcoordY( lid ) < strat_yo .and. boundcond( 1 ) == 2 )then
                creepz( lid ) = 1.e6_8
            endif
            ! North
            if( lcoordY( lid ) > strat_ym .and. boundcond( 2 ) == 2 )then
                creepz( lid ) = 1.e6_8
            endif
            ! West
            if( lcoordX( lid ) < strat_xo .and. boundcond( 3 ) == 2 )then
                creepz( lid ) = 1.e6_8
            endif
            ! East
            if( lcoordX( lid ) > strat_xm .and. boundcond( 4 ) == 2 )then
                creepz( lid ) = 1.e6_8
            endif
        enddo

        return

    end subroutine build_current_layer_diffusion_grid
    ! ============================================================================
    !> Subroutine find_maximum_elevation
    !! Subroutine find_maximum_elevation determines what the maximum topographic height
    !!  that each grid node can support whilst being stable and deposit the corresponding elevation.
    !! \param ks
    !<
    ! ============================================================================
    subroutine find_maximum_elevation( ks )

        integer :: k, locid, ks, i, lid

        real( tkind ) :: dist, elev, topnew, topmax
        real( tkind ), dimension( lnbPts ) :: zcreep

        zcreep = creepz

        ! Get the maximum topographic elevation for each node
        do i = 1, nbPts
            lid = lvertexID( i )
            if( lid > 0 )then
                if( difo( lid ) > 0.0_8 )then

                    topmax = toplimit

                    ! Loop over neighboring nodes
                    do k = 1, 8
                        locid = lngbID( lid, k )

                        if( locid > 0 )then

                            elev = zcreep( locid )
                            dist = ( lcoordY( locid ) - lcoordY( lid ) )**2 + &
                                ( lcoordX( locid ) - lcoordX( lid ) )**2
                            dist = sqrt( dist )

                            ! In case node is below sea level
                            if( elev < gsea%actual_sea )then
                                topnew = elev + sediment( ks )%slp_marine * dist
                                if( topnew > gsea%actual_sea + tor .and. sediment( ks )%slp_marine > 0.0_8 )&
                                    topnew = gsea%actual_sea + sediment( ks )%slp_aerial * &
                                    ( dist + ( elev - gsea%actual_sea ) / sediment( ks )%slp_marine )

                            ! Otherwise node is above water
                            else
                                topnew = elev + sediment( ks )%slp_aerial * dist
                            endif

                            topmax = min( topmax, topnew )
                        endif
                    enddo

                    if( topmax - zcreep( lid ) > 0.0_8 )then
                        topnew = zcreep( lid ) + difo( lid )
                        if( topnew <= topmax )then
                            cdif( lid ) = 1.0_8
                            creepz( lid ) = topnew
                        else
                            cdif( lid ) = ( topmax - zcreep( lid ) )
                            cdif( lid ) = cdif( lid ) / ( topnew - zcreep( lid ) )
                            creepz( lid ) = topmax
                        endif
                    endif
                    depo( lid, ks ) = depo( lid, ks ) +  difo( lid ) * cdif( lid )
                endif
            endif
        enddo

        return

    end subroutine find_maximum_elevation
    ! ============================================================================
    !> Subroutine parallel_find_maximum_elevation
    !! Subroutine find_maximum_elevation determines what the maximum topographic height
    !!  that each grid node can support whilst being stable and deposit the corresponding elevation.
    !! \param ks
    !<
    ! ============================================================================
    subroutine parallel_find_maximum_elevation( ks )

        integer :: k, locid, ks, i, lid

        integer :: IDlft, IDrgt
        integer :: rq1, rq2

        integer, dimension( mpi_status_size ) :: stt1, stt2

        real( tkind ) :: dist, elev, topnew, topmax

        real( tkind ), dimension( lstrat_X ) :: zdif
        real( tkind ), dimension( lnbPts ) :: zcreep

        zcreep = creepz

        ! Get the maximum topographic elevation for each node
        do i = 1, local_lnbPts
            lid = global_lnid( i )
            if( lid > 0 )then
                if( difo( lid ) > 0.0_8 )then

                    topmax = toplimit
                    ! Loop over neighboring nodes
                    do k = 1, 8
                        locid = lngbID( lid, k )

                        if( locid > 0 )then

                            elev = zcreep( locid )
                            dist = ( lcoordY( locid ) - lcoordY( lid ) )**2 + &
                                ( lcoordX( locid ) - lcoordX( lid ) )**2
                            dist = sqrt( dist )

                            ! In case node is below sea level
                            if( elev < gsea%actual_sea )then
                                topnew = elev + sediment( ks )%slp_marine * dist
                                if( topnew > gsea%actual_sea + tor .and. sediment( ks )%slp_marine > 0.0_8 )&
                                    topnew = gsea%actual_sea + sediment( ks )%slp_aerial * &
                                    ( dist + ( elev - gsea%actual_sea ) / sediment( ks )%slp_marine )

                            ! Otherwise node is above water
                            else
                                topnew = elev + sediment( ks )%slp_aerial * dist
                            endif

                            topmax = min( topmax, topnew )
                        endif
                    enddo

                    if( topmax - zcreep( lid ) > 0.0_8 )then
                        topnew = zcreep( lid ) + difo( lid )
                        if( topnew <= topmax )then
                            cdif( lid ) = 1.0_8
                            creepz( lid ) = topnew
                        else
                            cdif( lid ) = ( topmax - zcreep( lid ) )
                            cdif( lid ) = cdif( lid ) / ( topnew - zcreep( lid ) )
                            creepz( lid ) = topmax
                        endif
                    endif
                    depdif( lid ) = depdif( lid ) +  difo( lid ) * cdif( lid )
                endif
            endif
        enddo

        ! Send information to partition
        if( iam < gproc - 1 )then

            ! Send top row to other partition
            IDlft = ( part_rows( iam + 1 ) - 1 ) * lstrat_X + 1
            IDrgt = IDlft + lstrat_X - 1

            call mpi_isend( creepz( IDlft:IDrgt), lstrat_X, dbl_type, iam+1, 29, &
                badlands_comm_world, rq1, ierr )
            call mpi_request_free( rq1, ierr )

            ! Reveive top row + 1 from other partition
            zdif = 0.0_8
            call mpi_irecv( zdif, lstrat_X, dbl_type, iam+1, 28, &
                badlands_comm_world,  rq2, ierr )
            call mpi_wait( rq2, stt2, ierr )

            IDlft = IDlft + lstrat_X - 1
            do k = 1, lstrat_X
                creepz( IDlft + k ) = max( creepz( IDlft + k ), zdif( k ) )
            enddo
        endif

        if( iam > 0 )then

            !  Send row 1 to other partition
            IDlft = part_rows( iam ) * lstrat_X + 1
            IDrgt = IDlft + lstrat_X - 1

            call mpi_isend( creepz( IDlft:IDrgt), lstrat_X, dbl_type, iam-1, 28, &
                badlands_comm_world, rq2, ierr )
            call mpi_request_free( rq2, ierr )

            ! Receive row 1 - 1 from other partition
            zdif = 0.0_8
            call mpi_irecv( zdif, lstrat_X, dbl_type, iam-1, 29, &
                badlands_comm_world,  rq1, ierr )
            call mpi_wait( rq1, stt1, ierr )

            IDlft = IDlft - lstrat_X - 1
            do k = 1, lstrat_X
                creepz( IDlft + k ) = max( creepz( IDlft + k ), zdif( k ) )
            enddo
        endif

        return

    end subroutine parallel_find_maximum_elevation
    ! ============================================================================
    !> Subroutine get_available_space_remaining
    !! Subroutine get_available_space_remaining determines the space available for diffusion
    !! around a specific node.
    !<
    ! ============================================================================
    subroutine get_available_space_remaining( vd, spacsum, spac, down, ndown )

        integer :: k, locid, ndown, vd
        integer, dimension( 8 ) :: down

        real( tkind ) :: sp, dist, spacsum, slop
        real( tkind ), dimension( 8 ) :: spac

        ! Loop over neighboring cell
        do k = 1, 8
            locid = lngbID( vd, k )
            sp = creepz( vd ) - creepz( locid )
            dist = ( lcoordY( locid ) - lcoordY( vd ) )**2 + &
             ( lcoordX( locid ) - lcoordX( vd ) )**2
            dist = sqrt( dist )

            ! If there is some space
            if( sp > 0 )then
                slop = sp / dist
                if( slop >= minimum_slp )then
                    ndown = ndown + 1
                    down( k ) = 1
                    spac( k ) = sp
                    spacsum = spacsum + sp
                endif
            endif

        enddo

        return

    end subroutine get_available_space_remaining
    ! ============================================================================
    !> Subroutine update_sediment_layer_diffusion
    !! Subroutine update_sediment_layer_diffusion updates top sediment layers information
    !! as well as the global stratal on the strata grid.
    !<
    ! ============================================================================
    subroutine update_sediment_layer_diffusion

        integer :: k, ks, laynb, toplay, id

        real( tkind ) :: th, old_thick

        ! Update the stratigraphic layer
        do id = 1, nbPts
            k = lvertexID( id )
            laynb = v_nbLays( k )
            toplay = v_LaysID(  k, laynb )
            th = 0.0_8

            ! Get the updated thickness
            old_thick = 0.0_8
            do ks = 1, totgrn
                if( layID == toplay ) &
                    old_thick = old_thick + stratLays( laynb, k, ks )
                top_sedh( k ,ks ) = top_sedprev( k ,ks ) + depo( k, ks )
                top_sedprev( k ,ks ) = top_sedh( k ,ks )
                th =  th + top_sedh( k ,ks )
            enddo

            ! In case there is something to deposit
            if( th > tor )then

                ! If the top layer ID corresponds to the current layer ID
                if( layID == toplay )then
                    th = 0.0_8
                    ! Update layer sediment class thicknesses
                    if( gporo%compaction ) porosityLays( laynb, k, 1:totgrn ) = 0.0_8
                    do ks = 1, totgrn
                        th = th + top_sedh( k, ks )
                        stratLays( laynb, k, ks ) = top_sedh( k, ks )
                        if( gporo%compaction .and. top_sedh( k, ks ) > 0.0_8 ) &
                            porosityLays( laynb, k, ks ) = porosity( ks, 1 )
                    enddo
                    lcoordZ(  k ) = lcoordZ( k ) + ( th - old_thick )
                    soil_thick( k ) = soil_thick( k ) + ( th - old_thick )
                    if( soil_thick( k ) > dtb_marine .and. lcoordZ( k ) < gsea%actual_sea ) &
                        soil_thick( k ) = dtb_marine
                    if( soil_thick( k ) > dtb_aerial .and. lcoordZ( k ) >= gsea%actual_sea ) &
                        soil_thick( k ) = dtb_aerial
                    if( soil_thick( k ) < 0.0_8 ) soil_thick( k ) = 0.0_8

                ! Otherwise create a new layer
                else
                    v_nbLays( k ) = v_nbLays( k ) + 1
                    laynb = v_nbLays( k )
                    if( laynb > nbLays )then
                        print*,'Something went wrong went adding layer after diffusion'
                        call mpi_finalize( ierr )
                        stop
                    endif
                    v_LaysID( k, laynb ) = layID
                    ! Update layer sediment class thicknesses
                    th = 0.0_8
                    if( gporo%compaction ) porosityLays( laynb, k, 1:totgrn ) = 0.0_8
                    do ks = 1, totgrn
                        th = th + top_sedh( k, ks )
                        stratLays( laynb, k, ks ) = top_sedh( k, ks )
                        if( gporo%compaction .and. top_sedh( k, ks ) > 0.0_8 ) &
                            porosityLays( laynb, k, ks ) = porosity( ks, 1 )
                    enddo
                    lcoordZ(  k ) = lcoordZ( k ) + th
                    soil_thick( k ) = soil_thick( k ) + th
                    if( soil_thick( k ) > dtb_marine .and. lcoordZ( k ) < gsea%actual_sea ) &
                        soil_thick( k ) = dtb_marine
                    if( soil_thick( k ) > dtb_aerial .and. lcoordZ( k ) >= gsea%actual_sea ) &
                        soil_thick( k ) = dtb_aerial
                    if( soil_thick( k ) < 0.0_8 ) soil_thick( k ) = 0.0_8
                endif

            ! In case the layer is empty
            else
                ! If a layer was there delete it
                if( layID == toplay .and. old_thick > 0.0_8 )then
                    ! Update layer sediment class thicknesses
                    stratLays( laynb, k, 1:totgrn ) = 0.0_8
                    if( gporo%compaction ) &
                        porosityLays( laynb, k, 1:totgrn ) = 0.0_8
                    lcoordZ(  k ) = lcoordZ( k ) - old_thick
                    v_nbLays( k ) = v_nbLays( k ) - 1
                    soil_thick( k ) = soil_thick( k ) - old_thick
                    if( soil_thick( k ) > dtb_marine .and. lcoordZ( k ) < gsea%actual_sea ) &
                        soil_thick( k ) = dtb_marine
                    if( soil_thick( k ) > dtb_aerial .and. lcoordZ( k ) >= gsea%actual_sea ) &
                        soil_thick( k ) = dtb_aerial
                    if( soil_thick( k ) < 0.0_8 ) soil_thick( k ) = 0.0_8
                endif
            endif
        enddo

        ! Update DEM grid topography and soil thickness
        do k = 1, nbPts
            coordZ( k ) = lcoordZ( lvertexID( k ) )
        enddo

        return

    end subroutine update_sediment_layer_diffusion
    ! ============================================================================

end module diffusionfct

