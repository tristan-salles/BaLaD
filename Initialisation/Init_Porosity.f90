! ============================================================================
! Name        : Init_Porosity.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Init_Porosity.f90
!
! Description: it is used to read porosity look-up table and allocate porosity values over the
! mesh vertices. It is also used to modify mesh coordinates based on the induced compaction.
!
!<
! ============================================================================
module compaction

    use parallel
    use file_data
    use mesh_data
    use forces_data

    implicit none

contains

    ! ============================================================================
    !> Subroutine cmpt_compaction()
    !! calculates new porosities for each sediment cell caused by
    !! compaction due to an additional deposited layer.
    !! This subroutine calculates the lithostatic pressure of the overlying column
    !! and calls porosity_function() which calculates the porosity as a function of grain
    !! size distribution and effective pressure.
    !! Then it moves the topographic elevation by the compactional subsidence
    !! caused by compaction of all layers below the surface.
    !<
    ! ============================================================================
    subroutine cmpt_compaction

        integer :: kn, kl, nbl, ks, layerID
      
        real( tkind ) :: th, Plith
        real( tkind ), dimension( lnbPts ) :: newz

        newz = -1e6_8

        ! Loop through surface faces and calculate porosity for the new top deposit layer.
        call compute_top_layer_porosity

        ! Update porosities of underlying column due to compaction of newly deposited sediment
        do kn = 1, lnbPts
            nbl = v_nblays( kn )
            layerID = v_LaysID( kn, nbl )
            if( layerID == layID )then
                th = 0.0_8
                do ks = 1, totgrn
                    th = th + stratLays( nbl, kn, ks )
                enddo
                if( th > 0.0_8 )then
                    ! Lithostatic pressure
                    Plith = 0.0_8
                    ! Get porosity from top to bottom layers
                    do kl = nbl - 1, 1, -1
                        call compute_underlying_porosity( kl, kn, Plith )
                    enddo
                endif
            endif
            newz( kn ) = lcoordZ( kn )
        enddo

        do kn = 1, nbPts
            coordZ( kn ) = lcoordZ( lvertexID( kn ) )
        enddo

        return

    end subroutine cmpt_compaction
    ! ============================================================================
    !> Subroutine update_top_layer_porosity()
    !! update top layer porosities.
    !<
    ! ============================================================================
    subroutine  update_top_layer_porosity( k, kdem )

        integer :: kdem, k, laynb, ks, nbl

        real( tkind ) :: in_poro, oldh, newh

        nbl = v_nblays( k )
        laynb = v_LaysID( k, nbl )
        oldh = 0.0_8
        newh = 0.0_8

        ! Get porosity value and update sediment class and overall layer thicknesses
        if( laynb >= 1 )then

            do ks = 1, totgrn
                oldh = oldh + stratLays( nbl, k, ks )
                in_poro =  porosityLays(  nbl, k, ks  )
                porosityLays( nbl, k, ks  ) = porosity_function( ks, 0.0_8, in_poro )
                ! In case the porosity was not set for the considered layer virtually increase
                ! the sediment thickness to take into account the void space based on the
                ! newly calculated porosity value
                if( in_poro == 0.0_8 )then
                    stratLays( nbl, k, ks ) = stratLays( nbl, k, ks ) * &
                        ( 1 + porosityLays( nbl, k, ks ) )
                ! In case the porosity has decreased then decrease the sediment thickness
                ! accordingly to take into account the compaction of the layer.
                elseif( in_poro > porosityLays( nbl, k, ks ) )then
                    stratLays( nbl, k, ks ) = stratLays( nbl, k, ks ) * &
                        ( 1 + porosityLays( nbl, k, ks ) - in_poro )
                endif
                newh = newh + stratLays( nbl, k, ks )
            enddo

            ! Update top layer elevation
            lcoordZ( k ) = lcoordZ( k ) + newh - oldh
            coordZ( kdem ) = lcoordZ( k )

        endif

        return

    end subroutine  update_top_layer_porosity
    ! ============================================================================
    !> Subroutine compute_top_layer_porosity()
    !! computes new layer porosities.
    !<
    ! ============================================================================
    subroutine  compute_top_layer_porosity

        integer :: k, laynb, ks, nbl

        real( tkind ) :: in_poro, oldh, newh

        do k = 1, lnbPts
            nbl = v_nblays( k )
            laynb = v_LaysID( k, nbl )

            if( laynb == layID )then

                ! Get porosity value and update sediment class and overall layer thicknesses
                oldh = 0.0_8
                newh = 0.0_8
                do ks = 1, totgrn
                    oldh = oldh + stratLays( nbl, k, ks )
                    in_poro =  porosityLays( nbl, k, ks )
                    porosityLays(  nbl, k, ks  ) = porosity_function( ks, 0.0_8, in_poro )
                    ! In case the porosity was not set for the considered layer virtually increase
                    ! the sediment thickness to take into account the void space based on the
                    ! newly calculated porosity value
                    if( in_poro == 0.0_8 )then
                        stratLays( nbl, k, ks ) = stratLays( nbl, k, ks ) * &
                            ( 1 + porosityLays( nbl, k, ks ) )
                    ! In case the porosity has decreased then decrease the sediment thickness
                    ! accordingly to take into account the compaction of the layer.
                    elseif( in_poro > porosityLays( nbl, k, ks ) )then
                        stratLays( nbl, k, ks ) = stratLays( nbl, k, ks ) * &
                            ( 1 + porosityLays( nbl, k, ks ) - in_poro )
                    endif
                    newh = newh + stratLays( nbl, k, ks )
                enddo

                ! Update top layer elevation
                lcoordZ( k ) = lcoordZ( k ) + newh - oldh

            endif

        enddo

        return

    end subroutine  compute_top_layer_porosity
    ! ============================================================================
    !> Subroutine compute_underlying_porosity()
    !! calculates porosities and compaction values for vertices which coincided with top ones.
    !! \param kn, kl, nid, Plith, mass, nbl
    !<
    ! ============================================================================
    subroutine  compute_underlying_porosity( kl, nid, Plith )

        integer :: ks, kl, nid

        real( tkind ) :: in_poro, mass, Plith, oldh, newh

        mass = 0.0_8
        do ks = 1, totgrn
            ! Compute upper layer sediment mass
            mass = mass + stratLays(  kl+1, nid,  ks  ) * &
                sediment( ks )%density * ( 1.0_8 - porosityLays( kl+1, nid, ks ) )
            ! Compute mass due to water in sediment pore space
            mass = mass + stratLays(  kl+1, nid, ks ) * &
                sea_density * porosityLays( kl+1, nid, ks )
        enddo

        ! If cell has sediment, increment cumulative lithostatic pressure
        if( mass > tor )then

            ! Compute recursively lithostatic pressure for current sediment layer depth
            Plith = Plith + 1.0e-6_8 * mass * gravity

            oldh = 0.0_8
            newh = 0.0_8
            do ks = 1, totgrn
                oldh = oldh + stratLays( kl, nid, ks )
                in_poro =  porosityLays( kl, nid, ks )
                porosityLays( kl, nid, ks ) = porosity_function( ks, Plith, in_poro )
                ! In case the porosity was not set for the considered layer virtually increase
                ! the sediment thickness to take into account the void space based on the
                ! newly calculated porosity value
                if( in_poro == 0.0_8 )then
                    stratLays( kl, nid, ks ) = stratLays( kl, nid, ks ) * &
                        ( 1 + porosityLays( kl, nid, ks ) )
                ! In case the porosity has decreased then decrease the sediment thickness
                ! accordingly to take into account the compaction of the layer.
                elseif( in_poro > porosityLays( kl, nid, ks ) )then
                    stratLays( kl, nid, ks ) = stratLays( kl, nid, ks ) * &
                        ( 1 + porosityLays( kl, nid, ks ) - in_poro )
                endif
                newh = newh + stratLays( kl, nid, ks )
            enddo

            ! Update top layer elevation
            lcoordZ( kl ) = lcoordZ( kl ) + newh - oldh

        endif

        return

    end subroutine  compute_underlying_porosity
    ! ============================================================================
    !> Subroutine porosity_init()
    !! Calculates initial deposit porosities for each sediment cell.
    !<
    ! ============================================================================
    subroutine porosity_init

        integer :: k, kl, ks, nbl, laynb

        real( tkind ) :: in_poro, Plith, mass

        ! Loop through surface faces and calculate porosity for the initial deposit layers.
        do k = 1, lnbPts
            nbl = v_nblays( k )
            laynb = v_LaysID( k, nbl )
            Plith = 0.0_8
            do kl = nbl, 1, -1
                mass = 0.0_8
                laynb = v_LaysID( k, kl )
                if( laynb <= layID .and. laynb >= 1  )then
                    if( layID > laynb )then
                        do ks = 1, totgrn
                            mass = mass + stratLays( kl+1, k, ks ) * &
                                sediment( ks )%density * ( 1 - porosityLays( kl+1, k, ks  ) )
                            ! Mass due to water in pore space
                            mass = mass + stratLays( kl+1, k, ks ) * &
                                sea_density * porosityLays( kl+1, k, ks  )
                        enddo
                    endif
                    Plith = Plith + 1.0e-6_8 * mass * gravity
                    do ks = 1, totgrn
                        in_poro =  porosityLays( kl, k, ks )
                        porosityLays( kl, k, ks  ) = porosity_function( ks, Plith, in_poro )
                        if( porosityLays( kl, k, ks  ) > initdepporo )then
                            porosityLays( kl, k, ks  ) = initdepporo
                        endif
                    enddo
                endif
            enddo
        enddo

        return

    end subroutine porosity_init
    ! ============================================================================
    !> Subroutine porosity_function()
    !! Computes a variable porosity based upon effective subsurface pressure.
    !! \param dz_poro, prate, in_poro, out_poro
    !<
    ! ============================================================================
    function porosity_function( grn, LithPressure, previous_porosity ) result( new_porosity )

        integer :: grn, k

        real( tkind ) :: LithPressure, previous_porosity, new_porosity, interpolate_porosity

        new_porosity = previous_porosity

        loop_over_pressure: do k = 2, gporo%ePnb
            if( effPressure( k - 1 ) <= LithPressure .and. &
                effPressure( k ) > LithPressure ) exit loop_over_pressure
        enddo loop_over_pressure

        ! Calculate interpolated porosity based on effective subsurface pressure.
        interpolate_porosity = ( porosity( grn, k ) - porosity( grn, k - 1 ) ) / ( effPressure( k ) - effPressure( k - 1 ) )
        interpolate_porosity = interpolate_porosity * ( LithPressure - porosity( grn, k - 1 ) )
        interpolate_porosity = interpolate_porosity + porosity( grn, k - 1 )

        if( interpolate_porosity < previous_porosity .or. previous_porosity == 0.0_8 )&
            new_porosity = interpolate_porosity

    end function porosity_function
    ! ============================================================================
    !> Subroutine assign_porosity_table()
    !! Assigns the porosity look-up table values from the input file.
    !<
    ! ============================================================================
    subroutine assign_porosity_table

        logical :: found
        integer :: i, maxl, l, d
        integer :: iu, ios
        character( len=200 ) :: line

        ! Open porosity look-up table file
        iu = 80 + iam
        inquire( file=fporosity, exist=found )
        if(found)then
            open(iu,file=fporosity,status="old",action="read",iostat=ios)
            rewind(iu)
        else
            write(*,*)'Warning: the input file for porosity looku-up table cannot be found'
            write(*,*)' The code is looking for the following name:',trim(fporosity)
            stop
        endif

        ! Number of lines in the porosity file
        maxl = 1 + totgrn
        l = 0
        d = 0

        ! Read parameters
        do while ( l < maxl )
            read(iu,'(a200)') line
            i = len_trim( line )
            if( line(i:i) == char(13) ) i = i-1
            if( line(1:2) /= '* ')then
                ! Read effective pressure
                if( l == 0 )then
                    read( line(2:i),* ) effPressure( 1:gporo%ePnb )
                ! Read porosity values
                elseif( l < maxl - 1)then
                    d = d + 1
                    read( line(2:i),* ) porosity( d,1:gporo%ePnb )
                endif
                l = l + 1
            endif
        enddo

        ! Close file
        close( iu )

        if( d /= totgrn ) print*,'Something went wrong when reading porosity look-up table',d,totgrn

        return

    end subroutine assign_porosity_table
  ! ============================================================================

end module compaction
