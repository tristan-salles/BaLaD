! ============================================================================
! Name        : RegolithFormation.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file RegolithFormation.f90
!
! Description :  This module constrains the rates of regolith production, colluvial
! transport, and deposition over geologic time scales.
!
!<
! ============================================================================
!> Module regolith
!<
module regolith

    use fwpath
    use parallel
    use morpho
    use interpol
    use file_data
    use ode_data
    use flow_data
    use time_data
    use mesh_data
    use forces_data

    implicit none

    save

    public :: quick_sort

contains

    ! ============================================================================
    !> Subroutine quick_sort
    !! This subroutine computes landscape evolution through time
    !! grabbed from A millers web site http://users.bigpond.net.au/amiller/
    !! Quick sort routine from:
    !! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    !! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    !! Modified by Alan Miller to include an associated integer array which gives
    !! the positions of the elements in the original order.
    !! pjr added module declaration
    !! mvr modified integer array to intent inout - may now be any integer
    !!     array that gets sorted along with associated real array
    !<
    ! ============================================================================
    recursive subroutine quick_sort(list, order)

        real( tkind ), dimension(:), intent( inout ) :: list
        integer, dimension(:), intent( out ) :: order
        integer :: i

        do i = 1, size(list)
            order( i ) = i
        enddo

        call quick_sort_1(1, size(list))

    contains

        ! ============================================================================
        recursive subroutine quick_sort_1(left_end, right_end)

            integer, intent( in ) :: left_end, right_end

            ! Local variables
            integer :: i, j, itemp
            real( tkind ) :: reference, temp
            integer, parameter :: max_simple_sort_size = 6

            if (right_end < left_end + max_simple_sort_size) then
                ! Use interchange sort for small lists
                call interchange_sort(left_end, right_end)
            else
                ! Use partition ("quick") sort
                reference = list((left_end + right_end)/2)
                i = left_end - 1; j = right_end + 1

                do
                    ! Scan list from left end until element >= reference is found
                    do
                        i = i + 1
                        if (list(i) >= reference) exit
                    end do
                    ! Scan list from right end until element <= reference is found
                    do
                        j = j - 1
                        if (list(j) <= reference) exit
                    end do

                    if (i < j) then
                        ! Swap two out-of-order elements
                        temp = list(i); list(i) = list(j); list(j) = temp
                        itemp = order(i); order(i) = order(j); order(j) = itemp
                    else if (i == j) then
                        i = i + 1
                        exit
                    else
                        exit
                    end if
                end do

                if (left_end < j) call quick_sort_1(left_end, j)
                if (i < right_end) call quick_sort_1(i, right_end)
            end if

        end subroutine quick_sort_1
        ! ============================================================================
        subroutine interchange_sort(left_end, right_end)

            integer, intent( in ) :: left_end, right_end

            !  Local variables
            integer :: i, j, itemp
            real( tkind ) :: temp

            do i = left_end, right_end - 1
                do j = i+1, right_end
                    if (list(i) > list(j)) then
                        temp = list(i); list(i) = list(j); list(j) = temp
                        itemp = order(i); order(i) = order(j); order(j) = itemp
                    end if
                end do
            end do

        end subroutine interchange_sort
       ! ============================================================================

    end subroutine quick_sort
    ! ============================================================================
    !> Subroutine getFluxWeight
    !! This subroutine computes the matrix of flow fraction
    !<
    ! ============================================================================
    subroutine getFluxWeight

        integer :: k, k2, p

        real( tkind ) :: sum, maxslope
        real( tkind ), dimension( 8 ) :: slopeArray


        if( .not. allocated( fluxweight ) ) allocate( fluxweight( lnbPts, 8 ) )

        do k = 1, lnbPts

            sum = 0.0_8
            maxslope = 0.0_8
            slopeArray = 0.0_8

            do p = 1, 8
                fluxweight( k, p ) = 0.0_8
                k2 = lngbID( k, p )
                if( k2 > 0 )then
                    if( lfilldem( k ) > lfilldem( k2 ) )then
                        slopeArray( p ) = ( lfilldem( k ) - lfilldem( k2 ) ) / strat_dx
                        if( mod( p, 2 ) == 0 ) slopeArray( p ) = slopeArray( p ) / sqrt( 2.0_8 )
                    endif
                endif
                sum = sum + slopeArray( p )
            enddo

            do p = 1, 8
                if( slopeArray( p ) > 0.0_8 )then
                    fluxweight( k, p ) =  slopeArray( p ) / sum
                endif
            enddo

        enddo

        return

    end subroutine getFluxWeight
    ! ============================================================================
    !> Subroutine getFluxWeightTopo
    !! This subroutine computes the matrix of flow fraction
    !<
    ! ============================================================================
    subroutine getFluxWeightTopo

        integer :: k, k2, p

        real( tkind ) :: sum, maxslope
        real( tkind ), dimension( 8 ) :: slopeArray


        if( .not. allocated( fluxweightdem ) ) allocate( fluxweightdem( lnbPts, 8 ) )

        do k = 1, lnbPts

            sum = 0.0_8
            maxslope = 0.0_8
            slopeArray = 0.0_8

            do p = 1, 8
                fluxweightdem( k, p ) = 0.0_8
                k2 = lngbID( k, p )
                if( k2 > 0 )then
                    if( lcoordZ( k ) > lcoordZ( k2 ) )then
                        slopeArray( p ) = ( lcoordZ( k ) - lcoordZ( k2 ) ) / strat_dx
                        if( mod( p, 2 ) == 0 ) slopeArray( p ) = slopeArray( p ) / sqrt( 2.0_8 )
                    endif
                endif
                sum = sum + slopeArray( p )
            enddo

            do p = 1, 8
                if( slopeArray( p ) > 0.0_8 )then
                    fluxweightdem( k, p ) =  slopeArray( p ) / sum
                endif
            enddo

        enddo

        return

    end subroutine getFluxWeightTopo
    ! ============================================================================
    !> Subroutine soil_production
    !!
    !! This subroutine computes soil production which is assumed to declined exponentially
    !! with increasing soil thickness (Heimsath et al., 2000).
    !!
    !! Regolith production during each time step and each point is computed using Euler’s method,
    !! i.e., during each time step the regolith thickness at the beginning of the time step is used to
    !! compute the increase in regolith thickness during that time step due to regolith production.
    !<
    ! ============================================================================
    subroutine soil_production

        integer :: k, id, nlay, layerid, p, ks

        real( tkind ) :: secant, soil_prod, newsoil
        real( tkind ) :: total_thick, remain_thick, frac

        real( tkind ), dimension( totgrn ) :: sedh, lay_thick

        if( .not. allocated( prop ) ) allocate( prop( lnbPts, totgrn ) )
        prop = 0.0_8

        if( tnow < elapsed_time )then
            elapsed_time = tnow
            return
        endif

        do k = 1, nbPts

            id = lvertexID( k )

            ! Compute the secant
            secant = sqrt( ( slp( k ) )**2.0_8 + 1.0_8 )

            ! Manning's equation and the observed power law width-area relationship
            ! of bedrock channels (Whipple, 2004) enable P0 to be recast in terms of
            ! contributing area A via a power law relationship with an exponent of 3/8:
            ! P0s to be in the range of 0.01–0.1 m/kyr and As to be 1–10 m2.
            !P0 = P0s * ( facc( k ) * strat_dx**2 / As )**(3.0_8/8.0_8)

            ! Here we choose a simplified version where P0 is user defined as well as h0
!            if( facc( k ) > 1.0_8 )then
            ! Braun formulation
            soil_prod = P0 * secant * exp( -soil_thick( id ) / ( h0 * secant ) )
            ! Heimsath et al. formulation
!            soil_prod = P0 * exp( -soil_thick( id ) / h0 )
            ! Furbish and Fagherazzi formulation
!            soil_prod = P0 * ( soil_thick( id ) / h0 ) * exp( -soil_thick( id ) / h0 )
!            else
!                soil_prod = 0.0_8
!            endif

            ! Euler's method
            newsoil = soil_thick( id ) + ( tnow - elapsed_time ) * soil_prod

            ! From soil thickness defined nodes material composition allowed to be
            ! transported
            id = lvertexID( k )
            if( newsoil > 0.0_8 )then
                nlay = v_nbLays( id )
                layerid = v_LaysID( id, nlay )

                ! In case the top layer has enough sediment to diffuse
                if( layerid == layID )then
                    total_thick = 0.0_8
                    do ks = 1, totgrn
                        total_thick = total_thick + stratLays( nlay, id, ks )
                    enddo
                    if( total_thick >= newsoil )then
                        do ks = 1, totgrn
                            prop( id, ks ) = stratLays( nlay, id, ks ) / total_thick
                        enddo
                        goto 21
                    endif
                endif

                lay_thick = 0.0_8
                total_thick = 0.0_8
                remain_thick = newsoil
                soil_loop: do p = v_nbLays( id ), 1, -1

                    if( remain_thick == 0.0_8 ) exit soil_loop

                    ! Get the amount to erode for the considered layer
                    sedh( 1:totgrn ) = 0.0_8
                    do ks = 1, totgrn
                        sedh( ks ) = stratLays( p, id, ks )
                        total_thick = total_thick + stratLays( p, id, ks )
                    enddo

                    ! In case the layer contains enough materials
                    if( remain_thick < total_thick )then
                        frac = remain_thick / total_thick
                        do ks = 1, totgrn
                            sedh( ks ) = stratLays( p, id, ks ) * frac
                            stratLays( p, id, ks ) = stratLays( p, id, ks ) - sedh( ks )
                            lay_thick( ks )  = lay_thick( ks ) + sedh( ks )
                        enddo
                        exit soil_loop

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

                enddo soil_loop

                ! Update the soil layer
                if( nlay < v_nbLays( id ) .or. layerid < layID ) nlay = nlay + 1
                v_nbLays( id ) = nlay
                v_LaysID( id, nlay ) = layID
                total_thick = 0.0_8
                do ks = 1, totgrn
                    stratLays( nlay, id, ks ) = lay_thick( ks )
                    total_thick = total_thick + lay_thick( ks )
                    if( gporo%compaction ) porosityLays( nlay, id, ks ) = porosity( ks, 1 )
                enddo
                do ks = 1, totgrn
                    prop( id, ks ) = stratLays( nlay, id, ks ) / total_thick
                enddo
            else
                newsoil = 0.0_8
                prop( id, 1:totgrn ) = 0.0_8
            endif

21 continue

            ! Here we have the new soil thickness and the required sediments has been
            ! moved to the top layer for possible transport
            soil_thick( id ) = newsoil

        enddo


        return

    end subroutine soil_production
    ! ============================================================================


end module regolith
! ============================================================================
