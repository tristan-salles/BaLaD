! ============================================================================
! Name        : Morphometrics.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Morphometrics.f90
!
! Description :   Set of functions to perform morphometrics analysis of the
!!  topography.
!
!<
! ============================================================================
!> Module morpho
!<
module morpho

    use parallel
    use interpol
    use file_data
    use time_data
    use flow_data
    use mesh_data

    implicit none

    real( tkind ), dimension( 8 ) :: angleP = (/ 0.0_8, 45.0_8, 90.0_8, 135.0_8, 180.0_8, &
        225.0_8, 270.0_8, 315.0_8 /)

    integer, dimension( 8 ) :: direct8x = (/ 1, 1, 1, 0, -1, -1, -1, 0 /)
    integer, dimension( 8 ) :: direct8y = (/ -1, 0, 1, 1, 1, 0, -1, -1 /)

    real( tkind ), dimension( 9 ) :: circ9x = (/  -1.0_8/sqrt( 2.0_8 ), 0.0_8, 1.0_8/sqrt( 2.0_8 ), &
        -1.0_8, 0.0_8, 1.0_8, -1.0_8/sqrt( 2.0_8 ), 0.0_8, 1.0_8/sqrt( 2.0_8 ) /)
    real( tkind ), dimension( 9 ) :: circ9y = (/  1.0_8/sqrt( 2.0_8 ), 1.0_8, 1.0_8/sqrt( 2.0_8 ), &
        0.0_8, 0.0_8, 0.0_8, -1.0_8/sqrt( 2.0_8 ), -1.0_8, -1.0_8/sqrt( 2.0_8 ) /)

contains

    ! ============================================================================
    !> Subroutine planchon_dem_fill_algorithm
    !! No sink algorithm.
    !<
    ! ============================================================================
    subroutine planchon_dem_fill_algorithm

        logical :: flag
        integer :: p, k, k2, i, j

        flag = .true.

        if( .not. allocated( lfilldem ) ) allocate( lfilldem( lnbPts ) )
        if( .not. allocated( filldem ) ) allocate( filldem( nbPts ) )

        ! In case we need to fill holes within the DEM
        if( step_fill > 0.0_8 )then

            ! Update DEM borders elevation values
            p = 1
            do k = 1, lnbPts
                ! On East West border fix DEM elevation
                if( p == 1 .or. p == lstrat_X )then
                    lfilldem( k ) = lcoordZ( k )
                    if( p == 1 ) p = p + 1
                    if( p == lstrat_X ) p = 1
                ! On South border
                elseif( k <= lstrat_X )then
                    lfilldem( k ) = lcoordZ( k )
                    p = p + 1
                ! On the North border
                elseif( k > lnbPts - lstrat_X )then
                    lfilldem( k ) = lcoordZ( k )
                    p = p + 1
                ! Else give high values
                else
                    lfilldem( k ) = 1.e6_8
                    p = p + 1
                endif

            enddo

            ! Now find the sinks and fill them usng Planchon's method
            do while( flag )
                flag = .false.

                do j = 2, lstrat_Y - 1
                    do i = 2, lstrat_X - 1
                        k = ( j - 1 ) * lstrat_X + i
                        ! In case DEM elevation is greater than top stratal grid elevation
                        if( lfilldem( k ) > lcoordZ( k ) )then
                            ! Look at the neighbors
                            do p = 1, 8

                                k2 = lngbID( k, p )
                                ! In case stratal grid elevation greater than neighbors DEM elevation
                                if( lcoordZ( k ) >= lfilldem( k2 ) + step_fill )then
                                    lfilldem( k ) = lcoordZ( k )
                                    flag = .true.
                                ! Otherwise it is a sink and we perform sink filling
                                else
                                    if( lfilldem( k ) > lfilldem( k2 ) + step_fill )then
                                        lfilldem( k ) = lfilldem( k2 ) + step_fill
                                        flag = .true.
                                    endif
                                endif

                            enddo

                        endif

                    enddo
                enddo
            enddo

        ! In case the sink filling algorithm is not used just update the DEM elevation to
        ! the top stratal elevation
        else
            lfilldem = lcoordZ
        endif

        return

    end subroutine planchon_dem_fill_algorithm
    ! ============================================================================
    !> Subroutine circular_neighborhood_slp_aspect
    !!
    !! This subroutine interpolate relief value using quadratic polynomial interpolation.
    !<
    ! ============================================================================
    subroutine circular_neighborhood_slp_aspect( z, ort )

        real( tkind )  :: aspect, clls, z( 9 ), ort

        real( tkind ) :: p, q

        ! Quadratic polynomial coefficient
        p = ( sqrt( 2.0_8 ) * z(3) + 2.0_8 * z(6) + sqrt( 2.0_8 ) * z(9) - &
            sqrt( 2.0_8 ) * z(1) - 2.0_8 * z(4) - sqrt( 2.0_8 ) * z(7) ) &
            / ( 8.0_8 * strat_dx )
        q = ( sqrt( 2.0_8 ) * z(1) + 2.0_8 * z(2) + sqrt( 2.0_8 ) * z(3) - &
            sqrt( 2.0_8 ) * z(7) - 2.0_8 * z(8) - sqrt( 2.0_8 ) * z(9) ) &
            / ( 8.0_8 * strat_dx )

        ! Get the aspect
        if( q == 0.0_8 .and. p == 0.0_8 )then
            aspect = 0.0_8
        else
            aspect = atan2( q, p )
        endif
        aspect = 180.0_8 * aspect / pi + 180.0_8
        if( aspect < 0.0_8 )then
            clls = 90.0_8 - aspect
        elseif( aspect > 90.0_8 )then
            clls = 360.0 - aspect + 90.0_8
        else
            clls = 90.0_8 - aspect
        endif
        if( clls > 360.0_8 ) clls = clls - 360.0_8
        ort = pi * clls / 180.0_8

        return

    end subroutine circular_neighborhood_slp_aspect
    ! ============================================================================
    !> Subroutine circular_neighborhood_slope
    !! Compute DEM aspect based on circular neighborhood statistic.
    !<
    ! ============================================================================
    subroutine circular_neighborhood_aspect( pt, aspect )

        integer :: p

        real( tkind )  :: pt( 2 ), xi( 2 ), z( 9 ), zval, aspect

        do p = 1, 9
            xi( 1 ) = pt( 1 ) + circ9x( p ) * strat_dx
            xi( 2 ) = pt( 2 ) + circ9y( p ) * strat_dx
            call get_lagrangian_interpolated_relief_value( xi, zval )
            z( p ) = zval
        enddo

        ! Get morphometrics
        call circular_neighborhood_slp_aspect( z, aspect )

        return

    end subroutine circular_neighborhood_aspect
    ! ============================================================================
    !> Subroutine circular_neighborhood_statistic
    !! Compute DEM statistics based on circular neighborhood statistic.
    !<
    ! ============================================================================
    subroutine circular_neighborhood_statistic

        integer :: k, p, n

        real( tkind )  :: xi( 2 ), z( 9 ), zval

        ! Check that arrays for curvature and slope have been initialised
        if( .not. allocated( slp ) ) allocate( slp( nbPts ) )
        if( .not. allocated( vcurv ) ) allocate( vcurv( nbPts ) )
        if( .not. allocated( hcurv ) ) allocate( hcurv( nbPts ) )
        if( .not. allocated( orient ) ) allocate( orient( nbPts ) )

        ! Initialise
        slp = 0.0_8
        hcurv = -1.e6_8
        vcurv = -1.e6_8
        orient = 0.0_8

        do n = 1, local_nbPts
            k = global_nid( n )
!        do k = 1, nbPts
            ! Check the point is not on the border
            do p = 1, 9
                xi( 1 ) = lcoordX( lvertexID( k ) ) + circ9x( p ) * strat_dx
                xi( 2 ) = lcoordY( lvertexID( k ) ) + circ9y( p ) * strat_dx
                call get_lagrangian_interpolated_relief_value( xi, zval )
                z( p ) = zval
           enddo

           ! Get morphometrics
           call quadratic_circular_neighborhood( k, z )

        enddo

        return

    end subroutine circular_neighborhood_statistic
    ! ============================================================================
    !> Subroutine quadratic_circular_neighborhood
    !!
    !! This subroutine interpolate relief value using quadratic polynomial interpolation.
    !<
    ! ============================================================================
    subroutine quadratic_circular_neighborhood( k, z )

        integer :: k

        real( tkind )  :: aspect, clls, z( 9 )

        real( tkind ) :: r, t, s, p, q, u

        ! Quadratic polynomial coefficient
        r = ( z(1) + z(3) + 3.0_8 * z(4) + 3.0_8 * z(6) + z(7) + z(9) - z(2) - 8.0_8 * z(5) - z(8) ) &
            / ( 4.0_8 * strat_dx**2.0_8 )
        t = ( z(1) + z(3) + 3.0_8 * z(2) + 3.0_8 * z(8) + z(7) + z(9) - z(4) - 8.0_8 * z(5) - z(6) ) &
            / ( 4.0_8 * strat_dx**2.0_8 )
        s =  ( z(3) + z(7) - z(1) - z(9) ) / ( 2.0_8 * strat_dx**2.0_8 )
        p = ( sqrt( 2.0_8 ) * z(3) + 2.0_8 * z(6) + sqrt( 2.0_8 ) * z(9) - &
            sqrt( 2.0_8 ) * z(1) - 2.0_8 * z(4) - sqrt( 2.0_8 ) * z(7) ) &
            / ( 8.0_8 * strat_dx )
        q = ( sqrt( 2.0_8 ) * z(1) + 2.0_8 * z(2) + sqrt( 2.0_8 ) * z(3) - &
            sqrt( 2.0_8 ) * z(7) - 2.0_8 * z(8) - sqrt( 2.0_8 ) * z(9) ) &
            / ( 8.0_8 * strat_dx )
        u = z(5)

        ! Get the gradient
        slp( k ) = atan( sqrt( p*p + q*q ) )

        ! Get the aspect
        if( q == 0.0_8 .and. p == 0.0_8 )then
            aspect = 0.0_8
        else
            aspect = atan2( q, p )
        endif

        aspect = 180.0_8 * aspect / pi + 180.0_8
        if( aspect < 0.0_8 )then
            clls = 90.0_8 - aspect
        elseif( aspect > 90.0_8 )then
            clls = 360.0 - aspect + 90.0_8
        else
            clls = 90.0_8 - aspect
        endif
        if( clls > 360.0_8 ) clls = clls - 360.0_8
        orient( k ) = pi * clls / 180.0_8

        if( q == 0.0_8 .and. p == 0.0_8 )then
            hcurv( k ) = 0.0_8
            vcurv( k ) = 0.0_8
        else
            ! Get the horizontal curvature
            hcurv( k ) = - ( q*q*r - 2.0_8*p*q*s + p*p*t ) / ( ( p*p + q*q ) * sqrt( 1 + p*p + q*q ) )
            ! Get the vertical curvature
            vcurv( k ) = - ( p*p*r + 2.0_8*p*q*s + q*q*t ) / ( ( p*p + q*q ) * sqrt( ( 1 + p*p + q*q )**3.0_8 ) )
        endif

        return

    end subroutine quadratic_circular_neighborhood
    ! ============================================================================

end module morpho
! ============================================================================
