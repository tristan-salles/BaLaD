! ============================================================================
! Name        : IsostaticFlexure.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file IsostaticFlexure.f90
!
! Description :  This module performs a multigrid and high-order solution of isostasy equation.
!
! Reference: 3D modelling of flexural isostatic deformation
! authors: Li, F., Dyt, C., Griffiths, C.M.
! source: Computers & Geosciences v. 30, no.9/10, pp. 1105-1115.
!
!<
! ============================================================================
!> Module IsotaticFlexure 
!<
module isoflex

    use parallel
    use file_data
    use flow_data
    use time_data
    use mesh_data
    use forces_data

    implicit none

    integer :: nrow, ncol

    real( tkind ) :: constant1, constant2, constant3

    ! Fourth order schema parameters
    real( tkind ), parameter :: q1=216.0_8,q2=-18_8,q3=-84.0_8
    real( tkind ), parameter :: q4=-3.0_8,q5=792.0_8

    real( tkind ), parameter :: r1=144.0_8,r2=18.0_8,r3=-48.0_8
    real( tkind ), parameter :: r4=-6.0_8,r5=1.0_8/240.0_8

    real( tkind ), dimension( :,: ), allocatable :: w, w1, wx, wy, flex_load, flex_z
    real( tkind ), dimension( : ), allocatable :: load, prevload

contains

    ! ============================================================================
    !> Subroutine compute_isostatic_flexure
    !! This subroutine performs a multigrid and high-order solution of isostasy equation.
    !<
    ! ============================================================================
    subroutine compute_isostatic_flexure

        integer :: i, j, p, n2r, n2c, n4r, n4c, iter, m, nr1, nc1

        real( tkind ) :: diffload, flexdiff, flextot, flex_tor

        flex_tor = ( 10.0_8 )**( -3.0_8 + 2.0_8 * ( 23.0_8 - log10( flex_rigidity ) ) / 3.0_8 )  &
            * ( strat_dx / 27000.0_8 )**2

        ! Find the loading values
        call initialise_isostatic_flexure( diffload )

        ! If added load is small ignore isostatic changes
        if( abs( diffload ) < constant1 *  sediment( 1 )%density )then
            ! Put back the previous load as it was
            prevload = prevload - load
            return
        endif

        ! Initialise flexure arrays
        if( .not. allocated( wx ) ) allocate( wx( lstrat_Y, lstrat_X ) )
        if( .not. allocated( wy ) ) allocate( wy( lstrat_Y, lstrat_X ) )
        if( .not. allocated( w1 ) ) allocate( w1( lstrat_Y, lstrat_X ) )
        if( .not. allocated( flex_load ) ) allocate( flex_load( lstrat_Y, lstrat_X ) )
        if( .not. allocated( flex_z ) ) allocate( flex_z( lstrat_Y, lstrat_X ) )

        w = 0.0_8
        wx = 0.0_8
        wy = 0.0_8
        w1 = 0.0_8

        nrow = lstrat_Y
        ncol = lstrat_X
        nr1 = nrow - 1
        nc1 = ncol - 1

        p = 0
        do i = 1, nrow
            do j = 1, ncol
                p = p + 1
                flex_load( i, j ) = load( p )
                flex_z( i, j ) = lcoordZ( p )
            enddo
        enddo

        n2c = ( ncol - 2 ) / 2
        n2c = n2c * 2 + 2
        n2r = ( nrow - 2 ) / 2
        n2r = n2r * 2 + 2

        n4c = ( ncol - 2 ) / 4
        n4c = n4c * 4 + 2
        n4r = ( nrow - 2 ) / 4
        n4r = n4r * 4 + 2

        ! Iterate until a solution is reached
        iter = 0
        flexloop: do
            iter = iter + 1

            ! Finest grid first
            m = 1
            call fourth_order_schema( m, nr1, nc1, 4 )

            flexdiff = 0.0_8
            flextot = 0.0_8

            ! Compare with the last finest grid results
            do j = 2, ncol - 1
                do i = 2, nrow - 1
                    flexdiff = flexdiff + abs( w1( i, j ) - w( i, j ) )
                    flextot = flextot + abs( w( i, j ) )
                enddo
            enddo

            w1 = w
            if( flexdiff < flextot * flex_tor .or. iter > 1.e5 ) exit flexloop
            if( flextot == 0.0_8 ) exit flexloop

            ! Coarse grid
            m=2
            do j = 2 + m, n2c - m, m
                do i = 2 + m, n2r - m, m
                    w( i, j ) = 0.0625_8 * ( w( i - 1, j - 1 ) + w( i + 1, j - 1 ) + w( i - 1, j + 1 ) + w( i + 1, j + 1 ) + &
                        2.0_8 * ( w( i, j - 1 ) + w( i, j + 1 ) + w( i - 1, j ) + w( i + 1, j ) ) + 4.0_8 * w( i, j ) )

                    wx( i, j ) = 0.0625_8 * ( wx( i - 1, j - 1 ) + wx( i + 1, j - 1 ) + wx( i - 1, j + 1 ) + wx( i + 1, j + 1 ) + &
                        2.0_8 * ( wx( i, j - 1 ) + wx( i, j + 1 ) + wx( i - 1, j ) + wx( i + 1, j ) ) + 4.0_8 * wx( i, j ) )

                    wy( i, j ) = 0.0625_8 * ( wy( i - 1, j - 1 ) + wy( i + 1, j - 1 ) + wy( i - 1, j + 1 ) + wy( i + 1, j + 1 ) + &
                        2.0_8 * ( wy( i, j - 1 ) + wy( i, j + 1 ) + wy( i - 1, j ) + wy( i + 1, j ) ) + 4.0_8*wy( i, j ) )
                enddo
            enddo
            call fourth_order_schema( m, n2r, n2c, 4 )

            ! Coarsest grid
            m = 4
            do j = 2 + m, n4c - m, m
                do i = 2 + m, n4r - m, m
                    w( i, j ) = 0.0625_8 * ( w( i - 2, j - 2 ) + w( i + 2, j - 2 ) + w( i - 2, j + 2 ) + w( i + 2, j + 2 ) + &
                        2.0_8 * ( w( i, j - 2 ) + w( i, j + 2 ) + w( i - 2, j ) + w( i + 2, j ) ) + 4.0_8 * w( i, j ) )

                    wx( i, j ) = 0.0625_8 * ( wx( i - 2, j - 2 ) + wx( i + 2, j - 2 ) + wx( i - 2, j + 2 ) + wx( i + 2, j + 2 ) + &
                        2.0_8 * ( wx( i, j - 2 ) + wx( i, j + 2 ) + wx( i - 2, j ) + wx( i + 2, j ) ) + 4.0_8 * wx( i, j ) )

                    wy( i, j ) = 0.0625_8 * ( wy( i - 2, j - 2 ) + wy( i + 2, j - 2 ) + wy( i - 2, j + 2 ) + wy( i + 2, j + 2 ) + &
                        2.0_8 * ( wy( i, j - 2 ) + wy( i, j + 2 ) + wy( i -2, j ) + w( i + 2, j ) ) + 4.0_8 * wy( i, j ) )

                enddo
            enddo
            call update_boundary_values( m, n4r - m, n4c - m, w )
            call update_boundary_values( m, n4r - m, n4c - m, wx )
            call update_boundary_values( m, n4r - m, n4c - m, wy )

            call fourth_order_schema( m, n4r, n4c, 16 )

            ! Interpolate to the coarse grid
            do j = 2, n4c - m, 4
                do i = 2, n4r - m, 4
                    w( i, j + 2 ) = 0.5_8 * ( w( i, j ) + w( i, j + m ) )
                    w( i + 2, j ) = 0.5_8 * ( w( i, j ) + w( i + m, j ) )
                    w( i + 2, j + 2 ) = 0.25_8 * ( w( i, j ) + w( i, j + m ) + w( i + m, j ) + w( i + m, j + m ) )

                    wx( i, j + 2 ) = 0.5_8 * ( wx( i, j ) + wx( i, j + m ) )
                    wx( i + 2, j ) = 0.5_8 * ( wx( i, j ) + wx( i + m, j ) )
                    wx( i + 2, j + 2 ) = 0.25_8 * ( wx( i, j ) + wx( i, j + m ) + wx( i + m, j ) + wx( i + m, j + m ) )

                    wy( i, j + 2 ) = 0.5_8 * ( wy( i, j ) + wy( i, j + m ) )
                    wy( i + 2, j ) = 0.5_8 * ( wy( i, j ) + wy( i + m, j ) )
                    wy( i + 2, j + 2 ) = 0.25_8 * ( wy( i, j ) + wy( i, j + m ) + wy( i + m, j ) + wy( i + m, j + m ) )

                enddo
            enddo

            ! At the boundary of the coarse grid
            if( n2c > n4c )then
                do i = 2, n4r - m, 4
                    w( i, n2c ) = w( i, n4c )
                    w( i + 2, n2c ) = 0.5_8 * ( w( i, n4c ) + w( i + m, n4c ) )

                    wx( i, n2c ) = wx( i, n4c )
                    wx( i + 2, n2c ) = 0.5_8 * ( wx( i, n4c ) + wx( i + m, n4c ) )

                    wy( i, n2c ) = wy( i, n4c )
                    wy( i + 2, n2c ) = 0.5_8 * ( wy( i, n4c ) + wy( i + m, n4c ) )
                enddo
            endif

            if( n2r > n4r )then
                do j = 2, n4c - m, 4
                    w( n2r, j ) = w( n4r, j )
                    w( n2r, j+2 ) = 0.5_8 * ( w( n4r, j ) + w( n4r, j + m ) )

                    wx( n2r, j ) = wx( n4r, j )
                    wx( n2r, j+2 ) = 0.5_8 * ( wx( n4r, j ) + wx( n4r, j + m ) )

                    wy( n2r, j) = wy( n4r, j )
                    wy( n2r, j+2 ) =0.5_8 * ( wy( n4r, j ) + wy( n4r, j + m ) )
                enddo
            endif
            w( n2r, n2c ) = w( n4r, n4c )
            wx( n2r, n2c ) = wx( n4r, n4c )
            wy( n2r, n2c ) = wy( n4r, n4c )

            m = 2
            call fourth_order_schema( m, n2r, n2c, 4 )

            ! Interpolate to the fine grid
            do j = 2, n2c - m, 2
                do i = 2, n2r - m, 2
                    w( i, j + 1 ) = 0.5_8 * ( w( i, j ) + w( i, j + m ) )
                    w( i + 1, j ) = 0.5_8 *( w( i, j ) + w( i + m, j ) )
                    w( i + 1, j + 1 ) = 0.25_8 * ( w( i, j ) + w( i, j + m ) + w( i + m, j ) + w( i + m, j + m ) )

                    wx( i, j + 1 ) = 0.5_8 * ( wx( i, j ) + wx( i, j + m ) )
                    wx( i + 1, j ) = 0.5_8 *( wx( i, j ) + wx( i + m, j ) )
                    wx( i + 1, j + 1 ) = 0.25_8 * ( wx( i, j ) + wx( i, j + m ) + wx( i + m, j ) + wx( i + m, j + m ) )

                    wy( i, j + 1 ) = 0.5_8 * ( wy( i, j ) + wy( i, j + m ) )
                    wy( i + 1, j ) = 0.5_8 *( wy( i, j ) + wy( i + m, j ) )
                    wy( i + 1, j + 1 ) = 0.25_8 * ( wy( i, j ) + wy( i, j + m ) + wy( i + m, j ) + wy( i + m, j + m ) )
                enddo
            enddo

            ! At the boundary of the fine grid
            if( nc1 > n2c )then
                do m = 2, n2r - m, m
                    w( i, nc1 ) = w( i, n2c )
                    w( i + 1, nc1 ) = 0.5_8 * ( w( i, n2c ) + w( i + m, n2c ) )

                    wx( i, nc1 ) = wx( i, n2c )
                    wx( i + 1, nc1 ) = 0.5_8 * ( wx( i, n2c ) + wx( i + m, n2c ) )

                    wy( i, nc1 ) = w( i, n2c )
                    wy( i + 1, nc1 ) = 0.5_8 * ( wy( i, n2c ) + wy( i + m, n2c ) )
                enddo
            endif

            if( nr1 > n2r )then
                do j = 2, n2c - m, 2
                    w( nr1, j ) = w( n2r, j )
                    w( nr1, j + 1 ) = 0.5_8 * ( w( n2r, j ) + w( n2r, j + m ) )

                    wx( nr1, j ) = wx( n2r, j )
                    wx( nr1, j + 1 ) = 0.5_8 * ( wx( n2r, j ) + wx( n2r, j + m ) )

                    wy( nr1, j ) = wy( n2r, j )
                    wy( nr1, j + 1 ) = 0.5_8 * ( wy( n2r, j ) + wy( n2r, j + m ) )
                enddo
            endif

            w( nr1, nc1 ) = w( n2r, n2c )
            wx( nr1, nc1 ) = wx( n2r, n2c )
            wy( nr1, nc1 ) = wy( n2r, n2c )

            do i = 1, nrow
                w( i, 1 ) = w( i, 2 )
                w( i, ncol ) = w( i, ncol - 1 )

                wx( i, 1 ) = wx( i, 2 )
                wx( i, ncol ) = wx( i, ncol - 1 )

                wy( i, 1 ) = wy( i, 2 )
                wy( i, ncol ) = wy( i, ncol - 1 )
            enddo

            do j = 1, ncol
                w( 1, j ) = w( 2, j )
                w( nrow, j ) = w( nr1, j )

                wx( 1, j ) = wx( 2, j )
                wx( nrow, j ) = wx( nr1, j )

                wy( 1, j ) = wy( 2, j )
                wy( nrow, j )=wy( nr1, j )
            enddo


        enddo flexloop

        if( iter > 1.e5 .and. iam == 0 )print*,'Flexural isostasy did not converge.',&
            flexdiff,1.0e-18_8 * flextot

        ! Apply the flexure to the deposit layers
        p = 0
        do i = 1, nrow
            do j = 1, ncol
                p = p + 1
                lcoordZ( p ) = lcoordZ( p ) - w( i, j )
                lbase( p ) = lbase( p ) - w( i, j )
                flexure( p ) = -w( i, j )
            enddo
        enddo

        ! Broadcast top elevation to DEM grid
        do i = 1, nbPts
            coordZ( i ) = lcoordZ( lvertexID( i ) )
        enddo

        return

    end subroutine compute_isostatic_flexure
    ! ============================================================================
    !> Subroutine initialise_isostatic_flexure
    !! This subroutine initialise isostatic flexure conditions by computing the load
    !! from basement.
    !<
    ! ============================================================================
    subroutine initialise_isostatic_flexure( difftot )

        integer :: k, p, id, left, right, ks

        real( tkind ) :: temp, temp1, temp2, difftot, oldload

        if( .not. allocated( load ) ) allocate( load( lnbPts ) )
        load = 0.0_8

        if( .not. allocated( prevload ) )then
            allocate( prevload( lnbPts ) )
            prevload = 0.0_8
        endif

        if( .not. allocated( w ) ) allocate( w( lstrat_Y, lstrat_X ) )

        ! Declare parameters
        constant1 = gravity / flex_rigidity
        constant2 = constant1 * ( crust_density - sea_density )
        constant3 = constant1 * crust_density

        left = 1
        right = strat_X
        difftot = 0.0_8
        do k = 1, nbPts

            id = lvertexID( k )
            temp = 0.0_8
            do p = 1, v_nbLays( id )
                temp1 = 0.0_8
                temp2 = 0.0_8
                do ks = 1, totgrn
                    temp1 =  temp1 + stratLays( p, id, ks ) * sediment( ks )%density
                    ! Add load due to water in pore space
                    if( gporo%compaction )then
                        temp2 = temp2 + stratLays( p, id, ks ) * porosityLays( p, id, ks )
                    endif
                enddo
                temp = temp + temp1 + temp2 * sea_density
            enddo

            ! Add water load if required
            if( lcoordZ( id ) < gsea%actual_sea )then
                temp = temp + ( gsea%actual_sea - lcoordZ( id ) ) * sea_density
            endif


            ! Define the borders load
            if( k == left )then
                load( id - 1 ) = load( id )
                left = left + strat_X
            elseif( k == right )then
                load( id ) = load( id + 1 )
                right = right + strat_X
            endif

            ! Define the loading value
            oldload = prevload( id )
            prevload( id ) = constant1 * temp
            load( id ) = prevload( id ) - oldload
            difftot = difftot + load( id )
        enddo

        ! Define the South / North border loads
        prevload( 1 : lstrat_X ) = prevload( lstrat_X + 1 : 2 * lstrat_X )
        prevload( lnbPts - lstrat_X + 1 : lnbPts ) = prevload( lnbPts - 2 * lstrat_X + 1 : lnbPts - lstrat_X )
        load( 1 : lstrat_X ) = load( lstrat_X + 1 : 2 * lstrat_X )
        load( lnbPts - lstrat_X + 1 : lnbPts ) = load( lnbPts - 2 * lstrat_X + 1 : lnbPts - lstrat_X )

        return

    end subroutine initialise_isostatic_flexure
    ! ============================================================================
    !> Subroutine update_boundary_values
    !! This subroutine updates the boundary values.
    !<
    ! ============================================================================
    subroutine update_boundary_values( ks, nrw, ncl, temp )

        integer :: i, j, nrw, ncl, k1, ks
        real( tkind ), dimension( nrow, ncol ) :: temp

        ! We know the values from 2+val up to n2x or n2y
        ! find the values for the borders
        k1 = 2 + ks
        do i = 1, nrow
            do j = 1, k1 - 1
                temp( i, j ) = temp( i, k1 )
            enddo

            do j = ncl + 1, ncol
                temp( i, j ) = temp( i, ncl )
            enddo
        enddo

        do j = 1, ncol
            do i = 1, k1 - 1
                temp( i, j ) = temp( k1, j )
            enddo
            do i = nrw + 1, nrow
                temp( i, j ) = temp( nrw, j )
            enddo
        enddo

        return

    end subroutine update_boundary_values
    ! ============================================================================
    !> Subroutine fourth_order_schema
    !! This subroutine performs a fourth-order compact (nine grid points) finite difference
    !! approximation to solve the isotatic response.
    !<
    ! ============================================================================
    subroutine fourth_order_schema( m, nrw, ncl, step )

        integer :: i, j, nrw, ncl, step, m, iter, pass, isw, jsw, ks

        real( tkind ) :: dxm4, elev, flexdiff, flextot, piv, resw, w3
        real( tkind ), dimension( nrow, ncol ) :: ld

        ks=0
        if( m > 1) ks = m
        dxm4 = ( strat_dx * m )**4.0_8

        do j= 2, ncl, m
            do i = 2, nrw, m
                elev = flex_z( i, j )
                if( elev < gsea%actual_sea )then
                    ld( i, j ) = ( flex_load( i, j ) - constant2 * w( i, j ) ) * dxm4
                else
                    ld( i, j ) = ( flex_load( i, j ) - constant3 * w( i, j ) ) * dxm4
                endif
            enddo
        enddo
        call update_boundary_values( 0, nrw, ncl, ld )

        iter = 0
        schema_loop: do
            iter = iter + 1
            flexdiff = 0.0_8
            flextot = 0.0_8

            jsw = 1
            do pass = 1, 2
                isw = jsw
                do j = 2 + ks, ncl - ks, m
                    do i = isw + 1 + ks, nrw - ks, 2 * m
                        elev = flex_z( i, j )
                        if( elev < gsea%actual_sea )then
                            piv = q5 + 11.0_8 * constant2 * dxm4
                        else
                            piv = q5 + 11.0_8 * constant3 * dxm4
                        endif

                        wx( i, j )=(  r1 * ( w( i + m, j ) - w( i - m, j ) ) &
                            + r2 * ( w( i + m, j + m ) - w( i - m, j + m ) + w( i + m, j - m ) - w( i - m, j - m ) )  &
                            + r3 * ( wx( i + m, j ) + wx( i - m, j ) ) &
                            + r4 * ( wx( i + m, j + m ) + wx( i - m, j + m ) + wx( i + m, j - m ) + wx( i - m, j - m ) ) &
                            + r4 * ( wy( i + m, j + m ) - wy( i - m, j + m ) - wy( i + m, j - m ) + wy( i - m, j - m ) ) &
                            + ( ld( i + m, j ) - ld( i - m, j ) ) ) * r5

                        wy( i, j ) = ( r1 * ( w( i, j + m ) - w( i, j - m ) ) &
                            + r2 * ( w( i + m, j + m ) + w( i - m, j + m ) -w( i + m, j - m ) - w( i - m, j - m ) ) &
                            + r3 * ( wy( i, j + m ) + wy( i, j - m ) ) &
                            + r4 * ( wy( i + m, j + m ) + wy( i - m, j + m ) + wy( i + m, j - m ) + wy( i - m, j - m ) ) &
                            + r4 * ( wx( i + m, j + m ) -wx( i - m, j + m ) - wx( i + m, j - m ) + wx( i - m, j - m ) ) &
                            + ( ld( i, j + m ) - ld( i, j - m ) ) )*r5

                        w3 = ( q1 * ( w( i + m, j ) + w( i - m, j ) + w( i, j + m ) + w( i, j - m ) ) &
                            + q2 * ( w( i + m, j + m ) + w( i - m, j + m ) + w( i + m, j - m ) + w( i - m, j - m ) ) &
                            + q3 * ( wx( i + m, j ) - wx( i - m, j ) + wy( i, j + m ) - wy( i, j - m ) ) &
                            + q4 * ( wx( i + m, j + m ) - wx( i - m, j + m ) + wx( i + m, j - m ) - wx( i - m, j - m ) ) &
                            + q4 * ( wy( i + m, j + m ) + wy( i - m, j + m ) - wy( i + m, j - m ) - wy( i - m, j - m ) ) &
                            + ( ld( i + m, j ) + ld( i - m, j ) + ld( i, j + m ) + ld( i, j - m ) ) &
                            + 11.0_8 * flex_load( i, j ) * dxm4 ) / piv

                        if( elev < gsea%actual_sea )then
                            ld( i, j ) = ( flex_load( i, j ) - constant2 * w3 ) * dxm4
                        else
                            ld( i, j ) = ( flex_load( i, j ) - constant3 * w3 ) * dxm4
                        endif

                        resw = w3 - w( i, j )
                        flexdiff = flexdiff + abs( w3 - w( i, j ) )
                        flextot = flextot + abs( w3 )
                        w( i, j ) = w3

                    enddo
                    isw = m + 2 - isw

                enddo
                jsw = m + 2 - jsw

                call update_boundary_values( ks, nrw - ks, ncl - ks, w )
                call update_boundary_values( ks, nrw - ks, ncl - ks, wx )
                call update_boundary_values( ks, nrw - ks, ncl - ks, wy )
                call update_boundary_values( ks, nrw - ks, ncl - ks, ld )

            enddo

            if( flexdiff < flextot * 1.e-10_8 .or. iter > step ) exit schema_loop
            if( flextot == 0.0_8 ) exit schema_loop

        enddo schema_loop

        return

    end subroutine fourth_order_schema
    ! ============================================================================

end module isoflex
! ============================================================================
