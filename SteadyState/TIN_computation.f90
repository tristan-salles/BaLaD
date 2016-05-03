! ============================================================================
! Name        : TIN_computation.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file TIN_computation.f90
!
! Description :  This module performs a range of calculations and transformations on the TIN.
!
!<
! ============================================================================
!> Module tin_data
!<
module tin_data

    use morpho
    use parallel
    use mesh_data
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

contains

    ! ============================================================================
    !> Subroutine define_elevation_interpolation
    !! Subroutine define_elevation_interpolation is used to compute interpolated elevations.
    !<
    ! ============================================================================
    subroutine define_elevation_interpolation( xi, fcS, zint )

        integer :: fcS, id, k, n
        real( tkind ) :: xi( 2 ), dist( 4 ), maxdist, z( 4 ), zint, wgth( 4 ), sumwgth

        ! Compute elevation based on inverse weighted averaged distance
        ! Shepard's method adapted by Franke and Nielson, 1980
        n = 0
        z = 0.0_8
        zint = 0.0_8
        maxdist = 0.0_8
        sumwgth = 0.0_8
        do k = 1, 4
            id = lfptIDs( fcS, k )
            ! Compute distance
            dist( k ) = sqrt( ( xi( 1 ) - lcoordX( id ) )**2.0_8 + ( xi( 2 ) - lcoordY( id ) )**2.0_8 )
            maxdist = max( maxdist, dist( k ) )
            z( k ) = lcoordZ( id )
            if( dist( k ) == 0.0_8 )then
                zint = z( k )
                n = 1
            endif
        enddo
        if( n == 1 ) return

        do k = 1, 4
            wgth( k ) = ( ( maxdist - dist( k ) ) / ( maxdist * dist( k ) ) )**2.0_8
            sumwgth = wgth( k ) + sumwgth
            zint = zint + wgth( k ) * z( k )
        enddo
        zint = zint / sumwgth

        return

    end subroutine define_elevation_interpolation
    ! ============================================================================
    !> Subroutine find_closest_slope_value
    !! Subroutine find_closest_slope_value is used to find the closest slope for a given point.
    !<
    ! ============================================================================
    subroutine find_closest_slope_value( pt, slope )

        integer :: fcS, id, k, pickID

        real( tkind ) :: pt( 2 ), mindist, dist, slope

        ! Find DEM face
        call check_DEMface_containing_point( pt( 1 ), pt( 2 ), fcS )

        ! Get closest point ID
        pickID = -1
        mindist = 2.0_8 * strat_dx
        do k = 1, 4
            id = fptIDs( fcS, k )
            ! Compute distance
            dist = sqrt( ( pt( 1 ) - coordX( id ) )**2.0_8 + ( pt( 2 ) - coordY( id ) )**2.0_8 )
            if( dist < mindist )then
                mindist = dist
                pickID = id
            endif
        enddo

        ! Get slope
        slope = slp( pickID )

        return

    end subroutine find_closest_slope_value
    ! ============================================================================
    !> Subroutine check_DEMface_containing_point
    !! Subroutine check_DEMface_containing_point is used to determine the face ID
    !! on the DEM grid containing containing a specific point.
    !<
    ! ============================================================================
    subroutine check_DEMface_containing_point( x, y, fcS )

        integer :: Xid, Yid, fcS, id
        real( tkind ) :: x, y

        ! Define the propable row and column number
        Xid = int( ( x - strat_xo ) / strat_dx ) + 1
        Yid = int( ( y - strat_yo ) / strat_dx ) + 1

        ! Define the face id accordingly
        if( Yid == 1 .and. Xid < strat_X )then
            fcS = Xid
        elseif( Yid == 1 .and. Xid == strat_X )then
            fcS = Xid - 1
        elseif( Xid == 1 .and. Yid < strat_Y )then
            fcS = ( Yid - 1 )*( strat_X - 1 ) + 1
        elseif( Yid == strat_Y .and. Xid < strat_X )then
            fcS = ( Yid - 2 )*( strat_X - 1 ) + Xid
        elseif( Yid == strat_Y .and. Xid == strat_X )then
            fcS = ( Yid - 1 )*( Xid - 1 )
        elseif( Yid < strat_Y .and. Xid == strat_X )then
            fcS = Yid * ( Xid - 1 )
        else
            fcS = ( strat_X - 1 )* ( Yid - 1 ) + Xid
        endif

        ! Check that picked face is containing the point
        id = fptIDs( fcS, 1 )
        if( x - coordX( id ) > strat_dx .or. y - coordY( id ) > strat_dx ) then
            print*,'Something went wrong when looking for the DEM face containing the TIN node'
            print*,x,y,id,fcS,Xid,Yid,strat_xo,strat_xm,strat_yo,strat_ym
            stop
        endif

        return

    end subroutine check_DEMface_containing_point
    ! ============================================================================
    !> Subroutine check_Meshface_containing_point
    !! Subroutine check_Meshface_containing_point is used to determine the face ID
    !! on the mesh grid containing containing a specific point.
    !<
    ! ============================================================================
    subroutine check_Meshface_containing_point( x, y, fcS )

        integer :: Xid, Yid, fcS, id
        real( tkind ) :: x, y

        ! Define the propable row and column number
        Xid = int( ( x - lstrat_xo ) / strat_dx ) + 1
        Yid = int( ( y - lstrat_yo ) / strat_dx ) + 1

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
        id = lfptIDs( fcS, 1 )
        if( x - lcoordX( id ) > strat_dx .or. y - lcoordY( id ) > strat_dx ) then
            print*,'Something went wrong when looking for the stratal face containing the TIN node'
            print*,x,y,id,fcS,Xid,Yid,lstrat_xo,lstrat_xm,lstrat_yo,lstrat_ym
            stop
        endif

        return

    end subroutine check_Meshface_containing_point
    ! ============================================================================
    !> Subroutine find_TIN_SF_direction
    !! Subroutine find_TIN_SF_direction performs single flow direction over the TIN.
    !<
    ! ============================================================================
    subroutine find_TIN_SF_direction

        integer :: i, j, k, it, i0, i1, ilast
        integer, dimension( TIN_nbPts ) :: SFdir

        real( tkind ) :: maxel, dist, elast, gradmax, grad

        ! Allocate direction path
        if( .not. allocated( SFdirection ) ) allocate( SFdirection( TIN_nbPts ) )

        ! Mark boundaries as termination points
        SFdir = -10
        SFdirection = -10

        do
            ! Set the initial value to very deep
            maxel = -1.0e8_8

            ! Scan through and find the highest remaining point on the surface
            do j = 1, tinprocNb( iam + 1 )
                i = tinprocID( j )
                if( tcoordZ( i ) >= maxel .and. SFdir( i ) == -10 )then
                    it = i
                    maxel = tcoordZ( i )
                endif
            enddo

            ! Maxel should now contain the highest remaining point or
            ! if no points are left, then it should contain the initial value
            if( maxel == -1.0e8_8 ) exit
            ilast = it
            elast = maxel
            i1 = -1

            SF_loop: do

                ! From the last position, locate the path of steepest decent
                gradmax = -1.0e16_8
                grad = -1.0e16_8
                it = ilast

                do k = 1, tinv_Ngb( it )
                    i0 = tinv_NbfIDs( it, k )
                    if( i0 < 1 .or. i0 > TIN_nbPts )print*,'Something went wrong in the TIN single flow algorithm.'
                    grad = 0.0_8
                    dist = compute_distance( tcoordX( ilast ), tcoordY( ilast ), tcoordX( i0 ), tcoordY( i0 ) )
                    if( dist > 0 ) grad = ( elast - tcoordZ( i0 ) ) / dist
                    if( gradmax < grad )then
                        gradmax = grad
                        i1 = i0
                    endif
                enddo

                ! Depression or flat area.
                if( gradmax <= 0.0_8 )then
                    SFdir( ilast ) = -9
                    exit SF_loop
                else
                    SFdir( ilast ) = i1
                    if( i1 < 0 )then
                        SFdir( ilast ) = -9
                        exit SF_loop
                    endif
                    if( SFdir( i1 ) /= -10 ) exit SF_loop
                    ilast = i1
                    elast = tcoordZ( i1 )
                endif

            enddo SF_loop

        enddo

        call mpi_allreduce( SFdir, SFdirection, TIN_nbPts, int_type, &
            max_type, badlands_comm_world, ierr )

        return

    end subroutine find_TIN_SF_direction
    ! ============================================================================
    !> Subroutine search_all_TIN_face
    !! Subroutine search_all_TIN_face determines over the entire number of faces the one
    !! containing a given point.
    !<
    ! ============================================================================
    subroutine search_all_TIN_face( pt, fce )

        logical :: inside

        integer :: k, fID, fce

        real( tkind ), dimension( 2 ) :: pt
        type(kdtree2_result), dimension( mxnghb ) :: TFRslt

        call kdtree2_n_nearest( TFtree, pt, nn=mxnghb, results=TFRslt)
        fce = -1
        loop_TIN_face: do k = 1, mxnghb
            fID = TFRslt( k )%idx
            if( fID > 0 .and. fID <= TIN_faces )then
                ! Check is the point is in the picked face number
                call check_point_TIN_face( fID, pt, inside )
                if( inside )then
                    fce = fID
                    exit loop_TIN_face
                endif
            else
                print*,'Something went wrong when locating TIN face for flow walker.'
                print*,'FW position',iam,pt,strat_xo,strat_xm,strat_yo,strat_ym
                stop
            endif
        enddo loop_TIN_face

        if( fce < 1 .or. fce > TIN_faces )then
            print*,'Something went wrong when updating TIN face for flow walker.'
            stop
        endif

        return

    end subroutine search_all_TIN_face
    ! ============================================================================
    !> Subroutine check_point_TIN_face()
    !!
    !! Subroutine check_point_TIN_face checks if a given TIN face contained a point.
    !<
    ! ============================================================================
    subroutine check_point_TIN_face( fce, pt, inside )

        logical :: inside

        integer :: fce, kk, nid( 3 )
        real( tkind ) :: x( 3 ), y( 3 ), pt( 2 )

        inside = .false.

        nid = tinf_PtID( fce, 1:3 )

        ! Get face points
        x( 1 ) = tcoordX( nid( 1 ) )
        y( 1 ) = tcoordY( nid( 1 ) )
        x( 2 ) = tcoordX( nid( 2 ) )
        y( 2 ) = tcoordY( nid( 2 ) )
        x( 3 ) = tcoordX( nid( 3 ) )
        y( 3 ) = tcoordY( nid( 3 ) )

        call is_point_in_triangle( pt, x, y, kk )

        if( kk >= 0 ) inside = .true.

        return

    end subroutine check_point_TIN_face
    ! ============================================================================
    !> Subroutine is_point_in_triangle
    !! Subroutine is_point_in_triangle given a triangle check whether or not the point is
    !! inside the polygon.
    !<
    ! ============================================================================
    subroutine is_point_in_triangle( pt, xb, yb, kl )

        integer :: kl

        real( tkind ) :: pt( 2 ), xb( 3 ), yb( 3 )
        real( tkind ) :: det0, det1, det2

        kl = -1

        det0 = ( xb( 2 ) - xb( 1 ) )*( pt( 2 ) - yb( 1 ) ) - ( yb( 2 ) - yb( 1 ) )*( pt( 1 ) - xb( 1 ) )
        det1 = ( xb( 3 ) - xb( 2 ) )*( pt( 2 ) - yb( 2 ) ) - ( yb( 3 ) - yb( 2 ) )*( pt( 1 ) - xb( 2 ) )
        det2 = ( xb( 1 ) - xb( 3 ) )*( pt( 2 ) - yb( 3 ) ) - ( yb( 1 ) - yb( 3 ) )*( pt( 1 ) - xb( 3 ) )

        if( det0 >= 0 .and. det1 >= 0 .and. det2 >= 0 )then
            kl = 1
        elseif( det0 <= 0 .and. det1 <= 0 .and. det2 <= 0 )then
            kl = 1
        endif

        return

    end subroutine is_point_in_triangle
    ! ============================================================================

end module tin_data
! ============================================================================
