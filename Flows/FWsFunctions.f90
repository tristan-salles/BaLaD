! ============================================================================
! Name        : FWFunctions.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file FWFunctions.f90
!
! Description :  This module encapsulates functions to compute channel flows evolution
! over the TIN.
!
!<
! ============================================================================
!> Module StreamFunctions 
!<
module stream

    use tin
    use parallel
    use tin_data
    use flow_data
    use mesh_data
    use forces_data
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    real( tkind ) :: glb_slp

contains

    ! ============================================================================
    !> Subroutine calculate_local_slope_elevation
    !! Subroutine calculate_local_slope_elevation determines flow walkers slope vector x/y
    !! components and elevation.
    !<
    ! ============================================================================
    subroutine calculate_local_slope_elevation

        logical :: inside, vicinity_inside

        integer :: face, nface, i

        real( tkind ) :: pt( 2 )

        inside = .false.
        if( rainfws )then
            face = fw_face
            ! Set flow walker point coordinates
            pt( 1 ) = fw_xpos
            pt( 2 ) = fw_ypos
        else
            face = fr_face
            ! Set flow walker point coordinates
            pt( 1 ) = fr_xpos
            pt( 2 ) = fr_ypos
        endif

        ! Check if the flow walker is still in the same face
        if( face > 1 )then
            call check_point_TIN_face( face, pt, inside )

            if( inside )then
                nface = face
            else
                ! If not, find the TIN face containing the flow walker in the vicinity of the previous face
                nface = -1

                vicinity_inside = .false.
                vicinitySearch: do i = 1, tinf_Ngb( face )
                    if( tinf_NbfIDs( face, i ) > 0 )then
                        nface = tinf_NbfIDs( face, i )

                        call check_point_TIN_face( nface, pt, vicinity_inside )
                        if( vicinity_inside ) then
                            exit vicinitySearch
                        else
                            nface = -1
                        endif
                    endif
                enddo vicinitySearch

                ! In case the face is not found in the vicinity do a global search using
                ! the TIN kdtree
                if( .not. vicinity_inside )then
                    nface = -1
                   call search_all_TIN_face( pt, nface )
                endif
            endif

            if( nface > 0 .and. nface <= TIN_faces )then
                if( rainfws )then
                    fw_face = nface
                else
                    fr_face = nface
                endif
            else
                print *,iam,'Something went wrong locating TIN face during local slope calculation 1.'
                stop
            endif
        else
            print *,iam,'Something went wrong locating TIN face during local slope calculation 2.'
            stop
        endif

        ! Calculate the elevation and local slope of the flow walker
        call find_elevation_slope_TIN

        return

    end subroutine calculate_local_slope_elevation
    ! ============================================================================
    !> Subroutine find_elevation_slope_TIN
    !! Subroutine find_elevation_slope_TIN finds the z-coordinate of the flow walker
    !! inside a triangle based on the equation of the plane using the 3 points
    !! which defined the triangular face and compute the local slope of the plane.
    !<
    ! ============================================================================
    subroutine find_elevation_slope_TIN

        logical :: coincide

        integer :: face, fid( 3 )

        real( tkind ) :: x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, zpos
        real( tkind ) :: A, B, C, D

        coincide = .false.
        
        if( rainfws )then
            face = fw_face
            ! Set flow walker point coordinates
            x = fw_xpos
            y = fw_ypos
        else
            face = fr_face
            ! Set flow walker point coordinates
            x = fr_xpos
            y = fr_ypos
        endif

        ! Get the points defining the face
        fid = tinf_PtID( face, 1:3 )

        ! Allocate point coordinates
        x1 = tcoordX( fid( 1 ) )
        y1 = tcoordY( fid( 1 ) )
        z1 = tcoordZ( fid( 1 ) )

        x2 = tcoordX( fid( 2 ) )
        y2 = tcoordY( fid( 2 ) )
        z2 = tcoordZ( fid( 2 ) )

        x3 = tcoordX( fid( 3 ) )
        y3 = tcoordY( fid( 3 ) )
        z3 = tcoordZ( fid( 3 ) )

        ! Check if the flow walker coincide with one of the
        ! face node
        if( x == x1 .and. y == y1 )then
            z = z1
            zpos = z
            coincide = .true.
        elseif( x == x2 .and. y == y2 )then
            z = z2
            zpos = z
            coincide = .true.
        elseif(  x == x3 .and. y == y3 )then
            z = z3
            zpos = z
            coincide = .true.
        endif

        ! If the three vertices of the face have the same z-coordinate
        ! then the face is horizontal and there is no slope.
        if( z2 == z1 .and. z2 == z3 )then
            if( rainfws )then
                fw_zpos = z2
                fw_xslp = 0.0_8
                fw_yslp = 0.0_8
            else
                fr_zpos = z2
                fr_xslp = 0.0_8
                fr_yslp = 0.0_8
            endif
            return
        endif

        ! Calculate plane coefficients
        A = y1 * ( z2 - z3 ) + y2 * ( z3 - z1 ) + y3 * ( z1 - z2 )
        B = z1 * ( x2 - x3 ) + z2 * ( x3 - x1 ) + z3 * ( x1 - x2 )
        C = x1 * ( y2 - y3 ) + x2 * ( y3 - y1 ) + x3 * ( y1 - y2 )
        D = x1 * ( y2 * z3 - y3 * z2 ) + x2 * ( y3 * z1 - y1 * z3 ) + x3 * ( y1 * z2 - y2 * z1 )

        ! Plane equation: A*x + B*y + C*z = D
        z = D - A*x - B*y
        if( C /= 0.0_8 )then
            if( .not. coincide )then
                if( rainfws )then
                    fw_zpos = dble( z / C )
                else
                    fr_zpos = dble( z / C )
                endif
            else
                if( rainfws )then
                    fw_zpos = zpos
                else
                    fr_zpos = zpos
                endif
            endif
            ! Compute the slope in each direction
            if( rainfws )then
                fw_xslp = -A/C
                fw_yslp = -B/C
             else
                fr_xslp = -A/C
                fr_yslp = -B/C
             endif
        else
            print*,'Problem computing FW elevation and local slope'
            print*,'flow element position',x,y
            print*,'TIN face pt1',fid( 1 ),tcoordX( fid( 1 )),tcoordY( fid( 1 )),tcoordZ( fid( 1 ))
            print*,'TIN face pt2',fid( 2 ),tcoordX( fid( 2 )),tcoordY( fid( 2 )),tcoordZ( fid( 2 ))
            print*,'TIN face pt3',fid( 3 ),tcoordX( fid( 3 )),tcoordY( fid( 3 )),tcoordZ( fid( 3 ))
            stop
        endif

        ! If the FW enters the sea as a plume consider that there is no slope
        if( fw_zpos < gsea%actual_sea .and. rainfws ) then
            fw_xslp =  0.0_8
            fw_yslp =  0.0_8
        endif

        if( .not. rainfws .and. fr_zpos < gsea%actual_sea .and. fr_density < sea_density ) then
            fr_xslp =  0.0_8
            fr_yslp =  0.0_8
        endif

        return

    end subroutine find_elevation_slope_TIN
    ! ============================================================================
    !> Subroutine slope_SF_direction
    !! Subroutine slope_SF_direction calculates the slope experienced by the FW based on
    !! on single flow direction.
    !<
    ! ============================================================================
    subroutine slope_SF_direction

        integer :: face, i, k
        integer :: ptDir( 3 ), kpt, kpt2

        real( tkind ) :: dist, maxtravel
        real( tkind ) :: deltx, delty, deltz, delts2

        glb_slp = minimum_slp
        if( rainfws )then
            face = fw_face
        else
            face = fr_face
        endif

        ! Look at the direction of travel and determines which of the 3 surrounding face points
        ! is in the direction of movement
        call travel_distance_vector( ptDir )

        kpt = -1
        kpt2 = -1
        maxtravel = strat_dx
        i = ptDir( 1 )

        if( i > 0 )then
            kpt = i
            ! Compute maximal travel distance
            do k = 1, tinv_Ngb( kpt )
                kpt2 = tinv_NbfIDs( kpt, k )
                dist = compute_distance( tcoordX( kpt ), tcoordY( kpt ), tcoordX( kpt2 ), tcoordY( kpt2 ) )
                maxtravel = min( maxtravel, dist )
            enddo
            maxtravel = max( strat_dx, maxtravel )
        else
            print*,iam,'Something went wrong when finding point in the direction of travel.'
            stop
        endif

        ! If the SF direction is on the edge or on a flat area use local slope as best estimation
        if( SFdirection( kpt ) <= 0 ) return

        ! Otherwise use the single flow direction to estimate slope experienced by the flow walker
        kpt2 = SFdirection( kpt )

        ! If we are pointing towards a vertex from the face containing the flow walker
        ! get the flow direction from this point instead
        if( SFdirection( kpt2 ) > 0 )then
            if( kpt2 == ptDir( 1 ) .or. kpt2 == ptDir( 2 ) .or. kpt2 == ptDir( 3 ) ) &
                kpt2 = SFdirection( kpt2 )
        ! If the SF direction is on the edge or on a flat area use first point direction
        else
            kpt2 = kpt
        endif

        ! On a plain use local slope as best estimation
        if( rainfws .and. fw_zpos == tcoordZ( kpt ) .and. &
            tcoordZ( kpt2 ) == fw_zpos ) return
        if( .not. rainfws .and. fr_zpos == tcoordZ( kpt ) .and. &
            tcoordZ( kpt2 ) == fr_zpos ) return

        ! At the interface with ocean: use first pointing direction
        if( ( rainfws .and. fw_zpos >=  gsea%actual_sea .and. &
            tcoordZ( kpt ) < gsea%actual_sea ) .or. ( .not. rainfws .and. &
            fr_zpos >=  gsea%actual_sea .and. tcoordZ( kpt ) < gsea%actual_sea ) )then
            if( rainfws )then
                deltx = fw_xpos - tcoordX( kpt )
                delty = fw_ypos - tcoordY( kpt )
            else
                deltx = fr_xpos - tcoordX( kpt )
                delty = fr_ypos - tcoordY( kpt )
            endif
            delts2 = deltx**2.0_8 + delty**2.0_8
            if( delts2 /= 0.0_8 )then
                if( rainfws )then
                    deltz = fw_zpos + fw_h - gsea%actual_sea
                    fw_xslp = deltx * deltz / delts2
                    fw_yslp = delty * deltz / delts2
                else
                    deltz = fr_zpos + fr_h - gsea%actual_sea
                    fr_xslp = deltx * deltz / delts2
                    fr_yslp = delty * deltz / delts2
                    glb_slp= max( glb_slp, deltz / sqrt( delts2 ) )
                endif
            endif

        ! Near the ocean interface on the second direction
        elseif( ( rainfws .and. fw_zpos >=  gsea%actual_sea .and. &
            tcoordZ( kpt2 ) < gsea%actual_sea ) .or. ( .not. rainfws .and. &
            fr_zpos >=  gsea%actual_sea .and. tcoordZ( kpt2 ) < gsea%actual_sea ) )then
            if( rainfws )then
                deltx = fw_xpos - tcoordX( kpt2 )
                delty = fw_ypos - tcoordY( kpt2 )
                delts2 = deltx**2.0_8 + delty**2.0_8
                if( delts2 /= 0.0_8 )then
                    deltz = fw_zpos + fw_h - gsea%actual_sea
                    fw_xslp = deltx * deltz / delts2
                    fw_yslp = delty * deltz / delts2
                endif
            else
                deltx = fr_xpos - tcoordX( kpt2 )
                delty = fr_ypos - tcoordY( kpt2 )
                delts2 = deltx**2.0_8 + delty**2.0_8
                if( delts2 /= 0.0_8 )then
                    deltz = fr_zpos + fr_h - gsea%actual_sea
                    fr_xslp = deltx * deltz / delts2
                    fr_yslp = delty * deltz / delts2
                	glb_slp= max( glb_slp, deltz / sqrt( delts2 ) )
                endif
             endif
        else
            if( rainfws )then
                deltx = fw_xpos - tcoordX( kpt2 )
                delty = fw_ypos - tcoordY( kpt2 )
                delts2 = deltx**2 + delty**2
                if( delts2 /= 0.0_8 )then
                    deltz = fw_zpos + fw_h - tcoordZ( kpt2 )
                    fw_xslp = deltx * deltz / delts2
                    fw_yslp = delty * deltz / delts2
                endif
            else
                deltx = fr_xpos - tcoordX( kpt2 )
                delty = fr_ypos - tcoordY( kpt2 )
                delts2 = deltx**2 + delty**2
                if( delts2 /= 0.0_8 )then
                    deltz = fr_zpos + fr_h - tcoordZ( kpt2 )
                    fr_xslp = deltx * deltz / delts2
                    fr_yslp = delty * deltz / delts2
                	glb_slp= max( glb_slp, deltz / sqrt( delts2 ) )
                endif
            endif

        endif

        return

    end subroutine slope_SF_direction
    ! ============================================================================
    !> Subroutine travel_distance_vector
    !! Subroutine travel_distance_vector determines which of the 3 surrounding face points
    !! is in the direction of movement.
    !<
    ! ============================================================================
    subroutine travel_distance_vector( ptDir )

        logical :: intersect, in

        integer :: k, p1, p2, p3, face
        integer, dimension( 3 ) :: ptDir, fid

        real( tkind ) :: side, lside, d1, d2, d0
        real( tkind ) :: ratio1, ratio2, rate, fv, amplifier

        real( tkind ), dimension( 2, 4 ) :: vertex
        real( tkind ), dimension( 2 ) :: travel, pt1, pt2, pt3, ptintersect

        type(kdtree2), pointer :: tree
        type(kdtree2_result), dimension( 1 ) :: results

        ptDir = -1
        if( rainfws )then
            face = fw_face
        else
            face = fr_face
        endif

        ! Get the points defining the face
        fid = tinf_PtID( face, 1:3 )

        ! Allocate point coordinates
        vertex( 1, 1 ) = tcoordX( fid( 1 ) )
        vertex( 2, 1 ) = tcoordY( fid( 1 ) )

        vertex( 1, 2 ) = tcoordX( fid( 2 ) )
        vertex( 2, 2 ) = tcoordY( fid( 2 ) )

        vertex( 1, 3 ) = tcoordX( fid( 3 ) )
        vertex( 2, 3 ) = tcoordY( fid( 3 ) )

        if( rainfws )then
            ! Add the flow walker coordinates
            vertex( 1, 4 ) = fw_xpos
            vertex( 2, 4 ) = fw_ypos
            ! Flow walker velocity
            fv = sqrt( fw_xv**2 + fw_yv**2 )
        else
            ! Add the flow walker coordinates
            vertex( 1, 4 ) = fr_xpos
            vertex( 2, 4 ) = fr_ypos
            ! Flow walker velocity
            fv = sqrt( fr_xv**2 + fr_yv**2 )
        endif

        ! Initialise intersection point
        ptintersect = -9999999999.99_8
        rate = 1.0_8


        if( fv > 0.0_8 )then

            ! Get the maximum side lenght
            lside = 0.0_8
            side = sqrt( ( vertex( 1, 1 ) - vertex( 1, 2 ) )**2 + ( vertex( 2, 1 ) - vertex( 2, 2 ) )**2 )
            lside = max( lside, side )
            side = sqrt( ( vertex( 1, 3 ) - vertex( 1, 1 ) )**2 + ( vertex( 2, 3 ) - vertex( 2, 1 ) )**2 )
            lside = max( lside, side )
            side = sqrt( ( vertex( 1, 2 ) - vertex( 1, 3 ) )**2 + ( vertex( 2, 2 ) - vertex( 2, 3 ) )**2 )
            lside = max( lside, side )

            ! Use the maximum side lenght to enable intersection with TIN face edges
            ! in the direction of travel
            lside = lside * 1000.0_8
            amplifier = dble( 100.0_8 * lside / fv )

20          continue

            ! Based on flow walker velocity and position defined
            ! point in the direction of travel
            if( rainfws )then
                travel( 1 ) = vertex( 1, 4 ) + amplifier * fw_xv * rate
                travel( 2 ) = vertex( 2, 4 ) + amplifier * fw_yv * rate
            else
                travel( 1 ) = vertex( 1, 4 ) + amplifier * fr_xv * rate
                travel( 2 ) = vertex( 2, 4 ) + amplifier * fr_yv * rate
            endif
            pt3( 1 ) = vertex( 1, 4 )
            pt3( 2 ) = vertex( 2, 4 )
            in = .false.

            intersectionSegments: do k = 1, 3
                ! Allocate points
                p1= k
                pt1( 1 ) = vertex( 1, k ) * rate
                pt1( 2 ) = vertex( 2, k ) * rate
                if( k == 3 )then
                    p2 = k - 2
                    pt2( 1 ) = vertex( 1, k - 2 ) * rate
                    pt2( 2 ) = vertex( 2, k - 2 ) * rate
                else
                    p2 = k + 1
                    pt2( 1 ) = vertex( 1, k + 1 ) * rate
                    pt2( 2 ) = vertex( 2, k + 1 ) * rate
                endif
                if( k == 1 ) p3 = 3
                if( k == 2 ) p3 = 1
                if( k == 3 ) p3 = 2

                ! Check if segments intersect
                call find_segments_intersection( pt1, pt2, pt3, travel, intersect, ptintersect )

                ! If it intersects find the closest point to the intersection one
                if( intersect )then
                    in = .true.
                    d0 = sqrt( ( pt2( 1 ) - pt1( 1 ) )**2 + ( pt2( 2 ) - pt1( 2 ) )**2 )
                    d1 = sqrt( ( ptintersect( 1 ) - pt1( 1 ) )**2 + ( ptintersect( 2 ) - pt1( 2 ) )**2 )
                    d2 = sqrt( ( ptintersect( 1 ) - pt2( 1 ) )**2 + ( ptintersect( 2 ) - pt2( 2 ) )**2 )
                    ratio1 = d1 / d0
                    ratio2 = d2 / d0
                    if( ratio1 <= ratio2 )then
                        ptDir( 1 ) = fid( p1 )
                        ptDir( 2 ) = fid( p2 )
                        ptDir( 3 ) = fid( p3 )
                    else
                        ptDir( 1 ) = fid( p2 )
                        ptDir( 2 ) = fid( p1 )
                        ptDir( 3 ) = fid( p3 )
                    endif
                    exit intersectionSegments
                endif
            enddo intersectionSegments

            if( .not. in )then
                if( rate == 1.0_8 )then
                    rate = 100000.0_8
                    goto 20
                elseif( rate == 100000.0_8 )then
                    goto 21
                endif
            else
                goto 22
            endif
        endif

21      continue

        ! Cannot find the point in direction of travel then find the closest point
        tree => kdtree2_create( vertex, sort = .true., rearrange = .true. )
        call kdtree2_n_nearest_around_point(tree, idxin=4, nn=1, correltime=1, results=results)
        call kdtree2_destroy(tree)
        if( results( 1 )%idx == 1 )then
            ptDir( 1 ) = fid( 1 )
            ptDir( 2 ) = fid( 2 )
            ptDir( 3 ) = fid( 3 )
        elseif( results( 1 )%idx == 2 )then
            ptDir( 1 ) = fid( 2 )
            ptDir( 2 ) = fid( 1 )
            ptDir( 3 ) = fid( 3 )
        else
            ptDir( 1 ) = fid( 3 )
            ptDir( 2 ) = fid( 2 )
            ptDir( 3 ) = fid( 1 )
        endif

22      continue

        return

   end subroutine travel_distance_vector
   ! ============================================================================
   !> subroutine find_segments_intersection()
   !!  return false if no solution exists.
   !! \param pt1, pt2, pt3, pt4, solution, intersect
   !<
   ! ============================================================================
   subroutine find_segments_intersection( pt1, pt2, pt3, pt4, solution, intersect )

       ! Variables declaration
       logical :: solution

       integer :: sol, headingpt

       real( tkind ) :: num_a, num_b , denom, mua, mub
       real( tkind ), dimension( 2 ) :: pt1, pt2, pt3, pt4
       real( tkind ), dimension( 2 ) :: p13, p43, p21, intersect

       solution = .false.
       headingpt = 0
       intersect( : ) = -99999999.99_8
       p21( 1 ) = pt2( 1 ) - pt1( 1 )
       p21( 2 ) = pt2( 2 ) - pt1( 2 )
       p43( 1 ) = pt4( 1 ) - pt3( 1 )
       p43( 2 ) = pt4( 2 ) - pt3( 2 )
       p13( 1 ) = pt1( 1 ) - pt3( 1 )
       p13( 2 ) = pt1( 2 ) - pt3( 2 )
       denom = p43( 2 ) * p21( 1 ) - p43( 1 ) * p21( 2 )
       num_a = p43( 1 ) * p13( 2 ) - p43( 2 ) * p13( 1 )
       num_b = p21( 1 ) * p13( 2 ) - p21( 2 ) * p13( 1 )
       if ( denom == 0.0_8 )then
           ! Coincident lines
           if( num_a == 0.0_8 .and. num_b == 0.0_8 )then
               sol = 1
           else
               ! Parallel lines
               sol = 2
           endif
       else
           mua = num_a / denom
           mub = num_b / denom
           ! Intersecting
           if( mua >= 0.0_8 .and. mua <= 1.0_8 .and. mub >= 0.0_8 .and. mub <= 1.0_8 )then
               intersect( 1 ) = pt1( 1 ) + mua * ( pt2( 1 ) - pt1( 1 ) )
               intersect( 2 ) = pt1( 2 ) + mua * ( pt2( 2 ) - pt1( 2 ) )
               sol = 3
              ! Not intersecting
           else
               sol = 4
           endif
       endif
       if( sol == 3 ) solution = .true.

       return

   end subroutine find_segments_intersection
   ! ============================================================================
   !> Subroutine record_flow_walkers
   !! RecordFW records flow walkers position and velocity.
   !<
   ! ============================================================================
   subroutine record_flow_walkers( fwID )

      integer :: k, fwID

      ! Record flow walkers for output
      if( rec_nb < max_rec )then
        rec_nb = rec_nb + 1
        if( rainfws )then
            rfw_id( rec_nb ) = fwID
            rfw_x( rec_nb ) = fw_xpos
            rfw_y( rec_nb ) = fw_ypos
            rfw_z( rec_nb ) = fw_zpos
            rfw_v( rec_nb ) = sqrt( fw_xv**2 + fw_yv**2 )
            rfw_h( rec_nb ) = fw_h
        else
            rfw_id( rec_nb ) = -fwID
            rfw_x( rec_nb ) = fr_xpos
            rfw_y( rec_nb ) = fr_ypos
            rfw_z( rec_nb ) = fr_zpos
            rfw_v( rec_nb ) = sqrt( fr_xv**2 + fr_yv**2 )
            rfw_h( rec_nb ) = fr_h
        endif
     else
        rec_nb = max_rec
        do k = 1, max_rec - 1
             rfw_id( k ) = rfw_id( k + 1 )
             rfw_x( k ) = rfw_x( k + 1 )
             rfw_y( k ) = rfw_y( k + 1 )
             rfw_z( k ) = rfw_z( k + 1 )
             rfw_v( k ) = rfw_v( k + 1 )
             rfw_h( k ) = rfw_h( k + 1 )
        enddo
        if( rainfws )then
            rfw_id( rec_nb ) = fwID
            rfw_x( rec_nb ) = fw_xpos
            rfw_y( rec_nb ) = fw_ypos
            rfw_z( rec_nb ) = fw_zpos
            rfw_v( rec_nb ) = sqrt( fw_xv**2 + fw_yv**2 )
            rfw_h( rec_nb ) = fw_h
        else
            rfw_id( rec_nb ) = -fwID
            rfw_x( rec_nb ) = fr_xpos
            rfw_y( rec_nb ) = fr_ypos
            rfw_z( rec_nb ) = fr_zpos
            rfw_v( rec_nb ) = sqrt( fr_xv**2 + fr_yv**2 )
            rfw_h( rec_nb ) = fr_h
        endif
      endif


      return

   end subroutine record_flow_walkers
   ! ============================================================================

end module stream
! ============================================================================
