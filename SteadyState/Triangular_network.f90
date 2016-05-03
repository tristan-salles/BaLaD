! ============================================================================
! Name        : Triangular_network.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Triangular_network.f90
!
! Description :  This module is used to create the triangular surface used for flow paths
!                      computation.
!
!<
! ============================================================================
!> Module Triangular_network 
!<
module tin

    use morpho
    use parallel
    use interpol
    use tin_data
    use file_data
    use flow_data
    use time_data
    use mesh_data

    implicit none

    integer, dimension( : ), allocatable :: pick_refine

    real( tkind ), dimension( : ), allocatable :: tinfill

contains

	! ============================================================================
    !> Subroutine generate_TIN_surface
    !! Subroutine generate_TIN_surface generates the TIN surface.
    !<
    ! ============================================================================
    subroutine generate_TIN_surface

        if( .not. allocated( tcoordX ) )then
            allocate( tcoordX( lnbPts ) )
            allocate( tcoordY( lnbPts ) )
            allocate( tcoordZ( lnbPts ) )
        endif

        ! Define the TIN grid
        call create_TIN_grid

        ! Call Triangle to built the TIN based on cloud points
        if( iam == 0 ) call trianglegen( 0, ftincpp )
        call mpi_barrier( badlands_comm_world,ierr )

        ! Allocate the TIN grid
        call build_TIN_grid_allocation

        return

    end subroutine generate_TIN_surface
    ! ============================================================================
    !> Subroutine update_TIN_surface
    !! Subroutine update_TIN_surface updates the TIN surface to match with top surface.
    !<
    ! ============================================================================
    subroutine update_TIN_surface

        integer :: id, n, face

        real( tkind ) :: xi( 2 )

        id = 0
        do n = 1, lnbPts
            if( pick_refine( n ) >= 1 )then
                id = id + 1
                if( pick_refine( n ) == 1 )then
                    tcoordZ( id ) = lcoordZ( n )
                else
                    xi( 1 ) = tcoordX( id )
                    xi( 2 ) = tcoordY( id )
                    ! Get global face id containing the point
                    call check_Meshface_containing_point( xi( 1 ), xi( 2 ), face )
                    call define_elevation_interpolation( xi, face, tcoordZ( id ) )
                endif
            endif
        enddo

        ! Apply TIN depression-less algorithm
        call TIN_fill_algorithm

        ! Build single flow direction for flow
        call find_TIN_SF_direction

        return

    end subroutine update_TIN_surface
    ! ============================================================================
    !> Subroutine create_TIN_grid
    !! Subroutine create_TIN_grid create the TIN.
    !<
    ! ============================================================================
    subroutine create_TIN_grid

        ! Parameters Declaration
        integer :: iunit, ios, n, nbpt, id, iseed, face, pid

        real( tkind ) :: ran, ran1
        real( tkind ) :: xhigh, yhigh, x, y

        allocate( pick_refine( lnbPts ) )
        allocate( tinprocNb( nproc ) )
        allocate( tinprocID( lnbPts ) )

        pick_refine = 0
        tinprocNb = 0
        tinprocID = -1

        ! Get the vertex from the stratal mesh
        nbpt = 0
        xhigh = strat_xo + strat_dx
        yhigh = strat_yo + strat_dx
        pick_refine = 0

        do n = 1, lnbPts

            ! Border & low resolution nodes
            ! South side
            if( lcoordY( n ) <= strat_yo )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
               ! North side
            elseif( lcoordY( n ) >= strat_yo + strat_dx * ( strat_Y- 1 ) )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
               ! West side
            elseif( lcoordX( n ) <= strat_xo )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
               ! East side
            elseif( lcoordX( n ) >= strat_xo + strat_dx * ( strat_X- 1 ) )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
            elseif( lcoordX( n ) == xhigh .and. lcoordY( n ) == yhigh )then
                pick_refine( n ) = 3
                nbpt = nbpt + 1
                if( xhigh + strat_dx < strat_xo + strat_dx * ( strat_X - 1 ) )then
                    xhigh = xhigh + strat_dx
                    if( strat_xo + strat_dx * ( strat_X - 1 ) - xhigh < strat_dx )then
                        xhigh = strat_xo + strat_dx
                        yhigh = yhigh + strat_dx
                        if( strat_yo + strat_dx * ( strat_Y - 1 ) - yhigh < strat_dx )&
                                yhigh = strat_yo + strat_dx * lstrat_Y
                    endif
                else
                    xhigh = strat_xo + strat_dx
                    yhigh = yhigh + strat_dx
                    if( strat_yo + strat_dx * ( strat_Y - 1 ) - yhigh < strat_dx )&
                        yhigh = strat_yo + strat_dx * strat_Y
                endif
            endif
        enddo

        ! Create the TIN refined file
        TIN_nbPts = nbpt
        if( iam == 0 )then
            iunit = 17
            open(iunit,file=ftin,status="replace",action="write",iostat=ios)
            rewind(iunit)
            write(iunit, '(I10,1X,3(I2,1X))') TIN_nbPts, 2, 1, 0

            id = 0
            do n = 1, lnbPts
                if( pick_refine( n ) >= 1 )then
                    id = id + 1
                    if( pick_refine( n ) == 1 )then
                        x = lcoordX( n )
                        y = lcoordY( n )
                    else
                        call random_seed( size= iseed )
                        call random_number( harvest= ran )
                        x = lcoordX( n ) + strat_dx * 0.2_8 * ( ran * 2.0_8 - 1.0_8 )
                        call random_seed( size= iseed )
                        call random_number( harvest= ran1 )
                        y = lcoordY( n ) + strat_dx * 0.2_8 * ( ran1 * 2.0_8 - 1.0_8 )
                    endif
                    write(iunit,'(I10,1X,2(F16.3,1X),F4.1,1X,I2)')id,x,y,0.0,0
                    tcoordX( id ) = x
                    tcoordY( id ) = y
                endif
            enddo

            if( id /= TIN_nbPts )&
                print*,'Something went wrong: values mismatch in building TIN grid'

            write(iunit,*)' '
            close( iunit )

        endif

        ! Broadcast value to processors
        call mpi_bcast( tcoordX, lnbPts, dbl_type, 0, badlands_comm_world, ierr )
        call mpi_bcast( tcoordY, lnbPts, dbl_type, 0, badlands_comm_world, ierr )

        id = 0
        pid = 0
        do n = 1, lnbPts
            if( pick_refine( n ) >= 1 )then
                id = id + 1
                if( pick_refine( n ) == 1 )then
                    pid = pid + 1
                    tinprocNb( iam + 1 ) = tinprocNb( iam + 1 ) + 1
                    tinprocID( pid ) = id
                else
                    ! Distribute vertex based on DEM partitioning
                    call check_DEMface_containing_point( tcoordX( id ), tcoordY( id ), face )
                    if( fPid( face ) == iam )then
                        pid = pid + 1
                        tinprocNb( iam + 1 ) = tinprocNb( iam + 1 ) + 1
                        tinprocID( pid ) = id
                    endif
                endif
            endif
        enddo

        return

    end subroutine create_TIN_grid
    ! ============================================================================
    !> Subroutine build_TIN_grid_allocation
    !! Subroutine build_TIN_grid_allocation allocates TIN parameters: nodes coordinates,
    !! faces and grid.
    !<
    ! ============================================================================
    subroutine build_TIN_grid_allocation

        ! Parameters Declaration
        logical :: found, inside

        integer :: iunit, ios, n, nb, kn, n1, n2, fid, k, kp, nb1, km, p1
        integer :: ks, kw, ke, dimension, skip, skip1, skip2

        integer, dimension( 3 ) :: nids, nids2, id
        integer, dimension( mxnghb ) :: nid

        type(kdtree2_result), dimension( mxnghb ) :: TFRslt

        character(len=128) :: tElem

        ! Name of Triangle output files
        iunit = 444
        tElem = 'TIN.1.ele'
        call addpath2( tElem )

        ! Read elements
        inquire(FILE=tElem, EXIST=found)
        if(found)then
            open(iunit,file=tElem,status="old",action="read",iostat=ios)
            rewind(iunit)
        else
            write(6,*)'TIN grid cannot be found.',iam
            call mpi_finalize( ierr )
            stop
        endif

        read(iunit, *) TIN_faces, dimension, skip

        ! Read the TIN file face points IDs
        if( allocated( tinf_PtID ) ) deallocate( tinf_PtID )
        if( allocated( tinf_centroid ) ) deallocate( tinf_centroid )
        if( allocated( tinf_Ngb ) ) deallocate( tinf_Ngb )
        if( allocated( tinf_NbfIDs ) ) deallocate( tinf_NbfIDs )

        if( allocated( TFdata ) )then
            deallocate( TFdata )
            call kdtree2_destroy( TFtree )
        endif
        allocate( TFdata( 2, TIN_faces ) )

        allocate( tinf_PtID( TIN_faces, 3 ) )
        allocate( tinf_centroid( TIN_faces, 2 ) )
        allocate( tinf_Ngb( TIN_faces ) )
        allocate( tinf_NbfIDs( TIN_faces, mxnghb ) )

        do n = 1, TIN_faces
            read(iunit, *,end=123) nb, skip, skip1, skip2
            tinf_PtID( n, 1 ) = skip
            tinf_PtID( n, 2 ) = skip1
            tinf_PtID( n, 3 ) = skip2
            tinf_centroid( n, 1:2 ) = cmp_centroid( n )
            TFdata( 1, n ) =  tinf_centroid( n, 1 )
            TFdata( 2, n ) =  tinf_centroid( n, 2 )
        enddo
123     continue

        ! Close element file
        close( iunit )

        ! Create the kd-tree
        TFtree => kdtree2_create( TFdata, sort = .true., rearrange = .true. )

        ! For each TIN faces find the neighboring faces IDs
        nids = -1
        tinf_Ngb = 0
        tinf_NbfIDs = -1
        do n = 1,  TIN_faces
            nids = tinf_PtID( n, 1:3 )
            nb = 0
            call kdtree2_n_nearest_around_point(TFtree, idxin=n, nn=mxnghb, correltime=1, results=TFRslt)
            do kn = 1, mxnghb
                fid = TFRslt( kn )%idx
                nids2 = tinf_PtID( fid, 1:3 )
                inside = .false.

                same_points_face: do n1 = 1, 3
                    do n2 = 1, 3
                        if( nids2( n1 ) == nids( n2 ) )then
                            inside = .true.
                            exit same_points_face
                        endif
                    enddo
                enddo same_points_face

                if( inside )then
                    nb = nb + 1
                    if( nb > mxnghb )&
                        print*,'Something went worng the number of TIN face neighbors is greater than the maximum allowed.'
                    nid( nb ) = fid
                    if( nid( nb ) < 1 .or. nid( nb ) > TIN_faces )then
                        print*,'Something went worng when looking for face neighbors'
                        stop
                    endif
                endif
            enddo
            tinf_Ngb( n ) = nb
            tinf_NbfIDs( n, 1:nb ) = nid( 1:nb )
        enddo

        ! Built TIN node neighbors
        if( allocated( tinv_Ngb ) ) deallocate( tinv_Ngb )
        if( allocated( tinv_NbfIDs ) ) deallocate( tinv_NbfIDs )
        allocate( tinv_Ngb( TIN_nbPts ) )
        allocate( tinv_NbfIDs( TIN_nbPts, mxnghb ) )
        tinv_Ngb = 0
        tinv_NbfIDs = -1
        do n = 1, TIN_faces
            nb1 = tinf_Ngb( n )
            nids(1:3) = tinf_PtID( n, 1:3)
            do k = 1, 3
                if( tinv_Ngb( nids( k ) ) /= 0 ) goto 21
                if( k == 1 )then
                    tinv_NbfIDs( nids( k ), 1 ) = nids( 2 )
                    tinv_NbfIDs( nids( k ), 2 ) = nids( 3 )
                elseif( k == 2 )then
                    tinv_NbfIDs( nids( k ), 1 ) = nids( 1 )
                    tinv_NbfIDs( nids( k ), 2 ) = nids( 3 )
                else
                    tinv_NbfIDs( nids( k ), 1 ) = nids( 1 )
                    tinv_NbfIDs( nids( k ), 2 ) = nids( 2 )
                endif
                kp = 2

                do p1 = 1, nb1
                    fid = tinf_NbfIDs( n, p1 )
                    found = .false.
                    if( fid < 1 .or. fid > TIN_faces )then
                        print*,'Something went wrong when looking for TIN faces.'
                        stop
                    endif
                    id = tinf_PtID( fid, 1:3 )
                    do km = 1, 3
                        if( id( km ) == nids( k ) ) found = .true.
                    enddo

                    if( found )then
                        do km = 1, 3
                            if( id( km ) /= nids( k ) ) kp = kp + 1
                            if( id( km ) /= nids( k ) ) tinv_NbfIDs( nids( k ), kp ) = id( km )
                        enddo
                        ke = kp
                        do km = 1, kp - 1
                            do ks = km + 1, kp
                                if( tinv_NbfIDs( nids( k ), ks ) == tinv_NbfIDs( nids( k ), km ) )then
                                    tinv_NbfIDs( nids( k ), ks ) = -1
                                    ke = ke - 1
                                endif
                            enddo
                        enddo
                        do km = 1, kp - 1
                            if( tinv_NbfIDs( nids( k ), km ) == -1 )then
                                do kw = km + 1, kp
                                    tinv_NbfIDs( nids( k ), kw - 1 ) = tinv_NbfIDs( nids( k ), kw )
                                enddo
                            endif
                        enddo
                        if( ke < kp )then
                            do km = ke + 1, kp
                                tinv_NbfIDs( nids( k ), km ) = -1
                            enddo
                        endif
                        kp = ke
                    endif
                enddo
                tinv_Ngb( nids( k ) ) = kp
21              continue
            enddo
        enddo

        ! Apply TIN depression-less algorithm
        call TIN_fill_algorithm

        ! Build single flow direction for flow
        call find_TIN_SF_direction

        return

    end subroutine build_TIN_grid_allocation
    ! ============================================================================
    !> Subroutine TIN_fill_algorithm
    !! No sink algorithm.
    !<
    ! ============================================================================
    subroutine TIN_fill_algorithm

        logical :: flag

        integer :: p, k, k2

        flag = .true.

        if( allocated( tinfill ) ) deallocate( tinfill )
        allocate( tinfill( TIN_nbPts ) )

        ! In case we need to fill holes within the DEM
        if( step_fill > 0.0_8 )then

            ! Update DEM borders elevation values
            do k = 1, TIN_nbPts
                ! On East West border fix DEM elevation
                if( tcoordX( k ) == lstrat_xo .or. tcoordX( k ) == lstrat_xm )then
                    tinfill( k ) = tcoordZ( k )
                ! On South North border
                elseif( tcoordY( k ) == lstrat_yo .or. tcoordY( k ) == lstrat_ym )then
                    tinfill( k ) = tcoordZ( k )
                else
                    tinfill( k ) = 1.e6_8
                endif

            enddo

            ! Now find the sinks and fill them usng Planchon's method
            do while( flag )
                flag = .false.
                do k = 1, TIN_nbPts
                    ! In case DEM elevation is greater than top stratal grid elevation
                    if( tinfill( k ) > tcoordZ( k ) )then
                        ! Look at the neighbors
                        do p = 1, tinv_Ngb( k )
                            k2 = tinv_NbfIDs( k, p )
                            ! In case stratal grid elevation greater than neighbors DEM elevation
                            if( tcoordZ( k ) >= tinfill( k2 ) + step_fill )then
                                tinfill( k ) = tcoordZ( k )
                                flag = .true.
                            ! Otherwise it is a sink and we perform sink filling
                            else
                                if( tinfill( k ) > tinfill( k2 ) + step_fill )then
                                    tinfill( k ) = tinfill( k2 ) + step_fill
                                    flag = .true.
                                endif
                            endif
                        enddo
                    endif
                enddo
            enddo
        endif

        ! Change the elevation
        tcoordZ = tinfill

        return

    end subroutine TIN_fill_algorithm
    ! ============================================================================

end module tin
! ============================================================================
