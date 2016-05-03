! ============================================================================
! Name        : DiffusionTransport.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file DiffusionTransport.f90
!!
!! DiffusionTransport simulates sediment repartition over topography according to slope.
!!
!<
! ============================================================================
module diffusion

    use parallel
    use mesh_data
    use diffusionfct

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine soil_diffusion_transport
    !! Subroutine soil_diffusion_transport mimics creep by geometric redistribution of sediment
    !! that does not meet a slope criterion which is function of grain size.
    !<
    ! ============================================================================
    subroutine soil_diffusion_transport

        dstart = 0.0_8

        ! Update grid and face values
        call build_current_layer_diffusion_grid

        ! Distribute sediment according to slope criteria
        if( parallel_diff )then
            call parallel_redistribute_sediment_by_diffusion
        else
            call redistribute_sediment_by_diffusion
        endif

        ! Redistribute sediment into active cells
        call update_sediment_layer_diffusion

        return

    end subroutine soil_diffusion_transport
    ! ============================================================================
    !> Subroutine redistribute_sediment_by_diffusion
    !!
    !! Subroutine redistribute_sediment_by_diffusion redistributes sediment geometrically using slope criteria.
    !! Slope criteria are functions of grain size and sediment is diffused to downstream nodes.
    !!
    !! Diffusion of sediment is subdivided into diffusion cycles.
    !! In each cycle all excess material is diffused iteratively until everything is deposited.
    !!
    !! Diffusion only takes place to connected nodes.
    !! All material that has been deposited at unconnected nodes is not checked for slope.
    !!
    !! List of variables:
    !! \arg \c idown - size of kernel
    !! \arg \c depo - redistributed sediment
    !! \arg \c difo - sediment to be deposited at current node
    !! \arg \c difp - sediment to be deposited to downstream nodes
    !! \arg \c dstart - fraction of sediment used when starting diffusion
    !! \arg \c cdif - fraction of sediment to be deposited at current node
    !! \arg \c movout - sediment moved out of simulation area
    !! \arg \c spac - accomodation space at downstream nodes
    !! \arg \c icrow - row position change of downstream nodes in kernel
    !! \arg \c iccol - col position change of downstream nodes in kernel
    !!
    !<
    ! ============================================================================
    subroutine redistribute_sediment_by_diffusion

        integer :: k, kd, ks, vp, iter, locid, ndown, vd, vn

        integer, dimension( 8 ) :: down

        real( tkind ) :: difmax, fc
        real( tkind ) :: spacsum, cc

        real( tkind ), dimension( 8 ) :: spac

        ! Initialize deposition arrays
        depo = 0.0_8

        ! Deposit sediment in diff_nb increments
        do kd = 1, diff_nb

            ! Deposit grain sizes from maximum to minimum sediment grain weight
            do ks = 1, totgrn

                ! Initialise sediments
                difp = 0.0_8
                difo = 0.0_8

                ! Update diffusion based on initial sediment class thickness
                do vn = 1, nbPts
                    vp = lvertexID( vn )
                    difo( vp ) = dstart( vp, ks )
                enddo

                ! Move sediment to downstream nodes as long as there is sediment in the buffer.
                difmax = 1.0e6_8
                iter = 0
                sediment_diffusion: do

                    if( difmax <= diff_res .or. iter >= max_it_cyc ) exit sediment_diffusion
                    difmax = 0.0_8
                    iter = iter + 1

                    ! Determine the fraction of buffer that will be deposited at current node
                    cdif = 0.0_8

                    ! Get maximum elevation at current nodes to ensure sediment stability
                    call find_maximum_elevation( ks )

                    ! Calculate individual and cumulative accomodation space for each downstream node
                    do vd = 1, nbPts
                        vp = lvertexID( vd )

                        if( vp > 0 )then
                            if( difo( vp ) > 0.0_8 )then
                                if( cdif( vp ) < 1.0_8 )then
                                    spacsum = 0.0_8
                                    ndown = 0
                                    down(:) = 0
                                    spac(:) = 0.0_8
                                    call get_available_space_remaining( vp, spacsum, spac, down, ndown )
                                    cc = difo( vp ) * ( 1.0_8 - cdif( vp ) )
                                    difmax = difmax + cc

                                    ! If enough downstream space is available weight sediment by
                                    ! available space and diffuse sediment to downstream nodes
                                    if( spacsum >= cc )then

                                        ! Loop over neighboring nodes
                                        do k = 1, 8
                                            locid = lngbID( vp, k )
                                            if( locid > 0 ) &
                                                difp( locid ) = difp( locid ) + cc * spac( k ) / spacsum
                                        enddo

                                    ! If not enough downstream space is available weight sediment
                                    ! by available space and distribute rest equally among
                                    ! connected downstream nodes or if no downstream nodes are
                                    ! connected assume that enough downstream space is available
                                    ! and weight sediment by available space
                                    else
                                        if( ndown > 0 )then
                                            fc = ( cc - spacsum )/dble( ndown )
                                            ! Loop over neighboring nodes
                                            do k = 1, 8
                                                locid = lngbID( vp, k )
                                                if( locid > 0 )then
                                                    if( down( k ) == 1 )then
                                                        difp( locid ) = difp( locid ) + spac( k ) + fc
                                                    else
                                                        difp( locid ) = difp( locid ) + spac( k )
                                                    endif
                                                endif
                                            enddo
                                        ! If there are no downtream nodes, just distribute it all around
                                        ! and hope the next iteration looks after it
                                        else
                                            fc = cc / 8.0_8

                                            ! Loop over neighboring nodes
                                            do k = 1, 8
                                                locid = lngbID( vp, k )
                                                if( locid > 0 ) difp( locid ) = difp( locid ) + fc
                                            enddo
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo

                    ! Store sediment still to be diffused in difo for next iteration
                    difo = difp
                    difp = 0.0_8

                enddo sediment_diffusion
            enddo

        enddo

        return

    end subroutine redistribute_sediment_by_diffusion
    ! ============================================================================
    !> Subroutine parallel_redistribute_sediment_by_diffusion
    !!
    !! Subroutine parallel_redistribute_sediment_by_diffusion redistributes sediment geometrically
    !! using slope criteria as above but in parallel.
    !!
    !<
    ! ============================================================================
    subroutine parallel_redistribute_sediment_by_diffusion

        integer :: k, kd, ks, vp, iter, locid, ndown, vd, vn
        integer :: IDlft, IDrgt
        integer :: rq1, rq2

        integer, dimension( 8 ) :: down
        integer, dimension( mpi_status_size ) :: stt1, stt2

        real( tkind ) :: difmax, fc
        real( tkind ) :: spacsum, cc

        real( tkind ), dimension( 8 ) :: spac
        real( tkind ), dimension( lstrat_X ) :: sdif

        ! Initialize deposition arrays
        depo = 0.0_8

        ! Deposit sediment in diff_nb increments
        do kd = 1, diff_nb

            ! Deposit grain sizes from maximum to minimum sediment grain weight
            do ks = 1, totgrn

                ! Initialise sediments
                difp = 0.0_8
                difo = 0.0_8
                depdif = 0.0_8

                ! Update diffusion based on initial sediment class thickness
                do vn = 1, nbPts
                    vp = lvertexID( vn )
                    difo( vp ) = dstart( vp, ks )
                enddo

                ! Move sediment to downstream nodes as long as there is sediment in the buffer.
                difmax = 1.0e6_8
                iter = 0
                sediment_diffusion: do

                    if( difmax <= diff_res .or. iter >= max_it_cyc ) exit sediment_diffusion
                    difmax = 0.0_8
                    iter = iter + 1

                    ! Determine the fraction of buffer that will be deposited at current node
                    cdif = 0.0_8

                    ! Get maximum elevation at current nodes to ensure sediment stability
                    call parallel_find_maximum_elevation( ks )
                    
                    ! Calculate individual and cumulative accomodation space for each downstream node
                    do vd = 1, local_lnbPts
                        vp = global_lnid( vd )

                        if( vp > 0 )then
                            if( difo( vp ) > 0.0_8 )then
                                if( cdif( vp ) < 1.0_8 )then
                                    spacsum = 0.0_8
                                    ndown = 0
                                    down(:) = 0
                                    spac(:) = 0.0_8
                                    call get_available_space_remaining( vp, spacsum, spac, down, ndown )
                                    cc = difo( vp ) * ( 1.0_8 - cdif( vp ) )
                                    difmax = difmax + cc

                                    ! If enough downstream space is available weight sediment by
                                    ! available space and diffuse sediment to downstream nodes
                                    if( spacsum >= cc )then

                                        ! Loop over neighboring nodes
                                        do k = 1, 8
                                            locid = lngbID( vp, k )
                                            if( locid > 0 ) &
                                                difp( locid ) = difp( locid ) + cc * spac( k ) / spacsum
                                        enddo

                                    ! If not enough downstream space is available weight sediment
                                    ! by available space and distribute rest equally among
                                    ! connected downstream nodes or if no downstream nodes are
                                    ! connected assume that enough downstream space is available
                                    ! and weight sediment by available space
                                    else
                                        if( ndown > 0 )then
                                            fc = ( cc - spacsum )/dble( ndown )
                                            ! Loop over neighboring nodes
                                            do k = 1, 8
                                                locid = lngbID( vp, k )
                                                if( locid > 0 )then
                                                    if( down( k ) == 1 )then
                                                        difp( locid ) = difp( locid ) + spac( k ) + fc
                                                    else
                                                        difp( locid ) = difp( locid ) + spac( k )
                                                    endif
                                                endif
                                            enddo
                                        ! If there are no downtream nodes, just distribute it all around
                                        ! and hope the next iteration looks after it
                                        else
                                            fc = cc / 8.0_8

                                            ! Loop over neighboring nodes
                                            do k = 1, 8
                                                locid = lngbID( vp, k )
                                                if( locid > 0 ) difp( locid ) = difp( locid ) + fc
                                            enddo
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo

                    ! Send information to partition
                    if( iam < gproc - 1 )then

                        ! Send top row to other partition
                        IDlft = ( part_rows( iam + 1 ) - 1 ) * lstrat_X + 1
                        IDrgt = IDlft + lstrat_X - 1

                        call mpi_isend( difp( IDlft:IDrgt), lstrat_X, dbl_type, iam+1, 9, &
                            badlands_comm_world, rq1, ierr )
                        call mpi_request_free( rq1, ierr )

                        ! Reveive top row + 1 from other partition
                        sdif = 0.0_8
                        call mpi_irecv( sdif, lstrat_X, dbl_type, iam+1, 8, &
                            badlands_comm_world,  rq2, ierr )
                        call mpi_wait( rq2, stt2, ierr )

                        IDlft = IDlft + lstrat_X - 1
                        do k = 1, lstrat_X
                            difp( IDlft + k ) = difp( IDlft + k ) + sdif( k )
                        enddo
                    endif

                    if( iam > 0 )then

                        !  Send row 1 to other partition
                        IDlft = part_rows( iam ) * lstrat_X + 1
                        IDrgt = IDlft + lstrat_X - 1

                        call mpi_isend( difp( IDlft:IDrgt), lstrat_X, dbl_type, iam-1, 8, &
                            badlands_comm_world, rq2, ierr )
                        call mpi_request_free( rq2, ierr )

                        ! Receive row 1 - 1 from other partition
                        sdif = 0.0_8
                        call mpi_irecv( sdif, lstrat_X, dbl_type, iam-1, 9, &
                            badlands_comm_world,  rq1, ierr )
                        call mpi_wait( rq1, stt1, ierr )

                        IDlft = IDlft - lstrat_X - 1
                        do k = 1, lstrat_X
                            difp( IDlft + k ) = difp( IDlft + k ) + sdif( k )
                        enddo
                    endif

                    ! Store sediment still to be diffused in difo for next iteration
                    difo = difp
                    difp = 0.0_8

                    ! Get the maximum residual
                    call mpi_allreduce( mpi_in_place, difmax, 1, dbl_type, max_type, badlands_comm_world, ierr )

                enddo sediment_diffusion

                ! Get the deposit thickness for the considered sediment type
                call mpi_allreduce( mpi_in_place, depdif, lnbPts, dbl_type, max_type, badlands_comm_world, ierr )
                depo( 1:lnbPts, ks ) = depdif( 1:lnbPts )

            enddo

        enddo

        return

    end subroutine parallel_redistribute_sediment_by_diffusion
    ! ============================================================================

end module diffusion
! ============================================================================
