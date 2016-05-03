! ============================================================================
! Name        : Facc_algorithm.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Facc_algorithm.f90
!
! Description :  This module performs a multiple-flow direction computation to get the
! flow accumulation.
!
!<
! ============================================================================
!> Module Flow_accumulation 
!<
module flow_algo

    use parallel
    use flow_data
    use mesh_data

    implicit none

    integer, dimension( : ), allocatable :: degreeMatrix
    integer, dimension( :,: ), allocatable :: directMatrix, reverseMatrix

    real( tkind ), dimension( : ), allocatable :: oldfacc, newfacc
    real( tkind ), dimension( : ), allocatable :: oldfacr, newfacr

    real( tkind ), dimension( :,: ), allocatable :: weight

contains

    ! ============================================================================
    !> Subroutine compute_MFD
    !! This subroutine computes the multiple-flow direction accumulation
    !<
    ! ============================================================================
    subroutine compute_MFD

        call getCodeDirectMatrix
        call getReverseMatrix
        call getWeightMatrix
        call getDegreeFA
        !call parallel_getDegreeFA

        return

    end subroutine compute_MFD
    ! ============================================================================
    !> Subroutine getWeightMatrix
    !! This subroutine computes the matrix of flow fraction
    !<
    ! ============================================================================
    subroutine getWeightMatrix

        integer :: k, k2, p

        real( tkind ) :: sum, maxslope, fe
        real( tkind ), dimension( 8 ) :: slopeArray


        if( .not. allocated( weight ) ) allocate( weight( lnbPts, 8 ) )

        do k = 1, lnbPts

            fe = 0.0_8
            sum = 0.0_8
            maxslope = 0.0_8
            slopeArray = 0.0_8

            do p = 1, 8
                weight( k, p ) = 0.0_8
                k2 = lngbID( k, p )
                if( k2 > 0 )then
                    if( lfilldem( k ) > lfilldem( k2 ) )then
                        slopeArray( p ) = ( lfilldem( k ) - lfilldem( k2 ) ) / strat_dx
                        if( mod( p, 2 ) == 0 ) slopeArray( p ) = slopeArray( p ) / sqrt( 2.0_8 )
                    endif
                endif
                maxslope = max( slopeArray( p ), maxslope )
            enddo

            ! Flow partition exponent
            fe = 8.9_8 * min( maxslope, 1.0_8 ) + 1.1_8

            do p = 1, 8
                if( slopeArray( p ) > 0.0_8 )then
                    if( mod( p, 2 ) == 0 )then
                        sum = sum + ( (slopeArray( p ))**fe ) * sqrt( 2.0_8 ) / 4.0_8
                    else
                        sum = sum + ( (slopeArray( p ))**fe ) / 2.0_8
                    endif
                endif
            enddo

            do p = 1, 8
                if( slopeArray( p ) > 0.0_8 )then
                    if( mod( p, 2 ) == 0 )then
                        if( sum > 0.0_8 )then
                            weight( k, p ) =  ( (slopeArray( p ))**fe ) * sqrt( 2.0_8 ) / ( sum * 4.0_8 )
                        endif
                    else
                        if( sum > 0.0_8 )then
                            weight( k, p ) = ( (slopeArray( p ))**fe ) / ( sum * 2.0_8 )
                        endif
                    endif
                endif
            enddo

        enddo

        return

    end subroutine getWeightMatrix
    ! ============================================================================
    !> Subroutine getCodeDirectMatrix
    !! This subroutine determines the matrix of the flow direction
    !<
    ! ============================================================================
    subroutine getCodeDirectMatrix

        integer :: k, k2, p

        real( tkind ), dimension( 8 ) :: slopeArray

        if( .not. allocated( directMatrix ) ) allocate( directMatrix( lnbPts, 8 ) )

        do k = 1, lnbPts

            slopeArray = 0.0_8

            do p = 1, 8
                k2 = lngbID( k, p )
                if( k2 > 0 )then
                    if( lfilldem( k ) > lfilldem( k2 ) )then
                        slopeArray( p ) = lfilldem( k ) - lfilldem( k2 )
                        if( mod( p, 2 ) == 0 ) slopeArray( p ) = slopeArray( p ) / sqrt( 2.0_8 )
                    endif
                endif
            enddo

            directMatrix( k, 1:8 ) = 0
            do p = 1, 8
                if( slopeArray( p ) > 0.0_8 )then
                    directMatrix( k, p ) = 1
                endif
            enddo

        enddo

        return

    end subroutine getCodeDirectMatrix
    ! ============================================================================
    !> Subroutine getReverseMatrix
    !! This subroutine gets the reverse flow direction matrix
    !<
    ! ============================================================================
    subroutine getReverseMatrix

        integer :: k, mask, p

        if( .not. allocated( reverseMatrix ) ) allocate( reverseMatrix( lnbPts, 8 ) )

        reverseMatrix = 0

        do k = 1, lnbPts

            do p = 1, 8
                mask = directMatrix( k, p )
                if( mask == 1 )then
                    if( p == 1 )then
                        reverseMatrix( k + lstrat_X, 5 ) = 1
                    elseif( p == 2 )then
                        reverseMatrix( k + lstrat_X + 1, 6 ) = 1
                    elseif( p == 3 )then
                        reverseMatrix( k + 1, 7 ) = 1
                    elseif( p == 4 )then
                        reverseMatrix( k - lstrat_X + 1, 8 ) = 1
                    elseif( p == 5 )then
                        reverseMatrix( k - lstrat_X, 1 ) = 1
                    elseif( p == 6 )then
                        reverseMatrix( k - lstrat_X - 1, 2 ) = 1
                    elseif( p == 7 )then
                        reverseMatrix( k - 1, 3 ) = 1
                    elseif( p == 8 )then
                        reverseMatrix( k + lstrat_X - 1, 4 ) = 1
                    endif
                endif
            enddo

        enddo

        return

    end subroutine getReverseMatrix
    ! ============================================================================
    !> Subroutine getDegreeMatrix
    !! This subroutine gets the degree flow direction matrix
    !<
    ! ============================================================================
    subroutine getDegreeMatrix

        integer :: k, mask, p

        if( .not. allocated( degreeMatrix ) ) allocate( degreeMatrix( lnbPts ) )

        do k = 1, lnbPts
            degreeMatrix( k ) = 0
            do p = 1, 8
                mask = reverseMatrix( k, p )
                if( mask == 1 ) degreeMatrix( k ) = degreeMatrix( k ) + 1
            enddo
        enddo

        return

    end subroutine getDegreeMatrix
    ! ============================================================================
    !> Subroutine getDegreeFA
    !! This subroutine computes the flow accumulation
    !<
    ! ============================================================================
    subroutine getDegreeFA

        logical :: faflag

        integer :: k, dir, p

        real( tkind ) :: accu, prep

        call getDegreeMatrix

        if( .not. allocated( lfacc ) ) allocate( lfacc( lnbPts ) )
        if( .not. allocated( facc ) ) allocate( facc( nbPts ) )
        if( .not. allocated( lfacr ) ) allocate( lfacr( lnbPts ) )
        if( .not. allocated( facr ) ) allocate( facr( nbPts ) )

        lfacc = 1.0_8
        lfacr = precipitation

        faflag = .true.

        do while( faflag )

            faflag = .false.

            do k = 1, lnbPts

                if( degreeMatrix( k ) == 0 )then

                    degreeMatrix( k ) = -1
                    accu = 0.0_8
                    prep = 0.0_8
                    do p = 1, 8
                        dir = reverseMatrix( k, p )
                        if( dir == 1 )then
                            if( p == 1 )then
                                accu = accu + lfacc( k + lstrat_X )*weight( k + lstrat_X, 5 )
                                prep = prep + lfacr( k + lstrat_X )*weight( k + lstrat_X, 5 )
                            endif
                            if( p == 2 )then
                                accu = accu + lfacc( k + lstrat_X + 1 )*weight( k + lstrat_X + 1, 6 )
                                prep = prep + lfacr( k + lstrat_X + 1 )*weight( k + lstrat_X + 1, 6 )
                            endif
                            if( p == 3 )then
                                accu = accu + lfacc( k + 1 )*weight( k + 1, 7 )
                                prep = prep + lfacr( k + 1 )*weight( k + 1, 7 )
                            endif
                            if( p == 4 )then
                                accu = accu + lfacc( k - lstrat_X + 1 )*weight( k - lstrat_X + 1, 8 )
                                prep = prep + lfacr( k - lstrat_X + 1 )*weight( k - lstrat_X + 1, 8 )
                            endif
                            if( p == 5 )then
                                accu = accu + lfacc( k - lstrat_X )*weight(  k - lstrat_X, 1 )
                                prep = prep + lfacr( k - lstrat_X )*weight(  k - lstrat_X, 1 )
                            endif
                            if( p == 6 )then
                                accu = accu + lfacc( k - lstrat_X - 1 )*weight( k - lstrat_X - 1, 2 )
                                prep = prep + lfacr( k - lstrat_X - 1 )*weight( k - lstrat_X - 1, 2 )
                            endif
                            if( p == 7 )then
                                accu = accu + lfacc( k - 1 )*weight( k - 1, 3 )
                                prep = prep + lfacr( k - 1 )*weight( k - 1, 3 )
                            endif
                            if( p == 8 )then
                                accu = accu + lfacc( k + lstrat_X - 1 )*weight( k + lstrat_X - 1, 4 )
                                prep = prep + lfacr( k + lstrat_X - 1 )*weight( k + lstrat_X - 1, 4 )
                            endif
                        endif
                    enddo
                    lfacc( k ) = lfacc( k ) + accu
                    lfacr( k ) = lfacr( k ) + prep

                    do p = 1, 8
                        dir = directMatrix( k, p )
                        if( dir == 1 )then
                            if( p == 1 )&
                                degreeMatrix( k + lstrat_X ) = degreeMatrix( k + lstrat_X ) - 1
                            if( p == 2 )&
                                degreeMatrix( k + lstrat_X + 1 ) = degreeMatrix( k + lstrat_X + 1 ) - 1
                            if( p == 3 )&
                                degreeMatrix( k + 1 ) = degreeMatrix( k + 1 ) - 1
                            if( p == 4 )&
                                degreeMatrix( k - lstrat_X + 1 ) = degreeMatrix( k - lstrat_X + 1 ) - 1
                            if( p == 5 )&
                                degreeMatrix( k - lstrat_X ) = degreeMatrix( k - lstrat_X ) - 1
                            if( p == 6 )&
                                degreeMatrix( k - lstrat_X - 1 ) = degreeMatrix( k - lstrat_X - 1 ) - 1
                            if( p == 7 )&
                                degreeMatrix( k - 1 ) = degreeMatrix( k - 1 ) - 1
                            if( p == 8 )&
                                degreeMatrix( k + lstrat_X - 1 ) = degreeMatrix( k + lstrat_X - 1 ) - 1
                        endif
                    enddo
                    faflag = .true.

                endif

            enddo

        enddo

        return

    end subroutine getDegreeFA
    ! ============================================================================
    !> Subroutine parallel_getDegreeFA
    !! This subroutine computes the flow accumulation in parallel
    !<
    ! ============================================================================
    subroutine parallel_getDegreeFA

        logical :: faflag

        integer :: k, dir, p, n, intflag

        integer :: IDlft, IDrgt
        integer :: rq1, rq2, rq3, rq4

        integer, dimension( mpi_status_size ) :: stt1, stt2, stt3, stt4

        real( tkind ), dimension( lstrat_X ) :: pfa, pfr
        real( tkind ) :: accu, prep

        if( .not. allocated( lfacc ) ) allocate( lfacc( lnbPts ) )
        if( .not. allocated( facc ) ) allocate( facc( nbPts ) )
        if( .not. allocated( lfacr ) ) allocate( lfacr( lnbPts ) )
        if( .not. allocated( facr ) ) allocate( facr( nbPts ) )

        if( .not. allocated( oldfacc ) ) allocate( oldfacc( lnbPts ) )
        if( .not. allocated( newfacc ) ) allocate( newfacc( lnbPts ) )
        if( .not. allocated( oldfacr ) ) allocate( oldfacr( lnbPts ) )
        if( .not. allocated( newfacr ) ) allocate( newfacr( lnbPts ) )

        lfacc = 1.0_8
        lfacr = precipitation

        oldfacc = 1.0_8
        newfacc = 0.0_8
        oldfacr = precipitation
        newfacr = 0.0_8

        faflag = .true.

        do while( faflag )

            faflag = .false.
            intflag = 0

            ! For each processors
            do n = 1, local_lnbPts

                k = global_lnid( n )

                if( oldfacc( k ) > 0.0_8 )then

                    accu = 0.0_8
                    prep = 0.0_8

                    do p = 1, 8
                        dir = reverseMatrix( k, p )
                        if( dir == 1 )then
                            if( p == 1 )then
                                accu = accu + oldfacc( k + lstrat_X )*weight( k + lstrat_X, 5 )
                                prep = prep + oldfacr( k + lstrat_X )*weight( k + lstrat_X, 5 )
                            endif
                            if( p == 2 )then
                                accu = accu + oldfacc( k + lstrat_X + 1 )*weight( k + lstrat_X + 1, 6 )
                                prep = prep + oldfacr( k + lstrat_X + 1 )*weight( k + lstrat_X + 1, 6 )
                            endif
                            if( p == 3 )then
                                accu = accu + oldfacc( k + 1 )*weight( k + 1, 7 )
                                prep = prep + oldfacr( k + 1 )*weight( k + 1, 7 )
                            endif
                            if( p == 4 )then
                                accu = accu + oldfacc( k - lstrat_X + 1 )*weight( k - lstrat_X + 1, 8 )
                                prep = prep + oldfacr( k - lstrat_X + 1 )*weight( k - lstrat_X + 1, 8 )
                            endif
                            if( p == 5 )then
                                accu = accu + oldfacc( k - lstrat_X )*weight(  k - lstrat_X, 1 )
                                prep = prep + oldfacr( k - lstrat_X )*weight(  k - lstrat_X, 1 )
                            endif
                            if( p == 6 )then
                                accu = accu + oldfacc( k - lstrat_X - 1 )*weight( k - lstrat_X - 1, 2 )
                                prep = prep + oldfacr( k - lstrat_X - 1 )*weight( k - lstrat_X - 1, 2 )
                            endif
                            if( p == 7 )then
                                accu = accu + oldfacc( k - 1 )*weight( k - 1, 3 )
                                prep = prep + oldfacr( k - 1 )*weight( k - 1, 3 )
                            endif
                            if( p == 8 )then
                                accu = accu + oldfacc( k + lstrat_X - 1 )*weight( k + lstrat_X - 1, 4 )
                                prep = prep + oldfacr( k + lstrat_X - 1 )*weight( k + lstrat_X - 1, 4 )
                            endif
                        endif
                    enddo

                    if( accu > 0 )then
                        lfacc( k ) = lfacc( k ) + accu
                        lfacr( k ) = lfacr( k ) + prep
                        newfacc( k ) = accu
                        newfacr( k ) = prep
                        intflag = 1
                    endif

                endif

            enddo

            call mpi_allreduce( mpi_in_place, intflag, 1, int_type, max_type, badlands_comm_world, ierr )
            if( intflag == 1 ) faflag = .true.

            ! Send information to partition
            if( iam < gproc - 1 )then

                ! Send top row to other partition
                IDlft = ( part_rows( iam + 1 ) - 1 ) * lstrat_X + 1
                IDrgt = IDlft + lstrat_X - 1

                call mpi_isend( newfacc( IDlft:IDrgt), lstrat_X, dbl_type, iam+1, 13, &
                    badlands_comm_world, rq1, ierr )
                call mpi_request_free( rq1, ierr )

                call mpi_isend( newfacr( IDlft:IDrgt), lstrat_X, dbl_type, iam+1, 15, &
                    badlands_comm_world, rq3, ierr )
                call mpi_request_free( rq3, ierr )

                ! Reveive top row + 1 from other partition
                pfa = 0.0_8
                call mpi_irecv( pfa, lstrat_X, dbl_type, iam+1, 14, &
                    badlands_comm_world,  rq2, ierr )
                call mpi_wait( rq2, stt2, ierr )

                pfr = 0.0_8
                call mpi_irecv( pfr, lstrat_X, dbl_type, iam+1, 16, &
                    badlands_comm_world,  rq4, ierr )
                call mpi_wait( rq4, stt4, ierr )

                IDlft = IDlft + lstrat_X - 1
                do p = 1, lstrat_X
                    newfacc( IDlft + p ) = pfa( p )
                    newfacr( IDlft + p ) = pfr( p )
                enddo
            endif

            if( iam > 0 )then

                !  Send row 1 to other partition
                IDlft = part_rows( iam ) * lstrat_X + 1
                IDrgt = IDlft + lstrat_X - 1

                call mpi_isend( newfacc( IDlft:IDrgt), lstrat_X, dbl_type, iam-1, 14, &
                    badlands_comm_world, rq2, ierr )
                call mpi_request_free( rq2, ierr )

                call mpi_isend( newfacr( IDlft:IDrgt), lstrat_X, dbl_type, iam-1, 16, &
                    badlands_comm_world, rq4, ierr )
                call mpi_request_free( rq4, ierr )

                ! Receive row 1 - 1 from other partition

                pfa = 0.0_8
                call mpi_irecv( pfa, lstrat_X, dbl_type, iam-1, 13, &
                    badlands_comm_world,  rq1, ierr )
                call mpi_wait( rq1, stt1, ierr )

                pfr = 0.0_8
                call mpi_irecv( pfr, lstrat_X, dbl_type, iam-1, 15, &
                    badlands_comm_world,  rq3, ierr )
                call mpi_wait( rq3, stt3, ierr )

                IDlft = IDlft - lstrat_X - 1
                do p = 1, lstrat_X
                    newfacc( IDlft + p ) = pfa( p )
                    newfacr( IDlft + p ) = pfr( p )
                enddo
            endif

            ! Update flow
            oldfacc = newfacc
            oldfacr = newfacr
            newfacc = 0.0_8
            newfacr = 0.0_8

        enddo

        ! Broadcast flow and rain accumulation
        call mpi_allreduce( mpi_in_place, lfacc, lnbPts, dbl_type, max_type, badlands_comm_world, ierr )
        call mpi_allreduce( mpi_in_place, lfacr, lnbPts, dbl_type, max_type, badlands_comm_world, ierr )

        return

    end subroutine parallel_getDegreeFA
    ! ============================================================================

end module flow_algo
! ============================================================================
