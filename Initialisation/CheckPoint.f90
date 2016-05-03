! ============================================================================
! Name        : CheckPoint.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file CheckPoint.f90
!
! Description : CheckPoint is used to write the information usefull for continuing a SPModel
! experiment.
!
!<
! ============================================================================
module checkpoint

    use hdf5
    use isoflex
    use parallel
    use file_data
    use FoX_wxml
    use time_data
    use mesh_data
    use forces_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine checkpt_layers
    !! Subroutine checkpt_layers gather sediment layers information to the processors.
    !<
    ! ============================================================================
    subroutine checkpt_layers

        integer :: i, k

        if( .not. allocated( check_NbLay ) )then
            allocate( check_NbLay( lnbPts ) )
            allocate( check_base( lnbPts ) )
            allocate( soil_thickness( lnbPts ) )
            allocate( check_layID( lnbPts, nbLays ) )
            if( gporo%compaction ) &
                allocate( check_porosity( lnbPts, nbLays, totgrn ) )
            allocate( check_sed( lnbPts, nbLays, totgrn ) )
        endif

        check_NbLay = -1
        do i = 1, lnbPts
            check_NbLay( i ) = v_nbLays( i )
            check_base( i ) = lbase( i )
            soil_thickness( i ) = soil_thick( i )
            do k = 1, v_nbLays( i )
                check_sed( i , k , 1:totgrn ) = 0.0_8
                check_layID( i , k ) = v_LaysID( i, k )
                check_sed( i, k, 1:totgrn ) = stratLays( k, i, 1:totgrn )
                if( gporo%compaction ) &
                    check_porosity( i, k, 1:totgrn ) = porosityLays( k, i, 1:totgrn )
            enddo
        enddo

        return

    end subroutine checkpt_layers
    ! ============================================================================
    !> Subroutine write_checkpoint
    !!
    !! Subroutine write_checkpoint records layers space evolution through time.
    !! \param iter
    !<
    ! ============================================================================
    subroutine write_checkpoint

        logical :: compression

        integer :: k, p, q, n, rank
        integer :: hdferr
        real(tkind), dimension( : ),allocatable :: nodes, lays

        integer(hid_t) :: file_id, plist_id
        integer(hid_t) :: filespace, dset_id
        integer(hsize_t),dimension(2) :: dims

        character(len=128) :: text, checkf, stg

        checkf = ''
        checkf = 'checkpoint'
        call noblnk( checkf )
        text = '.'
        call append_str( checkf,text )
        call append_zero( checkf,layID )
        stg = '.h5'
        call append_str( checkf,stg )
        call noblnk(checkf)
        call addpath2( checkf )
        
        ! Gather all information from the processors
        call checkpt_layers
        
        ! Initialize predefined datatypes
        call h5open_f( hdferr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, hdferr )
        
		! Setup file access property list for MPI-IO access.
        call h5pcreate_f( h5p_file_access_f, plist_id, hdferr )
      
	    ! Create the file collectively.
        call h5fcreate_f( checkf, h5f_acc_excl_f, file_id, hdferr )
     
        ! ========================
        ! Stratal coordinates
        ! ========================
        dims( 1 ) = 3
        dims( 2 ) = lnbPts
        rank = 2
        ! Create dataspace and opens it for access
        call h5screate_simple_f( rank, dims, filespace, hdferr )
        allocate( nodes( 3 * lnbPts ) )
        n = 0
        do p = 1, lnbPts
            nodes( n + 1 ) = real( check_NbLay( p ) )
            nodes( n + 2 ) = check_base( p )
            nodes( n + 3 ) = soil_thickness( p )
            n = n + 3
        enddo
        text = ''
        text = "/DataVertices"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, hdferr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, hdferr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, hdferr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, hdferr )
            call h5pset_deflate_f( plist_id, 9, hdferr )
            call h5pset_chunk_f( plist_id, rank, dims, hdferr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, hdferr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, hdferr )
            call h5pclose_f( plist_id, hdferr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, hdferr )
        call h5sclose_f( filespace, hdferr )
  
        ! ========================
        ! Deposition space parameters
        ! ========================
        dims( 1 ) = 1 + totgrn
        if( gporo%compaction ) dims( 1 ) = dims( 1 ) + totgrn
        dims( 2 ) = lnbPts * layID
        rank = 2
        allocate( lays( dims( 1 ) * dims( 2 ) ) )
        n = 0
        lays = 0.0_8
        do p = 1, lnbPts
            do k = 1, check_NbLay( p )
                n = n + 1
                lays( n ) = real( check_layID( p , k ) )
                do q = 1, totgrn
                    n = n + 1
                    lays( n ) = check_sed( p , k , q )
                enddo
                if( gporo%compaction )then
                    do q = 1, totgrn
                        n = n + 1
                        lays( n ) = check_porosity( p , k , q )
                    enddo
                endif
            enddo
        enddo
     
        ! Create dataspace and opens it for access
        dims( 2 ) = n / dims( 1 )
        call h5screate_simple_f( rank, dims, filespace, hdferr )
        text = ''
        text = "/DepositsLayer"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, hdferr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, hdferr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, lays( 1:dims( 1 )*dims( 2) ), dims, hdferr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, hdferr )
            call h5pset_deflate_f( plist_id, 9, hdferr )
            call h5pset_chunk_f( plist_id, rank, dims, hdferr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, hdferr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, lays, dims, hdferr )
            call h5pclose_f( plist_id, hdferr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, hdferr )
        call h5sclose_f( filespace, hdferr )

        ! Close the file.
        call h5fclose_f( file_id, hdferr )
        ! Close interface
        call h5close_f( hdferr )

        deallocate( lays, nodes )

        return

    end subroutine write_checkpoint
    ! ============================================================================
    !> Subroutine read_checkpoint_simulation
    !! Subroutine read_checkpoint_simulation reads and creates a new simulation based on previous run.
    !<
    ! ============================================================================
    subroutine read_checkpoint_simulation

        logical :: simple, compression

        integer :: k, p, q, n, nbnodes, rank
        integer :: hdferr, filter, filter_all

        real(tkind) :: th
        real(tkind), dimension( : ),allocatable :: nodes, lays

        integer(hid_t) :: file_id, d_spc
        integer(hid_t) :: dset_id, dtype_id
        integer(hsize_t),dimension( 2 ) :: dims, maxdims

        character(len=128) :: text, stg

        ! Previous files for restart
        frestart = restartfolder
        call noblnk( frestart )
        stg = '/runfiles/checkpoint.'
        call append_str( frestart, stg )
        call append_zero( frestart, restart_iter )
        stg = '.h5'
        call append_str( frestart, stg )
        call noblnk( frestart )

        ! Initialize predefined datatypes
        call h5open_f( hdferr )

        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, hdferr )
        call h5zget_filter_info_f( h5z_filter_deflate_f, filter, hdferr )
        filter_all = ior( h5z_filter_encode_enabled_f, h5z_filter_decode_enabled_f )
        if( filter_all .ne. filter ) compression = .false.

        ! Open the file collectively.
        call h5fopen_f( frestart, h5f_acc_rdonly_f, file_id, hdferr )

        !-----------------------------------------------
        ! Layer coordinates
        text = ''
        text = "/DataVertices"
        call h5dopen_f( file_id, trim(text), dset_id, hdferr )
        call h5dget_type_f( dset_id, dtype_id, hdferr )
        call h5dget_space_f( dset_id, d_spc, hdferr )
        call h5sis_simple_f( d_spc, simple, hdferr )
        call h5sget_simple_extent_ndims_f( d_spc, rank, hdferr )
        call h5sget_simple_extent_dims_f( d_spc, dims, maxdims, hdferr )
        allocate( nodes( dims( 1 ) * dims( 2 ) ) )
        call h5dread_f( dset_id, h5t_native_double, nodes, dims, hdferr )

        nbnodes = int( dims( 2 ), 4 )

        if( nbnodes /= lnbPts .or. dims( 1 ) /= 3 )then
            print*,'Something is wrong in checkpoint file declaration.',iam
            stop
        endif

        n = 0
        do k = 1, lnbPts
            v_nbLays( k ) =  int( nodes( n + 1 ) )
            lbase( k ) = nodes( n + 2 )
            soil_thick( k ) = nodes( n + 3 )
            n = n + 3
        enddo

        ! Close the dataset
        call h5sclose_f( d_spc, hdferr )
        call h5dclose_f( dset_id, hdferr )

        !-----------------------------------------------
        ! Deposition space parameters
        text = ''
        text = "/DepositsLayer"
        call h5dopen_f( file_id, trim(text), dset_id, hdferr )
        call h5dget_type_f( dset_id, dtype_id, hdferr )
        call h5dget_space_f( dset_id, d_spc, hdferr )
        call h5sis_simple_f( d_spc, simple, hdferr )
        call h5sget_simple_extent_ndims_f( d_spc, rank, hdferr )
        call h5sget_simple_extent_dims_f( d_spc, dims, maxdims, hdferr )
        allocate( lays( dims( 1 ) * dims( 2 ) ) )
        call h5dread_f( dset_id, h5t_native_double, lays, dims, hdferr )

        n = 0
        do k = 1, lnbPts
            th = 0.0_8
            do p = 1, v_nbLays( k )
                n = n + 1
                v_LaysID( k, p ) = int( lays( n ) )
                do q = 1, totgrn
                    n = n + 1
                    stratLays( p, k, q ) = lays( n )
                    th = th + lays( n )
                enddo
                if( gporo%compaction )then
                    do q = 1, totgrn
                        n = n + 1
                        porosityLays( p, k, q ) = lays( n )
                    enddo
                endif
            enddo
            lcoordZ( k ) = lbase( k ) + th
        enddo

        ! Close the dataset
        call h5sclose_f( d_spc, hdferr )
        call h5dclose_f( dset_id, hdferr )

        ! Close the file.
        call h5fclose_f( file_id, hdferr )
        ! Close interface
        call h5close_f( hdferr )

        deallocate( nodes, lays )

        if( flexureon == 1 ) call initialise_isostatic_flexure( th )

        return

    end subroutine read_checkpoint_simulation
    ! ============================================================================

end module checkpoint
