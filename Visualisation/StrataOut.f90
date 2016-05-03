! ============================================================================
! Name        : StrataOut.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file StrataOut.f90
!
! Description : StrataOut is used to generate Hdf5 output of the stratal mesh.
!
!<
! ============================================================================
module strata_out

    use hdf5
    use tin_out
    use sea_out
    use parallel
    use fws_out
    use file_data
    use FoX_wxml
    use hydro_out
    use flow_data
    use time_data
    use mesh_data
    use forces_data

    implicit none

    public

    ! Output integer
    integer, dimension( : ),allocatable :: check_int

contains

    ! ============================================================================
    !> Subroutine xdmf_output
    !! Subroutine xdmf_output controls SP Model output frequency.
    !! \param iter
    !<
    ! ============================================================================
    subroutine xdmf_output( iter )

        integer :: iter

!        real( tkind ) :: ts_st3, ts_ed3

        ! TIN surface
!        ts_st3 = mpi_wtime( )
!        if( iam == 0 )then
!            call TIN_hdf5( iter )
!            call TIN_xmf( iter )
!            call TIN_series( iter )
!        endif
!        ts_ed3 = mpi_wtime( )
!        if( iam == 0 ) write(6,*)'TIN output completed ...'
!        if( iam == 0 ) write(6,100),'Time elapse (s):',ts_ed3 - ts_st3

        ! Stratigraphy
!        ts_st3 = mpi_wtime( )
        if( local_nbPts > 0 )then
            if( .not. allocated( check_int ) ) allocate( check_int( local_nbPts * 2 ) )
!            if( rain_event > 0 .or. masson == 1 ) call Hydro_hdf5( iter )
            call Hydro_hdf5( iter )
            if( transport_mode /= 1 ) call Strata_hdf5( iter )
        endif
!        ts_ed3 = mpi_wtime( )
!        if( iam == 0 ) write(6,*)'Hdf5 output completed ...'
!        if( iam == 0 ) write(6,100),'Time elapse (s):',ts_ed3 - ts_st3

!        ts_st3 = mpi_wtime( )
!        if( rain_event > 0 .or. masson == 1 )then
            call Hydro_xmf( iter )
            if( iam == 0 ) call Hydro_series( iter )
!        endif
!        ts_ed3 = mpi_wtime( )
!        if( iam == 0 ) write(6,*)'Hydro xmf output completed ...'
!        if( iam == 0 ) write(6,100),'Time elapse (s):',ts_ed3 - ts_st3

        ! Flow walkers
        if( fwsvis == 1 .and. transport_mode == 0 )then
            call FWs_hdf5( iter )
            call FW_xmf( iter )
            if( iam == 0 ) call FW_series( iter )
        endif

!        ts_st3 = mpi_wtime( )
        if( transport_mode /= 1 )then
            call Strata_xmf( iter )
            if( iam == 0 ) call Strata_series( iter )
        endif
!        ts_ed3 = mpi_wtime( )
!        if( iam == 0 ) write(6,*)'Strata xmf output completed ...'
!        if( iam == 0 ) write(6,100),'Time elapse (s):',ts_ed3 - ts_st3
!        call mpi_barrier( badlands_comm_world,ierr )

!        ts_st3 = mpi_wtime( )
        if( iam == 0 .and. seavis == 1 )then
            call Sea_xmf( iter )
            call Sea_series( iter )
        endif

!        ts_ed3 = mpi_wtime( )
!        if( iam == 0 ) write(6,*)'Sea xmf output completed ...'
!        if( iam == 0 ) write(6,100),'Time elapse (s):',ts_ed3 - ts_st3

!100   format(A18,F12.4)

        return

    end subroutine xdmf_output
    ! ============================================================================
    !> Subroutine Strata_hdf5
    !! Subroutine Strata_hdf5 uses Hdf5 to output SP Model stratal evolution.
    !! It records both connectivity and vertices.
    !<
    ! ============================================================================
    subroutine Strata_hdf5( iter )

        logical :: compression

        ! Parameters Declaration
        integer :: totnodes, totelems, id, id0, i, rank, iter, k, pid, idt, ks, gid, c, nnb
        integer, dimension( : ),allocatable :: lay_id, nlay, connect

        real(tkind) :: mgz, sedh, por
        real(tkind), dimension( : ),allocatable :: nodes, thick, mz, Ethick, Emz, Nporo
        real(tkind), dimension( :,: ),allocatable :: prop, Eprop

        character(len=128) :: text, stg, file

        integer(hid_t) :: file_id, plist_id
        integer(hid_t) :: filespace, dset_id
        integer(hsize_t),dimension(2) :: dims

        file = ''
        file = 'SMesh'
        call noblnk( file )
        stg = '.'
        call append_str( file,stg )
        call append_zero( file,iter )
        stg = '.p'
        call append_str( file,stg )
        call append_nb( file,iam )
        stg = '.h5'
        call append_str( file,stg )
        call addpath1( file )
        totnodes = local_nbPts * ( layID + 1  )
        totelems = local_Faces * ( layID )
        allocate( mz(totnodes) )
        allocate( nlay(totnodes) )
        allocate( thick(totnodes) )
        allocate( prop(totnodes,totgrn) )
        allocate( nodes(3*totnodes) )
        allocate( connect(8*totelems) )
        allocate( lay_id(local_nbPts) )
        allocate( Emz(totelems) )
        allocate( Ethick(totelems) )
        allocate( Eprop(totelems,totgrn) )
        if( gporo%compaction ) allocate( Nporo(totnodes) )

        check_int = 0

        ! Create nodes arrays
        id = 1
        id0 = 1
        idt = 1
        lay_id = 0
        mz = 0.0_8
        thick = 0.0_8
        prop = 0.0_8
        Eprop = 0.0_8
        Ethick = 0.0_8
        if( gporo%compaction ) Nporo = 0.0_8

        ! Create the basement nodes
        do i = 1, local_nbPts
            gid = global_nid( i )
            check_int( id0 ) = gid
            check_int( id0 + 1 ) = v_nbLays( lvertexID( gid ) )
            nodes( id ) = real( coordX( gid ) )
            nodes( id + 1 ) = real( coordY( gid ) )
            nodes( id + 2 ) = real( lbase( lvertexID( gid ) ) )
            thick( idt ) = 0.0_8
            prop( idt,1:totgrn ) = 0.0_8
            mz( idt ) = 0.0_8
            if( gporo%compaction )then
                Nporo( idt ) = 0.0_8
            endif
            nlay( idt ) = 0
            idt = idt + 1
            id0 = id0 + 2
            id = id + 3
        enddo

        ! Create the deposit layer nodes
        do k = 1, layID
            do i = 1, local_nbPts
                gid = global_nid( i )
                nodes( id ) = real( coordX( gid ) )
                nodes( id + 1 ) = real( coordY( gid ) )

                ! Check if there is a layer corresponding to the current layer ID
                if( lay_id( i ) < v_nbLays( lvertexID( gid ) ) )then
                    if( v_LaysID( lvertexID( gid ), lay_id( i ) + 1 ) == k )then
                        lay_id( i ) = lay_id( i ) + 1
                    endif
                endif
                if( v_LaysID(  lvertexID( gid ), lay_id( i ) ) /= k )then
                    ! Get the elevation from previous layer
                    pid = 3 * ( local_nbPts * ( k - 1 ) + i )
                    nodes( id + 2 ) = nodes( pid )
                    ! Define previous layer mean grain size in the current one
                    pid = local_nbPts * ( k - 1 ) + i
                    thick( idt ) =  0.0_8
                    prop( idt, 1:totgrn ) = prop( pid, 1:totgrn )
                    mz( idt ) = mz( pid )
                    if( gporo%compaction ) Nporo( idt ) = Nporo( pid )
!                    por = 0.0_8
!                    mgz = 0.0_8
!                    sedh = 0.0_8
!                    do ks = 1, totgrn
!                        sedh = sedh + real( stratLays( lay_id( i ), lvertexID( gid ), ks ) )
!                        prop( idt, ks ) = real( stratLays( lay_id( i ), lvertexID( gid ), ks ) )
!                        mgz = real( mgz + sediment( ks )%diameter *  stratLays( lay_id( i ), lvertexID( gid ), ks ) )
!                        if( gporo%compaction )then
!                            por = por + porosityLays( lay_id( i ), lvertexID( gid ), ks ) * prop( idt, ks )
!                        endif
!                    enddo
!                    if( sedh > 0.0_8 )then
!                        mz( idt ) = real( mgz / sedh )
!                        if( gporo%compaction ) Nporo( idt ) = ( por / sedh )
!                    else
!                        mz( idt ) = 0.0_8
!                        if( gporo%compaction ) Nporo( idt ) = 0.0_8
!                    endif
!                    thick( idt ) =  0.0_8
!                    mz( idt ) = mz( idt )
!                    mz( idt ) = mz( idt ) * 1000.0_8
                else
                    ! There is a layer take it for defining the nodes parameters
                    pid = 3 * ( local_nbPts * ( k - 1 ) + ( i  ) )
                    por = 0.0_8
                    mgz = 0.0_8
                    sedh = 0.0_8
                    do ks = 1, totgrn
                        sedh = sedh + real( stratLays( lay_id( i ), lvertexID( gid ), ks ) )
                        prop( idt, ks ) = real( stratLays( lay_id( i ), lvertexID( gid ), ks ) )
                        mgz = real( mgz + sediment( ks )%diameter *  stratLays( lay_id( i ), lvertexID( gid ), ks ) )
                        if( gporo%compaction )then
                            por = por + porosityLays( lay_id( i ), lvertexID( gid ), ks ) * prop( idt, ks )
                        endif
                    enddo
                    nodes( id + 2 ) = nodes( pid ) + sedh
                    thick( idt ) = sedh
                    if( sedh > 0.0_8 )then
                        mz( idt ) = real( mgz / sedh )
                        if( gporo%compaction ) Nporo( idt ) = ( por / sedh )
                    else
                        mz( idt ) = 0.0_8
                        if( gporo%compaction ) Nporo( idt ) = 0.0_8
                    endif
                    mz( idt ) = mz( idt )
                    mz( idt ) = mz( idt ) * 1000.0_8
                endif
                id = id + 3
                nlay( idt ) = k
                idt = idt + 1
            enddo
        enddo

        ! Create cell data
        id = 1
        idt = 1
        do i = 1, local_Faces
            connect( id : id + 3 ) = locf_points( i, 1:4 )
            id = id + 4
            connect( id : id + 3 ) = locf_points( i, 1:4 ) + local_nbPts
            Emz( idt ) = 0.0_8
            Eprop( idt, : ) = 0.0_8
            Ethick( idt ) = 0.0_8
            do c = 1, 4
                Ethick( idt ) = real( Ethick( idt ) + thick( connect( id - 1 + c ) ) )
                Emz( idt ) = real( Emz( idt ) + mz( connect( id - 1 + c ) ) )
                do ks = 1, totgrn
                    Eprop( idt, ks ) = real( Eprop( idt, ks ) + prop( connect( id - 1 + c ), ks ) )
                enddo
            enddo
            Emz( idt ) = Emz( idt ) * 0.25_8
            Ethick( idt ) = Ethick( idt ) * 0.25_8
            if( Ethick( idt ) > 0.0_8 )then
                do ks = 1, totgrn
                    Eprop( idt, ks ) = Eprop( idt, ks ) * 0.25_8
                    Eprop( idt, ks ) = Eprop( idt, ks ) / Ethick( idt )
                    Eprop( idt, ks ) = Eprop( idt, ks ) * 100.0_8
                    if( Eprop( idt, ks ) > 100.0_8 ) Eprop( idt, ks ) = 100.0_8
                enddo
            else
                Eprop( idt, : ) = 0.0_8
            endif
            idt = idt + 1
            id = id + 4
        enddo
        do k = 1, layID  - 1
            do i = 1, local_Faces
                Ethick( idt ) = 0.0_8
                connect( id : id + 3 ) = locf_points( i, 1:4 ) + k * local_nbPts
                id = id + 4
                connect( id : id + 3 ) = locf_points( i, 1:4 ) + ( k + 1 ) * local_nbPts
                nnb = 0
                Emz( idt ) = 0.0_8
                Eprop( idt, : ) = 0.0_8
                do c = 1, 4
                    Ethick( idt ) = real( Ethick( idt ) + thick( connect( id - 1 + c ) ) )
                    if( thick( connect( id - 1 + c ) ) > 0.0_8 )then
                        if( mz( connect( id - 1 + c ) ) > 0.0_8 ) nnb = nnb + 1
                        Emz( idt ) = real( Emz( idt ) + mz( connect( id - 1 + c ) ) )
                        do ks = 1, totgrn
                            Eprop( idt, ks ) = real( Eprop( idt, ks ) + prop( connect( id - 1 + c ), ks ) )
                        enddo
                    endif
                enddo
                Ethick( idt ) = Ethick( idt ) * 0.25_8
                if( nnb == 0 ) Emz( idt ) = 0.0_8
                if( nnb > 0 ) Emz( idt ) = Emz( idt ) / real( nnb )
                if( Ethick( idt ) > 0.0_8 )then
                    do ks = 1, totgrn
                        Eprop( idt, ks ) = Eprop( idt, ks ) * 0.25_8
                        Eprop( idt, ks ) = Eprop( idt, ks ) / Ethick( idt )
                        Eprop( idt, ks ) = Eprop( idt, ks ) * 100.0_8
                        if( Eprop( idt, ks ) > 100.0_8 ) Eprop( idt, ks ) = 100.0_8
                    enddo
                else
                    Eprop( idt, : ) = 0.0_8
                endif
                id = id + 4
                idt = idt + 1
            enddo
        enddo

        ! Initialize predefined datatypes
        call h5open_f( ierr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, ierr )

        ! Setup file access property list for MPI-IO access.
        call h5pcreate_f( h5p_file_access_f, plist_id, ierr )
        ! Create the file collectively.
        call h5fcreate_f( file, h5f_acc_trunc_f, file_id, ierr, access_prp = plist_id )

        ! ========================
        ! The Coordinates - vertices
        ! ========================
        dims( 1 ) = 3
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/vertices"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Thickness attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/thick"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, thick, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, thick, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Layer number type attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/layerID"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_integer, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_integer, nlay, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_integer , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_integer, nlay, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Mean grainsize attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/mgz"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, mz, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, mz, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        if( gporo%compaction )then
            ! ========================
            ! Porosity attribute
            ! ========================
            dims( 1 ) = 1
            dims( 2 ) = totnodes
            rank = 2
            call h5screate_simple_f( rank, dims, filespace, ierr )
            text = ''
            text = "/porosity"
            if( .not. compression )then
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
                ! Write the dataset
                call h5dwrite_f( dset_id, h5t_native_double, Nporo, dims, ierr, &
                    file_space_id = filespace, xfer_prp = plist_id )
            else
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
                call h5pset_deflate_f( plist_id, 9, ierr )
                call h5pset_chunk_f( plist_id, rank, dims, ierr )
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
                ! Write the dataset
                call h5dwrite_f( dset_id, h5t_native_double, Nporo, dims, ierr )
                call h5pclose_f( plist_id, ierr )
            endif
            ! Close the dataset
            call h5dclose_f( dset_id, ierr )
            call h5sclose_f( filespace, ierr )
        endif

        ! ========================
        ! Element thickness
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totelems
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/ethick"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Ethick, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Ethick, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Element mean grainsize
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totelems
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/emgz"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Emz, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Emz, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Element proportion
        ! ========================
        do ks = 1, totgrn
            dims( 1 ) = 1
            dims( 2 ) = totelems
            rank = 2
            call h5screate_simple_f( rank, dims, filespace, ierr )
            text = ''
            text = "/eperc"
            call append_nb( text,ks )
            if( .not. compression )then
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
                ! Write the datase
                call h5dwrite_f( dset_id, h5t_native_double, Eprop( 1:totelems, ks ), dims, ierr, &
                    file_space_id = filespace, xfer_prp = plist_id )
            else
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
                call h5pset_deflate_f( plist_id, 9, ierr )
                call h5pset_chunk_f( plist_id, rank, dims, ierr )
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
                ! Write the datase
                call h5dwrite_f( dset_id, h5t_native_double, Eprop( 1:totelems, ks ), dims, ierr )
                call h5pclose_f( plist_id, ierr )
            endif
            ! Close the dataset
            call h5dclose_f( dset_id, ierr )
            call h5sclose_f( filespace, ierr )
        enddo

        ! Close the file.
        call h5fclose_f( file_id, ierr )
        ! Close interface
        call h5close_f( ierr )

        deallocate( nodes, lay_id, thick, nlay, mz, prop, Emz, Eprop, Ethick, connect )

        return

    end subroutine Strata_hdf5
    ! ============================================================================
    !> Subroutine Strata_xmf
    !! Subroutine Strata_xmf generates the Strata XdmF file for the considered time simulation.
    !! It builds the file which calls all the hdf5 created.
    !<
    ! ============================================================================
    subroutine Strata_xmf( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: iter, totnodes, k, ks, nX, nY, nbZ
        integer,dimension( nproc ) :: nbX, nbY

        character(len=128) :: str, filename, file, h5file, stg
        character(len=128) :: filename1, filename2, filename3, filename4, filename5, filename0
        character(len=128) :: filename6, filename7, filename8, filename9, filename10, filename11, txt

        totnodes = local_nbPts * ( layID + 1 )

        ! Define decomposition component
        if( local_nbPts > 0 )then
            nX = int( ( coordX( global_nid( local_nbPts ) )  - coordX( global_nid( 1 ) ) ) / strat_dx + 1 )
            nY = int( ( coordY( global_nid( local_nbPts ) )  - coordY( global_nid( 1 ) ) ) / strat_dx + 1 )
        else
            nX = 0
            nY = 0
        endif
        nbZ = layID + 1
        call mpi_gather( nX,1,int_type, nbX,1,int_type,0,badlands_comm_world,ierr)
        call mpi_gather( nY,1,int_type, nbY,1,int_type,0,badlands_comm_world,ierr)

        if( iam == 0 )then
            file = ''
            file = 'Mesh'
            call noblnk( file )
            str = '.'
            call append_str( file,str )
            call append_zero( file,iter )
            str = '.xmf'
            call append_str( file,str )
            call addpath1( file )
            call xml_OpenFile( file, xf )
            ! Header
            call xml_AddDOCTYPE(xf, "Xdmf", "Xdmf.dtd")
            call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XInclude", "xi")
            call xml_NewElement(xf, "Xdmf" )
            call xml_AddAttribute(xf, "Version", "2.0" )
            call xml_NewElement( xf, "Domain" )
            call xml_NewElement( xf, "Grid" )
            call xml_AddAttribute( xf, "GridType", "Collection" )
            call xml_AddAttribute( xf, "CollectionType", "Spatial" )
            call xml_NewElement( xf, "Time" )
            call xml_AddAttribute( xf, "Type", "Single" )
            call xml_AddAttribute( xf, "Value", real( tnow ) )
            call xml_EndElement( xf, "Time" )
            ! Loop over processors
            do k = 1, gproc
                h5file = ''
                h5file = 'SMesh'
                call noblnk( h5file )
                stg = '.'
                call append_str( h5file,stg )
                call append_zero( h5file,iter )
                stg = '.p'
                call append_str( h5file,stg )
                call append_nb( h5file,k-1 )
                stg = '.h5'
                call append_str( h5file,stg )
                filename = h5file
                call noblnk( filename )
                filename1 = filename
                filename2 = filename
                filename3 = filename
                filename4 = filename
                filename5 = filename
                filename6 = filename
                filename7 = filename
                filename8 = filename
                filename9 = filename
                filename10 = filename
                if( gporo%compaction )then
                    filename0 = filename
                    str = ':/porosity'
                    call append_str( filename0,str )
                endif
                str = ':/vertices'
                call append_str( filename1,str )
                str = ':/thick'
                call append_str( filename2,str )
                str = ':/layerID'
                call append_str( filename3,str )
                str = ':/mgz'
                call append_str( filename7,str )
                str = ':/ethick'
                call append_str( filename8,str )
                str = ':/emgz'
                call append_str( filename9,str )
                str = ':/eperc'
                call append_str( filename10,str )
                ! Block begin
                call xml_NewElement( xf, "Grid" )
                str = 'SBlock.'
                call append_zero( str,iter )
                stg = '.p'
                call append_str( str,stg )
                call append_nb( str,k - 1 )
                call xml_AddAttribute( xf, "Name", trim( str ) )
                call xml_NewElement( xf, "Topology" )
                call xml_AddAttribute( xf, "TopologyType", "3DSMesh" )
                str = ' '
                call append_nb2( str, nbZ )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, nbX( k ) )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_EndElement( xf, "Topology" )
                ! Geometry
                call xml_NewElement( xf, "Geometry" )
                call xml_AddAttribute( xf, "Type", "XYZ" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, nbZ )
                call append_nb2( str, 3 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename1 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Geometry" )
                ! Layer thickness
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Layer Thickness" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, nbZ )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename2 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Layer number
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Layer Index" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Int" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename3 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Porosity
                if( gporo%compaction )then
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "Porosity" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    call xml_AddAttribute( xf, "Precision", "4" )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename0 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                endif
                ! Mean grain size
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Mean Grainsize" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename7 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Cell thickness
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Cell" )
                call xml_AddAttribute( xf, "Name", "Cell Thickness" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) - 1 )
                call append_nb2( str, nbY( k ) - 1 )
                call append_nb2( str, nbZ  - 1 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename8 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Cell mean grainsize
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Cell" )
                call xml_AddAttribute( xf, "Name", "Cell Mean Grain Size" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename9 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Cell proportion
                do ks = 1, totgrn
                    filename11 = filename10
                    call append_nb( filename11, ks )
                    txt = 'Prop_'
                    call append_str( txt,material_name( ks ) )
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Cell" )
                    call xml_AddAttribute( xf, "Name", trim( txt ) )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    call xml_AddAttribute( xf, "Precision", "8" )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename11 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                enddo
                call xml_EndElement( xf, "Grid" )
            enddo
            ! Footer
            call xml_EndElement( xf, "Grid" )
            call xml_EndElement( xf, "Domain" )
            call xml_EndElement( xf, "Xdmf" )
            call xml_Close( xf )
        endif

        return

    end subroutine Strata_xmf
    ! ============================================================================
    !> Subroutine Strata_series
    !! Subroutine Strata_series generates the XmL file for TIN surface.
    !<
    ! ============================================================================
    subroutine Strata_series( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: i, iter, it0
        character(len=128) :: filename, str, fname

        filename = 'Stratal_series.xdmf'
        call addpath1(filename)

        ! Header
        call xml_OpenFile( filename, xf )
        call xml_AddDOCTYPE(xf, "Xdmf", "Xdmf.dtd")
        call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XInclude", "xi")
        call xml_NewElement(xf, "Xdmf" )
        call xml_AddAttribute(xf, "Version", "2.0" )
        call xml_NewElement( xf, "Domain" )
        call xml_NewElement( xf, "Grid" )
        call xml_AddAttribute( xf, "GridType", "Collection" )
        call xml_AddAttribute( xf, "CollectionType", "Temporal" )
        it0 = 1
        ! Loop over time step
        do i = it0, iter
            ! Grid name
            fname = ''
            fname = 'Mesh'
            call noblnk( fname )
            str = '.'
            call append_str( fname,str )
            call append_zero( fname,i)
            str = '.xmf'
            call append_str( fname,str )
            call xml_NewElement( xf, "xi:include" )
            call xml_AddAttribute( xf, "href", trim( fname ) )
            call xml_AddAttribute( xf, "xpointer", "xpointer(//Xdmf/Domain/Grid)" )
            call xml_EndElement( xf, "xi:include" )
        enddo
        ! Footer
        call xml_EndElement( xf, "Grid" )
        call xml_EndElement( xf, "Domain" )
        call xml_EndElement( xf, "Xdmf" )
        call xml_Close( xf )

        return

    end subroutine Strata_series
  ! ============================================================================

end module strata_out
