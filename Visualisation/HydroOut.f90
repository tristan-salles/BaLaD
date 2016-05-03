! ============================================================================
! Name        : HydroOut.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file HydroOut.f90
!
! Description : HydroOut is used to generate Hdf5 output of the hydraulic parameter surface.
!
!<
! ============================================================================
module hydro_out

    use hdf5
    use parallel
    use file_data
    use time_data
    use FoX_wxml
    use flow_data
    use mesh_data

    implicit none

    public

    integer ::  tnodes
    real(tkind), dimension( : ),allocatable :: nodes, thflowo, vflowo, qflowo, wflowo, soil, flex

contains

    ! ============================================================================
    !> Subroutine Hydro_hdf5
    !! Subroutine Hydro_hdf5 uses Hdf5 to output SP Model hydraulics evolution.
    !! It records both connectivity and vertices.
    !<
    ! ============================================================================
    subroutine Hydro_hdf5( iter )

        logical :: compression

        ! Parameters Declaration
        integer :: id, i, rank, iter, kk, kkk
        character(len=128) :: text, stg, file

        integer(hid_t) :: file_id, plist_id
        integer(hid_t) :: filespace, dset_id
        integer(hsize_t),dimension(2) :: dims

        file = ''
        file = 'HMesh'
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
        if( iter == 1 )then
            tnodes =  local_nbPts
            allocate( nodes(3*tnodes) )
            allocate( thflowo(tnodes) )
            allocate( qflowo(tnodes) )
            allocate( vflowo(tnodes) )
            allocate( wflowo(tnodes) )
            allocate( soil(tnodes) )
            if( flexureon == 1 ) allocate( flex( tnodes ) )
        endif

        ! Create nodes arrays
        id = 1
        do i = 1, local_nbPts
            kk = global_nid( i )
            kkk = lvertexID( kk )
            nodes( id ) = real( coordX( kk ) )
            nodes( id + 1 ) = real( coordY( kk ) )
            nodes( id + 2 ) = real( coordZ( kk ) )
            id = id + 3
            thflowo( i ) = real( hflow( kk ) )
            qflowo( i ) = real( qflow( kk ) )
            if( num_src == 0 .and. rain_event == 0 .and. masson == 1 ) &
                qflowo( i ) = real( facc( kk ) )
            vflowo( i ) = real( vflow( kk ) )
            wflowo( i ) = real( wflow( kk ) )
            ! Soil thickness
            soil( i ) = real( soil_thick( kkk ) )
            ! Flexure
            if( flexureon == 1 ) flex( i ) = real( flexure( kkk ) )
            ! Slope
!            vflowo( i ) = real( slp( kk ) )
            ! Curvature
!            wflowo( i ) = real( hcurv( kk ) )
            ! FACC
!            qflowo( i ) = real( hcurv( kk ) + vcurv( kk ) )
        enddo

        ! Initialize predefined datatypes
        call h5open_f( ierr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, ierr )

        ! Setup file access property list for MPI-IO access.
        call h5pcreate_f( h5p_file_access_f, plist_id, ierr )

        ! Create the file collectively.
        call h5fcreate_f( file, h5f_acc_trunc_f, file_id, ierr, access_prp = plist_id )

        ! The Coordinates - vertices
        dims( 1 ) = 3
        dims( 2 ) = tnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/vertices"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
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

        ! Flow height
        dims( 1 ) = 1
        dims( 2 ) = tnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/hflow"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, thflowo, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, thflowo, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! Flow rate
        dims( 1 ) = 1
        dims( 2 ) = tnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/qflow"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, qflowo, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, qflowo, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! Flow velocity
        dims( 1 ) = 1
        dims( 2 ) = tnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/vflow"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, vflowo, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, vflowo, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! Flow width
        dims( 1 ) = 1
        dims( 2 ) = tnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/wflow"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, wflowo, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, wflowo, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! Soil thickness
        dims( 1 ) = 1
        dims( 2 ) = tnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/hsoil"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, soil, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, soil, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! Flexure
        if( flexureon == 1 )then
            dims( 1 ) = 1
            dims( 2 ) = tnodes
            rank = 2
            call h5screate_simple_f( rank, dims, filespace, ierr )
            text = ''
            text = "/flex"
            if( .not. compression )then
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
                ! Write the dataset
                call h5dwrite_f( dset_id, h5t_native_double, flex, dims, ierr, &
                    file_space_id = filespace, xfer_prp = plist_id )
            else
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
                call h5pset_deflate_f( plist_id, 9, ierr )
                call h5pset_chunk_f( plist_id, rank, dims, ierr )
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
                ! Write the dataset collectively
                call h5dwrite_f( dset_id, h5t_native_double, flex, dims, ierr )
                call h5pclose_f( plist_id, ierr )
            endif
            ! Close the dataset
            call h5dclose_f( dset_id, ierr )
            call h5sclose_f( filespace, ierr )
        endif

        ! Close the file.
        call h5fclose_f( file_id, ierr )
        ! Close interface
        call h5close_f( ierr )

        return

    end subroutine Hydro_hdf5
    ! ============================================================================
    !> Subroutine Hydro_xmf
    !! Subroutine Hydro_xmf generates the Hydraulics XdmF file for the considered time simulation.
    !! It builds the file which calls all the hdf5 created.
    !<
    ! ============================================================================
    subroutine Hydro_xmf( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: iter, nX, nY, k
        integer,dimension( nproc ) :: nbX, nbY

        character(len=128) :: str, filename, filename1, filename2, h5file, filename7
        character(len=128) :: filename3, filename4, filename5, filename6, file, stg

        ! Define decomposition component
        if( local_nbPts > 0 )then
            nX = int( ( coordX( global_nid( local_nbPts ) )  - coordX( global_nid( 1 ) ) ) / strat_dx + 1 )
            nY = int( ( coordY( global_nid( local_nbPts ) )  - coordY( global_nid( 1 ) ) ) / strat_dx + 1 )
        else
            nX = 0
            nY = 0
        endif
        call mpi_gather( nX,1,int_type, nbX,1,int_type,0,badlands_comm_world,ierr)
        call mpi_gather( nY,1,int_type, nbY,1,int_type,0,badlands_comm_world,ierr)

        if( iam == 0 )then
            file = ''
            file = 'HMesh'
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
            call xml_AddAttribute( xf, "Value",real( tnow ) )
            call xml_EndElement( xf, "Time" )
            ! Loop over processors
            do k = 1, gproc
                h5file = ''
                h5file = 'HMesh'
                call noblnk( h5file )
                str = '.'
                call append_str( h5file,str )
                call append_zero( h5file,iter )
                stg = '.p'
                call append_str( h5file,stg )
                call append_nb( h5file,k-1 )
                str = '.h5'
                call append_str( h5file,str )
                filename = h5file
                filename1 = filename
                filename2 = filename
                filename3 = filename
                filename4 = filename
                filename5 = filename
                filename6 = filename
                filename7 = filename
                str = ':/vertices'
                call append_str( filename1,str )
                str = ':/hflow'
                call append_str( filename2,str )
                str = ':/qflow'
                call append_str( filename3,str )
                str = ':/vflow'
                call append_str( filename4,str )
                str = ':/wflow'
                call append_str( filename5,str )
                str = ':/hsoil'
                call append_str( filename6,str )
                str = ':/flex'
                call append_str( filename7,str )
                ! Block begin
                call xml_NewElement( xf, "Grid" )
                str = 'GBlock.'
                call append_zero( str,iter )
                stg = '.p'
                call append_str( str,stg )
                call append_nb( str,k - 1 )
                call xml_AddAttribute( xf, "Name", trim(str) )
                call xml_NewElement( xf, "Topology" )
                call xml_AddAttribute( xf, "TopologyType", "3DSMesh" )
                str = ' '
                call append_nb2( str, 1 )
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
                call append_nb2( str, 1 )
                call append_nb2( str, 3 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename1 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Geometry" )
                ! Flow height
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Flow height" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, 1 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename2 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Flow rate
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Flow discharge" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, 1 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename3 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Flow velocity
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Flow velocity" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, 1 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename4 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Soil depth
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Soil thickness" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, 1 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename6 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Flexure
                if( flexureon == 1 )then
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "Flexure" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    call xml_AddAttribute( xf, "Precision", "4" )
                    str = ' '
                    call append_nb2( str, nbX( k ) )
                    call append_nb2( str, nbY( k ) )
                    call append_nb2( str, 1 )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename7 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                endif
                ! Flow width
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Flow width" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, 1 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename5 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                call xml_EndElement( xf, "Grid" )
            enddo
            ! Footer
            call xml_EndElement( xf, "Grid" )
            call xml_EndElement( xf, "Domain" )
            call xml_EndElement( xf, "Xdmf" )
            call xml_Close( xf )
        endif

        return

    end subroutine Hydro_xmf
    ! ============================================================================
    !> Subroutine Hydro_series
    !! Subroutine Hydro_series generates the XmL file for hydraulics surface.
    !<
    ! ============================================================================
    subroutine Hydro_series( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: i, iter, it0
        character(len=128) :: filename, str, fname

        filename = 'HMesh_series.xdmf'
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
            fname = 'HMesh'
            call noblnk( fname )
            str = '.'
            call append_str( fname,str )
            call append_zero( fname,i )
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

  end subroutine Hydro_series
  ! ============================================================================

end module hydro_out
! ============================================================================
