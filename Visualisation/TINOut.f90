! ============================================================================
! Name        : TINOut.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file TINOut.f90
!
! Description : TINOut is used to generate Hdf5 output of the TIN surface.
!
!<
! ============================================================================
module tin_out

    use hdf5
    use parallel
    use file_data
    use time_data
    use FoX_wxml
    use flow_data
    use mesh_data

    implicit none

    public

    integer ::  totnodes, totelems

contains

    ! ============================================================================
    !> Subroutine TIN_hdf5
    !! Subroutine TIN_hdf5 uses Hdf5 to output SP Model TIN evolution.
    !! It records both connectivity and vertices.
    !<
    ! ============================================================================
    subroutine TIN_hdf5( iter )

        logical :: compression

        ! Parameters Declaration
        integer :: id, i, rank, iter, kk
        integer, dimension( : ),allocatable :: connect
        real(tkind), dimension( : ),allocatable :: nodes

        character(len=128) :: text, stg, file

        integer(hid_t) :: file_id, plist_id
        integer(hid_t) :: filespace, dset_id
        integer(hsize_t),dimension(2) :: dims

        file = ''
        file = ftin_out
        call noblnk( file )
        stg = '.'
        call append_str( file,stg )
        call append_zero( file,iter )
        stg = '.h5'
        call append_str( file,stg )
        call addpath1( file )
        totnodes =  TIN_nbPts
        totelems =  TIN_faces

        allocate( nodes(3*totnodes) )
        allocate( connect(3*totelems) )

        ! Create nodes arrays
        id=1
        kk = 0
        do i = 1, TIN_nbPts
            nodes( id ) = real( tcoordX( i ) )
            nodes( id + 1 ) = real( tcoordY( i ) )
            nodes( id + 2 ) = real( tcoordZ( i ) )
            id = id + 3
            kk = kk +1
        enddo

        ! Initialize predefined datatypes
        call h5open_f( ierr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, ierr )

        ! Setup file access property list for MPI-IO access.
        call h5pcreate_f( h5p_file_access_f, plist_id, ierr )

        ! Create the file collectively.
        call h5fcreate_f( file, h5f_acc_trunc_f, file_id, ierr, access_prp = plist_id )

        ! Surface - connectivity
        id = 1
        kk = 0
        do i = 1, TIN_faces
            kk = kk + 1
            connect( id ) = tinf_PtID( i, 1 )
            connect( id + 1 ) = tinf_PtID( i, 2 )
            connect( id + 2 ) = tinf_PtID( i, 3 )
            id = id  + 3
        enddo

        dims( 1 ) = 3
        dims( 2 ) = totelems
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/connectivity"
        if( .not. compression )then
            dims( 1 ) = 1
            dims( 2 ) = totelems*3
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_integer, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_integer, connect, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            dims( 1 ) = 1
            dims( 2 ) = totelems*3
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_integer , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_integer, connect, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! The Coordinates - vertices
        dims( 1 ) = 3
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/vertices"
        if( .not. compression )then
            dims( 1 ) = 1
            dims( 2 ) = totnodes*3
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
            dims( 1 ) = 1
            dims( 2 ) = totnodes*3
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! Close the file.
        call h5fclose_f( file_id, ierr )
        ! Close interface
        call h5close_f( ierr )

        deallocate( nodes, connect )

        return

    end subroutine TIN_hdf5
    ! ============================================================================
    !> Subroutine TIN_xmf
    !! Subroutine TIN_xmf generates the TIN XdmF file for the considered time simulation.
    !! It builds the file which calls all the hdf5 created.
    !<
    ! ============================================================================
    subroutine TIN_xmf( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: iter
        character(len=128) :: str, filename, filename1, file

        file = ''
        file = ftin_out
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
        filename = ''
        filename = ftin_out
        call noblnk( filename )
        str = '.'
        call append_str( filename,str )
        call append_zero( filename,iter )
        str = '.h5'
        call append_str( filename,str )
        filename1 = filename
        str = ':/connectivity'
        call append_str( filename, str )
        str = ':/vertices'
        call append_str( filename1,str )
        ! Block begin
        call xml_NewElement( xf, "Grid" )
        str = 'Block_TIN.'
        call append_zero( str,iter )
        call xml_AddAttribute( xf, "Name", trim(str) )
        call xml_NewElement( xf, "Topology" )
        call xml_AddAttribute( xf, "Type", "Triangle" )
        call xml_AddAttribute( xf, "NumberOfElements", totelems )
        call xml_AddAttribute( xf, "BaseOffset", "1" )
        call xml_NewElement( xf, "DataItem" )
        call xml_AddAttribute( xf, "Format", "HDF" )
        call xml_AddAttribute( xf, "DataType", "Int" )
        str = ' '
        call append_nb2( str, totelems )
        call append_nb2( str, 3 )
        call xml_AddAttribute( xf, "Dimensions", trim( str ) )
        call xml_AddCharacters( xf,trim( filename ) )
        call xml_EndElement( xf, "DataItem" )
        call xml_EndElement( xf, "Topology" )
        ! Geometry
        call xml_NewElement( xf, "Geometry" )
        call xml_AddAttribute( xf, "Type", "XYZ" )
        call xml_NewElement( xf, "DataItem" )
        call xml_AddAttribute( xf, "Format", "HDF" )
        call xml_AddAttribute( xf, "NumberType", "Float" )
        call xml_AddAttribute( xf, "Precision", "8" )
        str = ' '
        call append_nb2( str, totnodes )
        call append_nb2( str, 3 )
        call xml_AddAttribute( xf, "Dimensions", trim( str ) )
        call xml_AddCharacters( xf, trim( filename1 ) )
        call xml_EndElement( xf, "DataItem" )
        call xml_EndElement( xf, "Geometry" )
        call xml_EndElement( xf, "Grid" )
        ! Footer
        call xml_EndElement( xf, "Grid" )
        call xml_EndElement( xf, "Domain" )
        call xml_EndElement( xf, "Xdmf" )
        call xml_Close( xf )

        return

    end subroutine TIN_xmf
    ! ============================================================================
    !> Subroutine TIN_series
    !! Subroutine TIN_series generates the XmL file for TIN surface.
    !<
    ! ============================================================================
    subroutine TIN_series( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: i, iter, it0
        character(len=128) :: filename, str, fname

        filename = 'TIN_series.xdmf'
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
            fname = ftin_out
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

    end subroutine TIN_series
  ! ============================================================================

end module TIN_out
