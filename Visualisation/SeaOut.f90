! ============================================================================
! Name        : SeaOut.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file SeaO.f90
!
! Description : SeaO is used to generate Hdf5 output of the water surface during SPModel
! simulation.
!
!<
! ============================================================================
module sea_out

   use hdf5
   use parallel
   use file_data
   use time_data
   use FoX_wxml
   use mesh_data
   use forces_data

   implicit none

   public

contains

   ! ============================================================================
   !> Subroutine Sea_xmf
   !! Subroutine Sea_xmf generates the water surface XdmF file for the considered time simulation.
   !<
   ! ============================================================================
   subroutine Sea_xmf( iter )

      ! Parameters Declaration
      type(xmlf_t) :: xf

      integer :: iter, i, j
      real( tkind ) :: minx,miny,maxx,maxy
      real( tkind ) :: pt(4, 3 )

      character(len=128) :: str, file

      ! Create nodes arrays
      minx = strat_xo
      miny = strat_yo
      maxx = strat_xm
      maxy = strat_ym
      pt( 1, 1 ) = minx
      pt( 2, 1 ) = minx
      pt( 3, 1 ) = maxx
      pt( 4, 1 ) = maxx
      pt( 1, 2 ) = miny
      pt( 1, 3 ) = gsea%actual_sea
      pt( 2, 2 ) = maxy
      pt( 2, 3 ) = gsea%actual_sea
      pt( 3, 2 ) = miny
      pt( 3, 3 ) = gsea%actual_sea
      pt( 4, 2 ) = maxy
      pt( 4, 3 ) = gsea%actual_sea
      file = ''
      file = 'Sea'
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
      call xml_AddAttribute( xf, "GridType", "Uniform" )
      call xml_NewElement( xf, "Time" )
      call xml_AddAttribute( xf, "Type", "Single" )
      call xml_AddAttribute( xf, "Value",real( tnow ) )
      call xml_EndElement( xf, "Time" )
      ! Topology
      call xml_NewElement( xf, "Topology" )
      call xml_AddAttribute( xf, "TopologyType", "Quadrilateral" )
      call xml_AddAttribute( xf, "Dimensions", "1" )
      call xml_NewElement( xf, "DataItem" )
      call xml_AddAttribute( xf, "NumberType", "Int" )
      call xml_AddAttribute( xf, "Dimensions", "1 4" )
      call xml_AddAttribute( xf, "Precision", "4" )
      call xml_AddAttribute( xf, "Format", "XML" )
      call xml_AddNewline( xf )
      call xml_AddCharacters(xf, "0 1 3 2")
      call xml_AddNewline( xf )
      call xml_EndElement( xf, "DataItem" )
      call xml_EndElement( xf, "Topology" )
      ! Geometry
      call xml_NewElement( xf, "Geometry" )
      call xml_AddAttribute( xf, "GeometryType", "XYZ" )
      call xml_NewElement( xf, "DataItem" )
      call xml_AddAttribute( xf, "NumberType", "Float" )
      call xml_AddAttribute( xf, "Dimensions", "4 3" )
      call xml_AddAttribute( xf, "Precision", "4" )
      call xml_AddAttribute( xf, "Format", "XML" )
      call xml_AddNewline( xf )
      do i = 1, 4
         str = ' '
         do j = 1, 3
            call append_nbreal( str, pt( i, j ) )
         enddo
         call xml_AddCharacters( xf, trim( str ) )
         call xml_AddNewline( xf )
      enddo
      call xml_EndElement( xf, "DataItem" )
      call xml_EndElement( xf, "Geometry" )
      ! Footer
      call xml_EndElement( xf, "Grid" )
      call xml_EndElement( xf, "Domain" )
      call xml_EndElement( xf, "Xdmf" )
      call xml_Close( xf )

      return

   end subroutine Sea_xmf
   ! ============================================================================
   !> Subroutine Sea_series
   !! Subroutine Sea_series generates the XmL file for sea surface.
   !<
   ! ============================================================================
   subroutine Sea_series( iter )

      ! Parameters Declaration
      type(xmlf_t) :: xf

      integer :: i, iter, it0
      character(len=128) :: filename, str, fname

      filename = 'Sea_series.xdmf'
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
         fname = 'Sea'
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

   end subroutine Sea_series
  ! ============================================================================

end module sea_out
