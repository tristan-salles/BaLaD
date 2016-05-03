! ============================================================================
! Name        : Mesh_parameters.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Mesh_parameters.f90
!
! Description :  It encapsulates the sediment and stratal module used by SPModel during compilation.
!
!<
! ============================================================================
!> Module mesh_data
!<
module mesh_data

    use parallel
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    ! Transport mode flag: 0 for hydro transport, 1 for erosion only, 2 for diffusion only
    integer :: transport_mode

    ! Visualisation flags
    integer :: seavis, fwsvis, masson, flexureon

    ! Borders size
    integer, parameter :: border = 2

    ! Initial deposit layer number
    integer :: InitDep

    ! Number of points on the large cartesian grid
    integer :: lstrat_X, lstrat_Y

    ! Xo & Yo coordinates of the SW in the large cartesian grid
    real( tkind ) :: lstrat_xo, lstrat_yo

    ! Xm & Ym coordinates of the NE in the large cartesian grid
    real( tkind ) :: lstrat_xm, lstrat_ym

    ! Number of points on the initial cartesian grid
    integer :: strat_X, strat_Y

    ! Xo & Yo coordinates of the SW in the initial cartesian grid
    real( tkind ) :: strat_xo, strat_yo

    ! Xo & Yo coordinates of the NE in the initial cartesian grid
    real( tkind ) :: strat_xm, strat_ym

    ! Grid spacing
    real( tkind ) :: strat_dx

    ! Generate random noise (in metre)
    real( tkind ) :: noisegen

    ! Number of points of the large grid
    integer :: lnbPts

    ! Number of points of the grid
    integer :: nbPts

    ! Number of points of a considered processor
    integer :: local_nbPts, local_lnbPts

    ! Number of faces of a considered processor
    integer :: local_Faces

    ! Number of faces of the DEM grid
    integer :: nbFcs

    ! Number of faces of the grid
    integer :: lnbFcs

    ! Number of total grains
    integer :: totgrn

    ! Number of siliciclastics
    integer :: silgrn

    ! Number of carbonates
    integer :: carbgrn

    ! Number of organics
    integer :: orggrn

    ! Minimum slope value (dz/dx)
    real( tkind ) :: minimum_slp

    ! Aerial maximal depth to bedrock
    real( tkind ) :: dtb_aerial

    ! Marine maximal depth to bedrock
    real( tkind ) :: dtb_marine

    ! Material name
    character(len=128), dimension( : ), allocatable :: material_name

    !> Sediment parameters type
    type sediment_parameters

        ! Diameter in millimetres
        real( tkind ) :: diameter

        ! Density in kg.m-3
        real( tkind ) :: density

        ! Fall velocity
        real( tkind ) :: vfall

        ! Marine stability slope
        real( tkind ) :: slp_marine

        ! Aerial stability slope
        real( tkind ) :: slp_aerial

        ! Erodibility coefficient (kg m-2 s-1)
        real( tkind ) :: ero_coeff

        ! Expected transport time (percentage of a year)
        real( tkind ) :: transport

        ! Diffusivity coefficient for aerial soil
        real( tkind ) :: creep_aerial

        ! Diffusivity coefficient for marine soil
        real( tkind ) :: creep_marine

    end type sediment_parameters
    type(sediment_parameters),dimension(:),allocatable :: sediment

    ! Define nodes coordinates for the large grid
    real( tkind ), dimension(:), allocatable :: lcoordX, lcoordY, lcoordZ, lbase

    ! Define neighbourhood array for the large grid
    integer, dimension( 4 ) :: boundcond

    ! Define neighbourhood array for the large grid
    integer, dimension( :,: ), allocatable :: lngbID

    ! Define nodes coordinates for the DEM grid
    real( tkind ), dimension(:), allocatable :: coordX, coordY, coordZ

    ! Define faces coordinates for the DEM grid
    real( tkind ), dimension(:), allocatable :: fcoordX, fcoordY

    ! Define neighbourhood array for the DEM grid
    integer, dimension( :,: ), allocatable :: ngbID

    ! Define face nodes ID array for the DEM grid
    integer, dimension( :,: ), allocatable :: fptIDs

    ! Define face nodes ID array for the grid
    integer, dimension( :,: ), allocatable :: lfptIDs

    ! Define large grid vertex ID for the DEM nodes
    integer, dimension( : ), allocatable :: lvertexID

    ! Define processor ID containing the faces of the DEM grid
    integer, dimension( : ), allocatable :: fPid

    ! Define isostatic flexure deformation
    real( tkind ),dimension( : ),allocatable :: flexure

    ! Global node ID of a considered processor
    integer, dimension( : ), allocatable :: global_nid, global_lnid

    ! Total number of deposit layers
    integer :: nbLays

    ! Current layer ID
    integer :: layID

    ! Define the number of stratigraphic layers deposited on each vertex
    integer, dimension(:), allocatable :: v_nbLays

    ! Define the stratigraphic layers ID on each vertex
    integer, dimension(:,:), allocatable :: v_LaysID

    ! Define stratigraphic layers composition
    real( tkind ), dimension(:,:,:), allocatable :: stratLays

    ! Define stratigraphic layers porosity
    real( tkind ), dimension(:,:,:), allocatable :: porosityLays

    ! Local stratal face points ID
    integer,dimension(:,:),allocatable :: locf_points

    ! Number of rows on each partition
    integer, dimension( : ), allocatable :: part_rows

    ! DEM min / max slopes
    real( tkind ) :: slpmax, slpmin

    ! KdTree data
    real( tkind ),dimension( :, : ), allocatable :: Fdata
    type(kdtree2), pointer :: Ftree

    !----------------------------
    ! Morphometric arrays
    !----------------------------

    ! Define depression-less DEM elevation
    real( tkind ), dimension(:), allocatable :: lfilldem, filldem

    ! Define flow accumulation for the considered DEM
    real( tkind ), dimension(:), allocatable :: lfacc, facc

    ! Define rain accumulation for the considered DEM
    real( tkind ), dimension(:), allocatable :: lfacr, facr

    ! Morphometrics variables
    real( tkind ), dimension(:), allocatable :: slp, vcurv, hcurv, orient

    !----------------------------
    ! TIN arrays
    !----------------------------

    ! Maximum number of neighbors for a given node
    integer, parameter :: mxnghb = 25

    ! TIN number of points
    integer :: TIN_nbPts

    ! Define vertex number per processor for the TIN
    integer, dimension(:), allocatable :: tinprocNb

    ! Define TIN vertex global ID present on each processor
    integer, dimension(:), allocatable :: tinprocID

    ! Define neigbor vertex number for the TIN grid vertex
    integer, dimension(:), allocatable :: tinv_Ngb

    ! Define neigbor vertex IDs for the TIN grid vertex
    integer, dimension(:,:), allocatable :: tinv_NbfIDs

    ! TIN single flow distance
    integer, dimension( : ), allocatable :: SFdirection

    ! TIN number of points
    integer :: TIN_faces

    ! Define nodes coordinates for the TIN grid
    real( tkind ), dimension(:), allocatable :: tcoordX, tcoordY, tcoordZ

    ! Define nodes ID for the TIN grid faces
    integer, dimension(:,:), allocatable :: tinf_PtID

    ! Define neigbor face number for the TIN grid faces
    integer, dimension(:), allocatable :: tinf_Ngb

    ! Define neigbor face IDs for the TIN grid faces
    integer, dimension(:,:), allocatable :: tinf_NbfIDs

    ! Define centroid for the TIN grid faces
    real( tkind ), dimension(:,:), allocatable :: tinf_centroid

    ! KdTree data
    real( tkind ),dimension( :, : ), allocatable :: TFdata
    type(kdtree2), pointer :: TFtree

    !----------------------------
    ! Creep grid parameters
    !----------------------------

    ! Initial parameters for soil production
    real( tkind ) :: h0, P0, initsoil

    ! Soil thickness
    real( tkind ), dimension( : ), allocatable :: soil_thick

    ! Influx and proportion of sediment in a given cell
    real( tkind ), dimension( :,: ), allocatable :: sedin, prop

    ! MFD weights
    real( tkind ), dimension( :,: ), allocatable :: fluxweight, fluxweightdem

    ! Number of diffusion cycle
    integer, parameter :: diff_nb = 2

    ! Maximum number of iterations for the diffusion
    integer, parameter :: max_it_cyc = 10000

    real( tkind ), parameter :: diff_res = 1.e-4_8
    real( tkind ), parameter :: toplimit = 1.e10_8

    ! Creep grid elevation
    real( tkind ), dimension( : ), allocatable :: creepz

    ! Stratal sediment composition
    real( tkind ), dimension(:,:),allocatable :: top_sedh

    ! Previous stratal sediment composition
    real( tkind ), dimension(:,:),allocatable :: top_sedprev

    real( tkind ), dimension(:,:), allocatable :: depo
    real( tkind ), dimension(:,:), allocatable :: dstart
    real( tkind ), dimension( : ), allocatable :: cdif
    real( tkind ), dimension( : ), allocatable :: difo
    real( tkind ), dimension( : ), allocatable :: difp

    !----------------------------
    ! Restart and checkpointing parameters
    !----------------------------
    logical :: restart, checkpointing

    ! Check pointing frequency
    integer :: checkfreq

    ! Restart iteration
    integer :: restart_iter

    integer, dimension(:), allocatable :: check_NbLay
    integer, dimension(:,:), allocatable :: check_layID
    real( tkind ), dimension( : ), allocatable :: check_base
    real( tkind ), dimension( : ), allocatable :: soil_thickness
    real( tkind ), dimension( : ), allocatable :: check_load
    real( tkind ), dimension(:,:,:), allocatable :: check_sed
    real( tkind ), dimension(:,:,:), allocatable :: check_porosity

contains

    ! ============================================================================
    !> Function cmp_centroid()
    !! Function cmp_centroid is used to compute the centroid of the face.
    !! \param fce
    !<
    ! ============================================================================
    function cmp_centroid( fce ) result( centroid )

        integer :: fce, k
        integer, dimension( 3 ) :: nids

        real( tkind ), dimension( 2 ) :: centroid

        nids = tinf_PtID( fce , 1:3 )
        centroid = 0.0_8
        do k = 1, 3
            centroid( 1 ) = centroid( 1 ) + tcoordX( nids( k ) )
            centroid( 2 ) = centroid( 2 ) + tcoordY( nids( k ) )
        enddo

        do k = 1 , 2
            centroid( k ) = centroid( k ) / 3
        enddo

    end function cmp_centroid
    ! ============================================================================
    !> Function compute_distance()
    !! Function distance is used to return the distance between two nodes of a considered face
    !! \param node1, node2
    !<
    ! ============================================================================
    function compute_distance( node1x, node1y, node2x, node2y  ) result( dist )

        real( tkind ), intent( in ) :: node1x, node1y, node2x, node2y
        real( tkind ) :: dist

        dist = ( node1x - node2x )**2
        dist = dist + ( node1y  - node2y )**2
        dist = sqrt( dble( dist ) )

    end function compute_distance
    ! ============================================================================

end module mesh_data
! ============================================================================
