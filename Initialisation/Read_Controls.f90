! ============================================================================
! Name        : Read_Controls.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Read_Controls.f90
!
! Description : Read_Controls is used to gather the information within the SPModel XmL
! input file to build the external control forces (tectonic & eustatism.
!
!<
! ============================================================================
module forces

    use parallel
    use FoX_sax
    use file_data
    use flow_data
    use time_data
    use mesh_data
    use forces_data
    use FoX_common

    implicit none

    public

    integer :: dispn

    ! Sea Parameters
    logical, save :: seasection = .false.
    logical, save :: in_Sea = .false.
    logical, save :: in_SeaFile = .false.
    character(len=128), save :: SeaFile

    ! Displacements Parameters
    logical, save :: vdispsection = .false.
    logical, save :: in_DispF = .false.
    logical, save :: in_DispNb = .false.
    logical, save :: in_DispFile = .false.
    logical, save :: in_DispET = .false.
    logical, save :: in_DispST = .false.
    logical, save :: in_disp = .false.
    character(len=128), save :: DispF,DispNb,DispFile,DispET,DispST

    ! Porosity Parameters
    logical, save :: in_Porosity = .false.
    logical, save :: in_PorVal = .false.
    logical, save :: in_EffP = .false.
    logical, save :: in_pFile = .false.
    character(len=128), save :: porval, effp, pfile

    ! Hemipelagic Parameters
    logical, save :: hemisection = .false.
    logical, save :: in_Hemi = .false.
    logical, save :: in_HemiFile = .false.
    character(len=128), save :: HemiFile

    ! Flexure Paramaters
    logical, save :: in_crust = .false.
    logical, save :: in_rigidity  = .false.
    logical, save :: in_flex = .false.
    logical, save :: in_flextime = .false.
    character(len=128), save :: flextime, crust, rigidity


    ! Underworld Parameters
    logical, save :: udwsection = .false.
    logical, save :: in_SynchT = .false.
    logical, save :: in_SynchF = .false.
    logical, save :: in_UDW = .false.
    character(len=128), save :: synchF,synchT

contains

    ! ============================================================================
    subroutine startDocument_handler

    end subroutine startDocument_handler
    ! ============================================================================
    subroutine startElement_handler(namespaceURI, localname, name, atts)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name
        type(dictionary_t), intent(in) :: atts

        ! Sea Element
        if (name=='SeaLevelFluctuations') in_Sea = .true.
        if (name=='SeaLevelFluctuations') seasection = .true.
        if(in_Sea) call SseaElement_handler(name)

        ! Vertical Displacements Parameters
        if (name=='Displacement') in_DispF = .true.
        if (in_DispF) vdispsection = .true.
        if(in_DispF) call SvdispElement_handler(name)
        if (name=='disp')then
            in_disp = .true.
            dispn = dispn + 1
        endif
        if(in_disp) call SdispElement_handler(name)

        ! Isostatic flexure
        if (name=='Flexure') in_flex = .true.
        if( in_flex ) flexureon = 1
        if( in_flex ) call SflexElement_handler( name )

        ! Porosity Element
        if (name=='Porosity') in_Porosity = .true.
        if(in_Porosity) call SporosityElement_handler(name)

        ! Hemipelagic Element
        if (name=='HemipelagicRates') in_Hemi = .true.
        if (name=='HemipelagicRates') hemisection = .true.
        if(in_Hemi) call SHemiElement_handler(name)

        ! Underworld Element
        if (name=='UnderworldPlugin') in_UDW = .true.
        if (in_UDW) udw_plug = .true.
        if (in_UDW) call SudwElement_handler(name)

    end subroutine startElement_handler
    ! ============================================================================
    subroutine endElement_handler(namespaceURI, localname, name)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name

        ! Sea Element
        call EseaElement_handler(name)

        ! Vertical Displacements Parameters
        call EvdispElement_handler(name)
        call EdispElement_handler(name)

        ! Isostatic flexure
        call EflexElement_handler( name )

        ! Porosity Element
        call EporosityElement_handler(name)

        ! Hemipelagic Element
        call EHemiElement_handler(name)

        ! UDW Element
        call EudwElement_handler(name)

    end subroutine endElement_handler
    ! ============================================================================
    subroutine characters_handler(chars)

        character(len=*), intent(in) :: chars

        ! Ocean Element
        if(in_Sea) call sea_characters_handler(chars)

        ! Vertical Displacements Parameters
        if(in_DispF) call vdisp_characters_handler(chars)
        if(in_disp) call disp_characters_handler(chars)

        ! Flexure Element
        if(in_flex) call flex_characters_handler(chars)

        ! Get Porosity Element
        if(in_Porosity) call porosity_characters_handler(chars)

        ! Hemipelagic Element
        if(in_Hemi) call hemi_characters_handler(chars)

        ! Get UDW parameter
        if(in_UDW) call udw_characters_handler(chars)

    end subroutine characters_handler
    ! ============================================================================
    subroutine SseaElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='seaLevelFile') in_SeaFile = .true.

    end subroutine SseaElement_handler
    ! ============================================================================
    subroutine SvdispElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='nbDispInterval') in_DispNb = .true.

    end subroutine SvdispElement_handler
    ! ============================================================================
    subroutine SflexElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='flexuralInterval') in_flextime = .true.
        if (name=='flexuralRigidity') in_rigidity = .true.
        if (name=='crustDensity') in_crust = .true.

    end subroutine SflexElement_handler
    ! ============================================================================
    subroutine SdispElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='dispFile') in_DispFile = .true.
        if (name=='startDispTime') in_DispST = .true.
        if (name=='endDispTime') in_DispET = .true.

    end subroutine SdispElement_handler
    ! ============================================================================
    subroutine SporosityElement_handler(name)
        character(len=*), intent(in) :: name

        if (name=='depositPorosity') in_PorVal = .true.
        if (name=='effPressureNb') in_EffP = .true.
        if (name=='porosityFile') in_pFile = .true.

    end subroutine SporosityElement_handler
    ! ============================================================================
    subroutine SudwElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='usyncFolder') in_synchF = .true.
        if (name=='usyncTime') in_synchT = .true.

    end subroutine SudwElement_handler
    ! ============================================================================
    subroutine SHemiElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='hemipelagicFile') in_HemiFile = .true.

    end subroutine SHemiElement_handler
    ! ============================================================================
    subroutine EseaElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='seaLevelFluctuations') in_Sea = .false.
        if (name=='seaLevelFile') in_SeaFile = .false.

    end subroutine EseaElement_handler
    ! ============================================================================
    subroutine EHemiElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='HemipelagicRates') in_Hemi = .false.
        if (name=='hemipelagicFile') in_HemiFile = .false.

    end subroutine EHemiElement_handler
    ! ============================================================================
    subroutine EvdispElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='nbDispInterval') in_DispNb = .false.

    end subroutine EvdispElement_handler
    ! ============================================================================
    subroutine EflexElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='flexuralInterval') in_flextime = .false.
        if (name=='flexuralRigidity') in_rigidity = .false.
        if (name=='crustDensity') in_crust = .false.

    end subroutine EflexElement_handler
    ! ============================================================================
    subroutine EdispElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='startDispTime') in_DispST = .false.
        if (name=='endDispTime') in_DispET = .false.
        if (name=='dispFile') in_DispFile = .false.

    end subroutine EdispElement_handler
    ! ============================================================================
    subroutine EporosityElement_handler(name)
        character(len=*), intent(in) :: name

        if (name=='depositPorosity') in_PorVal = .false.
        if (name=='effPressureNb') in_EffP = .false.
        if (name=='porosityFile') in_pFile = .false.

    end subroutine EporosityElement_handler
    ! ============================================================================
    subroutine EudwElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='usyncFolder') in_synchF = .false.
        if (name=='usyncTime') in_synchT = .false.

    end subroutine EudwElement_handler
    ! ============================================================================
    subroutine sea_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if( in_SeaFile)then
            fsea = chars
            gsea%sealevel = .true.
        endif

    end subroutine sea_characters_handler
    ! ============================================================================
    subroutine hemi_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if( in_HemiFile)then
            fhemi = chars
            hemi_flag = .true.
        endif

    end subroutine hemi_characters_handler
    ! ============================================================================
    subroutine vdisp_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if(in_DispNb ) then
            DispNb = chars
            call rts(DispNb, gdisp%event )
            allocate( fdisp( gdisp%event ) )
            allocate( gdisp_time( gdisp%event, 2 ) )
        endif

    end subroutine vdisp_characters_handler
    ! ============================================================================
    subroutine disp_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_DispFile ) fdisp( dispn ) = chars
        if(in_DispST ) then
            DispST = chars
            call rts(DispST, gdisp_time( dispn, 1 ) )
        endif
        if(in_DispET ) then
            DispET = chars
            call rts(DispET, gdisp_time( dispn, 2 ) )
        endif

    end subroutine disp_characters_handler
    ! ============================================================================
    subroutine flex_characters_handler(chars)

        character(len=*), intent(in) :: chars

         if( in_flextime ) then
            flextime = chars
            call rts(flextime, flex_int )
        endif
        if( in_crust ) then
            crust = chars
            call rts(crust, crust_density )
        endif
        if(in_rigidity ) then
            rigidity = chars
            call rts(rigidity, flex_rigidity )
        endif

    end subroutine flex_characters_handler
    ! ============================================================================
    subroutine porosity_characters_handler(chars)
        character(len=*), intent(in) :: chars

        if (in_EffP) then
            gporo%compaction = .true.
            effp = chars
            call rts(effp,gporo%ePnb)
        elseif( in_PorVal )then
            porval = chars
            call rts(porval,initdepporo)
        elseif (in_pFile) then
            fporosity = chars
        endif

    end subroutine porosity_characters_handler
    ! ============================================================================
    subroutine udw_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if(in_synchF ) then
            outdir3 = chars
        endif
        if(in_synchT ) then
            synchT = chars
            call rts(synchT, udw_time )
        endif

    end subroutine udw_characters_handler
    ! ============================================================================
    subroutine endDocument_handler

    end subroutine endDocument_handler
    ! ============================================================================
    !> Subroutine xml_reader_forces()
    !! reads input that defines the SP Model experiment (XmL)
    !<
    ! ============================================================================
    subroutine xml_reader_forces

        type(xml_t) :: xf
        integer :: ios, k
        character(len=128) :: file


        dispn = 0
        gdisp%event = 0
        gsea%sealevel = .false.
        gporo%compaction = .false.
        initdepporo = 0.0_8
        hemi_flag = .false.
        udw_plug = .false.
        flexureon = 0
        flex_int = 1.e5
        flex_rigidity = 1.e22_8
        crust_density = 3300.0_8

        ! Open file
        call open_xml_file(xf, finput , ios)
        if (ios/=0) then
            print*,'---------------------'
            print*, 'Error opening input file for parsing XmL'
            print*,'---------------------'
            stop
        endif

        ! Parse file at first and find data
        call parse(xf, &
            startDocument_handler = startDocument_handler, &
            startElement_handler = startElement_handler, &
            endElement_handler = endElement_handler, &
            characters_handler = characters_handler, &
            endDocument_handler = endDocument_handler &
            )

        ! Close file
        call close_xml_t(xf)

        ! UDW plugin files definition
        if( udw_plug )then
            fudw = 'topsurface.vtk'
            fudisp = 'uw_output.ascii'
            maestro = 'maestro'
            file = ''
            file = fudw
            call addpath3( file )
            fudw = file
            file = ''
            file = maestro
            call addpath3( file )
            maestro = file
            file = ''
            file = fudisp
            call addpath3( file )
            fudisp = file

            ! Assign displacements
            gdisp%event = int( ( time_end - time_start ) / udw_time )
            if( allocated( fdisp ) ) deallocate( fdisp )
            allocate( fdisp( gdisp%event ) )
            if( allocated( gdisp_time ) ) deallocate( gdisp_time )
            allocate( gdisp_time( gdisp%event, 2 ) )

            ! Vertical displacement parameter
            do k = 1, gdisp%event
                fdisp( k ) = fudisp
                gdisp_time( k, 1 ) = time_start + udw_time * ( k - 1 )
                gdisp_time( k, 2 ) = gdisp_time( k, 1 ) +  udw_time
            enddo

        endif

    end subroutine xml_reader_forces
  ! ============================================================================

end module forces
