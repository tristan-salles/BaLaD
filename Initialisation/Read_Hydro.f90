! ============================================================================
! Name        : Read_Hydro.f90
! Author      : Tristan Salles
! Copyright (C) 2014
! ============================================================================
!> \file Read_Hydro.f90
!
! Description : Read_Hydro is used to gather the information within the XmL input file to
! get the processes parameters used in the run.
!
!<
! ============================================================================
module xml_hydro

    use parallel
    use FoX_sax
    use file_data
    use flow_data
    use mesh_data
    use FoX_common

    implicit none

    public

    integer :: srcNb, rainn

    ! Erosion Deposition Parameters
    logical, save :: in_EDParam = .false.
    logical, save :: edparamsection = .false.
    logical, save :: in_Manning = .false.
    logical, save :: in_ManningHypo = .false.
    logical, save :: in_ManningHyper = .false.
    logical, save :: in_DTBaa = .false.
    logical, save :: in_DTBm = .false.
    logical, save :: in_TransM = .false.
    logical, save :: in_mTa = .false.
    logical, save :: in_nTa = .false.
    logical, save :: in_mDe = .false.
    logical, save :: in_nDe = .false.
    logical, save :: in_eTr = .false.
    character(len=128), save :: Manni,ManniHypo,ManniHyper,DTBaa,DTBm, TransM, nDe, mDe, nTa, mTa, eTr

    ! Hydraulic Parameters
    logical, save :: in_HYDParam = .false.
    logical, save :: hydparamsection = .false.
    logical, save :: in_StepFill = .false.
    logical, save :: in_kwd = .false.
    logical, save :: in_awd = .false.
    logical, save :: in_bwd = .false.
    character(len=128), save :: StepFill,kwd,awd,bwd

    ! Mass Wasting Parameters
    logical, save :: in_MWParam = .false.
    logical, save :: mwparamsection = .false.
    logical, save :: in_mon = .false.
    logical, save :: in_cma = .false.
    logical, save :: in_cmm = .false.
    logical, save :: in_faca = .false.
    logical, save :: in_facm = .false.
    logical, save :: in_ctsla = .false.
    logical, save :: in_ctslm = .false.
    character(len=128), save :: mon,cma,cmm,faca,facm, ctsla, ctslm

    ! Rain Parameters Grid
    logical, save :: rainmapsection = .false.
    logical, save :: in_RainGrid = .false.
    logical, save :: in_RainNb = .false.
    logical, save :: in_RainMFile = .false.
    logical, save :: in_RainET = .false.
    logical, save :: in_RainST = .false.
    logical, save :: in_RainTM = .false.
    logical, save :: in_rain = .false.
    character(len=128), save :: RainTM,RainGrid, RainNb,RainMFile,RainET,RainST

    ! Sources Parameters
    logical, save :: sourcesparamsection = .false.
    logical, save :: sourcesection = .false.
    logical, save :: in_Sources = .false.
    logical, save :: in_SourcesNb = .false.
    logical, save :: in_FWNb = .false.
    logical, save :: in_source = .false.
    logical, save :: in_t1 = .false.
    logical, save :: in_t2 = .false.
    logical, save :: in_x = .false.
    logical, save :: in_y = .false.
    logical, save :: in_xrange = .false.
    logical, save :: in_yrange = .false.
    logical, save :: in_Vx = .false.
    logical, save :: in_Vy = .false.
    logical, save :: in_Q = .false.
    logical, save :: in_Qs = .false.
    logical, save :: in_sedConcentration = .false.
    logical, save :: in_flowType = .false.
    character(len=128), save :: srcnum, srcel, srct1, srct2, srcx, srcy, srcxr, srcyr, srcvx, srcvy
    character(len=128), save :: srcQ, srcQs, srcsC, srcT, fwNb

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

        ! Rain Parameters Grid
        if (name=='RainfallGrid') in_RainGrid = .true.
        if (in_RainGrid) rainmapsection = .true.
        if(in_RainGrid) call SrainmapElement_handler(name)
        if (name=='rain')then
            in_rain = .true.
            rainn = rainn + 1
        endif
        if(in_rain) call SrainfieldElement_handler(name)

        ! Hydraulic Parameters Element
        if (name=='HydraulicsParameters') in_HYDParam = .true.
        if(in_HYDParam) hydparamsection = .true.
        if(in_HYDParam) call ShydroElement_handler(name)

        ! Mass Wasting Parameters Element
        if (name=='MassWastingParameters') in_MWParam = .true.
        if(in_MWParam) mwparamsection = .true.
        if(in_MWParam) call SmwElement_handler(name)

        ! Sources Element
        if (name=='Sources') in_Sources = .true.
        if(in_Sources) sourcesparamsection = .true.
        if(in_Sources) call SsourcesElement_handler(name)
        if (name=='sourceVal') in_source = .true.
        if(in_source) sourcesection = .true.
        if(in_source) call SsourceElement_handler(name)

        ! Erosion Deposition Parameters Element
        if (name=='EroDepParameters') in_EDParam = .true.
        if(in_EDParam) edparamsection = .true.
        if(in_EDParam) call SerodepElement_handler(name)

    end subroutine startElement_handler
    ! ============================================================================
    subroutine endElement_handler(namespaceURI, localname, name)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name

        ! Rain Parameters Grid
        call ErainmapElement_handler(name)
        call ErainfieldElement_handler(name)

        ! Hydraulic Parameters Element
        call EhydroElement_handler(name)

        ! Mass Wasting Parameters Element
        call EmwElement_handler(name)

        ! Sources Element
        call EsourcesElement_handler(name)
        call EsourceElement_handler(name)

        ! Erosion Deposition Parameters Element
        call EerodepElement_handler(name)

    end subroutine endElement_handler
    ! ============================================================================
    subroutine characters_handler(chars)

        character(len=*), intent(in) :: chars

        ! Rain Parameters Grid
        if(in_RainGrid) call rainmap_characters_handler(chars)
        if(in_rain) call rainfield_characters_handler(chars)

        ! Hydraulic Element
        if(in_HYDParam) call hydro_characters_handler(chars)

        ! Mass Wasting Element
        if(in_MWParam) call mw_characters_handler(chars)

        ! Get Sources Element
        if(in_Sources)  call sources_characters_handler(chars)
        if(in_source)   call source_characters_handler(chars)

        ! Erosion Deposition Element
        if(in_EDParam) call erodepo_characters_handler(chars)

    end subroutine characters_handler
    ! ============================================================================
    subroutine SsourcesElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='realeaseInterval') in_FWNb = .true.
        if (name=='totalSourcesNb') in_SourcesNb = .true.

    end subroutine SsourcesElement_handler
    ! ============================================================================
    subroutine SsourceElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='ts') in_t1 = .true.
        if (name=='te') in_t2 = .true.
        if (name=='x') in_x = .true.
        if (name=='y') in_y = .true.
        if (name=='xRange') in_xrange = .true.
        if (name=='yRange') in_yrange = .true.
        if (name=='Vx') in_Vx = .true.
        if (name=='Vy') in_Vy = .true.
        if (name=='Q') in_Q = .true.
        if (name=='Qs') in_Qs = .true.
        if(name=='sedConcentration') in_sedConcentration = .true.
        if (name=='flowType') in_flowType = .true.

    end subroutine SsourceElement_handler
    ! ============================================================================
    subroutine SrainmapElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='nbRainfallInterval') in_RainNb = .true.

    end subroutine SrainmapElement_handler
    ! ============================================================================
    subroutine SrainfieldElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='rainFile') in_RainMFile = .true.
        if (name=='startTime') in_RainST = .true.
        if (name=='endTime') in_RainET = .true.
        if (name=='rainEvent') in_RainTM = .true.

    end subroutine SrainfieldElement_handler
    ! ============================================================================
    subroutine ShydroElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='stepFill') in_StepFill = .true.
        if (name=='kwidth') in_kwd = .true.
        if (name=='awidth') in_awd = .true.
        if (name=='bwidth') in_bwd = .true.

    end subroutine ShydroElement_handler
    ! ============================================================================
    subroutine SmwElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='massMovementOn') in_mon = .true.
        if (name=='maxConcAerial') in_cma = .true.
        if (name=='maxConcMarine') in_cmm = .true.
        if (name=='factorAerial') in_faca = .true.
        if (name=='factorMarine') in_facm = .true.
        if (name=='criticalSlopeAerial') in_ctsla = .true.
        if (name=='criticalSlopeMarine') in_ctslm = .true.

    end subroutine SmwElement_handler
    ! ============================================================================
    subroutine SerodepElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='transportMode') in_TransM = .true.
        if (name=='ntransport') in_nTa = .true.
        if (name=='mtransport') in_mTa = .true.
        if (name=='ndetach') in_nDe = .true.
        if (name=='mdetach') in_mDe = .true.
        if (name=='fwErosionTrigger') in_eTr = .true.
        if (name=='manningOpenChannel') in_Manning = .true.
        if (name=='manningHyperpycnal') in_ManningHyper = .true.
        if (name=='manningHypopycnal') in_ManningHypo = .true.
        if (name=='maxSoilDepthAerial') in_DTBaa = .true.
        if (name=='maxSoilDepthMarine') in_DTBm = .true.

    end subroutine SerodepElement_handler
    ! ============================================================================
    subroutine EsourcesElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='realeaseInterval') in_FWNb = .false.
        if (name=='totalSourcesNb') in_SourcesNb = .false.

    end subroutine EsourcesElement_handler
    ! ============================================================================
    subroutine EsourceElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='ts') in_t1 = .false.
        if (name=='te') in_t2 = .false.
        if (name=='x') in_x = .false.
        if (name=='y') in_y = .false.
        if (name=='xRange') in_xrange = .false.
        if (name=='yRange') in_yrange = .false.
        if (name=='Vx') in_Vx = .false.
        if (name=='Vy') in_Vy = .false.
        if (name=='Q') in_Q = .false.
        if (name=='Qs') in_Qs = .false.
        if(name=='sedConcentration') in_sedConcentration = .false.
        if (name=='flowType') in_flowType = .false.

    end subroutine EsourceElement_handler
    ! ============================================================================
    subroutine ErainmapElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='nbRainfallInterval') in_RainNb = .false.

    end subroutine ErainmapElement_handler
    ! ============================================================================
    subroutine ErainfieldElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='startTime') in_RainST = .false.
        if (name=='endTime') in_RainET = .false.
        if (name=='rainFile') in_RainMFile = .false.
        if (name=='rainEvent') in_RainTM = .false.

    end subroutine ErainfieldElement_handler
    ! ============================================================================
    subroutine EhydroElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='stepFill') in_StepFill = .false.
        if (name=='kwidth') in_kwd = .false.
        if (name=='awidth') in_awd = .false.
        if (name=='bwidth') in_bwd = .false.

    end subroutine EhydroElement_handler
    ! ============================================================================
    subroutine EmwElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='massMovementOn') in_mon = .false.
        if (name=='maxConcAerial') in_cma = .false.
        if (name=='maxConcMarine') in_cmm = .false.
        if (name=='factorAerial') in_faca = .false.
        if (name=='factorMarine') in_facm = .false.
        if (name=='criticalSlopeAerial') in_ctsla = .false.
        if (name=='criticalSlopeMarine') in_ctslm = .false.

    end subroutine EmwElement_handler
    ! ============================================================================
    subroutine EerodepElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='transportMode') in_TransM = .false.
        if (name=='ntransport') in_nTa = .false.
        if (name=='mtransport') in_mTa = .false.
        if (name=='ndetach') in_nDe = .false.
        if (name=='mdetach') in_mDe = .false.
        if (name=='fwErosionTrigger') in_eTr = .false.
        if (name=='manningOpenChannel') in_Manning = .false.
        if (name=='manningHyperpycnal') in_ManningHyper = .false.
        if (name=='manningHypopycnal') in_ManningHypo = .false.
        if (name=='maxSoilDepthAerial') in_DTBaa = .false.
        if (name=='maxSoilDepthMarine') in_DTBm = .false.

    end subroutine EerodepElement_handler
    ! ============================================================================
    subroutine sources_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_FWNb) then
            fwNb = chars
            call rts(fwNb, nbfw_src )
        elseif (in_SourcesNb) then
            srcnum = chars
            call rts(srcnum, num_src )
            allocate( fws_num(num_src) )
            allocate( fws_type(num_src) )
            allocate( fws_tstrt(num_src) )
            allocate( fws_tend(num_src) )
            allocate( fws_xposition(num_src) )
            allocate( fws_yposition(num_src) )
            allocate( fws_xrange(num_src) )
            allocate( fws_yrange(num_src) )
            allocate( fws_xvel(num_src) )
            allocate( fws_yvel(num_src) )
            allocate( fws_qfl(num_src) )
            allocate( fws_volume(num_src) )
            allocate( fws_sedconc(num_src) )
            allocate( fws_sedperc(num_src, silgrn ) )
            allocate( fws_sedcharge(num_src, silgrn ) )
            fws_sedperc=0.0_8
            fws_sedcharge=0.0_8
            fws_xrange=0.0_8
            fws_yrange=0.0_8
        endif

    end subroutine sources_characters_handler
    ! ============================================================================
    subroutine source_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_t1) then
            srcNb = srcNb + 1
            srct1 = chars
            call rts(srct1,fws_tstrt(srcNb))
        elseif(in_t2) then
            srct2 = chars
            call rts(srct2,fws_tend(srcNb))
        elseif(in_x)then
            srcx = chars
            call rts(srcx,fws_xposition(srcNb))
            fws_xposition(srcNb) = fws_xposition(srcNb)
        elseif(in_y) then
            srcy = chars
            call rts(srcy,fws_yposition(srcNb))
            fws_yposition(srcNb) = fws_yposition(srcNb)
        elseif(in_xrange)then
            srcxr = chars
            call rts(srcxr,fws_xrange(srcNb))
            fws_xrange(srcNb) = fws_xrange(srcNb) * 0.5_8
        elseif(in_yrange) then
            srcyr = chars
            call rts(srcyr,fws_yrange(srcNb))
            fws_yrange(srcNb) = fws_yrange(srcNb) * 0.5_8
        elseif(in_Vx)then
            srcvx = chars
            call rts(srcvx,fws_xvel(srcNb))
        elseif(in_Vy) then
            srcvy = chars
            call rts(srcvy,fws_yvel(srcNb))
        elseif(in_Q) then
            srcQ = chars
            call rts(srcQ,fws_qfl(srcNb))
        elseif(in_Qs)then
            srcQs = chars
            call rts(srcQs,fws_sedconc(srcNb))
        elseif(in_sedConcentration) then
            srcsC = chars
            call rts(srcsC,fws_sedperc(srcNb,1:silgrn))
        elseif(in_flowType)then
            srcT = chars
            call rts(srcT,fws_type(srcNb))
        endif

    end subroutine source_characters_handler
    ! ============================================================================
    subroutine rainmap_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_RainNb) then
            RainNb= chars
            call rts(RainNb, rain_event )
            allocate( frainmap( rain_event ) )
            allocate( rain_tend( rain_event ) )
            allocate( rain_tstart( rain_event ) )
            allocate( rain_duration( rain_event ) )
        endif

    end subroutine rainmap_characters_handler
    ! ============================================================================
    subroutine rainfield_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_RainMFile) then
            frainmap( rainn ) = chars
        endif
        if(in_RainST ) then
            RainST = chars
            call rts(RainST, rain_tstart( rainn ) )
        endif
        if(in_RainET ) then
            RainET = chars
            call rts(RainET, rain_tend( rainn ) )
        endif
        if(in_RainTM ) then
            RainTM = chars
            call rts(RainTM, rain_duration( rainn ) )
        endif

    end subroutine rainfield_characters_handler
    ! ============================================================================
    subroutine hydro_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_StepFill) then
            StepFill = chars
            call rts(StepFill, step_fill )
        endif
        if(in_kwd ) then
            kwd = chars
            call rts(kwd, kwidth )
        endif
        if(in_awd ) then
            awd = chars
            call rts(awd, awidth )
        endif
        if(in_bwd ) then
            bwd = chars
            call rts(bwd, bwidth )
        endif

    end subroutine hydro_characters_handler
    ! ============================================================================
    subroutine mw_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_mon) then
            mon = chars
            call rts(mon, masson )
        endif
        if (in_cma) then
            cma = chars
            call rts(cma, debris_max_conc_aerial )
        endif
        if(in_cmm ) then
            cmm = chars
            call rts(cmm, debris_max_conc_marine )
        endif
        if(in_facm ) then
            facm = chars
            call rts(facm, debris_factor_marine )
        endif
        if(in_faca ) then
            faca = chars
            call rts(faca, debris_factor_aerial )
        endif
        if(in_ctsla ) then
            ctsla = chars
            call rts(ctsla, critical_slope_aerial )
        endif
        if(in_ctslm ) then
            ctslm = chars
            call rts(ctslm, critical_slope_marine )
        endif

    end subroutine mw_characters_handler
    ! ============================================================================
    subroutine erodepo_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_TransM) then
            TransM = chars
            call rts( TransM, transport_mode )
        endif
        if (in_mTa) then
            mTa = chars
            call rts( mTa, mtransport )
        endif
        if (in_nTa) then
            nTa = chars
            call rts( nTa, ntransport )
        endif
        if (in_eTr) then
            eTr = chars
            call rts( eTr, erosion_trigger )
        endif
        if (in_nDe) then
            nDe = chars
            call rts( nDe, ndetach )
        endif
        if (in_mDe) then
            mDe = chars
            call rts( mDe, mdetach )
        endif
        if (in_Manning) then
            Manni = chars
            call rts(Manni, manning_open )
        endif
        if (in_ManningHypo) then
            ManniHypo = chars
            call rts(ManniHypo, manning_hypo )
        endif
        if (in_ManningHyper) then
            ManniHyper = chars
            call rts(ManniHyper, manning_hyper )
        endif
        if( in_DTBaa ) then
            DTBaa = chars
            call rts(DTBaa, dtb_aerial )
        endif
        if( in_DTBm ) then
            DTBm = chars
            call rts(DTBm, dtb_marine )
        endif

    end subroutine erodepo_characters_handler
    ! ============================================================================
    subroutine endDocument_handler

    end subroutine endDocument_handler
    ! ============================================================================
    !> Subroutine xml_reader_hydraulic()
    !! reads input that defines the experiment (XmL)
    !<
    ! ============================================================================
    subroutine xml_reader_hydraulic

        type(xml_t) :: xf
        integer :: ios

        srcNb = 0
        nbfw_src = 1
        num_src = 0
        kwidth = 5._8
        awidth = 0.5_8
        bwidth = 0.0_8
        manning_open = 0.02_8
        manning_hyper = 0.01_8
        manning_hypo = 0.08_8
        dtb_marine = 2.5_8
        dtb_aerial = 5.0_8
        step_fill = 2.0_8
        rainID = 0
        masson = 0
        rainn = 0
        transport_mode = 0
        critical_slope_aerial = 90
        critical_slope_marine = 90
        mtransport = 1.5_8
        ntransport = 1.0_8
        erosion_trigger = 0.01_8
        mdetach = 1.0_8 / 3.0_8
        ndetach = 2.0_8 / 3.0_8

        ! Open file
        call open_xml_file(xf, finput , ios)
        if( ios/=0 )then
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

        critical_slope_aerial = tan( critical_slope_aerial * pi / 180.0_8 )
        critical_slope_marine = tan( critical_slope_marine * pi / 180.0_8 )

    end subroutine xml_reader_hydraulic
    ! ============================================================================

end module xml_hydro
