! ============================================================================
! Name        : Read_Strata.f90
! Author      : Tristan Salles
! Copyright (C) 2014
! ============================================================================
!> \file Read_Strata.f90
!
! Description : Read_Strata is used to gather the information within the XmL input file to build
! the stratigraphic layers used in badlands.
!
!<
! ============================================================================
module xml_strata

    use parallel
    use FoX_sax
    use file_data
    use flow_data
    use mesh_data
    use forces_data
    use FoX_common

    implicit none

    public

    integer :: depF, sedNb, siNb, caNb, orNb, deform, claNb

    ! Output Directory
    logical, save :: in_OutDir = .false.
    logical, save :: in_OutSea = .false.
    logical, save :: in_OutFWs = .false.
    character(len=128), save :: outsea, outfws

    ! Strata Parameters
    logical, save :: stratasection = .false.
    logical, save :: in_Strata = .false.
    logical, save :: in_StrataGrid = .false.
    logical, save :: in_GridXs = .false.
    logical, save :: in_GridYs = .false.
    logical, save :: in_ProcGrid = .false.
    logical, save :: in_Bounds = .false.
    logical, save :: in_Noise = .false.
    character(len=128), save :: Noise, StrataGrid,grid_dxs,grid_dys, procgrid, bounds

    ! Sediment Parameters
    logical, save :: sedparamsection = .false.
    logical, save :: siparamsection = .false.
    logical, save :: caparamsection = .false.
    logical, save :: orparamsection = .false.
    logical, save :: in_Sediments = .false.
    logical, save :: in_Minslp = .false.
    logical, save :: in_Nbsiliciclastics = .false.
    logical, save :: in_Nborganics = .false.
    logical, save :: in_Nbcarbonates = .false.
    logical, save :: in_si = .false.
    logical, save :: in_ca = .false.
    logical, save :: in_or = .false.
    logical, save :: in_Matname = .false.
    logical, save :: in_Diameter = .false.
    logical, save :: in_Density = .false.
    logical, save :: in_Erode = .false.
    logical, save :: in_Time = .false.
    logical, save :: in_CrDa = .false.
    logical, save :: in_CrDm = .false.
    logical, save :: in_Mar = .false.
    logical, save :: in_Aer = .false.
    character(len=128), save :: sinum, minslp, ornum, canum, CrDa, CrDm
    character(len=128), save :: SedDiameter,SedDensity,SedErode, SedTime, Mar, Aer

    ! Deposit Parameters
    logical, save :: depositsection = .false.
    logical, save :: in_Deposit = .false.
    logical, save :: in_DepoFile = .false.
    logical, save :: in_DepNb = .false.
    logical, save :: in_dep = .false.
    character(len=128), save :: Deposit, DepFile, DepNb, depN

    ! Deposit Parameters
    logical, save :: soilsection = .false.
    logical, save :: in_Soil = .false.
    logical, save :: in_ho = .false.
    logical, save :: in_po = .false.
    logical, save :: in_hi = .false.
    character(len=128), save :: ho, po, hi


    ! Check pointing Parameters
    logical, save :: checksection = .false.
    logical, save :: in_Check = .false.
    logical, save :: in_CheckFreq = .false.
    character(len=128), save :: checkfr

    ! Restart Parameters
    logical, save :: restartsection = .false.
    logical, save :: in_RestartFields = .false.
    logical, save :: in_RestartFile = .false.
    logical, save :: in_RestartIt = .false.
    character(len=128), save :: rstit

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

        ! Output Element
        if (name=='outputDir') in_OutDir = .true.
        if (name=='outputSea') in_OutSea = .true.
        if (name=='outputFWs') in_OutFWs = .true.

        ! Strata Element
        if (name=='Strata_grid') in_Strata = .true.
        if (name=='Strata_grid') stratasection = .true.
        if(in_Strata) call SstrataElement_handler(name)

        ! Sediment Element
        if (name=='Materials') in_Sediments = .true.
        if(in_Sediments) sedparamsection = .true.
        if(in_Sediments) call SsedElement_handler(name)
        if (name=='si') in_si = .true.
        if (name=='carb') in_ca = .true.
        if (name=='org') in_or = .true.

        if(in_si) siparamsection = .true.
        if(in_si) call SsiElement_handler(name)

        if(in_ca) caparamsection = .true.
        if(in_ca) call ScaElement_handler(name)

        if(in_or) orparamsection = .true.
        if(in_or) call SorElement_handler(name)

        ! Deposit Element
        if (name=='Init_deposit') in_Deposit = .true.
        if(in_Deposit) depositsection = .true.

        if(in_Deposit) call SdepositElement_handler(name)
        if (name=='dep') in_dep = .true.
        if(in_dep) call SdepElement_handler(name)

        ! Regolith Element
        if (name=='Soil_definition') in_Soil = .true.
        if(in_Soil) soilsection = .true.
        if(in_Soil) call SsoilElement_handler(name)

        ! Restart Parameters
        if (name=='Restart') in_RestartFields = .true.
        if (name=='Restart') restartsection = .true.
        if(in_RestartFields) call SrestartElement_handler(name)

        ! Check pointing Element
        if (name=='CheckPointing') in_Check = .true.
        if(in_Check) checksection = .true.
        if(in_Check) call ScheckElement_handler(name)


    end subroutine startElement_handler
    !============================================================================
    subroutine endElement_handler(namespaceURI, localname, name)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name

        ! Output Element
        if (name=='outputDir') in_OutDir = .false.
        if (name=='outputSea') in_OutSea = .false.
        if (name=='outputFWs') in_OutFWs = .false.

        ! Mesh Element
        call EstrataElement_handler(name)

        ! Sediment Element
        call EsedElement_handler(name)
        call EsiElement_handler(name)
        call EcaElement_handler(name)
        call EorElement_handler(name)

        ! Deposit Element
        call EdepositElement_handler(name)
        call EdepElement_handler(name)

        ! Regolith Element
        call EsoilElement_handler(name)

        ! Restart Parameters
        call ErestartElement_handler(name)

        ! Check pointing Element
        call EcheckElement_handler(name)

    end subroutine endElement_handler
    ! ============================================================================
    subroutine characters_handler(chars)

        character(len=*), intent(in) :: chars

        ! Get Output Dircetory Name
        if (in_OutDir) then
            outdir=''
            outdir = chars
        endif

        ! Get Output Sea flag
        if (in_OutSea) then
            outsea = chars
            call rts(outsea,seavis)
        endif

        ! Get Output FWs flag
        if (in_OutFWs) then
            outfws = chars
            call rts(outfws,fwsvis)
        endif

        ! Strata Element
        if(in_Strata) call strata_characters_handler(chars)

        ! Get Deposit File Name
        if (in_Deposit)   call deposit_characters_handler(chars)
        if (in_dep)   call dep_characters_handler(chars)

        ! Get Sediment Element
        if(in_Sediments) call sed_characters_handler(chars)
        if( in_si .or. in_ca .or. in_or )then
            if( .not. allocated( material_name ) ) allocate( material_name( totgrn ) )
            if( .not. allocated( sediment ) ) allocate(sediment(totgrn))
        endif
        if(in_si .and. silgrn > 0 ) call si_characters_handler(chars)
        if(in_ca .and. carbgrn > 0 ) call ca_characters_handler(chars)
        if(in_or .and. orggrn > 0 ) call or_characters_handler(chars)

        ! Restart Parameters
        if(in_Soil) call soil_characters_handler(chars)

        ! Restart Parameters
        if(in_RestartFields) call restart_characters_handler(chars)

        ! Get Check Pointing Element
        if (in_Check)   call check_characters_handler(chars)

    end subroutine characters_handler
    ! ============================================================================
    subroutine SstrataElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='gridName') in_StrataGrid = .true.
        if (name=='gridX') in_GridXs = .true.
        if (name=='gridY') in_GridYs = .true.
        if (name=='procGrid') in_ProcGrid = .true.
        if (name=='boundConds') in_Bounds = .true.
        if (name=='noiseDist') in_Noise = .true.

    end subroutine SstrataElement_handler
    ! ============================================================================
    subroutine SdepositElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='layersNb') in_DepNb = .true.

    end subroutine SdepositElement_handler
    ! ============================================================================
    subroutine SdepElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='fileDep') in_DepoFile = .true.

    end subroutine SdepElement_handler
    ! ============================================================================
    subroutine SsedElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='nbSiliciclastics') in_Nbsiliciclastics = .true.
        if (name=='nbOrganics') in_Nborganics = .true.
        if (name=='nbCarbonates') in_Nbcarbonates = .true.
        if (name=='minSlope') in_Minslp = .true.

    end subroutine SsedElement_handler
    ! ============================================================================
    subroutine SsiElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='materialName') in_Matname = .true.
        if (name=='diameter') in_Diameter = .true.
        if (name=='density') in_Density = .true.
        if (name=='slopeMarine') in_Mar = .true.
        if (name=='slopeAerial') in_Aer = .true.
        if (name=='erodibility') in_Erode = .true.
        if (name=='transportEff') in_Time = .true.
        if (name=='creepAerial') in_CrDa = .true.
        if (name=='creepMarine') in_CrDm = .true.

    end subroutine SsiElement_handler
    ! ============================================================================
    subroutine ScaElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='materialName') in_Matname = .true.
        if (name=='diameter') in_Diameter = .true.
        if (name=='density') in_Density = .true.
        if (name=='slopeMarine') in_Mar = .true.
        if (name=='slopeAerial') in_Aer = .true.
        if (name=='erodibility') in_Erode = .true.
        if (name=='transportEff') in_Time = .true.
        if (name=='creepAerial') in_CrDa = .true.
        if (name=='creepMarine') in_CrDm = .true.

    end subroutine ScaElement_handler
    ! ============================================================================
    subroutine SorElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='materialName') in_Matname = .true.
        if (name=='diameter') in_Diameter = .true.
        if (name=='density') in_Density = .true.
        if (name=='slopeMarine') in_Mar = .true.
        if (name=='slopeAerial') in_Aer = .true.
        if (name=='erodibility') in_Erode = .true.
        if (name=='transportEff') in_Time = .true.
        if (name=='creepAerial') in_CrDa = .true.
        if (name=='creepMarine') in_CrDm = .true.

    end subroutine SorElement_handler
    ! ============================================================================
    subroutine ScheckElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='frequency') in_CheckFreq = .true.

    end subroutine ScheckElement_handler
    ! ============================================================================
    subroutine SsoilElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='soilProductionho') in_ho = .true.
        if (name=='soilProductionPo') in_po = .true.
        if (name=='soilInitialThickness') in_hi = .true.

    end subroutine SsoilElement_handler
    ! ============================================================================
    subroutine SrestartElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='restartFolder') in_RestartFile = .true.
        if (name=='restartIt') in_RestartIt = .true.

    end subroutine SrestartElement_handler
    ! ============================================================================
    subroutine EstrataElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='Strata_grid') in_Strata = .false.
        if (name=='gridName') in_StrataGrid = .false.
        if (name=='gridX') in_GridXs = .false.
        if (name=='gridY') in_GridYs = .false.
        if (name=='procGrid') in_ProcGrid = .false.
        if (name=='boundConds') in_Bounds = .false.
        if (name=='noiseDist') in_Noise = .false.

    end subroutine EstrataElement_handler
    ! ============================================================================
    subroutine EsedElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='nbSiliciclastics') in_Nbsiliciclastics = .false.
        if (name=='nbOrganics') in_Nborganics = .false.
        if (name=='nbCarbonates') in_Nbcarbonates = .false.
        if (name=='minSlope') in_Minslp = .false.

    end subroutine EsedElement_handler
    ! ============================================================================
    subroutine EdepositElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='layersNb') in_DepNb = .false.

    end subroutine EdepositElement_handler
    ! ============================================================================
    subroutine EdepElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='fileDep') in_DepoFile = .false.

    end subroutine EdepElement_handler
    ! ============================================================================
    subroutine EsiElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='materialName') in_Matname = .false.
        if (name=='diameter') in_Diameter = .false.
        if (name=='density') in_Density = .false.
        if (name=='slopeMarine') in_Mar = .false.
        if (name=='slopeAerial') in_Aer = .false.
        if (name=='erodibility') in_Erode = .false.
        if (name=='creepAerial') in_CrDa = .false.
        if (name=='creepMarine') in_CrDm = .false.
        if (name=='transportEff')then
            in_Time = .false.
            in_si = .false.
        endif

    end subroutine EsiElement_handler
    ! ============================================================================
    subroutine EcaElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='materialName') in_Matname = .false.
        if (name=='diameter') in_Diameter = .false.
        if (name=='density') in_Density = .false.
        if (name=='slopeMarine') in_Mar = .false.
        if (name=='slopeAerial') in_Aer = .false.
        if (name=='erodibility') in_Erode = .false.
        if (name=='creepAerial') in_CrDa = .false.
        if (name=='creepMarine') in_CrDm = .false.
        if (name=='transportEff')then
            in_Time = .false.
            in_ca = .false.
        endif

    end subroutine EcaElement_handler
    ! ============================================================================
    subroutine EorElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='materialName') in_Matname = .false.
        if (name=='diameter') in_Diameter = .false.
        if (name=='density') in_Density = .false.
        if (name=='slopeMarine') in_Mar = .false.
        if (name=='slopeAerial') in_Aer = .false.
        if (name=='erodibility') in_Erode = .false.
        if (name=='creepAerial') in_CrDa = .false.
        if (name=='creepMarine') in_CrDm = .false.
        if (name=='transportEff')then
            in_Time = .false.
            in_or = .false.
        endif

    end subroutine EorElement_handler
    ! ============================================================================
    subroutine EcheckElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='frequency') in_CheckFreq = .false.

    end subroutine EcheckElement_handler
    ! ============================================================================
    subroutine ErestartElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='restartFolder') in_RestartFile = .false.
        if (name=='restartIt') in_RestartIt = .false.

    end subroutine ErestartElement_handler
    ! ============================================================================
    subroutine EsoilElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='soilProductionho') in_ho = .false.
        if (name=='soilProductionPo') in_po = .false.
        if (name=='soilInitialThickness') in_hi = .false.

    end subroutine EsoilElement_handler
    ! ============================================================================
    subroutine strata_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if( in_StrataGrid)then
            fstrata = chars
        elseif (in_GridXs) then
            grid_dxs = chars
            call rts(grid_dxs,lstrat_X)
        elseif (in_GridYs) then
            grid_dys = chars
            call rts(grid_dys,lstrat_Y)
        elseif( in_ProcGrid )then
            procgrid = chars
            call rts(procgrid,gproc)
            if( gproc > nproc ) gproc = nproc
        elseif( in_Bounds )then
            bounds = chars
            call rts(bounds,boundcond(1:4))
        elseif( in_Noise )then
            Noise = chars
            call rts(Noise,noisegen)
        endif

    end subroutine strata_characters_handler
    ! ============================================================================
    subroutine check_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_CheckFreq) then
            checkfr = chars
            checkpointing = .true.
            call rts(checkfr,checkfreq)
        endif

   end subroutine check_characters_handler
   ! ============================================================================
   subroutine restart_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_RestartFile) then
            restartfolder = chars
            restart = .true.
        elseif(in_RestartIt) then
            rstit = chars
            call rts(rstit, restart_iter )
        endif

   end subroutine restart_characters_handler
   ! ============================================================================
   subroutine soil_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_po) then
            po = chars
            call rts(po, P0 )
        elseif(in_ho) then
            ho = chars
            call rts(ho, h0 )
        elseif(in_hi) then
            hi = chars
            call rts(hi, initsoil )
        endif

   end subroutine soil_characters_handler
    ! ============================================================================
    subroutine deposit_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if(in_DepNb)then
            depN = chars
            call rts(depN,InitDep)
            allocate( fdep( InitDep ) )
        endif

    end subroutine deposit_characters_handler
    ! ============================================================================
    subroutine dep_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if(in_DepoFile)then
            depF = depF + 1
            if( depF <= InitDep ) fdep( depF ) = chars
        endif

    end subroutine dep_characters_handler
    !============================================================================
    subroutine sed_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if(in_Nbsiliciclastics)then
            sinum = chars
            call rts(sinum, silgrn )
            totgrn = totgrn + silgrn
        endif
        if(in_Nbcarbonates)then
            canum = chars
            call rts(canum, carbgrn )
            totgrn = totgrn + carbgrn
        endif
        if(in_Nborganics)then
            ornum = chars
            call rts(ornum, orggrn )
            totgrn = totgrn + orggrn
        endif
        if(in_Minslp)then
            minslp = chars
            call rts(minslp, minimum_slp )
        endif

    end subroutine sed_characters_handler
    !============================================================================
    subroutine si_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_Matname) then
            siNb = siNb + 1
            sedNb = siNb
            material_name(sedNb) = chars
        elseif (in_Diameter) then
            SedDiameter = chars
            call rts(SedDiameter,sediment(sedNb)%diameter)
        elseif(in_Density) then
            SedDensity = chars
            call rts(SedDensity,sediment(sedNb)%density)
        elseif(in_Mar)then
            Mar = chars
            call rts(Mar,sediment(sedNb)%slp_marine)
        elseif(in_Aer)then
            Aer = chars
            call rts(Aer,sediment(sedNb)%slp_aerial)
        elseif(in_Erode)then
            SedErode = chars
            call rts(SedErode,sediment(sedNb)%ero_coeff)
        elseif(in_Time)then
            SedTime = chars
            call rts(SedTime,sediment(sedNb)%transport)
        elseif(in_CrDa)then
            CrDa = chars
            call rts(CrDa,sediment(sedNb)%creep_aerial)
        elseif(in_CrDm)then
            CrDm = chars
            call rts(CrDm,sediment(sedNb)%creep_marine)
        endif

    end subroutine si_characters_handler
    !============================================================================
    subroutine ca_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_Matname) then
            caNb = caNb + 1
            sedNb = silgrn + caNb
            material_name(sedNb) = chars
        elseif (in_Diameter) then
            SedDiameter = chars
            call rts(SedDiameter,sediment(sedNb)%diameter)
        elseif(in_Density) then
            SedDensity = chars
            call rts(SedDensity,sediment(sedNb)%density)
        elseif(in_Mar)then
            Mar = chars
            call rts(Mar,sediment(sedNb)%slp_marine)
        elseif(in_Aer)then
            Aer = chars
            call rts(Aer,sediment(sedNb)%slp_aerial)
        elseif(in_Erode)then
            SedErode = chars
            call rts(SedErode,sediment(sedNb)%ero_coeff)
        elseif(in_Time)then
            SedTime = chars
            call rts(SedTime,sediment(sedNb)%transport)
        elseif(in_CrDa)then
            CrDa = chars
            call rts(CrDa,sediment(sedNb)%creep_aerial)
        elseif(in_CrDm)then
            CrDm = chars
            call rts(CrDm,sediment(sedNb)%creep_marine)
        endif

    end subroutine ca_characters_handler
    !============================================================================
    subroutine or_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if (in_Matname) then
            orNb = orNb + 1
            sedNb = silgrn + carbgrn + orNb
            material_name(sedNb) = chars
        elseif (in_Diameter) then
            SedDiameter = chars
            call rts(SedDiameter,sediment(sedNb)%diameter)
        elseif(in_Density) then
            SedDensity = chars
            call rts(SedDensity,sediment(sedNb)%density)
        elseif(in_Mar)then
            Mar = chars
            call rts(Mar,sediment(sedNb)%slp_marine)
        elseif(in_Aer)then
            Aer = chars
            call rts(Aer,sediment(sedNb)%slp_aerial)
        elseif(in_Erode)then
            SedErode = chars
            call rts(SedErode,sediment(sedNb)%ero_coeff)
        elseif(in_Time)then
            SedTime = chars
            call rts(SedTime,sediment(sedNb)%transport)
        elseif(in_CrDa)then
            CrDa = chars
            call rts(CrDa,sediment(sedNb)%creep_aerial)
        elseif(in_CrDm)then
            CrDm = chars
            call rts(CrDm,sediment(sedNb)%creep_marine)
        endif

    end subroutine or_characters_handler
    ! ============================================================================
    subroutine endDocument_handler

    end subroutine endDocument_handler
    ! ============================================================================
    !> Subroutine xml_reader_strata()
    !! reads input that defines the experiment (XmL)
    !<
    ! ============================================================================
    subroutine xml_reader_strata

        type(xml_t) :: xf
        integer :: ios, ks
!        real( tkind ) :: maxtrans, mintrans, maxdia

        depF = 0
        sedNb = 0
        siNb = 0
        caNb = 0
        orNb = 0
        claNb = 0
        InitDep = 0
        silgrn = 0
        orggrn = 0
        carbgrn = 0
        totgrn = 0
        gproc = 1
        seavis = 0
        fwsvis = 0
        noisegen = 0.0_8
        boundcond = 0
        restart = .false.
        checkpointing = .false.
        h0 = 0.5_8
        P0 = 50.0_8 *1.0e-6_8
        initsoil = 0.0_8

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

        ! Fall velocities computation based on Zhang
        do ks = 1, totgrn
            sediment( ks )%diameter = sediment( ks )%diameter * 0.001_8
            if( sediment(ks)%diameter /= 0.0_8 )then
                if( sediment(ks)%diameter >= 0.062 * 0.001_8 )then
                    sediment(ks)%vfall = settling_velocity_zhang(ks)
                else
                    sediment(ks)%vfall = 0.5_8 * 0.001_8
                endif
            else
                sediment(ks)%vfall = 0.0_8
            endif
        enddo

        ! Check pointing processor
        checkproc = nproc - 1
        if( gproc < nproc ) checkproc = gproc + 1

        ! Get expected transport time
!        maxtrans = 0.0_8
!        mintrans = 2.0_8
!        maxdia = 0.0_8
!        do ks = 1, totgrn
!            maxtrans = max( maxtrans, sediment(ks)%transport )
!            mintrans = min( mintrans, sediment(ks)%transport )
!            maxdia = max( maxdia, sediment(ks)%diameter )
!        enddo
!
!        do ks = 1, totgrn
!            sediment(ks)%transport = mintrans + maxtrans * ( 1 - &
!                ( sediment(sedNb)%diameter )**( 1.0_8 / maxdia ) &
!                / ( maxdia )**( 1.0_8 / maxdia ) )
!        enddo


    end subroutine xml_reader_strata
    ! ============================================================================
    !> Function settling_velocity_zhang()
    !! Function settling_velocity_zhang is used to compute settling velocity for a considered type of grain
    !! Based on Zhang 1993
    !! \param sed
    !<
    ! ============================================================================
    function settling_velocity_zhang( sed ) result( ws )

      integer, intent( in ) :: sed

      real( tkind ) :: Cvisc, delg, ws

      Cvisc = 1.0e-6_8 * 13.95_8 / sediment( sed )%diameter
      delg=( sediment( sed )%density - fluid_density ) * gravity / fluid_density

      ws = sqrt( Cvisc**2 + 1.09_8 * delg * sediment( sed )%diameter )
      ws = ws - Cvisc

    end function settling_velocity_zhang
    ! ============================================================================

end module xml_strata
