! ============================================================================
! Name        : Init_model.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Init_model.f90
!
! Description :  This module performs initialises the parameters used in the experiment based on input file.
!
!<
! ============================================================================
!> Module init_phase
!<
module init_phase

    use forces
    use parallel
    use file_data
    use xml_time
    use flow_data
    use time_data
    use xml_hydro
    use xml_strata
    use checkpoint
    use compaction
    use forces_data
    use hemipelagite

    implicit none

contains

    ! ============================================================================
    !> Subroutine initialisation_phase
    !! This subroutine performs
    !<
    ! ============================================================================
    subroutine initialisation_phase

        integer :: iunit, i, p, ks, pos

        real( tkind ) :: ratcheck, dens, flowtime

        character(len=128) :: stg

        ! Read simulation parameters from XmL
        call simulation_initialisation

        ! Delete flatten XmL file from the local directory
        if( iam == 0 )then
            iunit = 423
            open( iunit, file=finput, form='unformatted' )
            close( iunit, status='delete' )
        endif

        ! Read the basement file
        call read_basement

        if( .not. restart )then

            ! Allocate depositional arrays
            nbLays = InitDep + int( ( time_end - time_start ) / time_layer )
            if( allocated( stratLays ) ) deallocate( stratLays )
            allocate( stratLays( nbLays, lnbPts, totgrn ) )
            if( gporo%compaction )then
                if( allocated( porosityLays ) ) deallocate( porosityLays )
                allocate( porosityLays( nbLays, lnbPts, totgrn ) )
                porosityLays = 0.0_8
            endif
            if( allocated( v_nbLays ) ) deallocate( v_nbLays )
            allocate( v_nbLays( lnbPts ) )
            if( allocated( v_LaysID ) ) deallocate( v_LaysID )
            allocate( v_LaysID( lnbPts, nbLays ) )
            stratLays = 0.0_8
            v_LaysID = 0
            v_nbLays = 0

            ! Read initial deposit file
            if( InitDep > 0 ) call read_initial_deposit
            layID = InitDep

        else

            ! Allocate depositional arrays
            nbLays = restart_iter + int( ( time_end - time_start ) / time_layer )
            if( allocated( stratLays ) ) deallocate( stratLays )
            allocate( stratLays( nbLays, lnbPts, totgrn ) )
            if( gporo%compaction )then
                if( allocated( porosityLays ) ) deallocate( porosityLays )
                allocate( porosityLays( nbLays, lnbPts, totgrn ) )
                porosityLays = 0.0_8
            endif
            if( allocated( v_nbLays ) ) deallocate( v_nbLays )
            allocate( v_nbLays( lnbPts ) )
            if( allocated( v_LaysID ) ) deallocate( v_LaysID )
            allocate( v_LaysID( lnbPts, nbLays ) )
            stratLays = 0.0_8
            v_LaysID = 0
            v_nbLays = 0

            ! Read previous deposit data
            call read_checkpoint_simulation
            layID = restart_iter

        endif

        ! Create the DEM coordinates
        call define_DEM_grid

        ! Define band partitioning
        allocate( part_rows( nproc ) )
        call define_band_partitioning

        ! Setup simulation time
        tnow = time_start
        tsampling = tnow

        ! Define time tolerance accuracy
        i = 0
        if( time_start /= 0.0_8 ) i = int( log10( abs( time_start ) ) )
        if( time_end /= 0.0_8 ) i = max( i , int( log10( abs( time_end ) ) ) )
        i = i - 13
        time_tolerance = max( tor, 10.0_8**i )

        ! Setup computation file names
        call get_cpp_file_name_extensions

        ! Setup diffusion transport grid
        allocate( difo( lnbPts ) )
        allocate( difp( lnbPts ) )
        allocate( top_sedh( lnbPts, totgrn ) )
        allocate( top_sedprev( lnbPts, totgrn ) )
        top_sedh = 0.0_8
        top_sedprev = 0.0_8

        allocate( creepz( lnbPts ) )
        allocate( depo( lnbPts, totgrn ) )
        allocate( dstart( lnbPts, totgrn ) )
        allocate( cdif( lnbPts ) )

        ! Sea fluctuations
        if( gsea%sealevel ) call read_sealevel_file

        ! Allocate displacement fields
        gdisp%disptime = 0.0_8
        if( gdisp%event > 0 .and. .not. udw_plug ) then
             new_disp = .true.
            ! Read and assign vertical displacement rate
            call assign_vertical_displacements_rate( 1 )
            ! Find displacements for considered time
            call update_local_displacements
        endif

        if( num_src == 0 ) goto 48

        ! Sources parameters
        do p = 1, num_src
            ratcheck = 0.0_8
            do ks = 1, silgrn
                ratcheck = ratcheck + fws_sedperc( p , ks )
            enddo
            if( ratcheck /= 100.0_8 )then
                print*,'The source sediment percentages do not sum up to 100.',p
                stop
            endif

            if( fws_tstrt( p ) > fws_tend( p ) )then
                print*,'The source time declaration is wrong.',p
                stop
            endif

            pos = int( ( fws_xposition( p ) - strat_xo ) / strat_dx + 1 )
            if( pos < 1 .or. pos > strat_X )then
                print*,'The source position declaration is wrong.',p
                stop
            endif
            pos = int( ( fws_yposition( p ) - strat_yo ) / strat_dx + 1 )
            if( pos < 1 .or. pos > strat_Y )then
                print*,'The source position declaration is wrong2.',p
                stop
            endif

            if( fws_tend( p ) - fws_tstrt( p ) < flow_int ) then
                flowtime = fws_tend( p ) - fws_tstrt( p )
            else
                flowtime = flow_int
            endif

            ! Transform source concentration from Mt/year to cubic metres per second
            dens = 0.0_8
            do ks = 1, silgrn
                dens = dens + fws_sedperc( p, ks ) * sediment( ks )%density * 0.01_8
            enddo
            fws_sedconc( p ) = fws_sedconc( p ) * 1.e9_8 / ( dens * secyear )

            ! Volume of each sediment transported (cubic metres)
            do ks = 1, silgrn
                fws_sedcharge( p, ks ) = fws_sedperc( p, ks ) * 0.01_8 * fws_sedconc( p ) &
                     * secyear * flowtime / nbfw_src
            enddo

            ! Volume of water transported (cubic metres)
            fws_volume( p ) = fws_qfl( p ) * secyear * flowtime / nbfw_src

        enddo

48 continue

        ! Allocate and initialise porosity table
        if( gporo%compaction )then
            allocate( effPressure( gporo%ePnb ) )
            allocate( porosity( totgrn, gporo%ePnb ) )
            call assign_porosity_table
            if( .not. restart ) call porosity_init
        endif

        ! Hemipelagic section
        hemi_mat = 0
        if( hemi_flag )then
            last_hemi = tnow
            do i = 1, totgrn
                stg = material_name( i )
                if( stg( 1:4 ) == 'Hemi' .or. stg( 1:4 ) == 'hemi' ) hemi_mat = i
            enddo
            if( hemi_mat == 0 )then
                 print*,'WARNING: To use hemipelagic sedimentation it is mandatory to named Hemi one of the materials'
                 stop
            endif
            ! Read hemipelagic file
            call read_hemipelagite_file
        endif

        return

    end subroutine initialisation_phase
    ! ============================================================================
    !> Subroutine simulation_initialisation
    !! Allocates all parameters that defines an experiment
    !<
    ! ============================================================================
    subroutine simulation_initialisation

        ! Parameters Declaration
        logical :: found
        integer :: l1, l2

        character(len=128) :: command, stg, fildir

        gsea%actual_sea = 0.0_8

        ! Default output directory
        outdir=''
        outdir(1:7)='outputs'
        outputs=''
        runfiles=''
        outputs(1:7)='outputs'
        runfiles(1:8)='runfiles'
        ftin_out = 'GTIN'

        ! Parsing input file
        call xml_reader_time
        call xml_reader_strata
        call xml_reader_hydraulic
        call xml_reader_forces

        ! Create the output directories
        if( iam == 0 )then
            call noblnk(outdir)
            fildir = outdir
            if( rain_event == 0 )then
                stg = '/outputs/Stratal_series.xdmf'
            else
                stg = '/outputs/HMesh_series.xdmf'
            endif
            call noblnk(stg)
            call append_str( fildir,stg )
            call noblnk(fildir)
            inquire(file=fildir,exist=found)
            if( .not. found )then
                command=' '
                command(1:6)='mkdir '
                l1 = len_trim(outdir)
                command(7:l1+7)=outdir
                call term_command( command )
                command(l1+7:l1+7)='/'
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+15)=outputs(1:7)
                call term_command( stg )
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+16)=runfiles(1:8)
                call term_command( stg )
                outdir1(1:l1)=outdir(1:l1)
                outdir1(l1+1:l1+1)='/'
                outdir1(l1+2:l1+9)=outputs(1:7)
                outdir2(1:l1)=outdir(1:l1)
                outdir2(l1+1:l1+1)='/'
                outdir2(l1+2:l1+10)=runfiles(1:8)
            else
                command=' '
                command(1:6) = 'rm -r '
                l1 = len_trim(outdir)
                command(7:l1+7)=outdir
                call term_command( command )
                command=' '
                command(1:6)='mkdir '
                l1 = len_trim(outdir)
                command(7:l1+7)=outdir
                call term_command( command )
                command(l1+7:l1+7)='/'
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+15)=outputs(1:7)
                call term_command( stg )
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+16)=runfiles(1:8)
                call term_command( stg )
                outdir1(1:l1)=outdir(1:l1)
                outdir1(l1+1:l1+1)='/'
                outdir1(l1+2:l1+9)=outputs(1:7)
                outdir2(1:l1)=outdir(1:l1)
                outdir2(l1+1:l1+1)='/'
                outdir2(l1+2:l1+10)=runfiles(1:8)
            endif

            ! Underworld Lecode shared folder
            if( udw_plug )then
                call noblnk(outdir3)
                fildir = outdir3
                stg = '/topsurface.vtk'
                call noblnk(stg)
                call append_str( fildir,stg )
                call noblnk(fildir)
                inquire(file=fildir,exist=found)
                if( .not. found )then
                    command=' '
                    command(1:6)='mkdir '
                    l1 = len_trim(outdir3)
                    command(7:l1+7)=outdir3
                    call term_command( command )
                else
                    command=' '
                    command(1:6) = 'rm -r '
                    l1 = len_trim(outdir3)
                    command(7:l1+7)=outdir3
                    call term_command( command )
                    command=' '
                    command(1:6)='mkdir '
                    l1 = len_trim(outdir3)
                    command(7:l1+7)=outdir3
                    call term_command( command )
                endif
            endif

            ! Copy XML file in the output folder
            command = ''
            command(1:3) = 'cp '
            l1 = len_trim(finput)
            command(4:l1+4) = finput(1:l1)
            command(l1+4:l1+4)= ' '
            l2 = len_trim(outdir)
            command(l1+5:l1+5+l2) = outdir(1:l2)
            command(l1+5+l2:l1+5+l2) = '/'
            call term_command( command )
        endif

        ! Broadcast file name
        call mpi_bcast( outdir,128,mpi_character,0,badlands_comm_world,ierr )
        call mpi_bcast( outdir1,128,mpi_character,0,badlands_comm_world,ierr )
        call mpi_bcast( outdir2,128,mpi_character,0,badlands_comm_world,ierr )
        call mpi_bcast( outdir3,128,mpi_character,0,badlands_comm_world,ierr )

        return

    end subroutine simulation_initialisation
    ! ============================================================================
    !> Subroutine read_basement
    !! Reads the basement file
    !<
    ! ============================================================================
    subroutine read_basement

        integer :: n, iu, ios, seed

        real( tkind ) :: hval

        lnbPts =  lstrat_X * lstrat_Y

        if( allocated( lcoordX ) ) deallocate( lcoordX )
        allocate( lcoordX( lnbPts ) )
        if( allocated( lcoordY ) ) deallocate( lcoordY )
        allocate( lcoordY( lnbPts ) )
        if( allocated( lcoordZ ) ) deallocate( lcoordZ )
        allocate( lcoordZ( lnbPts ) )
        if( allocated( lbase ) ) deallocate( lbase )
        allocate( lbase( lnbPts ) )
        if( allocated( soil_thick ) ) deallocate( soil_thick )
        allocate( soil_thick( lnbPts ) )
        soil_thick = initsoil
        if( flexureon == 1 )then
            if( .not. allocated( flexure ) ) allocate( flexure( lnbPts ) )
            flexure = 0.0_8
        endif

        ! Read DEM input file
        iu = 1
        open(iu,file=fstrata, status="old", action="read", iostat=ios )
        do n = 1, lnbPts
            read( iu, * ) lcoordX( n ), lcoordY( n ), lbase( n )
            if( noisegen > 0.0_8 .and. iam == 0 )then
                call random_seed( size = seed )
                call random_number( harvest = hval )
                lbase( n ) = lbase( n ) + hval * noisegen
            endif
        enddo
        close( iu )

        if( noisegen > 0.0_8 )&
            call mpi_bcast( lbase, lnbPts, dbl_type, 0, badlands_comm_world, ierr )

        ! Define top elevation
        lcoordZ = lbase

        ! Define grid extent
        strat_dx = lcoordX( 2 ) - lcoordX( 1 )
        lstrat_xo = lcoordX( 1 )
        lstrat_yo = lcoordY( 1 )
        lstrat_xm = lcoordX( lnbPts )
        lstrat_ym = lcoordY( lnbPts )

        ! Allocate neighbourhood
        allocate( lngbID( lnbPts, 8 ) )
        lngbID = -1
        do n = 1, lnbPts
            ! Corners
            if( n == 1 )then
                lngbID( n, 1 ) = n + lstrat_X
                lngbID( n, 2 ) = n + lstrat_X + 1
                lngbID( n, 3 ) = n + 1
            elseif( n == lstrat_X )then
                lngbID( n, 1 ) = n + lstrat_X
                lngbID( n, 7 ) = n - 1
                lngbID( n, 8 ) = n + lstrat_X - 1
            elseif( n == lnbPts - lstrat_X + 1 )then
                lngbID( n, 3 ) = n + 1
                lngbID( n, 4 ) = n - lstrat_X + 1
                lngbID( n, 5 ) = n - lstrat_X
            elseif( n == lnbPts )then
                lngbID( n, 5 ) = n - lstrat_X
                lngbID( n, 6 ) = n - lstrat_X - 1
                lngbID( n, 7 ) = n - 1
            ! Borders
            elseif( lcoordY( n ) == lstrat_yo )then
                lngbID( n, 1 ) = n + lstrat_X
                lngbID( n, 2 ) = n + lstrat_X + 1
                lngbID( n, 3 ) = n + 1
                lngbID( n, 7 ) = n - 1
                lngbID( n, 8 ) = n + lstrat_X - 1
            elseif( lcoordX( n ) == lstrat_xo )then
                lngbID( n, 1 ) = n + lstrat_X
                lngbID( n, 2 ) = n + lstrat_X + 1
                lngbID( n, 3 ) = n + 1
                lngbID( n, 4 ) = n - lstrat_X + 1
                lngbID( n, 5 ) = n - lstrat_X
            elseif( lcoordX( n ) == lstrat_xo + ( lstrat_X - 1 ) * strat_dx )then
                lngbID( n, 1 ) = n + lstrat_X
                lngbID( n, 5 ) = n - lstrat_X
                lngbID( n, 6 ) = n - lstrat_X - 1
                lngbID( n, 7 ) = n - 1
                lngbID( n, 8 ) = n + lstrat_X - 1
            elseif( lcoordY( n ) == lstrat_yo + ( lstrat_Y - 1 ) * strat_dx )then
                lngbID( n, 3 ) = n + 1
                lngbID( n, 4 ) = n - lstrat_X + 1
                lngbID( n, 5 ) = n - lstrat_X
                lngbID( n, 6 ) = n - lstrat_X - 1
                lngbID( n, 7 ) = n - 1
            else
                lngbID( n, 1 ) = n + lstrat_X
                lngbID( n, 2 ) = n + lstrat_X + 1
                lngbID( n, 3 ) = n + 1
                lngbID( n, 4 ) = n - lstrat_X + 1
                lngbID( n, 5 ) = n - lstrat_X
                lngbID( n, 6 ) = n - lstrat_X - 1
                lngbID( n, 7 ) = n - 1
                lngbID( n, 8 ) = n + lstrat_X - 1
            endif
        enddo

        return

    end subroutine read_basement
    ! ============================================================================
    !> Subroutine read_initial_deposit
    !! Reads the initial deposit
    !<
    ! ============================================================================
    subroutine read_initial_deposit

        logical :: found

        integer :: p, iu, ios, n, ks

        real( tkind ) :: sedh

        ! Read through the initial deposit layers
        do p = 1, InitDep

            iu = 40 + iam
            ! Open deposit file
            inquire(FILE=fdep( p ), EXIST=found)
            if(found)then
                open(iu,file=fdep( p ),status="old",action="read",iostat=ios)
                rewind(iu)
            else
                write(6,*)'Cannot find the deposit file:',trim( fdep( p ) )
                stop
            endif

            ! Read the sediment thicknesses
            do n = 1, lnbPts
                read(iu, *) stratLays( p, n, 1:totgrn )
                sedh = 0.0_8
                do ks = 1, totgrn
                    sedh = stratLays( p, n, ks ) + sedh
                enddo
                lcoordZ( n ) = lcoordZ( n ) + sedh
                if( sedh > 0 )then
                    v_nbLays( n ) = v_nbLays( n ) + 1
                    v_LaysID( n, p ) = p
                endif
            enddo
            close( iu )

        enddo

        return

    end subroutine read_initial_deposit
    ! ============================================================================
    !> Subroutine define_large_DEM_grid
    !! Define the grid of the large DEM.
    !<
    ! ============================================================================
    subroutine define_DEM_grid

        integer :: n, k, p, m, id

        strat_X = lstrat_X - 2 * border
        strat_Y = lstrat_Y - 2 * border

        nbPts = strat_X * strat_Y
        nbFcs = ( strat_X - 1 ) * ( strat_Y - 1 )
        lnbFcs = ( lstrat_X - 1 ) * ( lstrat_Y - 1 )

        if( allocated( coordX ) ) deallocate( coordX )
        allocate( coordX( nbPts ) )
        if( allocated( coordY ) ) deallocate( coordY )
        allocate( coordY( nbPts ) )
        if( allocated( coordZ ) ) deallocate( coordZ )
        allocate( coordZ( nbPts ) )
        if( allocated( fcoordX ) ) deallocate( fcoordX )
        allocate( fcoordX( nbFcs ) )
        if( allocated( fcoordY ) ) deallocate( fcoordY )
        allocate( fcoordY( nbFcs ) )
        if( allocated( fptIDs ) ) deallocate( fptIDs )
        allocate( fptIDs( nbFcs, 4 ) )
        if( allocated( lfptIDs ) ) deallocate( lfptIDs )
        allocate( lfptIDs( lnbFcs, 4 ) )
        if( allocated( fPid ) ) deallocate( fPid )
        allocate( fPid( nbFcs ) )

        ! Define large grid extent
        strat_xo = lstrat_xo + strat_dx * border
        strat_yo = lstrat_yo + strat_dx * border
        n = 0
        do k = 1, strat_Y
            do p = 1, strat_X
                n = n + 1
                coordX( n ) = strat_xo + ( p - 1 ) * strat_dx
                coordY( n ) = strat_yo + ( k - 1 ) * strat_dx
            enddo
        enddo
        strat_xm = coordX( nbPts )
        strat_ym = coordY( nbPts )

        ! Allocate neighbourhood for the DEM
        allocate( ngbID( nbPts, 8 ) )
        ngbID = -1
        do n = 1, nbPts
            ! Corners
            if( n == 1 )then
                ngbID( n, 1 ) = n + strat_X
                ngbID( n, 2 ) = n + strat_X + 1
                ngbID( n, 3 ) = n + 1
            elseif( n == strat_X )then
                ngbID( n, 1 ) = n + strat_X
                ngbID( n, 7 ) = n - 1
                ngbID( n, 8 ) = n + strat_X - 1
            elseif( n == nbPts - strat_X + 1 )then
                ngbID( n, 3 ) = n + 1
                ngbID( n, 4 ) = n - strat_X + 1
                ngbID( n, 5 ) = n - strat_X
            elseif( n == nbPts )then
                ngbID( n, 5 ) = n - strat_X
                ngbID( n, 6 ) = n - strat_X - 1
                ngbID( n, 7 ) = n - 1
            ! Borders
            elseif( coordY( n ) == strat_yo )then
                ngbID( n, 1 ) = n + strat_X
                ngbID( n, 2 ) = n + strat_X + 1
                ngbID( n, 3 ) = n + 1
                ngbID( n, 7 ) = n - 1
                ngbID( n, 8 ) = n + strat_X - 1
            elseif( coordX( n ) == strat_xo )then
                ngbID( n, 1 ) = n + strat_X
                ngbID( n, 2 ) = n + strat_X + 1
                ngbID( n, 3 ) = n + 1
                ngbID( n, 4 ) = n - strat_X + 1
                ngbID( n, 5 ) = n - strat_X
            elseif( coordX( n ) == strat_xo + ( strat_X - 1 ) * strat_dx )then
                ngbID( n, 1 ) = n + strat_X
                ngbID( n, 5 ) = n - strat_X
                ngbID( n, 6 ) = n - strat_X - 1
                ngbID( n, 7 ) = n - 1
                ngbID( n, 8 ) = n + strat_X - 1
            elseif( coordY( n ) == strat_yo + ( strat_Y - 1 ) * strat_dx )then
                ngbID( n, 3 ) = n + 1
                ngbID( n, 4 ) = n - strat_X + 1
                ngbID( n, 5 ) = n - strat_X
                ngbID( n, 6 ) = n - strat_X - 1
                ngbID( n, 7 ) = n - 1
            else
                ngbID( n, 1 ) = n + strat_X
                ngbID( n, 2 ) = n + strat_X + 1
                ngbID( n, 3 ) = n + 1
                ngbID( n, 4 ) = n - strat_X + 1
                ngbID( n, 5 ) = n - strat_X
                ngbID( n, 6 ) = n - strat_X - 1
                ngbID( n, 7 ) = n - 1
                ngbID( n, 8 ) = n + strat_X - 1
            endif
        enddo

        allocate( lvertexID( nbPts ) )
        lvertexID = -1
        p = 0
        n = 0
        do k = 1, lstrat_Y
            do m = 1, lstrat_X
                p = p + 1
                if( lcoordX( p ) >= strat_xo .and. lcoordX( p ) <= strat_xm .and. &
                    lcoordY( p ) >= strat_yo .and. lcoordY( p ) <= strat_ym )then
                    n = n + 1
                    lvertexID( n ) = p
                endif
            enddo
        enddo

        ! Allocate the faces values for stratal grid
        id = 1
        k = 1
        do n = 1, lnbFcs
            ! Points
            lfptIDs( n, 1 ) = id
            lfptIDs( n, 2 ) = id + 1
            lfptIDs( n, 3 ) = lstrat_X + id + 1
            lfptIDs( n, 4 ) = lstrat_X + id
            id = id + 1
            if( id == k * lstrat_X )then
                id = k * lstrat_X + 1
                k = k + 1
            endif
        enddo

        ! Allocate the faces values for DEM
        id = 1
        k = 1
        do n = 1, nbFcs
            ! Points
            fptIDs( n, 1 ) = id
            fptIDs( n, 2 ) = id + 1
            fptIDs( n, 3 ) = strat_X + id + 1
            fptIDs( n, 4 ) = strat_X + id
            ! Coordinates
            fcoordX( n ) = real( coordX( id ) + strat_dx * 0.50_8 )
            fcoordY( n ) = real( coordY( id ) + strat_dx * 0.50_8 )
            id = id + 1
            if( id == k * strat_X )then
                id = k * strat_X + 1
                k = k + 1
            endif
        enddo

        return

    end subroutine define_DEM_grid
    ! ============================================================================
    !> Subroutine define_band_partitioning
    !! Subroutine define_band_partitioning used to send the partitioned grid values.
    !<
    ! ============================================================================
    subroutine define_band_partitioning

        integer :: n, avecol, extras, cols, k, nID, m, p, p1, strat_width

        integer, dimension( nbPts ) :: nodepShare
        integer, dimension( nbPts, 2 ) :: nodepShareID

        real( tkind ) :: maxP( gproc ), ptmax

        part_rows = 0

        ! Check higher cell number in X or Y direction
!        if( strat_Y >= strat_X )then
            avecol = int( ( strat_Y - 1 ) / gproc )
            extras = mod( ( strat_Y - 1 ), gproc )
!        else
!            avecol = int( ( strat_X - 1 ) / gproc )
!            extras = mod( ( strat_X - 1 ), gproc )
!        endif

        if( iam == 0 )then

            ! Get the number of columns for each processor
            do k = 0, gproc - 1
                cols = avecol
                if( k < extras ) cols = avecol + 1
                if( k == 0 )then
!                    if( strat_Y >= strat_X )then
                        maxP( k + 1 ) = strat_yo + cols * strat_dx
!                    else
!                        maxP( k + 1 ) = strat_xo + cols * strat_dx
!                    endif
                else
                    maxP( k + 1 ) = maxP( k ) + cols * strat_dx
                endif
                part_rows( k+1 ) = cols + 1
            enddo

            ! Loop over the faces to find the processor they belong to.
            do n = 1, nbFcs
                partition_face: do k = 1, gproc
!                    if( strat_Y >= strat_X )then
                        ptmax = fcoordY( n )
!                    else
!                        ptmax = fcoordX( n )
!                    endif
                    if( ptmax < maxP( k ) )then
                        fPid( n ) = k - 1
                        exit partition_face
                    endif
                enddo partition_face
            enddo

        endif

        ! Broadcast values globally
        call mpi_bcast( part_rows, nproc, int_type, 0, badlands_comm_world, ierr )
        call mpi_bcast( fPid, nbFcs, int_type, 0, badlands_comm_world, ierr )
!        if( strat_Y >= strat_X )then
            local_lnbPts = part_rows( iam + 1 ) * lstrat_X
!        else
!            local_lnbPts = part_rows( iam + 1 ) * lstrat_Y
!        endif

        if( iam == 0 .or. iam == gproc - 1 ) local_lnbPts = local_lnbPts + lstrat_X
        ! Get the number of columns for each processor on the large DEM
        part_rows( 1 ) = part_rows( 1 ) + 1
        do k = 2, gproc
            part_rows( k ) = part_rows( k - 1 ) + part_rows( k )
            if( k == gproc ) part_rows( k ) = part_rows( k ) + 1
        enddo

        ! Find the number of faces lying in each partition
        k = 0
        do n = 1, nbFcs
            if( fPid( n ) == iam )then
                k = k + 1
                local_Faces = k
            endif
        enddo

        ! Define processor id for each node based on faces
        nodepShare = 0
        do n = 1, nbFcs
            do k = 1, 4
                nID = fptIDs( n, k )
                if( nID > 0 .and. nID <= nbPts )then
                    if( nodepShare( niD ) > 0 )then
                        do m = 1, nodepShare( nID )
                            if( nodepShareID( nID, m ) == fPid( n ) )goto 20
                        enddo
                    endif
                    nodepShare( nID ) =  nodepShare( nID ) + 1
                    if( nodepShare( nID ) > 2 )print*,'Something went wrong when allocating processor ID to nodes'
                    nodepShareID( nID, nodepShare( nID ) ) = fPid( n )
                else
                    print*,'Something went wrong when looking for the points ID of stratal faces'
                    stop
                endif
20          continue
            enddo
        enddo

        ! Define local nodes number for each processors
        nID = 0
        do n = 1, nbPts
            do k = 1, nodepShare( n )
                if( nodepShareID( n, k ) == iam )then
                    nID = nID + 1
                    local_nbPts = nID
                endif
            enddo
        enddo

        if( local_nbPts > 0 )then

            ! For global points find local id
            if( allocated( global_nid ) )  deallocate( global_nid )
            allocate( global_nid( local_nbPts ) )
            nID = 0
            do n = 1, nbPts
                do k = 1, nodepShare( n )
                    if( nodepShareID( n, k ) == iam )then
                        nID = nID + 1
                        global_nid( nID ) = n
                    endif
                enddo
            enddo

            ! Compute faces points ID
            if( allocated( locf_points ) )  deallocate( locf_points )
            allocate( locf_points( local_Faces, 4 ) )
            p = 1
            p1 = 1
            strat_width = strat_X
            do k = 1, local_Faces
                locf_points( k, 1 ) = p
                locf_points( k, 2 ) = p + 1
                locf_points( k, 3 ) = p + strat_width
                locf_points( k, 4 ) = p + strat_width + 1
                if( p1 < strat_width - 1 )then
                    p1 = p1 + 1
                    p = p + 1
                else
                    p1 = 1
                    p = p + 2
                endif
            enddo

            if( allocated( global_lnid ) )  deallocate( global_lnid )
            allocate( global_lnid( local_lnbPts ) )
            p = 0
            p1 = 1
            do n = 1, lstrat_Y
                do k = 1, lstrat_X
                    p = p + 1
                    if( iam == 0 )then
                        if( n <= part_rows( iam + 1 ) )then
                            global_lnid( p1 ) = p
                            p1 = p1 + 1
                        endif
                    elseif( iam < gproc )then
                        if( n > part_rows( iam ) .and. n <= part_rows( iam + 1 ) )then
                            global_lnid( p1 ) = p
                            p1 = p1 + 1
                        endif
                    endif
                enddo
            enddo

        endif

        return

    end subroutine define_band_partitioning
    ! ============================================================================

end module init_phase
! ============================================================================
