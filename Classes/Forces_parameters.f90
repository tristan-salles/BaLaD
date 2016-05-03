! ============================================================================
! Name        : Forces_parameters.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Forces_parameters.f90
!
! Description :  This module encapsulates the external forces module.
!
!<
! ============================================================================
!> Module forces_data
!<
module forces_data

    use parallel
    use file_data
    use flow_data
    use time_data
    use mesh_data

    implicit none

   ! Water density
   real( tkind ), parameter :: fluid_density = 1015.0_8

   ! Sea density
   real( tkind ), parameter :: sea_density = 1027.0_8

    ! Flexural rigidity (Nm)
    real( tkind ) :: flex_rigidity

    ! Crust density
    real( tkind ) :: crust_density

   !> Sea-level parameters
   type sl_par
      ! Lower sea-level elevation
      real( tkind )  :: sea1
      ! Upper sea-level elevation
      real( tkind ) :: sea2
      ! Lower sea-level time
      real( tkind ) :: time1
      ! Upper sea-level time
      real( tkind ) :: time2
      ! Simulation current time
      real( tkind ) :: tsim
      ! Interpolated sea-level elevation
      real( tkind ) :: sl_int
   end type sl_par

    !> Sea level fluctuations type
    type sea_fluc
        ! Output sea level
        logical  :: output
        ! Flag for sea level fluctuations
        logical :: sealevel
        ! Number of sea fluctuation events
        integer  :: event
        ! Old sea level value
        real( tkind ) :: last_sea
        ! Actual sea level value
        real( tkind ) :: actual_sea
    end type sea_fluc
    type( sea_fluc ) :: gsea
    ! Sea-level elevation data
    real( tkind ),dimension(:,:),allocatable :: sealvl

    ! Is there a new displacement field to load
    logical :: new_disp

    !> Vertical displacement type
    type vert_disp
        ! Flag for displacements
        logical :: vdisp
        ! Actual displacement event
        integer :: actual
        ! Displacements event numbers
        integer  :: event
        ! Time from last displacement call
        real( tkind ) :: lastdisp
        ! Time elapsed between current and previous displacement time
        real( tkind ) :: disptime
    end type vert_disp
    type( vert_disp ) :: gdisp, vdisp
    ! Displacement file unit to get a continuous displacement field
    integer( tkind ),dimension(:),allocatable :: gdisp_fill, vdisp_fill

    ! Displacement time values in the input file
    real( tkind ),dimension(:,:),allocatable :: gdisp_time, vdisp_time

    ! Displacement rate
    real( tkind ),dimension(:),allocatable :: disp, disprate

    !> Porosity type
    real( tkind ) :: initdepporo
    type porosity_lt
        ! Flag for sea level fluctuations
        logical :: compaction
        ! Effective pressure number
        integer :: ePnb
    end type porosity_lt
    type( porosity_lt ) :: gporo
    ! Effective pressure
    real( tkind ),dimension(:),allocatable :: effPressure
    ! Porosity table values
    real( tkind ),dimension(:,:),allocatable :: porosity

    !> Hemipelagic parameters
    logical :: hemi_flag

    ! Hemipelagic material number
    integer :: hemi_mat

    integer :: contourhemi_nb

    ! Hemipelagites time
    real( tkind ) :: last_hemi

    real( tkind ),dimension(:,:),allocatable :: hemipelagic

contains

    ! ============================================================================
    !> Subroutine read_sealevel_file()
    !! Reads formatted data of sea level fluctuations.
    !<
    ! ============================================================================
    subroutine read_sealevel_file

        ! Parameters Declaration
        logical :: found
        integer :: i, k, i2
        integer :: iu, ios

        character(len=128) :: text

        ! Find and open the sea level file
        iu = 45 + iam
        inquire( file=fsea, exist=found )
        if(found)then
            open(iu,file=fsea,status="old",action="read",iostat=ios)
            rewind(iu)
        else
            print*,'Cannot find sea level file.',iam
            call mpi_finalize( ierr )
            stop
        endif

        ! Find the number of events
        k = 0
        do
            read(iu,'(a128)',iostat=ios,end=30) text
            k = k + 1
        enddo
30      continue
        rewind(iu)

        ! Define the number of sea fluctuation events & allocate array
        gsea%event = k
        allocate( sealvl( 2, gsea%event ) )

        ! Sea-level fluctuations
        do k = 1, gsea%event
            read(iu, '(a128)' ,iostat=ios, end=60 )text
            i2=len_trim(text)
            if( text( i2:i2 ) == char(13) )then
                i2=i2-1
            endif
            read(text(1:i2),*,iostat=ios)( sealvl( i, k ), i = 1, 2 )
            if( ios /= 0 )then
                print*,'Error reading sea level file.',iam
                call mpi_finalize( ierr )
                stop
            endif
            if( k > 1 )then
                if( sealvl( 1, k ) <= sealvl( 1, k - 1 ) )then
                    print*,'Problem in sea level file definition.',iam
                    call mpi_finalize( ierr )
                    stop
                endif
            endif
        enddo

60      close(iu)

        gsea%last_sea = 0.0_8

        ! Update the sea level value
        call eustatism

        return

    end subroutine read_sealevel_file
    ! ============================================================================
    !> Subroutine eustatism()
    !! Performs the change of water elevation according to sea-level fluctuations.
    !<
    ! ============================================================================
    subroutine eustatism

        ! Parameters Declaration
        integer :: i
        type( sl_par ) :: sl_param

        if( .not. gsea%sealevel ) return

        gsea%last_sea = gsea%actual_sea

        i = 1
        ! Find time for the considered event
        do
            if( i > gsea%event )exit
            if( tnow < sealvl( 1, i ) )exit
            sl_param%time1 = sealvl( 1, i )
            sl_param%sea1 = sealvl( 2, i )
            i = i + 1
        enddo

        ! Linear interpolation of sea level elevation based on current time
        if( i <= gsea%event )then
            sl_param%time2 = sealvl( 1, i )
            sl_param%sea2 = sealvl( 2, i )
            sl_param%tsim = tnow
            call sealvl_interpolation( sl_param )
            gsea%actual_sea = sl_param%sl_int
        ! If event greater than max defined allocate with last value
        else
            gsea%actual_sea = sl_param%sea1
            sl_param%time2 = sealvl( 1, i - 1 )

        endif

        ! Declare the next sea level fluctuation update
        if( ( sl_param%time2 < tsampling ) .and. ( sl_param%time2 - tnow ) > time_tolerance )&
            tsampling = sl_param%time2 + tor

        return

    end subroutine eustatism
    ! ============================================================================
    !> Subroutine sealvl_interpolation()
    !! Performs the sea-level interpolation function.
    !!
    !! \arg \c use precision_data
    !!
    !! \param sl_param
    !<
    ! ============================================================================
    subroutine sealvl_interpolation( sl_param )

        ! Parameters Declaration
        type( sl_par ) :: sl_param

        sl_param%sl_int = ( sl_param%tsim - sl_param%time1 ) / ( sl_param%time2 - sl_param%time1 )
        sl_param%sl_int = sl_param%sl_int * ( sl_param%sea2 - sl_param%sea1 )
        sl_param%sl_int = sl_param%sl_int + sl_param%sea1

        return

    end subroutine sealvl_interpolation
    ! ============================================================================
    !> Subroutine assign_vertical_displacements_rate()
    !! Performs a simple vertical displacement fields over the nodes of the stratal mesh.
    !! It generates the vertical_displacements field which contains the
    !! displacements in the Z direction for the simulated period.
    !! \param period
    !<
    ! ============================================================================
    subroutine assign_vertical_displacements_rate( period )

        ! Parameters Declaration
        logical :: found

        integer,intent( in ) :: period

        integer :: k, i, kn
        integer :: iu, ios

        if( .not. allocated( gdisp_fill ) ) allocate( gdisp_fill( gdisp%event ) )
        if( .not. allocated( disp ) ) allocate( disp( lnbPts ) )
        if( .not. allocated( disprate ) ) allocate( disprate( lnbPts ) )

        ! First displacement period allocated required arrays
        if( period == 1 )then
            if( gdisp_time( 1, 1 ) > time_start ) vdisp%event = 1

            ! Find the vertical displacement event number
            do k = 1, gdisp%event
                vdisp%event = vdisp%event + 1
                if( gdisp_time( k, 1 ) >= gdisp_time( k, 2 ) )then
                    print*,'Error in displacements time declaration in event.',k,iam
                    call mpi_finalize( ierr )
                    stop
                endif
                if( k < gdisp%event )then
                    if( gdisp_time( k, 2 ) > gdisp_time( k + 1, 1 ) )then
                        print*,'Error in displacements time declaration between events.',iam
                        call mpi_finalize( ierr )
                        stop
                    endif
                    if( gdisp_time( k, 2 ) /= gdisp_time( k + 1, 1 ) ) vdisp%event = vdisp%event + 1
                endif
            enddo

            ! If last displacement period stops before the end of the simulation
            ! add a new event during the last period
            if( gdisp_time( gdisp%event, 2 ) < time_end ) vdisp%event = vdisp%event + 1

            ! In case the number of vertical displacement is different from the number of
            ! of tectonic events defined in the XmL input file defines the new event parameters
            if( vdisp%event /= gdisp%event )then

                allocate( vdisp_time( vdisp%event, 2 ), vdisp_fill( vdisp%event ) )
                kn = 1
                if( gdisp_time( 1, 1 ) > time_start )then
                    vdisp_time( kn, 1 ) = time_start
                    vdisp_time( kn, 2 ) = gdisp_time( 1, 1 )
                    vdisp_fill( kn ) = 0
                    kn = kn + 1
                endif

                do k = 1, gdisp%event
                    vdisp_time( kn, 1 ) = gdisp_time( k, 1 )
                    vdisp_time( kn, 2 ) = gdisp_time( k, 2 )
                    vdisp_fill( kn ) = k
                    kn = kn + 1
                    if( k < gdisp%event )then
                        if( gdisp_time( k, 2 ) /= gdisp_time( k + 1, 1 ) )then
                                vdisp_time( kn, 1 ) = gdisp_time( k, 2 )
                                vdisp_time( kn, 2 ) = gdisp_time( k + 1, 1 )
                                vdisp_fill( kn ) = 0
                                kn = kn + 1
                        endif
                    endif
                enddo

                if( gdisp_time( gdisp%event, 2 ) < time_end )then
                    vdisp_time( kn, 2 ) = time_end
                    vdisp_time( kn, 1 ) = gdisp_time( gdisp%event, 2 )
                    vdisp_fill( kn ) = 0
                endif
                deallocate( gdisp_time, gdisp_fill )
                allocate( gdisp_time( vdisp%event, 2 ), gdisp_fill( vdisp%event ) )

                gdisp%event =  vdisp%event
                do k = 1,  vdisp%event
                    gdisp_time( k, 1 ) = vdisp_time( k, 1 )
                    gdisp_time( k, 2 ) = vdisp_time( k, 2 )
                    gdisp_fill( k ) = vdisp_fill( k )
                enddo
                deallocate( vdisp_time, vdisp_fill )

            ! Otherwise define the displacement event number directly
            else
                do k = 1, gdisp%event
                    gdisp_fill( k ) = k
                enddo
            endif

        endif

        ! Find actual event number
        gdisp%actual = 0
        do k = 1, gdisp%event
            if( tnow < gdisp_time( k, 2 ) .and. tnow >= gdisp_time( k, 1 ) ) gdisp%actual = k
        enddo

        ! Open displacements field file for the considered event number
        kn = gdisp%actual
        if( gdisp_fill( kn ) > 0 )then
            inquire( file=fdisp( gdisp_fill( kn )), exist=found )
            if( .not. found )then
                print*,'Cannot find displacement file.',iam
                call mpi_finalize( ierr )
                stop
            endif

            ! Allocate displacement values
            iu = 80 + iam
            open(iu,file=fdisp( gdisp_fill( kn )),status="old",action="read",iostat=ios)
            rewind( iu )

            ! Assign displacement array
            do i = 1, lnbPts
                read( iu, * ) disp( i )
            enddo

            close( iu )
        ! If this event doesn't correspond to any defined one initialise displacement values to 0
        else
            disp = 0.0_8
        endif

        return

    end subroutine assign_vertical_displacements_rate
    ! ============================================================================
    !> Subroutine update_local_displacements()
    !! Update displacement on the mesh node according to the displacement fields.
    !<
    ! ============================================================================
    subroutine update_local_displacements

        ! Parameters Declaration
        integer :: k
        integer :: period

        ! Take displacement event period number
        period = gdisp%actual

        ! Assign array
        do k =1, lnbPts
            ! Compute displacement in metres
            disprate( k ) = dble( disp( k ) / &
                ( (gdisp_time( period, 2 ) - gdisp_time( period, 1 ) ) * secyear ) )
        enddo

        return

    end subroutine update_local_displacements
    ! ============================================================================
    !> Subroutine apply_displacements()
    !! Apply displacement on the mesh node according to the displacement fields.
    !! This will read the previous SPModel mesh and create a new one based on
    !! displacements values and previous node and elements connectivity.
    !<
    ! ============================================================================
    subroutine apply_displacements

        ! Parameters Declaration
        integer :: k, n

        ! Get the elapsed time since previous displacement call
        gdisp%disptime = ( tnow - gdisp%lastdisp ) * secyear

        ! Apply the displacement to the deposit layers
        do k = 1, lnbPts
            lcoordZ( k ) = lcoordZ( k ) + disprate( k ) * gdisp%disptime
            lbase( k ) = lbase( k ) + disprate( k ) * gdisp%disptime
        enddo

        ! Broadcast top elevation to DEM grid
        do n = 1, nbPts
            coordZ( n ) = lcoordZ( lvertexID( n ) )
        enddo

        return

    end subroutine apply_displacements
    ! ============================================================================
    !> Subroutine load_displacements()
    !! load_displacements is used to load a new displacement field based on user defined
    !! displacement data.
    !! \param period
    !<
    ! ============================================================================
    subroutine load_displacements( period )

        ! Parameters Declaration
        integer,intent( in ) :: period

        ! Get new displacement fields
        call assign_vertical_displacements_rate( period )

        ! Update displacement locally
        call update_local_displacements

        ! Set displacement flag to false as we've just updated it
        new_disp = .false.

        return

    end subroutine load_displacements
    ! ============================================================================

end module forces_data
! ============================================================================
