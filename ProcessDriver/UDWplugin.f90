! ============================================================================
! Name        : UDWplugin.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file UDWplugin.f90
!
! Description :  This module is used to couple UDW and Badlands using the vertical only
! component.
!
!<
! ============================================================================
!> Module UDWplugin 
!<
module UDWplugin

    use parallel
    use file_data
    use time_data
    use forces_data

    implicit none

contains

    ! ============================================================================
    !> Subroutine create_vtk_topsurface
    !! Create the top elevation file for Underworld.
    !<
    ! ============================================================================
    subroutine create_vtk_topsurface( uorl )

        integer :: iunit, ios, k, p, n, uorl

        ! Open a file tunnel
        iunit = 20
        open(iunit,file=fudw,status="replace",action="write",iostat=ios)
        rewind(iunit)

        ! Conform vtk top file format
        write( iunit,'(a26)')'# vtk DataFile Version 2.0'
        write( iunit,'(a25)')'LECODE to Underworld Data'
        write( iunit,'(a5)')'ASCII'
        write( iunit,'(a25)')'DATASET STRUCTURED_POINTS'
        write( iunit,'(a10,1x,i10,1x,i10,1x,i2)')'DIMENSIONS',lstrat_X,lstrat_Y,1
        write(iunit,'(a6,1x,f12.3,1x,f12.3,1x,f12.3)')'ORIGIN',lstrat_xo,lstrat_yo,0.0
        write(iunit,'(a7,1x,f12.3,1x,f12.3,1x,f12.3)')'SPACING',strat_dx,strat_dx,1.0
        write(iunit,'(a10,1x,i10)')'POINT_DATA',lstrat_X*lstrat_Y
        write(iunit,'(a26)')'SCALARS elevation double 1'
        write(iunit,'(a20)')'LOOKUP_TABLE default'

        ! Write vtk file
        p = 0
        do k = 1, lstrat_Y
            do n = 1, lstrat_X - 1
                p = p + 1
                write(iunit,'(f12.3,1x)',advance='no') lCoordZ( p )
            enddo
            p = p + 1
            write(iunit,'(f12.3)') lCoordZ( p )
        enddo

        write(iunit,'(a21)')'SCALARS fieldID int 1'
        write(iunit,'(a20)')'LOOKUP_TABLE default'

        p = 0
        do k = 1, lstrat_Y
            do n = 1, lstrat_X - 1
                p = p + 1
                write(iunit,'(i6,1x)',advance='no') p
            enddo
            p = p + 1
            write(iunit,'(i6)') p
        enddo

        close( iunit )

        ! Update/Create the maestro file
        open(iunit,file=maestro,status="replace",action="write",iostat=ios)
        rewind(iunit)

        if( uorl == 0 )then
            write(iunit,'(a1)')'U'
            write(iunit,'(a1)')' '
        else
            write(iunit,'(a1)')'L'
            write(iunit,'(a1)')' '
        endif
        close( iunit )

        return

    end subroutine create_vtk_topsurface
    ! ============================================================================
    !> Subroutine udw_plugin_wait_function
    !! Wait for Underworld to compute the displacement field.
    !<
    ! ============================================================================
    subroutine udw_plugin_wait_function

        integer :: iunit, ios, k

        character(len=1) :: charac

        charac = 'U'

        if( iam == 0 )then

112     continue

            do while( charac /= 'L' )
                iunit = 12

                ! Read the maestro file
                open(iunit,file=maestro,status="old",action="read",iostat=ios)
                rewind(iunit)

                ! Read character
                read(iunit,'(a1)',err=112,end=112,iostat=k) charac
                if( k /= 0 ) goto 112
                close( iunit )

                call Sleep( 1 )

           enddo
       endif

       call mpi_barrier( badlands_comm_world, ierr )

       return

    end subroutine udw_plugin_wait_function
    ! ============================================================================

end module UDWplugin
! ============================================================================
