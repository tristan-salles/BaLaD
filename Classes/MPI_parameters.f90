! ============================================================================
! Name        : MPI_parameters.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file MPI_parameters.f90
!
! Description :  This module encapsules parameters used to setup MPI components.
!
!<
! ============================================================================
!> Module parallel
!<
module parallel

#include "mpif.h"

    ! INTEL COMPILER
!    use ifport

    ! Error code
    integer :: ierr

    ! Processor ID
    integer :: iam

    ! Number of processors for FW
    integer :: nproc

    ! Number of processors for grid
    integer :: gproc

    ! Number of processors for checkpointing
    integer :: checkproc

    ! MPI integer type communicator
    integer :: int_type

    ! MPI double type communicator
    integer :: dbl_type

    ! MPI logical type communicator
    integer :: lgc_type

    ! MPI max type communicator
    integer :: max_type

    ! MPI min type communicator
    integer :: min_type

    ! MPI sum type communicator
    integer :: sum_type

    ! Badlands communicator
    integer :: badlands_comm_world

    ! New communicator
    integer :: new_comm


    ! Single or double precision
    integer, parameter :: singlep = kind(0.0)
    integer, parameter :: doublep = kind(0.0d0)

    integer, parameter :: tkind = doublep
    integer, parameter :: sizesp = 4
    integer, parameter :: sizedp = 8

contains

    ! ============================================================================
    !> Subroutine term_command.
    !!
    !<
    ! ============================================================================
    subroutine term_command( cmds )

        logical( 4 ) :: result
        character(len=128) :: cmds

        result = .false.

        ! INTEL fortran compiler
        !result = systemqq( cmds )
        ! GNU fortran compiler
        call system( cmds )

        return

    end subroutine term_command
    ! ============================================================================
    !> Subroutine count_lines_file.
    !! Counts the number of line in a given file.
    !<
    ! ============================================================================
    subroutine count_lines_file( filetoread, nlines )

        integer :: lunit, stat, nlines
        character(len=128) :: filetoread

        lunit = 12098
        nlines = 0
        open( lunit, file = filetoread )
        lp: do
            read( lunit,*, iostat = stat )
            if( stat /= 0 ) exit lp
            nlines = nlines + 1
        enddo lp

        close (lunit)

        return

    end subroutine count_lines_file
    ! ============================================================================
    !> Subroutine mpi_specific_types_declaration
    !! This subroutine defines MPI types communicator
    !<
    ! ============================================================================
    subroutine mpi_specific_types_declaration

        badlands_comm_world = mpi_comm_world
        int_type = mpi_integer
        dbl_type = mpi_double_precision
        lgc_type = mpi_logical
        max_type = mpi_max
        min_type = mpi_min
        sum_type = mpi_sum

        return

    end subroutine mpi_specific_types_declaration
    ! ============================================================================

end module parallel
! ============================================================================
