! ============================================================================
! Name        : badlands.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file badlands.f90
!!
!! The is the main function of badlands.
!!
!<
! ============================================================================
program badlands

    use steady
    use fwpath
    use parallel
    use file_data
    use run_phase
    use init_phase
    use strata_out
    use time_data

    implicit none

    integer :: m

    real( tkind ) :: ts_st, ts_ed

    ! Initialise MPI components
    call mpi_init( ierr )
    call mpi_comm_size( mpi_comm_world, nproc, ierr )
    call mpi_comm_rank( mpi_comm_world, iam, ierr )
    ts_st = mpi_wtime( )
    call mpi_specific_types_declaration

    ! Get experiment file name
    m = iargc()
    fin_noflat = ' '
    if( m < 1 ) then
        if( iam == 0 ) write(6,*)'use: mpirun -np X ./lecode <input-file-name> '
        if( iam == 0 ) write(6,*)'where X refers to the number of processors for the run. '
        call mpi_finalize( ierr )
        stop
    elseif( m == 1 )then
        call getarg( 1, fin_noflat )
        if( iam == 0 ) call read_input_to_flatten
        call mpi_bcast( finput,128,mpi_character,0,badlands_comm_world,ierr )
    endif

    ! Initialisation phase
    call badlands_initialise

    ! Execution phase
    call badlands_run

    ! Finalise phase
    call badlands_finalise

    ts_ed = mpi_wtime( )
    if( iam == 0 ) write(6,*)'Simulation completed ...'
    if( iam == 0 ) write(6,100),'Time elapse (s):',ts_ed - ts_st

    ! Finalise MPI components
    call mpi_finalize( ierr )

100   format(A18,F12.4)

end program badlands
