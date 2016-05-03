! ============================================================================
! Name        : FIle_parameters.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file FIle_parameters.f90
!
! Description :  This module  encapsulates the file module used by badlands.
!
!<
! ============================================================================
!> Module file_data gives the file and directory names used for badlands simulations.
!<
module file_data

    use parallel

    implicit none

    ! Underworld time
    logical :: udw_plug

    ! Underworld time
    real( tkind ) :: udw_time

    ! Underworld top surface csv file
    character(len=128) :: fudw

    ! Underworld displacement ascii file
    character(len=128) :: fudisp

    ! Underworld maestro file
    character(len=128) :: maestro

    ! Output directory
    character(len=128) :: outdir

    ! Run files directory name
    character(len=128) :: outdir1

    ! Output directory name
    character(len=128) :: outdir2

    ! UDW directory name
    character(len=128) :: outdir3

    ! Ocean directory name
    character(len=128) :: outdir4

    ! Flatten XmL input file
    character(len=128) :: finput

    ! Output subdirectory name
    character(len=128) :: outputs

    ! Runfiles subdirectory name
    character(len=128) :: runfiles

    ! Sea level fluctuations file
    character(len=128) :: fsea

    ! XInclude XmL input file
    character(len=128) :: fin_noflat

    ! XIncluded input file
    character(len=128) :: included_file

    ! Temporary input files
    character(len=128) :: temp_file, temp_input

    ! Rain map files
    character(len=128), dimension( : ),allocatable :: frainmap

    ! Displacement files
    character(len=128), dimension( : ),allocatable :: fdisp

    ! Strata output files
    character(len=128) :: fstrato

    ! Strata nodes file
    character(len=128) :: fstrata

    ! Deposit files
    character(len=128), dimension( : ),allocatable :: fdep

    ! TIN nodes / elements file
    character(len=128) :: ftin, ftinCPP, ftin_out

    ! Porosity lookup table file
    character(len=128) :: fporosity

    ! Hemipelagic file
    character(len=128) :: fhemi

    ! Restart folder name
    character(len=128) :: restartfolder

    ! Restart file name
    character(len=128) :: frestart

contains

    ! ============================================================================
    !> Subroutine noblnk
    !! Subroutine noblnk is used to clear any blank from the end of a string.
    !! This subroutine is similar to the one implemented in Sedsim.
    !! \param string
    !<
    ! ============================================================================
    subroutine noblnk( string )

        integer:: i, j, lg
        character(len=128) :: string

        lg = len( string )
        do
            if( lg <= 0 .or. string( lg:lg ) /= ' ' ) exit
            lg = lg - 1
        enddo
        if ( lg > 0 ) then
            ! find first non-blank character
            i = 1
            do
                if(i > lg .or. ( string( i:i ) /= ' ' .and. string /= ' ' ) ) exit
                i = i + 1
            end do
            ! determine end of continuous (non-blank) string
            j = i
            do
                if( j > lg )then
                    exit
                elseif( string( j:j ) == ' ' )then
                    exit
                elseif( string == '  ' )then
                    exit
                elseif( j == 128 )then
                    exit
                endif
                j = j + 1
            end do
            ! j points to first blank position or past end of string; adjust to last
            ! non-blank position in string
            j = min( j - 1, lg )
            string = string( i:j )
            if( j < len( string ) ) string( j + 1:len( string ) ) = ' '
        else
            ! there were only blanks in string
            string = ' '
        end if

        return

    end subroutine noblnk
    ! ============================================================================
    !> Subroutine addpath()
    !! Subroutine addpath is used to add a path onto any filenames.
    !! \param fname
    !<
    ! ============================================================================
    subroutine addpath( fname )

        integer:: pathlen, flen
        character(len=128) :: fname, dname, dummy

        ! for files to be read, they'll be in the session path
        dname=' '
        call noblnk(outdir)
        pathlen = len_trim(outdir)
        dname(1:pathlen)=outdir
        dname(pathlen+1:pathlen+1) = '/'
        pathlen = pathlen+1
        call noblnk(fname)
        flen = len_trim(fname)
        dummy=' '
        dummy=fname
        fname=' '
        fname(1:pathlen)=dname(1:pathlen)
        fname(pathlen+1:pathlen+flen)=dummy(1:flen)

        return

    end subroutine addpath
    ! ============================================================================
    !> Subroutine addpath1()
    !! Subroutine addpath1 is used to add a path onto any filenames.
    !! \param fname
    !<
    ! ============================================================================
    subroutine addpath1( fname )

        integer:: pathlen, flen
        character(len=128):: fname, dname, dummy

        ! for files to be read, they'll be in the session path
        dname=' '
        call noblnk(outdir1)
        pathlen = len_trim(outdir1)
        dname(1:pathlen)=outdir1
        dname(pathlen+1:pathlen+1) = '/'
        pathlen = pathlen+1
        call noblnk(fname)
        flen = len_trim(fname)
        dummy=' '
        dummy=fname
        fname=' '
        fname(1:pathlen)=dname(1:pathlen)
        fname(pathlen+1:pathlen+flen)=dummy(1:flen)

        return

    end subroutine addpath1
    ! ============================================================================
    !> Subroutine addpath2()
    !! Subroutine addpath2 is used to add a path onto any filenames.
    !! \param fname
    !<
    ! ============================================================================
    subroutine addpath2( fname )

        integer:: pathlen,flen
        character(len=128):: fname, dname, dummy

        ! for files to be read, they'll be in the session path
        dname=' '
        call noblnk(outdir2)
        pathlen = len_trim(outdir2)
        dname(1:pathlen)=outdir2
        dname(pathlen+1:pathlen+1) = '/'
        pathlen = pathlen+1
        call noblnk(fname)
        flen = len_trim(fname)
        dummy=' '
        dummy=fname
        fname=' '
        fname(1:pathlen)=dname(1:pathlen)
        fname(pathlen+1:pathlen+flen)=dummy(1:flen)

        return

    end subroutine addpath2
    ! ============================================================================
    !> Subroutine addpath3()
    !! Subroutine addpath3 is used to add a path onto any filenames.
    !! \param fname
    !<
    ! ============================================================================
    subroutine addpath3( fname )

        integer:: pathlen,flen
        character(len=128):: fname, dname, dummy

        ! for files to be read, they'll be in the session path
        dname=' '
        call noblnk(outdir3)
        pathlen = len_trim(outdir3)
        dname(1:pathlen)=outdir3
        dname(pathlen+1:pathlen+1) = '/'
        pathlen = pathlen+1
        call noblnk(fname)
        flen = len_trim(fname)
        dummy=' '
        dummy=fname
        fname=' '
        fname(1:pathlen)=dname(1:pathlen)
        fname(pathlen+1:pathlen+flen)=dummy(1:flen)

        return

    end subroutine addpath3
    ! ============================================================================
    !> Subroutine addpath4()
    !! Subroutine addpath4 is used to add a path onto any filenames.
    !! \param fname
    !<
    ! ============================================================================
    subroutine addpath4( fname )

        integer:: pathlen,flen
        character(len=128):: fname, dname, dummy

        ! for files to be read, they'll be in the session path
        dname=' '
        call noblnk(outdir4)
        pathlen = len_trim(outdir4)
        dname(1:pathlen)=outdir4
        dname(pathlen+1:pathlen+1) = '/'
        pathlen = pathlen+1
        call noblnk(fname)
        flen = len_trim(fname)
        dummy=' '
        dummy=fname
        fname=' '
        fname(1:pathlen)=dname(1:pathlen)
        fname(pathlen+1:pathlen+flen)=dummy(1:flen)

        return

    end subroutine addpath4
    ! ============================================================================
    !> Function rmpath()
    !! Function rmpath is used to remove a path onto any filenames.
    !! \param pname
    !<
    ! ============================================================================
    function rmpath( pname ) result( fname )

        integer:: pathlen, plen
        character(len=128):: fname, pname, dname

        ! for files to be read, they'll be in the session path
        dname=' '
        call noblnk(outdir1)
        pathlen = len_trim(outdir1)
        dname(1:pathlen)=outdir1
        dname(pathlen+1:pathlen+1) = '/'
        pathlen = pathlen+1
        call noblnk(pname)
        plen = len_trim(pname)
        fname=' '
        fname(1:plen - pathlen)=pname(pathlen+1:plen)
        call noblnk( fname )

    end function rmpath
    ! ============================================================================
    !> Subroutine append_str()
    !! Subroutine append_str is used to append two strings.
    !! \param stg1,stg2
    !<
    ! ============================================================================
    subroutine append_str( stg1,stg2 )

        integer :: l1, l2
        character(len=128) :: stg1, stg2

        l1 = len_trim(stg1)
        l2 = len_trim(stg2)
        stg1( l1+1:l1+l2 ) = stg2
        call noblnk(stg1)

        return

    end subroutine append_str
    ! ============================================================================
    !> Subroutine append_nb()
    !! Subroutine append_nb is used to append a string with an integer
    !! \param stg1,i
    !<
    ! ============================================================================
    subroutine append_nb( stg1, i )

        integer :: l1, l2, i
        character(len=128) :: stg1, stg2

        l1 = len_trim(stg1)
        write(stg2,'(i10)')i
        call noblnk(stg2)
        l2 = len_trim(stg2)
        stg1(l1+1:l1+l2) = stg2
        call noblnk(stg1)

        return

    end subroutine append_nb
    ! ============================================================================
    !> Subroutine append_nbreal()
    !! Subroutine append_nbreal is used to append a string with a real number
    !! \param stg1,n
    !<
    ! ============================================================================
    subroutine append_nbreal( stg1, n )

        real( tkind ) :: n
        integer :: l1, l2
        character(len=128) :: stg1, stg2

        l1 = len_trim(stg1)
        write(stg2,'(f30.4)')n
        call noblnk(stg2)
        l2 = len_trim(stg2)
        stg1(l1+2:l1+l2+1) = stg2

        return

    end subroutine append_nbreal
    ! ============================================================================
    !> Subroutine append_nb2()
    !! Subroutine append_nb2 is used to append a string with an integer number followed by a space
    !! \param stg1,i
    !<
    ! ============================================================================
    subroutine append_nb2( stg1, i )

        integer :: l1, l2, i
        character(len=128) :: stg1, stg2

        l1 = len_trim(stg1)
        write(stg2,'(i10)')i
        call noblnk(stg2)
        l2 = len_trim(stg2)
        stg1(l1+2:l1+l2+1) = stg2

        return

    end subroutine append_nb2
    ! ============================================================================
    !> Subroutine append_zero()
    !! Subroutine append_zero is used to append a string with a number and fill it with 0
    !! in order to have a 4 digit integer
    !! \param stg1,i
    !<
    ! ============================================================================
    subroutine append_zero( stg1, i )

        integer :: l1, l2, i
        character(len=128) :: stg1, stg2, str

        l2 = len_trim(stg1)
        write(stg2,'(i10)')i
        call noblnk(stg2)
        l1 = len_trim(stg2)
        str = ''
        if( l1 == 1 )then
            str(1:3) = '000'
            call append_str(str,stg2)
        elseif( l1 == 2 )then
            str(1:2) = '00'
            call append_str(str,stg2)
        elseif( l1 == 3 )then
            str(1:1) = '0'
            call append_str(str,stg2)
        endif
        l1 = len_trim(str)
        stg1(l2+1:l2+l1) = str
        call noblnk(stg1)

        return

    end subroutine append_zero
    ! ============================================================================
    !> Subroutine read_input_to_flatten()
    !! Subroutine read_input_to_flatten is used to flatten the XmL input file for LECODE.
    !<
    ! ============================================================================
    subroutine read_input_to_flatten

        logical :: find
        integer :: len1, iu11, iu12

        ! Create flatten input file name
        len1 = len( trim( fin_noflat ) )

        finput = ' '
        finput = fin_noflat( 1:len1-4 )
        finput( len1 - 4 + 1: len1-4+12 ) = '-flatten.xml'

        ! Create temporary input file
        temp_file = 'temporary_input.xml'
        temp_input = 'temp_input.xml'

        call copy_file( fin_noflat, temp_input )

        find = .true.

        ! Check for XInclude statement
        do while( find )
            call check_xinclude_statement( find )
        enddo

        ! Create the flatten file
        call copy_file( temp_input, finput )

        ! Clean working directory from temporary file
        iu11 = 229
        open( iu11, file=temp_input, form='unformatted' )
        close( iu11, status='delete' )

        iu12 = 227
        open( iu12, file=temp_file, form='unformatted' )
        close( iu12, status='delete' )

        return

    end subroutine read_input_to_flatten
    ! ============================================================================
    !> Subroutine copy_file()
    !! Subroutine copy_file is used to copy a file to another one.
    !<
    ! ============================================================================
    subroutine copy_file( f1, f2 )

        integer :: iu1, iu2, io, io2

        character(len=128) :: f1, f2
        character(len=256) :: txt

        iu1 = 231
        open(iu1,file=f1,status="old",action="read",iostat=io)
        rewind( iu1 )
        iu2 = 241
        open(iu2, file=f2,status="replace",action="write",iostat=io2)
        rewind( iu2 )

110     read( iu1, '(a256)',iostat=io, end=120 ) txt
        write( iu2,'(A)') trim( txt )
        goto 110
120 continue

    close( iu1 )
    close( iu2 )

    return

end subroutine copy_file
! ============================================================================
!> Subroutine check_xinclude_statement()
!! Subroutine check_xinclude_statement is used to find within the input file the
!! XInclude statement.
!<
! ============================================================================
subroutine check_xinclude_statement( find )

    integer :: iut, iuf, iuc, ios, ios2, ios3, position, p
    logical :: find, found

    character( len = 256 ) :: line, line2

    find = .false.
    p = 0

    ! Open file
    inquire( file=temp_input, exist=found )
    if( .not. found )then
        print*,'Problem: the input file cannot be found: ',trim(temp_input)
        stop
    endif

    iut = 123
    open(iut,file=temp_input,status="old",action="read",iostat=ios)
    rewind( iut )
    iuf = 124
    open(iuf,file=temp_file,status="replace",action="write",iostat=ios2)
    rewind( iuf )

10  read( iut, '(a256)',iostat=ios, end=20 ) line
    p = p + 1
    position = index( line, 'xi:include' )
    if( position /= 0 )then

        find = .true.
        call get_pointed_file_name( line )
        inquire( file=included_file, exist=found )
        if( .not. found )then
            print*,'Problem: the include file cannot be found: ',trim(included_file)
            stop
        endif
        iuc = 125
        open(iuc,file=included_file,status="old",action="read",iostat=ios3)
        rewind( iuc )
11      read( iuc, '(a256)',iostat=ios3, end=21 ) line2
        write( iuf,'(A)') trim( line2 )
        goto 11
21  continue
    close( iuc )

else
    write( iuf,'(A)') trim( line )
endif
goto 10
20 continue

   close( iut )
   close( iuf )

   ! Copy new changes temporary
   call copy_file( temp_file, temp_input )

   return

   end subroutine check_xinclude_statement
   ! ============================================================================
   !> Subroutine get_pointed_file_name()
   !! Subroutine get_pointed_file_name extract from a given line the name of the XInclude
   !! file.
   !<
   ! ============================================================================
   subroutine get_pointed_file_name( line )

       integer :: position, position2, len

       character( len = 256 ) :: line
       position = index( line, 'href=' )
       position2 = index( line, '.xml' )
       if( position == 0 )then
           print*,'Problem: the picked line does not have any href element: ',trim( line )
           print*,'All included files should have a href element.'
           stop
       endif
       if( position2 == 0 )then
           print*,'Problem: the picked line does not have any xml file: ',trim( line )
           print*,'All included files should have a xml extension.'
           stop
       endif

       included_file = ' '
       len = ( position2 + 3 ) - ( position + 6 ) + 1
       included_file( 1 : len ) = line( position+6:position2+3 )

       return

   end subroutine get_pointed_file_name
   ! ============================================================================
   !> Subroutine get_cpp_file_name_extensions()
   !! Subroutine get_cpp_file_name_extensions creates CPP file names.
   !<
   ! ============================================================================
   subroutine get_cpp_file_name_extensions

       ftin = 'TIN.node'
       call noblnk( ftin )
       call addpath2( ftin )
       ftinCPP = ''
       ftinCPP = ftin
       call noblnk( ftinCPP )
       ftinCPP( len( ftinCPP ) : len( ftinCPP ) ) = CHAR(0)

       return

   end subroutine get_cpp_file_name_extensions
   ! ============================================================================

end module file_data
! ============================================================================
