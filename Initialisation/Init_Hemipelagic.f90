! ============================================================================
! Name        : Init_Hemipelagic.f90
! Author      : tristan salles
! Copyright (C) 2014
! ============================================================================
!> \file Init_Hemipelagic.f90
!
! Description :  This module defines hemipelagic growth rate.
!
!<
! ============================================================================
!> Module hemipelagite
!<
module hemipelagite

    use parallel
    use file_data
    use time_data
    use forces_data

    implicit none

contains

    ! ============================================================================
    !> Subroutine read_hemipelagite_file()
    !! Reads formatted data of hemipelagic sedimentation rates.
    !<
    ! ============================================================================
    subroutine read_hemipelagite_file

        ! Parameters Declaration
        logical :: found
        integer :: i, k, i2
        integer :: iu, ios

        character(len=128) :: text

        ! Find and open the sea level file
        iu = 45 + iam
        inquire( file=fhemi, exist=found )
        if(found)then
            open(iu,file=fhemi,status="old",action="read",iostat=ios)
            rewind(iu)
        else
            print*,'Cannot find hemipelagic file.',iam
            call mpi_finalize( ierr )
            stop
        endif

        ! Find the number of hemipelagic sedimentation rate
        k = 0
        do
            read(iu,'(a128)',iostat=ios,end=30) text
            k = k + 1
        enddo
30      continue
        rewind(iu)

        ! Allocate hemipelagic array
        contourhemi_nb = k
        if( allocated( hemipelagic ) ) deallocate( hemipelagic )
        allocate( hemipelagic( 2, contourhemi_nb ) )

        ! Read hemipelagic rate for each contour
        do k = 1, contourhemi_nb
            read(iu, '(a128)' ,iostat=ios, end=60 )text
            i2=len_trim(text)
            if( text( i2:i2 ) == char(13) )then
                i2=i2-1
            endif
            read(text(1:i2),*,iostat=ios)( hemipelagic( i, k ), i = 1, 2 )
            if( ios /= 0 )then
                print*,'hemipelagic file reading problem',iam
                call mpi_finalize( ierr )
                stop
            endif
            if( k > 1 )then
                if( hemipelagic( 1, k-1 ) <= hemipelagic( 1, k ) )then
                    print*,'hemipelagic should be sorted in decreasing contour value',iam
                    call mpi_finalize( ierr )
                    stop
                endif
            endif
        enddo

60      close(iu)

        return

    end subroutine read_hemipelagite_file
    ! ============================================================================
    !> Subroutine hemipelagites()
    !! Performs the sedimentation based on hemipelagic rate.
    !<
    ! ============================================================================
    subroutine hemipelagites

        integer :: k, i, lb, layerID

        real( tkind ) :: hemi_rate1, hemi_contour1, hemi_rate2, hemi_contour2
        real( tkind ) :: ztop, hemi_contour, hemi_rate, hemi_time, hemi_dep


        do k = 1, lnbPts

            ! Sea bed elevation
            ztop = tCoordZ( k )

            ! As long as we are under water
            if( ztop < gsea%actual_sea )then
                i = 1
                ! Find the corresponding hemipelagic rate for the
                ! considered bathymetry.
                pick_hemipelagic_contour: do
                    if( i > contourhemi_nb )exit
                    if( ztop - gsea%actual_sea > hemipelagic( 1, i ) ) &
                        exit pick_hemipelagic_contour
                    hemi_contour1 = hemipelagic( 1, i )
                    hemi_rate1 = hemipelagic( 2, i )
                    i = i + 1
                enddo pick_hemipelagic_contour

                ! Update hemipelagic rate based on contour number
                hemi_rate = 0.0_8
                if( i <= contourhemi_nb )then
                    hemi_contour2 = hemipelagic( 1, i )
                    hemi_rate2 = hemipelagic( 2, i )
                    hemi_contour = ztop - gsea%actual_sea
                    ! Interpolate hemipelagic rate based on distance between current
                    ! depth and contours
                    call hemi_interpolation( hemi_contour, hemi_rate1, hemi_contour1, &
                        hemi_rate2, hemi_contour2, hemi_rate )
                else
                    hemi_rate = hemipelagic( 2, contourhemi_nb )
                endif

                ! Elapsed time since last call to hemipelagite function
                hemi_time = tnow - last_hemi

                ! Get the deposit thickness from hemipelagic rate ( m/yr )
                hemi_dep = hemi_rate * hemi_time

                ! In case there is sufficient deposit, add a new layer
                if( hemi_dep > 0.001_8 )then

                    lb = v_nblays( k )
                    layerID = v_LaysID(  k, lb )

                    ! Create a new layer in case a current layer doesn't exist
                    if( layID > layerID )then
                        lb = lb + 1
                        stratLays( lb, k, 1:totgrn ) = 0.0_8
                        v_nblays( k ) = lb
                        v_LaysID(  k ,  lb  ) = layID
                        lCoordZ( k  ) = lCoordZ( k  ) + hemi_dep
                        stratLays(  lb, k, hemi_mat  ) = hemi_dep
                        if( gporo%compaction ) porosityLays( lb, k, 1:totgrn ) = 0.0_8
                        if( gporo%compaction ) &
                            porosityLays( lb, k, hemi_mat ) = porosity( hemi_mat, 1 )

                    ! Otherwise add on an existing layer
                    else
                        lCoordZ( k  ) = lCoordZ( k  ) + hemi_dep
                        stratLays(  lb, k, hemi_mat  ) = hemi_dep + stratLays(  lb, k, hemi_mat  )
                        if( gporo%compaction ) &
                            porosityLays( lb, k, hemi_mat ) = porosity( hemi_mat, 1 )
                    endif

                endif

            endif
        enddo

        ! Update last call to hemipelagite function to actual time
        last_hemi = tnow

        return

    end subroutine hemipelagites
    ! ============================================================================
    !> Subroutine hemi_interpolation()
    !! Performs the hemipelagic sedimentation rate interpolation function.
    !!
    !! \arg \c use precision_data
    !!
    !! \param hemi_contour, hemi_rate1, hemi_contour1, hemi_rate2, hemi_contour2, hemi_rate
    !<
    ! ============================================================================
    subroutine hemi_interpolation( hemi_contour, hemi_rate1, hemi_contour1, hemi_rate2, &
        hemi_contour2, hemi_rate )

        real( tkind ) :: hemi_contour, hemi_rate1, hemi_contour1, hemi_rate2, hemi_contour2, hemi_rate

        hemi_rate = ( hemi_contour - hemi_contour1 ) / ( hemi_contour2 - hemi_contour1 )
        hemi_rate = hemi_rate * ( hemi_rate2 - hemi_rate1 )
        hemi_rate = hemi_rate + hemi_rate1

        return

    end subroutine hemi_interpolation
    ! ============================================================================

end module hemipelagite
! ============================================================================
