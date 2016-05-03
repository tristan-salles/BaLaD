! ============================================================================
! Name        : FWsEroDep.f90
! Author      : tristan salles
! Copyright (C) 2014
!============================================================================
!> \file FWsEroDep.f90
!!
!! File FWsEroDep determines the new concentration and deposit layer caracteristics
!! due to river / turbidity current erosion and deposition.
!!
!<
! ============================================================================
module stream_erodep

    use stream
    use parallel
    use morpho
    use flow_data
    use time_data
    use mesh_data
    use forces_data

    implicit none

contains

    ! ============================================================================
    !> Subroutine apply_erosion_deposition_rules
    !! Subroutine apply_erosion_deposition_rules change the sediment load of the flow
    !! walkers by eroding or depositing to the grid.
    !!
    !! If sediment concentration divided by sediment transportability is greater than
    !! flow's transport capacity, deposition occurs.  If sediment concentration
    !! divided by sediment transportability is less than transport capacity, erosion
    !! occurs provided that critical shear stress for material at water-sediment
    !! interface is exceeded.
    !<
    ! ============================================================================
    subroutine apply_erosion_deposition_rules

        integer :: mFcS, dFcS, idM, idD, ks, k

        integer, dimension( 4 ) :: idDEM, idMESH

        real( tkind ) :: th, dtb, regolith_h, slsed, volsed
        real( tkind ) :: deposit_fraction, eromax, changeC, depomax

        real( tkind ), dimension( 4 ) :: wght
        real( tkind ), dimension( totgrn ) :: regolith, nregolith, sedcharge
        real( tkind ), dimension( totgrn ) :: conc, conceq, erodep, poro


        ! Find DEM face containing the flow walker
        call check_DEMface_containing_point( fr_xpos, fr_ypos, dFcS )
        idDEM( 1: 4 ) = fptIDs( dFcS, 1:4 )

        ! Find stratigraphic face containing the flow walker
        call check_Meshface_containing_point( fr_xpos, fr_ypos, mFcS )
        idMESH( 1: 4 ) = lfptIDs( mFcS, 1:4 )

        ! Flow concentration
         conc = 0.0_8
         do ks = 1, totgrn
            conc( ks ) = fr_sedcharge( ks ) / fr_vol
        enddo

        call inverse_distance_weighting( idMESH, wght )

        ! Loop over the points which form the face and get the erosion/deposition thicknesses
        do k = 1, 4

            ! Get ids
            idD = idDEM( k )
            idM = idMESH( k )

            if( idD > 0 .and. wght( k ) > 0.0_8 )then

                ! Get maximum erosion thickness to prevent formation of hole
                call get_maximum_erosion_thickness( idM, eromax )

                ! Get maximum deposition thickness
                call get_maximum_deposition_thickness( idM, depomax )

                ! Relation between soil depth and slope to compute depth to bedrock
                dtb = soil_thick( idM )

                ! In case the depth to bedrock is greater than the maximum erosion change it to match
                if( dtb > eromax .or. dtb == 0.0_8 ) dtb = eromax

                ! Get the thickness of the regolith layer and its composition
                call get_regolith_composition( idM, dtb, regolith_h, regolith, poro )
                if( v_nblays( idM ) == 0 ) dtb = 0.0_8

                ! Based on regolith thickness compute the equilibrium concentration
                call get_sediment_equilibrium_concentration( regolith_h, regolith, conceq )

                ! For each sediment determine the change in concentration for the considered fw
                th = 0.0_8
                do ks = 1, totgrn
                    ! Change in concentration
                    changeC = -( conc( ks ) - conceq( ks ) ) * fr_dt * sediment( ks )%vfall / fr_h

                    ! Get the thickness of sediment eroded or deposited by the flow
                    !-------------
                    ! Note: this is purely sediments thickness so the porosity needs to be taken into account
                    ! when defining stratigraphic changes
                    !-------------
                    erodep( ks ) = -changeC * fr_vol / ( strat_dx**2 )
                    if( erodep( ks ) < 0.0_8  .and. dtb == 0.0_8 ) erodep( ks ) = 0.0_8

                    ! Apply a weighting coefficient based on distance between flow walker and node
                    erodep( ks ) = erodep( ks ) * wght( k )

                    ! In case of deposition
                    if( erodep( ks ) > 0.0_8 )then
                        ! Compute sediment thickness based on sediment porosity
                        if( gporo%compaction )then
                            th = th + erodep( ks ) / ( 1.0_8 - porosity( ks, 1 ) )
                        else
                            th = th + erodep( ks )
                        endif
                    endif
                enddo
                deposit_fraction = 1.0_8
                if( th > depomax ) deposit_fraction = depomax / th

                ! Potential erosion/deposition thickness to reach equilibrium concentration
                nregolith = regolith
                sedcharge = 0.0_8
                do ks = 1, totgrn

                    ! Apply erosion limited by sediment availability in the regolith layer
                    if( erodep( ks ) < 0.0_8 )then
                        ! In case of erosion, the thickness takes into account the
                        ! the layer porosity
                        if( regolith( ks ) * ( 1.0_8 - poro( ks ) ) > -erodep( ks ) )then
                            sedcharge( ks ) = -erodep( ks ) * strat_dx**2
                            nregolith( ks ) = nregolith( ks ) + erodep( ks ) /  ( 1.0_8 - poro( ks ) )
                        else
                            sedcharge( ks ) = regolith( ks ) * ( 1.0_8 - poro( ks ) ) * strat_dx**2
                            nregolith( ks ) = 0.0_8
                        endif

                    ! Apply deposition limited by sediment availability in the flow
                    elseif( erodep( ks ) > 0.0_8 )then
                        erodep( ks ) = deposit_fraction * erodep( ks )
                        if( gporo%compaction )then
                            ! In case of deposition, the thickness takes into account the
                            ! the layer porosity
                            if( erodep( ks ) * strat_dx**2 >= fr_sedcharge( ks )  )then
                                sedcharge( ks ) = -fr_sedcharge( ks )
                                nregolith( ks ) = nregolith( ks ) + fr_sedcharge( ks ) &
                                    / ( ( 1 - porosity( ks , 1 ) ) * strat_dx**2 )
                            elseif( erodep( ks ) * strat_dx**2 < fr_sedcharge( ks ) )then
                                sedcharge( ks ) = -erodep( ks ) * strat_dx**2
                                nregolith( ks ) = nregolith( ks ) + erodep( ks ) / ( 1.0_8 - porosity( ks , 1 ) )
                            endif
                        else
                            if( erodep( ks ) * strat_dx**2 >= fr_sedcharge( ks )  )then
                                sedcharge( ks ) = -fr_sedcharge( ks )
                                nregolith( ks ) = nregolith( ks ) + fr_sedcharge( ks ) / strat_dx**2
                            elseif( erodep( ks ) * strat_dx**2 < fr_sedcharge( ks ) )then
                                sedcharge( ks ) = -erodep( ks ) * strat_dx**2
                                nregolith( ks ) = nregolith( ks ) + erodep( ks )
                            endif
                       endif
                    endif
                enddo

                ! Define flow walker new sediment properties and flow density
                slsed = 0.0_8
                volsed = 0.0_8
                do ks = 1, totgrn
                    fr_sedcharge( ks ) = fr_sedcharge( ks ) + sedcharge( ks )
                    volsed = volsed + fr_sedcharge( ks )
                    slsed = slsed + fr_sedcharge( ks ) * sediment( ks )%density
                enddo
                fr_vol = fr_vol0 + volsed
                fr_density = slsed / fr_vol + fluid_density

                ! Update the stratigraphic layer
                call update_stratigraphic_layer( idD, idM, regolith, nregolith )

            endif
        enddo


        return

    end subroutine apply_erosion_deposition_rules
    ! ============================================================================
    !> Subroutine find_Mesh_elevation
    !! Subroutine find_Mesh_elevation is used to determine the elevation of a point in the
    !! mesh grid.
    !<
    ! ============================================================================
    subroutine find_Mesh_elevation

        integer :: Xid, Yid, fcS, id( 4 ), k

        real( tkind ) :: wght( 4 )

        ! Define the propable row and column number
        Xid = int( ( fr_xpos - lstrat_xo ) / strat_dx ) + 1
        Yid = int( ( fr_ypos - lstrat_yo ) / strat_dx ) + 1

        ! Define the face id accordingly
        if( Yid == 1 .and. Xid < lstrat_X )then
            fcS = Xid
        elseif( Yid == 1 .and. Xid == lstrat_X )then
            fcS = Xid - 1
        elseif( Xid == 1 .and. Yid < lstrat_Y )then
            fcS = ( Yid - 1 )*( lstrat_X - 1 ) + 1
        elseif( Yid == lstrat_Y .and. Xid < lstrat_X )then
            fcS = ( Yid - 2 )*( lstrat_X - 1 ) + Xid
        elseif( Yid == lstrat_Y .and. Xid == lstrat_X )then
            fcS = ( Yid - 1 )*( Xid - 1 )
        elseif( Yid < lstrat_Y .and. Xid == lstrat_X )then
            fcS = Yid * ( Xid - 1 )
        else
            fcS = ( lstrat_X - 1 )* ( Yid - 1 ) + Xid
        endif

        ! Check that picked face is containing the point
        id( 1:4 ) = lfptIDs( fcS, 1:4 )

        ! Weigth
        call inverse_distance_weighting( id, wght )

        fr_zstrat = 0.0_8
        do k = 1, 4
            fr_zstrat = fr_zstrat  + wght( k ) * lcoordZ( id( k ) )
        enddo

        return

    end subroutine find_Mesh_elevation
    ! ============================================================================
    !> Subroutine inverse_distance_weighting
    !! Subroutine inverse_distance_weighting get weight of the fw on each points
    !<
    ! ============================================================================
    subroutine inverse_distance_weighting( id, weights )

        integer :: k
        integer, dimension( 4 ) :: id

        real( tkind ) :: dis, sumdist, sumwgth
        real( tkind ), dimension( 4 ) :: weights, dist

        weights = 0.0_8
        sumdist = 0.0_8
        face_loop: do k = 1, 4
            dis = ( fr_xpos - lcoordX( id( k ) ) )**2 + ( fr_ypos - lcoordY( id( k ) ) )**2
            if( dis == 0.0_8 )then
                weights = 0.0_8
                weights( k ) = 1.0_8
                return
            endif
            dist( k ) = 1 / dis
            sumdist = sumdist + dist( k )
        enddo face_loop

        sumwgth = 0.0_8
        do k = 1, 4
            weights( k ) = dist( k ) / sumdist
            sumwgth = weights( k ) + sumwgth
        enddo

        do k = 1, 4
            weights( k ) = weights( k ) / sumwgth
        enddo

        sumwgth = 0.0_8
        do k = 1, 4
            sumwgth = weights( k ) + sumwgth
        enddo

        if( ( sumwgth - 1.0_8 ) > tor )then
            print*,'Problem when computing inverse weighting values',iam,sumwgth
            call mpi_finalize( ierr )
            stop
        endif

        return

    end subroutine inverse_distance_weighting
    ! ============================================================================
    !> Subroutine get_maximum_erosion_thickness
    !! Subroutine get_maximum_erosion_thickness get the maximum erosion thickness
    !<
    ! ============================================================================
    subroutine get_maximum_erosion_thickness( idM, eromax )

        integer :: idM, n, k

        real( tkind ) :: elev, minz, eromax

        eromax = 0.0_8

        ! Current layer elevation
        elev = lcoordZ( idM )
        minz = elev

        ! Find neighborhood minimum elevation
        do n = 1, 8
            k = lngbID( idM, n )
            if(  k > 0 ) minz = min( minz, lcoordZ( k ) )
        enddo

        ! Maximum erosion to prevent formation of holes
        eromax = 0.8_8 * ( elev - minz )
        if( eromax < 0.0_8 ) eromax = 0.0_8

        if( elev > gsea%actual_sea .and. minz < gsea%actual_sea ) eromax = 0.8_8 * ( elev - gsea%actual_sea )

        return

    end subroutine get_maximum_erosion_thickness
    ! ============================================================================
    !> Subroutine get_maximum_deposition_thickness
    !! Subroutine get_maximum_deposition_thickness get the maximum deposition thickness
    !<
    ! ============================================================================
    subroutine get_maximum_deposition_thickness( idM, depomax )

        integer :: idM, n, k

        real( tkind ) :: elev, maxz, depomax

        depomax = 0.0_8

        ! Current layer elevation
        elev = lcoordZ( idM )
        maxz = -1e6_8

        ! Find neighborhood maximum elevation
        do n = 1, 8
            k = lngbID( idM, n )
            if(  k > 0 ) maxz = max( maxz, lcoordZ( k ) )
        enddo

        ! Maximum erosion to prevent formation of holes
        depomax = ( maxz - elev ) / 10.0_8
        if( depomax < 0.0_8 ) depomax = 0.0_8
        depomax = depomax + fr_h

        return

    end subroutine get_maximum_deposition_thickness
    ! ============================================================================
    !> Subroutine get_regolith_composition
    !! Subroutine get_regolith_composition gets regolith layer composition information.
    !<
    ! ============================================================================
    subroutine get_regolith_composition( idM, dbed, regolith_h, regolith, poro  )

        integer :: idM, nlay, p, ks

        real( tkind ) :: dbed, dtb, th, cumulative, regolith_h, portion_h, frac
        real( tkind ), dimension( totgrn ) :: sedh, regolith, poro

        nlay = v_nblays( idM )
        poro = 0.0_8
        regolith_h = 0.0_8
        regolith = 0.0_8
        if( dbed > 0.0_8 .and. nlay >= 1 )then
            dtb = max( dbed, 1.e-4_8 )
            cumulative = 0.0_8

            ! Get layer number
            regolith = 0.0_8
            regolith_h = 0.0_8

            layer_loop: do p = nlay, 1, -1
                th = 0.0_8
                sedh = 0.0_8
                do ks = 1, totgrn
                    th = th + stratLays( p, idM, ks )
                    sedh( ks ) = stratLays( p, idM, ks )
                enddo
                cumulative = cumulative + th

                ! If the cumulative thickness is above
                ! the regolith thickness, take a portion of it
                if( dtb < cumulative )then

                    portion_h = dtb - regolith_h
                    frac = portion_h / th
                    sedh = 0.0_8
                    th = 0.0_8
                    do ks = 1, totgrn
                        sedh( ks ) = stratLays( p, idM, ks ) * frac
                        regolith( ks ) = regolith( ks ) + sedh( ks )
                        th = th + sedh( ks )
                        if( gporo%compaction )&
                            poro( ks ) = poro( ks ) + porosityLays( p, idM, ks ) * sedh( ks )
                    enddo
                    regolith_h = regolith_h + th
                    if( gporo%compaction )then
                        do ks = 1, totgrn
                            if( regolith( ks ) > 0.0_8 ) poro( ks ) = poro( ks ) / regolith( ks )
                        enddo
                    endif
                    exit layer_loop
                ! In case the regolith thickness is not reached,
                ! records the layer composition
                else
                    th = 0.0_8
                    do ks = 1, totgrn
                        th = th + stratLays( p, idM, ks )
                        regolith( ks ) = regolith( ks ) + stratLays( p, idM, ks )
                        if( gporo%compaction )&
                            poro( ks ) = poro( ks ) + porosityLays( p, idM, ks ) * sedh( ks )
                    enddo
                    regolith_h = regolith_h + th
                endif
            enddo layer_loop

        ! Erosion not possible but sedimentary layer exists
        elseif( nlay >= 1 )then
            th = 0.0_8
            sedh = 0.0_8
            regolith = 0.0_8
            regolith_h = 0.0_8
            do ks = 1, totgrn
                th = th + stratLays( nlay, idM, ks )
                regolith( ks ) = stratLays( nlay, idM, ks )
                if( gporo%compaction )&
                    poro( ks ) = porosityLays( nlay, idM, ks )
             enddo
             regolith_h = regolith_h + th
        ! On the basement
        else
            regolith_h = 1000.0_8
            regolith = 0.0_8
            regolith( totgrn ) = 1000.0_8

        endif

        return

    end subroutine get_regolith_composition
    ! ============================================================================
    !> Subroutine get_sediment_equilibrium_concentration
    !! Subroutine get_sediment_equilibrium_concentration computes equilibrium sediment
    !! concentration of each size class of sediment load.
    !<
    ! ============================================================================
    subroutine get_sediment_equilibrium_concentration( regolith_h, regolith, conceq )

        integer :: ks, kk, sed_nb

        real( tkind ) :: ratio, sumprop, d50, fv, regolith_h
        real( tkind ) :: wg1, wg2, idis1, idis2, sdis, swg

        real( tkind ), dimension( totgrn ) :: regolith, perc_sed, cum_perc
        real( tkind ), dimension( totgrn ) :: sflux, bflux, conceq, prop, ph, pe

        conceq = 0.0_8

        ! Get the fraction of each grain present in the top layer
        sumprop = 0.0_8
        do ks = 1, totgrn
            prop( ks ) = regolith(  ks  )
            prop( ks ) = prop( ks ) /  regolith_h
            sumprop = sumprop + prop( ks )
        enddo

        ! Ensure the sediment proportions are suming to 1.0
        ! minor numerical precision
        if( abs( sumprop - 1.0_8 ) > tor )then
            ratio = 1.0_8 / sumprop
            do ks = 1, totgrn
                prop( ks ) = prop( ks ) * ratio
            enddo
        endif

        ! Get the particles size percentage and cumulative distribution
        sed_nb = 0
        perc_sed = 0.0_8
        cum_perc = 0.0_8
        d50_loop: do ks = 1, totgrn
            perc_sed( ks ) = prop( ks ) * 100.0_8
            if( ks == 1 )then
                cum_perc( ks ) = perc_sed( ks )
            else
                cum_perc( ks ) =  cum_perc( ks - 1 ) + perc_sed( ks )
            endif
            if( cum_perc( ks ) > 50.0_8 )then
                if( ks == 1 )then
                    d50 = sediment( ks )%diameter
                    goto 28
                else
                    if( cum_perc( ks - 1 ) < 50.0_8 )then
                        sed_nb = ks - 1
                        exit d50_loop
                    endif
                 endif
            endif
        enddo d50_loop

        ! Find d50 grain size for this particular distribution
        if( sed_nb > 0 )then
            idis1 = 50.0_8 - cum_perc( sed_nb )
            idis2 = cum_perc( sed_nb + 1 ) - 50.0_8
            if( idis1 < tor )then
                d50 = sediment( sed_nb )%diameter
                goto 28
            elseif( idis2 < tor )then
                d50 = sediment( sed_nb + 1 )%diameter
                goto 28
            endif
            idis1 = 1 / idis1
            idis2 = 1 / idis2
            sdis = idis1 + idis2
            wg1 = idis1 / sdis
            wg2 = idis2 / sdis
            swg = wg1 + wg2
            wg1 = wg1 / swg
            wg2 = wg2 / swg
            d50 = wg1 * sediment( sed_nb )%diameter + wg2 * sediment( sed_nb+1 )%diameter
        else
            find_sed: do ks = 1, totgrn
                if( perc_sed( ks ) > 0.0_8 )then
                    d50 = sediment( ks )%diameter
                    exit find_sed
                endif
            enddo find_sed
        endif
        if( totgrn == 1 ) d50 = sediment( 1 )%diameter
28 continue

        ! Exposed / hidden probabilities in the top layer
        do kk = 1, totgrn
            ph( kk ) = 0.0_8
            pe( kk ) = 0.0_8
            ! Get hidden probabilities
            do ks = 1, totgrn
                ph( kk ) = ph( kk ) + prop( ks ) * sediment( ks )%diameter / &
                    ( sediment( ks )%diameter + sediment( kk )%diameter )
            enddo
            if( ph( kk ) < tor )then
                print*,'Something went wrong computing hidden probabilities.'
                stop
            endif
            ! Deduce exposed probabilities
            pe( kk ) =  1.0_8 - ph( kk )
            if( pe( kk ) < 0.0_8 )then
                ph( kk ) = 1.0_8
                pe( kk ) = 0.0_8
            endif
        enddo

        ! Suspension mode
        sflux = 0.0_8
        call suspension_sediment_transport( pe, ph, prop, sflux )

        ! Bedload mode
        bflux = 0.0_8
        call bedload_sediment_transport( pe, ph, prop, d50, bflux )

        ! Equilibrium sediment concentration
        fv = sqrt( fr_xv**2 + fr_yv**2 )

        if( fv * fr_h > 0.0_8 )then
            do ks = 1, totgrn
                conceq( ks ) = ( bflux( ks ) + sflux( ks ) ) / ( fr_h * fv )
            enddo
        else
            conceq = 0.0_8
        endif


        return

    end subroutine get_sediment_equilibrium_concentration
    ! ============================================================================
    !> Subroutine bedload_sediment_transport
    !! Subroutine bedload_sediment_transport determines bedload sediment transport.
    !! Wu et al. (2000) proposed a formula to calculate the fractional bed load transport capacity.
    !<
    ! ============================================================================
    subroutine bedload_sediment_transport( pe, ph, prop, d50, bflux )

        integer :: ks

        real( tkind ) :: d50, ntilde
        real( tkind ) :: theta, m, fv, Sf, tau_b, xnn

        real( tkind ), dimension( totgrn ) :: bflux, prop, ph, pe, tau_c, psi

        theta = 0.03_8
        m = -0.6_8

        ! Set Manning's coefficient
        xnn = manning_fct( gsea%actual_sea - fr_zstrat, fr_density )

        ntilde = d50**(1.0_8/6.0_8) / 20.0_8

        ! Bottom shear stress
        fv = sqrt( fr_xv**2 + fr_yv**2 )
        Sf = xnn**2 * fv**2 / ( fr_h**( 4.0_8 / 3.0_8 ) )
        tau_b = fluid_density * gravity * fr_h * Sf

        ! Criterion for incipient motion and bedload rate
        psi( : ) = 0.0_8
        tau_c( : ) = 0.0_8
        do ks = 1, totgrn
            tau_c( ks ) = theta * sediment( ks )%diameter * ( pe( ks ) / ph( ks ) )**m
            tau_c( ks ) = tau_c( ks ) * ( sediment( ks )%density - fluid_density ) * gravity
            if( tau_c( ks ) <= 0.0_8 )then
                bflux( ks ) = 0.0_8
            else
                if( ( ntilde / xnn )**1.5_8 * tau_b / tau_c( ks ) - 1.0_8 < 0.0_8 )then
                    psi( ks ) = 0.0_8
                    bflux( ks ) = 0.0_8
                else
                    psi( ks ) = 0.0053_8 * (  ( ntilde / xnn )**1.5_8 * tau_b / tau_c( ks ) - 1.0_8 )**2.2_8
                    bflux( ks ) = psi( ks ) * prop( ks ) * sqrt( sediment( ks )%diameter**3 * gravity * &
                        ( sediment( ks )%density  / fr_density - 1.0_8 ) )
                endif
            endif
        enddo

        return

    end subroutine bedload_sediment_transport
    ! ============================================================================
    !> Subroutine suspension_sediment_transport
    !! Subroutine suspension_sediment_transport determines suspension sediment transport capacity.
    !! Wu et al (2000) proposed a formula to calculate the fractional suspended load transport capacity.
    !<
    ! ============================================================================
    subroutine suspension_sediment_transport( pe, ph, prop, sflux )

        integer :: ks

        real( tkind ) :: theta, m, fv, Sf, tau_b, xnn, fall

        real( tkind ), dimension( totgrn ) :: sflux, prop, ph, pe, tau_c, psi

        theta = 0.03_8
        m = -0.6_8

        ! Set Manning's coefficient
        xnn = manning_fct( gsea%actual_sea - fr_zstrat, fr_density )

        ! Bottom shear stress
        fv = sqrt( fr_xv**2 + fr_yv**2 )
        Sf = xnn**2 * fv**2 / ( fr_h**( 4.0_8 / 3.0_8 ) )
        tau_b = fluid_density * gravity * fr_h * Sf

        ! Criterion for incipient motion and suspended load rate
        psi( : ) = 0.0_8
        tau_c( : ) = 0.0_8
        do ks = 1, totgrn
            tau_c( ks ) = theta * sediment( ks )%diameter
            tau_c( ks ) = tau_c( ks ) * ( pe( ks ) / ph( ks ) )**m
            tau_c( ks ) = tau_c( ks ) * ( sediment( ks )%density - fluid_density ) * gravity

            if( tau_c( ks ) > 0.0_8 )then
                if( tau_b / tau_c( ks ) - 1 < 0.0_8 )then
                    psi( ks ) = 0.0_8
                    sflux( ks ) = 0.0_8
                else
                    fall = sediment( ks )%vfall
                    psi( ks ) = 2.62e-5_8 * ( ( tau_b / tau_c( ks ) - 1.0_8 ) * fv / fall )**1.74_8
                    sflux( ks ) = psi( ks ) * prop( ks ) * sqrt( sediment( ks )%diameter**3 * gravity * &
                            ( sediment( ks )%density  / fr_density - 1.0_8 ) )
                endif
            else
                sflux( ks ) = 0.0_8
            endif
        enddo

        return

    end subroutine suspension_sediment_transport
    ! ============================================================================
    !> Subroutine update_stratigraphic_layer
    !! Subroutine update_stratigraphic_layer updates layer composition information.
    !<
    ! ============================================================================
    subroutine update_stratigraphic_layer( idD, idM, regolith, nregolith  )

        integer :: idM, idD, nlay, p, ks, k, minlay, play, newlay

        real( tkind ) :: th, hero
        real( tkind ), dimension( totgrn ) :: nregolith, regolith

        ! In case of erosion
        nlay = v_nblays( idM )
        minlay = nlay
        hero = 0.0_8

        do ks = 1, totgrn
            th = nregolith( ks ) - regolith( ks )

            if( th < 0.0_8 )then

                stratal_loop: do p = nlay, 1, -1

                    minlay = min( p, minlay )

                    ! In case there is sufficient sediment to erode take only what is required
                    if( -th <= stratLays( p, idM, ks ) )then
                        stratLays( p, idM, ks ) = stratLays( p, idM, ks ) + th
                        hero = hero - th
                        if( stratLays( p, idM, ks ) < tor )then
                            stratLays( p, idM, ks ) = 0.0_8
                            if( gporo%compaction ) porosityLays( p, idM, ks ) = 0.0_8
                        endif
                        exit stratal_loop
                    ! Otherwise grab all sediment present in the layer
                    else
                        th = th + stratLays( p, idM, ks )
                        hero = hero + stratLays( p, idM, ks )
                        stratLays( p, idM, ks ) = 0.0_8
                        if( gporo%compaction ) porosityLays( p, idM, ks ) = 0.0_8
                    endif

                enddo stratal_loop

            endif
        enddo

        ! Update layers ID
        newlay = nlay
        do p = minlay, nlay
            th = 0.0_8
            do ks = 1, totgrn
                th = th + stratLays( p, idM, ks )
            enddo
            if( th <= tor )then
                hero = hero + th
                stratLays( p, idM, 1:totgrn ) = 0.0_8
                newlay = newlay - 1
            endif
        enddo

        ! In case some layers have been entirely eroded
        if( newlay < nlay )then
           v_nblays( idM ) = newlay
           play = nlay
           layer_loop: do p = minlay, nlay
                th = 0.0_8
                do ks = 1, totgrn
                    th = th + stratLays( p, idM, ks )
                enddo
                if( play < newlay )then
                    print*,'Something went wrong went updating sedimentary layer after erosion',iam
                    call mpi_finalize( ierr )
                    stop
                endif
                if( v_LaysID( idM, p ) > 0 )then
                    if( th == 0.0_8 )then
                        do k = p, nlay - 1
                            v_LaysID( idM, k ) = v_LaysID( idM, k + 1 )
                            do ks = 1, totgrn
                                stratLays( k, idM, ks ) = stratLays( k + 1, idM, ks )
                            enddo
                        enddo
                        v_LaysID( idM, play ) = 0
                        stratLays( play, idM, 1:totgrn ) = 0.0_8
                        play = play - 1
                    endif
                endif
            enddo layer_loop
        endif

        ! Update elevation
        lcoordZ( idM ) = lcoordZ( idM ) - hero
        coordZ( idD ) = lcoordZ( idM )

        ! Update depth to bedrock
        soil_thick( idM ) = soil_thick( idM ) - hero
        if( soil_thick( idM ) < 0.0_8 ) soil_thick( idM ) = 0.0_8
        if( soil_thick( idM ) > dtb_marine .and. lcoordZ( idM ) < gsea%actual_sea ) &
            soil_thick( idM ) = dtb_marine
        if( soil_thick( idM ) > dtb_aerial .and. lcoordZ( idM ) >= gsea%actual_sea ) &
            soil_thick( idM ) = dtb_aerial

        ! In case of deposition
        nlay = v_nblays( idM )
        th = 0.0_8
        do ks = 1, totgrn
            if( nregolith( ks ) > regolith( ks ) )then
                if( layID == v_LaysID( idM, nlay ) )then
                    stratLays( nlay, idM, ks )  = stratLays( nlay, idM, ks ) + nregolith( ks ) - regolith( ks )
                    if( gporo%compaction ) porosityLays( nlay, idM, ks ) = 0.0_8
                    if( gporo%compaction .and. stratLays( nlay, idM, ks ) > 0.0_8 ) &
                        porosityLays( nlay, idM, ks ) = porosity( ks, 1 )
                else
                    nlay = nlay + 1
                    v_nblays( idM ) = nlay
                    v_LaysID( idM, nlay ) = layID
                    stratLays( nlay, idM, 1:totgrn ) = 0.0_8
                    if( gporo%compaction ) porosityLays( nlay, idM, ks ) = 0.0_8
                    stratLays( nlay, idM, ks ) = nregolith( ks ) - regolith( ks )
                    if( gporo%compaction .and. stratLays( nlay, idM, ks ) > 0.0_8 ) &
                        porosityLays( nlay, idM, ks ) = porosity( ks, 1 )
                endif
                th = th + nregolith( ks ) - regolith( ks )
            endif
        enddo

        ! Update elevation
        lcoordZ( idM ) = lcoordZ( idM ) + th
        coordZ( idD ) = lcoordZ( idM )

        ! Update depth to bedrock
        soil_thick( idM ) = soil_thick( idM ) + th
        if( soil_thick( idM ) > dtb_marine .and. lcoordZ( idM ) < gsea%actual_sea ) &
            soil_thick( idM ) = dtb_marine
        if( soil_thick( idM ) > dtb_aerial .and. lcoordZ( idM ) >= gsea%actual_sea ) &
            soil_thick( idM ) = dtb_aerial

        return

    end subroutine update_stratigraphic_layer
    ! ============================================================================
    !> Function manning_fct
    !! Function manning sets manning's coefficient.
    !! When flow walker's location is in the sea and the density of fresh water
    !! including suspended sediment is less than density of sea water, hypopycnal
    !! flow occurs, and manning's coefficient for hypopycnal flow is applied.
    !! when flow walker is located in the sea and the density of fresh water,
    !! including suspended sediment, is greater than density of sea water, then
    !! hyperpycnal flow occurs, and manning's coefficient for hyperpycnal flow is
    !! employed.
    !! \param depth, den
    !<
    ! ============================================================================
    function manning_fct(depth, den) result( manning )

        real( tkind ) :: manning, depth, den

        if( depth > 0.0_8 )then
            if( den > sea_density )then
                if( depth > 0.5_8 )then
                    manning = manning_hyper
                else
                    manning = manning_open - ( manning_open - manning_hyper ) * depth * 2.0_8
                endif
            elseif( den <= sea_density )then
                manning = manning_hypo
            endif
        else
            manning = manning_open
        endif

        return

    end function manning_fct
    ! ============================================================================
    
end module stream_erodep
