! ============================================================================
! Name        : Lag_interpol.f90
! Author      : Tristan Salles
! Copyright (C) 2014
! ============================================================================
!> \file  Lag_interpol.f90
!!
!! Description : This module performs Lagragian Polynomial Interpolation. The formula is
!! f(x) = sum_i=1^M w_i fa_i
!! fa_i is the function to interpolate, assumed to be an array
!! sa_i are the collocation points where the function is known
!! w_i  are the interpolation weights which depend on the relative distance
!!      of the  point x to the collocation points sa_i.
!
!! The weights are such that if x coincides with one of the points sa,
!! sa_j say,then w_j = 1 and w_i=0 for i/=j. This guarantees that
!! f(x)=fa_j.
!!
!! The formula for the weights is
!!                               (x-sa(k))
!! w_i = product_{k=1,k/=i}^M --------------
!!                             (sa(i)-sa(k))
!! In expanded form the formula looks as follows
!!         (x-sa(1))     (x-sa(2))   ...   (x-sa(i-1))     (x-sa(i+1))   ...   (x-sa(M-1))     (x-sa(M))
!! w_i = ------------- ------------- --- --------------- --------------- --- --------------- --------------
!!       (sa(i)-sa(1)) (sa(i)-sa(2)) ... (sa(i)-sa(i-1)) (sa(i)-sa(i+1)) ... (sa(i)-sa(M-1)) (sa(i)-sa(M))
!!
!! Lagrangian polynomial interpolation works very well if the point x
!! is centered in the interpolation stencil.  If it is off-centered and
!! the number of collocation points is large, polynomial interpolation
!! on an equally spaced grid is bad. It is imperative then to center the
!! stencil with respect to the point of interest as much as possible. It
!! is also important to use a subset of all available points.Sequential DEM preprocessing algorithm
!!
!! The module here contains several functions and subroutines
!!
!! StencilLocator:
!!     Finds the interpolation stencil such that the interpolation is as
!!     centered as possible This stencil is identified by its beginning
!!     index and the number of points in it.
!! LagrangeInterp:
!!     Interpolates the array fa at a point x in a single step
!! LagrangeInterp1:
!!     Computes the weights of the interpolation. It is useful if several
!!     functions need to be interpolated at the same point The duplicate
!!     evaluation of w_i is then avoided.
!! LagrangeInterp2:
!!     Applies the weights to the function computed previously.
!! LagrangeInterp2D:
!!     interpolates in 2-dimensions.
!!
!! Reference:
!! http://www.rsmas.miami.edu/personal/miskandarani/Courses/MSC321
!!
!<
! ============================================================================
module interpol

  use parallel
  use mesh_data
  use flow_data

  implicit none

  ! max stencil width, polynomial degree (minterp-1)
  integer, parameter, private :: minterp=4

  ! right stencil-offset
  integer, parameter, private :: moffr=minterp/2

  ! left stencil-offset
  integer, parameter, private :: moffl=moffr-mod(minterp-1,2)

  ! number of stencil data
  integer, parameter :: nstencil = 4

  ! grid dimension
  integer, dimension( 2 ) :: grid_dim, grid_dim2
  ! grid dx
  real( tkind ), dimension( 2 ) :: grid_dx
  ! grid minimum
  real( tkind ), dimension( 2 ) :: grid_min, grid_min2
  ! grid fields
  real( tkind ), dimension( : ), allocatable  :: x_int, y_int
  ! grid fields
  real( tkind ), dimension( : ), allocatable  :: x_int2, y_int2
  ! interpolate fields
  real( tkind ), dimension( :, : ), allocatable  :: vx_int, vy_int, z_int, h_int, w_int, fq_int
  ! updated velocity
  real( tkind ), dimension( 2)  :: Veval

contains

  ! ============================================================================
  !> Subroutine StencilLocator
  !!
  !! This subroutine finds the subset of points that most effectively centers
  !! the interpolation stencil with respect to the interpolation
  !! point.
  !<
  ! ============================================================================
  subroutine StencilLocator( istencil, xi, xmin, dx, M )

    ! number of grid points
    integer, intent(in) :: M
    ! interpolation point
    real( tkind ), intent(in) :: xi
    ! leftmost end of grid and grid spacing
    real( tkind ), intent(in) :: xmin,dx
    ! stencil info
    integer, intent(out) :: istencil(nstencil)

    integer :: ic,ie,ib,ninterp

    ! Identify Cell Location
    ! cell where xa(ic) <= xi <= xa(ic+1)
    ic = int( (xi-xmin)/dx ) + 1

    ! Stencil identification version 1
    ! interpolation stencil shrinks near the grid boundaries
    ! index of interpolation stencil start, must be >= 1
    ! ib = max(1,ic-moffl)
    ! index of interpolation stencil end, must be <= M
    ! ie = min(M,ic+moffr)

    ! Stencil identification version 2
    ! interpolation stencil does not shrink near the grid boundaries, always maintained at minterp
    ! index of interpolation stencil start, must be >= 1
    ib = max(1,ic-moffl)
    ! index of interpolation stencil end, must be <= M
    ie = ib+minterp-1
    if (ie > M) then
       ie = M
       ib = ie - minterp + 1
    endif

    ! number of points involved in interpolation
    ninterp = (ie-ib+1)

    !  Bundle data in an array for efficiency
    istencil(1) = ib
    istencil(2) = ie
    istencil(3) = ninterp
    istencil(4) = ic

    return

  end subroutine StencilLocator
  ! ============================================================================
  !> Function LagrangeInterp
  !!
  !! Compute an approximation to the function passing through the
  !! M points (sa,fa) at x. Single step function.
  !<
  ! ============================================================================
  real( tkind ) function LagrangeInterp(x,sa,fa,M)

    ! number of points where function is known
    integer, intent(in) :: M
    ! collocation points
    real( tkind ), intent(in) :: sa(M)
    ! Function to be interpolated
    real( tkind ), intent(in) :: fa(M)
    ! point where we need to do the interpolation
    real( tkind ), intent(in) :: x

    integer :: k,i
    real( tkind )  :: weight,f

    f = 0.d0
    do i = 1,M
       weight = 1.d0;
       do k = 1,i-1
          weight = weight * (x-sa(k))/(sa(i)-sa(k));
       enddo
       do k = i+1,M
          weight = weight * (x-sa(k))/(sa(i)-sa(k));
       enddo
       f = f + weight*fa(i)
    enddo
    LagrangeInterp = f

    return

  end function LagrangeInterp
  ! ============================================================================
  !> Subroutine LagrangeInterp1
  !!
  !! This subroutine computes only the weights w_i needed for the interpolation.
  !! It is useful if we need to interpolate several functions fa at the
  !! same point x.
  !<
  ! ============================================================================
  subroutine LagrangeInterp1(weight,x,sa,M)

    ! number of points where function is known
    integer, intent(in) :: M
    ! collocation points
    real( tkind ), intent(in) :: sa(M)
    ! point where we need to do the interpolation
    real( tkind ), intent(in) :: x
    ! weights of interpolation
    real( tkind ), intent(out):: weight(M)

    integer :: k,i

    do i = 1,M
       weight(i) = 1.d0;
       do k = 1,i-1
          weight(i) = weight(i) * (x-sa(k))/(sa(i)-sa(k));
       enddo
       do k = i+1,M
          weight(i) = weight(i) * (x-sa(k))/(sa(i)-sa(k));
       enddo
    enddo

    return

  end subroutine LagrangeInterp1
  ! ============================================================================
  !> Function LagrangeInterp2
  !!
  !! Apply the pre-computed weights, obtained from a previous call to LagrangeInterp1,
  !! to the interpolation formula. This subroutine performs the sum only.
  !<
  ! ============================================================================
  real( tkind ) function LagrangeInterp2(fa,weight,M)

    ! number of points where function is known
    integer, intent(in) :: M
    ! weights of interpolation
    real( tkind ), intent(in) :: weight(M)
    ! Function to be interpolated
    real( tkind ), intent(in) :: fa(M)

    ! approximated output
    real( tkind )  :: f

    integer :: i

    f = 0.d0
    do i = 1,M
       f = f + weight(i)*fa(i)
    enddo
    LagrangeInterp2 = f

    return

  end function LagrangeInterp2
  ! ============================================================================
  !> Subroutine LagrangeInterp2D
  !!
  !! This subroutine computes 2D Lagrange interpolation by tensor products of 1D interpolation.
  !<
  ! ============================================================================
  subroutine LagrangeInterp2D(fi,f2d,xi,xa,ya,xmin,dx,M)

    integer, intent(in) :: M(2)
    real( tkind ), intent(in)  :: xmin(2),dx(2)
    real( tkind ), intent(in)  :: xa(M(1)),ya(M(2))
    real( tkind ), intent(in)  :: xi(2)
    real( tkind ), intent(in)  :: f2d(M(1),M(2))
    real( tkind ), intent(out)  :: fi

    integer, parameter :: ndim=2
    integer :: id,j,ji
    integer :: ist(nstencil,ndim), ib, ninterp
    real( tkind )  :: weight(minterp), f1d(minterp)

    ! locate multi-dimensional interpolation stencil
    do id = 1,ndim
       call StencilLocator(ist(1,id),xi(id),xmin(id),dx(id), M(id))
    enddo

    id = 1; ib=ist(1,id); ninterp=ist(3,id)

    ! 2-step interpolation: weight computations
    call LagrangeInterp1(weight,xi(id),xa(ib),ninterp)
    do j = ist(1,2),ist(2,2)
       ji = j - ist(1,2) + 1
       ! 2-step interpolation: weight application
       f1d(ji) = LagrangeInterp2(f2d(ib,j),weight,ninterp)
    enddo

    id = 2; ib=ist(1,id); ninterp=ist(3,id)
    ! 2-step interpolation: weight computations
    call LagrangeInterp1(weight,xi(id),ya(ib),ninterp)

    ! 2-step interpolation: b. weight application again
    fi = LagrangeInterp2(f1d,weight,ninterp)

    return

  end subroutine LagrangeInterp2D
  ! ============================================================================
  !> Subroutine lagrangian_topofield_declaration
  !!
  !! This subroutine initialise the lagrangian grid to interpolate the relief elevation.
  !<
  ! ============================================================================
  subroutine lagrangian_topofield_declaration

    integer :: i, j, p

    grid_min( 1 ) = lstrat_xo
    grid_min( 2 ) = lstrat_yo

    grid_dim( 1 ) = lstrat_X
    grid_dim( 2 ) = lstrat_Y

    grid_dx( 1 ) = strat_dx
    grid_dx( 2 ) = strat_dx

    if( .not. allocated( x_int ) )then
        allocate( x_int( lstrat_X ) )
        ! x-grid location
        do i = 1, grid_dim( 1 )
            x_int( i ) = ( i - 1 ) * grid_dx( 1 ) + lstrat_xo
        enddo
    endif

    if( .not. allocated( y_int ) )then
        allocate( y_int( lstrat_Y ) )
        ! y-grid location
        do i = 1, grid_dim( 2 )
            y_int( i ) = ( i - 1 ) * grid_dx( 2 ) + lstrat_yo
        enddo
    endif

    if( .not. allocated( z_int ) ) allocate( z_int( lstrat_X, lstrat_Y ) )

    p = 0
    do j = 1, grid_dim( 2 )
       do i = 1, grid_dim( 1 )
            p = p + 1
            z_int( i, j ) = lfilldem( p )
       enddo
    enddo

    return

  end subroutine lagrangian_topofield_declaration
  ! ============================================================================
  !> Subroutine get_lagrangian_interpolated_relief_value
  !!
  !! This subroutine interpolate relief value.
  !<
  ! ============================================================================
  subroutine get_lagrangian_interpolated_relief_value( xi, zval )

    real( tkind )  :: xi( 2 ), zval

    ! Interpolation of relief field
    call LagrangeInterp2D( zval, z_int, xi, x_int, y_int, grid_min, grid_dx, grid_dim )

    return

  end subroutine get_lagrangian_interpolated_relief_value
  ! ============================================================================
  !> Subroutine lagrangian_velocityfield_declaration
  !!
  !! This subroutine initialise the lagrangian grid to interpolate velocity fields.
  !<
  ! ============================================================================
  subroutine lagrangian_velocityfield_declaration

    integer :: i, j, p

    grid_min2( 1 ) = strat_xo
    grid_min2( 2 ) = strat_yo

    grid_dim2( 1 ) = strat_X
    grid_dim2( 2 ) = strat_Y

    grid_dx( 1 ) = strat_dx
    grid_dx( 2 ) = strat_dx

    if( .not. allocated( vx_int ) ) allocate( vx_int( strat_X, strat_Y ) )
!    if( .not. allocated( vy_int ) ) allocate( vy_int( strat_X, strat_Y ) )

    if( .not. allocated( x_int2 ) )then
        allocate( x_int2( strat_X ) )

        ! x-grid location
        do i = 1, grid_dim2( 1 )
            x_int2( i ) = ( i - 1 ) * grid_dx( 1 ) + strat_xo
        enddo

    endif

    if( .not. allocated( y_int2 ) )then
        allocate( y_int2( strat_Y ) )

        ! y-grid location
        do i = 1, grid_dim2( 2 )
            y_int2( i ) = ( i - 1 ) * grid_dx( 2 ) + strat_yo
        enddo

    endif

    p = 0
    do j = 1, grid_dim2( 2 )
       do i = 1, grid_dim2( 1 )
            p = p + 1
            vx_int( i, j ) = vflow( p ) ! * sin( orient( p ) )
!            vy_int( i, j ) = vflow( p ) * cos( orient( p ) )
       enddo
    enddo

    return

  end subroutine lagrangian_velocityfield_declaration
  ! ============================================================================
  !> Subroutine get_lagrangian_interpolated_velocityfield_values
  !!
  !! This subroutine interpolates velocity field for steady state.
  !<
  ! ============================================================================
  subroutine get_lagrangian_interpolated_velocityfield_values( xi, val )

    real( tkind )  :: xi( 2 ), val !( 2 )

    ! Interpolation for vx velocity field
    call LagrangeInterp2D( val, vx_int, xi, x_int2, y_int2, grid_min2, grid_dx, grid_dim2 )

    ! Interpolation for vy velocity field
    !call LagrangeInterp2D( val( 2 ), vy_int, xi, x_int2, y_int2, grid_min2, grid_dx, grid_dim2 )

    return

  end subroutine get_lagrangian_interpolated_velocityfield_values
  ! ============================================================================
  !> Subroutine lagrangian_qflowfield_declaration
  !!
  !! This subroutine initialise the lagrangian grid to interpolate flow rate fields.
  !<
  ! ============================================================================
  subroutine lagrangian_qflowfield_declaration

    integer :: i, j, p

    grid_min2( 1 ) = strat_xo
    grid_min2( 2 ) = strat_yo

    grid_dim2( 1 ) = strat_X
    grid_dim2( 2 ) = strat_Y

    grid_dx( 1 ) = strat_dx
    grid_dx( 2 ) = strat_dx

    if( .not. allocated( fq_int ) ) allocate( fq_int( strat_X, strat_Y ) )

    if( .not. allocated( x_int2 ) )then
        allocate( x_int2( strat_X ) )

        ! x-grid location
        do i = 1, grid_dim2( 1 )
            x_int2( i ) = ( i - 1 ) * grid_dx( 1 ) + strat_xo
        enddo

    endif

    if( .not. allocated( y_int2 ) )then
        allocate( y_int2( strat_Y ) )

        ! y-grid location
        do i = 1, grid_dim2( 2 )
            y_int2( i ) = ( i - 1 ) * grid_dx( 2 ) + strat_yo
        enddo

    endif

    p = 0
    do j = 1, grid_dim2( 2 )
       do i = 1, grid_dim2( 1 )
            p = p + 1
            fq_int( i, j ) = qflow( p )
       enddo
    enddo

    return

  end subroutine lagrangian_qflowfield_declaration
  ! ============================================================================
  !> Subroutine get_lagrangian_interpolated_qflowfield_values
  !!
  !! This subroutine interpolates flow rates for steady state.
  !<
  ! ============================================================================
  subroutine get_lagrangian_interpolated_qflowfield_values( xi, val )

    real( tkind )  :: xi( 2 ), val

    ! Interpolation for flow rate field
    call LagrangeInterp2D( val, fq_int, xi, x_int2, y_int2, grid_min2, grid_dx, grid_dim2 )

    return

  end subroutine get_lagrangian_interpolated_qflowfield_values
  ! ============================================================================
  !> Subroutine lagrangian_flowwidthfield_declaration
  !!
  !! This subroutine initialise the lagrangian grid to interpolate flow width fields.
  !<
  ! ============================================================================
  subroutine lagrangian_flowwidthfield_declaration

    integer :: i, j, p

    grid_min2( 1 ) = strat_xo
    grid_min2( 2 ) = strat_yo

    grid_dim2( 1 ) = strat_X
    grid_dim2( 2 ) = strat_Y

    grid_dx( 1 ) = strat_dx
    grid_dx( 2 ) = strat_dx

    if( .not. allocated( w_int ) ) allocate( w_int( strat_X, strat_Y ) )

    if( .not. allocated( x_int2 ) )then
        allocate( x_int2( strat_X ) )

        ! x-grid location
        do i = 1, grid_dim2( 1 )
            x_int2( i ) = ( i - 1 ) * grid_dx( 1 ) + strat_xo
        enddo

    endif

    if( .not. allocated( y_int2 ) )then
        allocate( y_int2( strat_Y ) )

        ! y-grid location
        do i = 1, grid_dim2( 2 )
            y_int2( i ) = ( i - 1 ) * grid_dx( 2 ) + strat_yo
        enddo

    endif

    p = 0
    do j = 1, grid_dim2( 2 )
       do i = 1, grid_dim2( 1 )
            p = p + 1
            w_int( i, j ) = wflow( p )
       enddo
    enddo

    return

  end subroutine lagrangian_flowwidthfield_declaration
  ! ============================================================================
  !> Subroutine get_lagrangian_interpolated_flowwidthfield_values
  !!
  !! This subroutine interpolates flow width for steady state.
  !<
  ! ============================================================================
  subroutine get_lagrangian_interpolated_flowwidthfield_values( xi, val )

    real( tkind )  :: xi( 2 ), val

    ! Interpolation for flow rate field
    call LagrangeInterp2D( val, w_int, xi, x_int2, y_int2, grid_min2, grid_dx, grid_dim2 )

    return

  end subroutine get_lagrangian_interpolated_flowwidthfield_values
  ! ============================================================================
  !> Subroutine lagrangian_flowheightfield_declaration
  !!
  !! This subroutine initialise the lagrangian grid to interpolate flow height fields.
  !<
  ! ============================================================================
  subroutine lagrangian_flowheightfield_declaration

    integer :: i, j, p

    grid_min2( 1 ) = strat_xo
    grid_min2( 2 ) = strat_yo

    grid_dim2( 1 ) = strat_X
    grid_dim2( 2 ) = strat_Y

    grid_dx( 1 ) = strat_dx
    grid_dx( 2 ) = strat_dx

    if( .not. allocated( h_int ) ) allocate( h_int( strat_X, strat_Y ) )

    if( .not. allocated( x_int2 ) )then
        allocate( x_int2( strat_X ) )

        ! x-grid location
        do i = 1, grid_dim2( 1 )
            x_int2( i ) = ( i - 1 ) * grid_dx( 1 ) + strat_xo
        enddo

    endif

    if( .not. allocated( y_int2 ) )then
        allocate( y_int2( strat_Y ) )

        ! y-grid location
        do i = 1, grid_dim2( 2 )
            y_int2( i ) = ( i - 1 ) * grid_dx( 2 ) + strat_yo
        enddo

    endif

    p = 0
    do j = 1, grid_dim2( 2 )
       do i = 1, grid_dim2( 1 )
            p = p + 1
            h_int( i, j ) = hflow( p )
       enddo
    enddo

    return

  end subroutine lagrangian_flowheightfield_declaration
  ! ============================================================================
  !> Subroutine get_lagrangian_interpolated_flowheightfield_values
  !!
  !! This subroutine interpolates flow heights for steady state.
  !<
  ! ============================================================================
  subroutine get_lagrangian_interpolated_flowheightfield_values( xi, val )

    real( tkind )  :: xi( 2 ), val

    ! Interpolation for vx velocity field
    call LagrangeInterp2D( val, h_int, xi, x_int2, y_int2, grid_min2, grid_dx, grid_dim2 )

    return

  end subroutine get_lagrangian_interpolated_flowheightfield_values
  ! ============================================================================

end module interpol
! ============================================================================
