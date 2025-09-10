subroutine bas_powers(iconv,i,nc,nx,ny,nz,cc)

  ! powers of x, y, z for a given type of basis function

  ! input:
  ! iconv = integer indicating which FILE47 basis set convention we're using
  !         (1 = ORCA, 2 = Gaussian, 3 = Molden2AIM)
  ! i = integer with nbo basis function label (1, 101, 102, ...)

  ! output:
  ! nc = number of components in case of spherical basis (nc=1 otherwise)
  !      (e.g., 3 Cartesian components for dz2 = 2zz-xx-yy)
  ! nx(:), ny(:), nz(:) = arrays with powers of x,y,z per component
  ! cc(:) = coefficients for components, prefactors depend on iconv

  ! For basis functions that are not simply the same as one Cartesian function,
  ! i.e., dz2, d(xx-yy), and most spherical functions with l>=3, the
  ! coefficients come from IOData (https://iodata.readthedocs.io) but 
  ! were tweaked. For example, IOData's
  !
  ! python tools/harmonics.py none python 2
  !
  ! gave a set of d-functions for which all coefficients had to be
  ! divided by sqrt(3) to give properly normalized contracted
  ! functions with an ORCA generated FILE47. The contraction
  ! coefficients from FILE47 must be multiplied into the basis set to
  ! give normalized functions either (i) from plain unnormalized
  ! Cartesians, or (ii) from the spherical functions constructed with
  ! IOData coefficients after the kind of adjustment by a global
  ! factor mentioned above.

  ! Gaussian convention is: primitives must be normalized. The
  ! coefficients for the spherical functions are then given by IOData's
  ! python tools/harmonics.py L2 python 3

  ! Example:
  ! bas_powers(1, 255,nc,nx,ny,nz,cc)
  ! returns nc=3, nx(:)=(0,2,0), ny(:)=(0,0,2), nz(:)=(2,0,0), and
  ! cc(:) = (1,-0.5,-0.5)/sqrt(3), which translates into the function
  ! 1/sqrt(3) (z**2 - 0.5x**2 - 0.5y**2) 

  use definitions
  use constants_parameters
  
  implicit none

  integer(KINT), intent(in) :: iconv ! basis set convention used
  integer(KINT), intent(in) :: i
  integer(KINT), intent(out) :: nc, nx(maxbcomp), ny(maxbcomp), nz(maxbcomp)
  real(KREAL), intent(out) :: cc(maxbcomp)
  real(KREAL) :: rtmp
  integer :: dbg, x,y,z, ic

  character(10) pname ! subroutine name for debug output

  ! ---------------------------------------------------------------------------
  
  dbg = 0 ! debug switch
  pname = 'bas_powers'

  if (iconv.le.0 .or. iconv.gt.3) &
    stop pname//': basis convention flag iconv out of range'

  nc = 0
  nx(:) = 0
  ny(:) = 0
  nz(:) = 0
  cc(:) = zero

  ! in the following, for each basis function type and for each
  ! Cartesian component, where applicable, of spherical functions, we
  ! assign:
  ! x, y, z = integer power of variable x, y, z
  ! ic = component number (>1 if there is more than one component
  !                        e.g. for dz2 = 2zz-xx-yy there are 3)
  ! powers for each component and normalization of each component
  ! from self-overlap of Cartesian primitive are stored in the arrays
  ! nx(:), ny(:), nz(:), cc(:) 

  select case(i)

    ! -----------
    ! s functions
    ! -----------
    
  case (1, 51) ! s
    nc = 1    
    ic = 1
    x = 0; y = 0; z = 0 
    cc(ic) = one
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': s function detected'

    ! ------------
    ! p functions:
    ! ------------
    
  case (101, 151) ! px
    nc = 1 
    ic = 1
    x = 1; y = 0; z = 0 
    cc(ic) = one
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': px function detected'

  case (102, 152) ! py
    nc = 1 
    ic = 1
    x = 0; y = 1; z = 0
    cc(ic) = one
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': py function detected'

  case (103, 153) ! pz
    nc = 1
    ic = 1
    x = 0; y = 0; z = 1
    cc(ic) = one
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': pz function detected'

    ! ------------
    ! d functions:
    ! ------------

  case (201) ! dxx
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xx function detected'

  case (202, 251) ! dxy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(three)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xy function detected'

  case (203, 252) ! dxz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(three)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 0; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xz function detected'

  case (204) ! dyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1 
    ic = 1
    x = 0; y = 2; z = 0
    cc(ic) = rtmp
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yy function detected'

  case (205, 253) ! dyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(three)
    end select
    nc = 1 
    ic = 1 
    x = 0; y = 1; z = 1
    cc(ic) = rtmp
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yz function detected'

  case (206) ! dzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1 
    x = 0; y = 0; z = 2
    cc(ic) = rtmp
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//':zz function detected'

  case (254) ! xx-yy
    select case(iconv)
    case (1)
      rtmp = one/sqrt(three)
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 2
    ic = 1 
    x = 2; y = 0; z = 0 
    cc(ic) = rtmp * 0.86602540378443865_KREAL
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2 
    x = 0; y = 2; z = 0 
    cc(ic) = -rtmp*0.86602540378443865_KREAL 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xx-yy function detected'

  case (255) ! 2zz-xx-yy
    select case(iconv)
    case (1)
      rtmp = one/sqrt(three)
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 3
    ic = 1 
    x = 0; y = 0; z = 2 
    cc(ic) = rtmp
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    x = 2; y = 0; z = 0
    cc(ic) = -half*rtmp
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3 
    x = 0; y = 2; z = 0
    cc(ic) = -half*rtmp
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': 2zz-xx-yy function detected'

    ! ------------
    ! f functions:
    ! ------------
    
  case (301) ! fxxx
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 3; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxx function detected'

  case (302) ! fxxy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxy function detected'

  case (303) ! fxxz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxz function detected'

  case (304) ! fxyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 2; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyy function detected'

  case (305) ! fxyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(15.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 1; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyz function detected'

  case (306) ! fxzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 0; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xzz function detected'

  case (307) ! fyyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 3; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyy function detected'

  case (308) ! fyyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 2; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyz function detected'

  case (309) ! fyzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 1; z = 2 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yzz function detected'

  case (310) ! fzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 0; z = 3 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': zzz function detected'

  case (351) ! f0
    ! -(3*x**2*z + 3*y**2*z - 2*z**3)
    select case(iconv)
    case (1)
      rtmp = one/sqrt(60.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = half
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.67082039324993691_KREAL ! 3*sqrt(5)/10
    x = 2; y = 0; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.67082039324993691_KREAL
    x = 0; y = 2; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.0_KREAL
    x = 0; y = 0; z = 3 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': f0 function detected'

  case (352) ! f1+
    ! x**3 + x*y**2 - 4*x*z**2
    select case(iconv)
    case (1)
      rtmp = one/sqrt(40.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(3.0_KREAL/8.0_KREAL)
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.61237243569579452_KREAL
    x = 3; y = 0; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.27386127875258306_KREAL
    x = 1; y = 2; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*4.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.0954451150103322_KREAL
    x = 1; y = 0; z = 2 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': f1+ function detected'

  case (353) ! f1-
    ! x**2*y + y**3 - 4*y*z**2
    select case(iconv)
    case (1)
      rtmp = one/sqrt(40.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(3.0_KREAL/8.0_KREAL)
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.27386127875258306_KREAL
    x = 2; y = 1; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.61237243569579452_KREAL
    x = 0; y = 3; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*4.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.0954451150103322_KREAL
    x = 0; y = 1; z = 2 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': f1- function detected'

  case (354) ! f2+
    ! x**2*z - y**2*z
    select case(iconv)
    case (1)
      rtmp = one/sqrt(4.0_KREAL)
    case(2)
      rtmp = sqrt(3.0_KREAL/4.0_KREAL)
    case(3)
      rtmp = sqrt(30.0_KREAL/8.0_KREAL)
    end select
    nc = 2
    ic = 1
    cc(ic) = rtmp*1.0_KREAL
    x = 2; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*1.0_KREAL
    x = 0; y = 2; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': f2+ function detected'

  case (355) ! f2-
    ! x*y*z
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(15.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 1; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': f2- function detected'

  case (356) ! f3+
    ! -x**3 + 3*x*y**2
    select case(iconv)
    case (1)
      rtmp = one/sqrt(24.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      ! note the sign change
      rtmp = -sqrt(5.0_KREAL/8.0_KREAL)
    end select
    nc = 2
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.79056941504209483_KREAL
    x = 3; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.0606601717798213_KREAL
    x = 1; y = 2; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': f3+ function detected'

  case (357) ! f3-
    ! -3*x**2*y + y**3
    select case(iconv)
    case (1)
      rtmp = one/sqrt(24.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      ! note the sign change
      rtmp = -sqrt(5.0_KREAL/8.0_KREAL)
    end select
    nc = 2
    ic = 1
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.0606601717798213_KREAL
    x = 2; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.79056941504209483_KREAL
    x = 0; y = 3; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': f3- function detected'

    ! ------------
    ! g-functions:
    ! ------------

  case (401) ! gxxxx
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 4; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxx function detected'

  case (402) ! gxxxy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(7.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 3; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxy function detected'

  case (403) ! gxxxz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(7.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 3; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxz function detected'

  case (404) ! gxxyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(35.0_KREAL/3.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 2; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxyy function detected'

  case (405) ! gxxyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(35.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 1; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxyz function detected'

  case (406) ! gxxzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(35.0_KREAL/3.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 0; z = 2 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxzz function detected'

  case (407) ! gxyyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(7.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 3; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyyy function detected'

  case (408) ! gxyyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(35.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 2; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyyz function detected'

  case (409) ! gxyzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(35.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 1; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyzz function detected'

  case (410) ! gxzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(7.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 0; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xzzz function detected'

  case (411) ! gyyyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 4; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyyy function detected'

  case (412) ! gyyyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(7.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 3; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyyz function detected'

  case (413) ! gyyzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(35.0_KREAL/3.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 2; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyzz function detected'

  case (414) ! gyzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(7.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 1; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yzzz function detected'

  case (415) ! gzzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 0; z = 4
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': zzzz function detected'

  case (451) ! g0
    !  3*x**4 + 6*x**2*y**2  - 24*x**2*z**2 + 3*y**4 - 24*y**2*z**2 + 8*z**4
    select case(iconv)
    case (1)
      rtmp = one/sqrt(2240.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = one/8.0_KREAL
    end select
    nc = 6
    ic = 1
    cc(ic) = rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.375_KREAL
    x = 4; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*6.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.21957751641341997_KREAL
    x = 2; y = 2; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*24.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.87831006565367986_KREAL
    x = 2; y = 0; z = 2 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.375_KREAL
    x = 0; y = 4; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 5
    cc(ic) = -rtmp*24.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.87831006565367986_KREAL
    x = 0; y = 2; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 6
    cc(ic) = rtmp*8.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.0_KREAL
    x = 0; y = 0; z = 4 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g0 function detected'

  case (452) ! g1+
    !  -3*x**3*z - 3*x*y**2*z + 4*x*z**3
    select case(iconv)
    case (1)
      rtmp = one/sqrt(56.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL/8.0_KREAL)
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.89642145700079523_KREAL
    x = 3; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.40089186286863658_KREAL
    x = 1; y = 2; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*4.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1952286093343936_KREAL
    x = 1; y = 0; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g1+ function detected'

  case (453) ! g1-
    !  -3*x**2*y*z - 3*y**3*z + 4*y*z**3
    select case(iconv)
    case (1)
      rtmp = one/sqrt(56.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL/8.0_KREAL)
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.40089186286863658_KREAL
    x = 2; y = 1; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.89642145700079523_KREAL
    x = 0; y = 3; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*4.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1952286093343936_KREAL
    x = 0; y = 1; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g1- function detected'

  case (454) ! g2+
    !  -(x**4 - 6*x**2*z**2 - y**4 + 6*y**2*z**2)
    select case(iconv)
    case (1)
      rtmp = one/sqrt(112.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL/16.0_KREAL)
    end select
    nc = 4
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.55901699437494742_KREAL
    x = 4; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*6.0_KREAL
    if (iconv.eq.2) cc(ic) =  0.98198050606196572_KREAL
    x = 2; y = 0; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.55901699437494742_KREAL
    x = 0; y = 4; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = -rtmp*6.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.98198050606196572_KREAL
    x = 0; y = 2; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g2+ function detected'

  case (455) ! g2-
    !  -(x**3*y + x*y**3 - 6*x*y*z**2)
    select case(iconv)
    case (1)
      rtmp = one/sqrt(28.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(5.0_KREAL/4.0_KREAL)
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.42257712736425829_KREAL
    x = 3; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.42257712736425829_KREAL
    x = 1; y = 3; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*6.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1338934190276817_KREAL
    x = 1; y = 1; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g2- function detected'

  case (456) ! g3+
    !  -(x**3*z - 3*x*y**2*z)
    select case(iconv)
    case (1)
      rtmp = one/sqrt(8.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(35.0_KREAL/8.0_KREAL)
    end select
    nc = 2
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.79056941504209483_KREAL
    x = 3; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.0606601717798213_KREAL
    x = 1; y = 2; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g3+ function detected'

  case (457) ! g3-
    !  -3*x**2*y*z + y**3*z
    select case(iconv)
    case (1)
      rtmp = one/sqrt(8.0_KREAL)
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(35.0_KREAL/8.0_KREAL)
    end select
    nc = 2
    ic = 1
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.0606601717798213_KREAL
    x = 2; y = 1; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.79056941504209483_KREAL
    x = 0; y = 3; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g3- function detected'

  case (458) ! g4+
    !  -(x**4 - 6*x**2*y**2 + y**4)
    select case(iconv)
    case (1)
      rtmp = one/8.0_KREAL
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(35.0_KREAL/64.0_KREAL)
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.73950997288745201_KREAL
    x = 4; y = 0; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*6.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.299038105676658_KREAL
    x = 2; y = 2; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.73950997288745201_KREAL
    x = 0; y = 4; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g4+ function detected'

  case (459) ! g4-
    !  -(x**3*y - x*y**3)
    select case(iconv)
    case (1)
      rtmp = one/2.0_KREAL
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(35.0_KREAL/4.0_KREAL)
    end select
    nc = 2
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1180339887498948_KREAL
    x = 3; y = 1; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.1180339887498948_KREAL
    x = 1; y = 3; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': g4- function detected'
    
    ! ------------
    ! h-functions:
    ! ------------

  case (501) ! hxxxxx
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 5; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxxx function detected'

  case (502) ! hxxxxy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = 3.0_KREAL
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 4; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxxy function detected'

  case (503) ! hxxxxz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = 3.0_KREAL
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 4; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxxz function detected'

  case (504) ! hxxxyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(21.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 3; y = 2; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxyy function detected'

  case (505) ! hxxxyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(63.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 3; y = 1; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxyz function detected'

  case (506) ! hxxxzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(21.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 3; y = 0; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxxzz function detected'

  case (507) ! hxxyyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(21.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 3; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxyyy function detected'

  case (508) ! hxxyyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(105.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 2; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxyyz function detected'

  case (509) ! hxxyzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(105.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 1; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxyzz function detected'

  case (510) ! hxxzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(21.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 2; y = 0; z = 3 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xxzzz function detected'

  case (511) ! hxyyyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = 3.0_KREAL
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 4; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyyyy function detected'

  case (512) ! hxyyyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(63.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 3; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyyyz function detected'

  case (513) ! hxyyzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(105.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 2; z = 2 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyyzz function detected'

  case (514) ! hxyzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(63.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 1; z = 3 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xyzzz function detected'

  case (515) ! hxzzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = 3.0_KREAL
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 1; y = 0; z = 4
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': xzzzz function detected'

  case (516) ! hyyyyy
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 5; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyyyy function detected'

  case (517) ! hyyyyz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = 3.0_KREAL
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 4; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyyyz function detected'

  case (518) ! hyyyzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(21.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 3; z = 2 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyyzz function detected'

  case (519) ! hyyzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(21.0_KREAL)
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 2; z = 3 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yyzzz function detected'

  case (520) ! hyzzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = 3.0_KREAL
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 1; z = 4
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': yzzzz function detected'

  case (521) ! hzzzzz
    select case(iconv)
    case (1)
      rtmp = one
    case(2)
      rtmp = one
    case(3)
      rtmp = one
    end select
    nc = 1
    ic = 1
    cc(ic) = rtmp
    x = 0; y = 0; z = 5
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': zzzzz function detected'

    ! order of functions in Cartesian basis should be:
    ! x**5 x**4*y x**4*z x**3*y**2 x**3*y*z x**3*z**2 x**2*y**3
    ! x**2*y**2*z x**2*y*z**2 x**2*z**3 x*y**4 x*y**3*z x*y**2*z**2
    ! x*y*z**3 x*z**4 y**5 y**4*z y**3*z**2 y**2*z**3 y*z**4 z**5

  case (551) ! h0
    ! 15*x**4*z +30*x**2*y**2*z -40*x**2*z**3 +15*y**4*z  -40*y**2*z**3 +8*z**5
    select case(iconv)
    case (1)
      rtmp = one/8.0_KREAL
    case(2)
      rtmp = one
    case(3)
      rtmp = one/8.0_KREAL
    end select
    nc = 6
    ic = 1
    cc(ic) = rtmp*15.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.625_KREAL
    x = 4; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*30.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.36596252735569994_KREAL
    x = 2; y = 2; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*40.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.0910894511799619_KREAL
    x = 2; y = 0; z = 3 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = rtmp*15.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.625_KREAL
    x = 0; y = 4; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 5
    cc(ic) = -rtmp*40.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.0910894511799619_KREAL
    x = 0; y = 2; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 6
    cc(ic) = rtmp*8.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.0_KREAL
    x = 0; y = 0; z = 5
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h0 function detected'

  case (552) ! h1+
    ! x**5 +2*x**3*y**2 -12*x**3*z**2 + x*y**4 -12*x*y**2*z**2 +8*x*z**4
    select case(iconv)
    case (1)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL)/8.0_KREAL
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL)/8.0_KREAL
    end select
    nc = 6
    ic = 1
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.48412291827592711_KREAL
    x = 5; y = 0; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.21128856368212914_KREAL
    x = 3; y = 2; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*12.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.2677313820927749_KREAL
    x = 3; y = 0; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.1613743060919757_KREAL
    x = 1; y = 4; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 5
    cc(ic) = -rtmp*12.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.56694670951384084_KREAL
    x = 1; y = 2; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 6
    cc(ic) = rtmp*8.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.2909944487358056_KREAL
    x = 1; y = 0; z = 4
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h1+ function detected'

  case (553) ! h1-
    ! x**4*y +2*x**2*y**3  -12*x**2*y*z**2 + y**5 -12*y**3*z**2 + 8*y*z**4
    select case(iconv)
    case (1)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL)/8.0_KREAL
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL)/8.0_KREAL
    end select
    nc = 6
    ic = 1
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.1613743060919757_KREAL
    x = 4; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.21128856368212914_KREAL
    x = 2; y = 3; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*12.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.56694670951384084_KREAL
    x = 2; y = 1; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.48412291827592711_KREAL
    x = 0; y = 5; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 5
    cc(ic) = -rtmp*12.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.2677313820927749_KREAL
    x = 0; y = 3; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 6
    cc(ic) = rtmp*8.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.2909944487358056_KREAL
    x = 0; y = 1; z = 4
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h1- function detected'

  case (554) ! h2+
    ! -(x**4*z) + 2*x**2*z**3 + y**4*z - 2*y**2*z**3
    select case(iconv)
    case (1)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL*7.0_KREAL)/4.0_KREAL
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL*7.0_KREAL)/4.0_KREAL
    end select
    nc = 4
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.85391256382996653_KREAL
    x = 4; y = 0; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1180339887498948_KREAL
    x = 2; y = 0; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.85391256382996653_KREAL
    x = 0; y = 4; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = -rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.1180339887498948_KREAL
    x = 0; y = 2; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h2+ function detected'

  case (555) ! h2-
    ! -(x**3*y*z) - x*y**3*z + 2*x*y*z**3
    select case(iconv)
    case (1)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL*7.0_KREAL)/2.0_KREAL
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(3.0_KREAL*5.0_KREAL*7.0_KREAL)/2.0_KREAL
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.64549722436790281_KREAL
    x = 3; y = 1; z = 1 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.64549722436790281_KREAL
    x = 1; y = 3; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.2909944487358056_KREAL
    x = 1; y = 1; z = 3
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h2- function detected'

  case (556) ! h3+
    ! x**5 - 2*x**3*y**2 - 8*x**3*z**2 - 3*x*y**4 + 24*x*y**2*z**2
    select case(iconv)
    case (1)
      rtmp = sqrt(2.0_KREAL*5.0_KREAL*7.0_KREAL)/16.0_KREAL
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(2.0_KREAL*5.0_KREAL*7.0_KREAL)/16.0_KREAL
    end select
    nc = 5
    ic = 1
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.52291251658379722_KREAL
    x = 5; y = 0; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.22821773229381921_KREAL
    x = 3; y = 2; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*8.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.91287092917527686_KREAL
    x = 3; y = 0; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = -rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.52291251658379722_KREAL
    x = 1; y = 4; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 5
    cc(ic) = rtmp*24.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.224744871391589_KREAL
    x = 1; y = 2; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h3+ function detected'

  case (557) ! h3-
    ! -(-3*x**4*y - 2*x**2*y**3 + 24*x**2*y*z**2 + y**5 - 8*y**3*z**2)
    select case(iconv)
    case (1)
      rtmp = sqrt(2.0_KREAL*5.0_KREAL*7.0_KREAL)/16.0_KREAL
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(2.0_KREAL*5.0_KREAL*7.0_KREAL)/16.0_KREAL
    end select
    nc = 5
    ic = 1
    cc(ic) = rtmp*3.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.52291251658379722_KREAL
    x = 4; y = 1; z = 0 
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*2.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.22821773229381921_KREAL
    x = 2; y = 3; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*24.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.224744871391589_KREAL
    x = 2; y = 1; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 4
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.52291251658379722_KREAL
    x = 0; y = 5; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 5
    cc(ic) = rtmp*8.0_KREAL
    if (iconv.eq.2) cc(ic) = -0.91287092917527686_KREAL
    x = 0; y = 3; z = 2
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h3- function detected'

  case (558) ! h4+
    ! -(x**4*z - 6*x**2*y**2*z + y**4*z)
    select case(iconv)
    case (1)
      rtmp = sqrt(5.0_KREAL*7.0_KREAL*9.0_KREAL)/8.0_KREAL
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(5.0_KREAL*7.0_KREAL*9.0_KREAL)/8.0_KREAL
    end select
    nc = 3
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.73950997288745201_KREAL
    x = 4; y = 0; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*6.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.299038105676658_KREAL
    x = 2; y = 2; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.73950997288745201_KREAL
    x = 0; y = 4; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h4+ function detected'

  case (559) ! h4-
    ! -(x**3*y*z - x*y**3*z)
    select case(iconv)
    case (1)
      rtmp = sqrt(5.0_KREAL*7.0_KREAL*9.0_KREAL)/2.0_KREAL
    case(2)
      rtmp = one
    case(3)
      ! note the minus sign
      rtmp = -sqrt(5.0_KREAL*7.0_KREAL*9.0_KREAL)/2.0_KREAL
    end select
    nc = 2
    ic = 1
    cc(ic) = -rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1180339887498948_KREAL
    x = 3; y = 1; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.1180339887498948_KREAL
    x = 1; y = 3; z = 1
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h4- function detected'

  case (560) ! h5+
    ! x**5 - 10*x**3*y**2 + 5*x*y**4
    select case(iconv)
    case (1)
      rtmp = sqrt(2.0_KREAL*7.0_KREAL*9.0_KREAL)/16.0_KREAL
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(2.0_KREAL*7.0_KREAL*9.0_KREAL)/16.0_KREAL
    end select
    nc = 3
    ic = 1
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.70156076002011401_KREAL
    x = 5; y = 0; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*10.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.5309310892394863_KREAL
    x = 3; y = 2; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*5.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1692679333668567_KREAL
    x = 1; y = 4; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h5+ function detected'

  case (561) ! h5-
    ! 5*x**4*y - 10*x**2*y**3 + y**5
    select case(iconv)
    case (1)
      rtmp = sqrt(2.0_KREAL*7.0_KREAL*9.0_KREAL)/16.0_KREAL
    case(2)
      rtmp = one
    case(3)
      rtmp = sqrt(2.0_KREAL*7.0_KREAL*9.0_KREAL)/16.0_KREAL
    end select
    nc = 3
    ic = 1
    cc(ic) = rtmp*5.0_KREAL
    if (iconv.eq.2) cc(ic) = 1.1692679333668567_KREAL
    x = 4; y = 1; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 2
    cc(ic) = -rtmp*10.0_KREAL
    if (iconv.eq.2) cc(ic) = -1.5309310892394863_KREAL
    x = 2; y = 3; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    ic = 3
    cc(ic) = rtmp*1.0_KREAL
    if (iconv.eq.2) cc(ic) = 0.70156076002011401_KREAL
    x = 0; y = 5; z = 0
    nx(ic) = x ; ny(ic) = y; nz(ic) = z
    if (dbg>0) write (out,*) pname//': h5- function detected'
    
    
    ! -----------------------------------------
    ! here is the end of the rope. Abort, if we
    ! haven't found a match yet:
    ! -----------------------------------------
    
  case default
    stop 'bas_powers: b.f. label not recognized/not implemented'
    
  end select

  return
end subroutine bas_powers
