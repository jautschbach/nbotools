
module constants_parameters

  ! define some constants and parameters used in various places
  ! throughout the code

  ! (c) 2024 Jochen Autschbach

  use definitions

  implicit none

  ! numerical constants:
  
  real(KREAL), parameter :: zero=0.0_KREAL
  real(KREAL), parameter :: one=1.0_KREAL
  real(KREAL), parameter :: two=2.0_KREAL
  real(KREAL), parameter :: three=3.0_KREAL
  real(KREAL), parameter :: four=4.0_KREAL
  real(KREAL), parameter :: half=0.5_KREAL
  real(KREAL), parameter :: third=one/three
  real(KREAL), parameter :: fourth=0.25_KREAL
  real(KREAL), parameter :: pi=four*atan(one)
  real(KREAL), parameter :: small=1E-5_KREAL
  real(KREAL), parameter :: tiny=1E-10_KREAL
  real(KREAL), parameter :: tnsy=1E-15_KREAL

  ! physical constants
  real(KREAL), parameter :: angstrom=0.529177210544_KREAL

  ! dimensioning for some fixed-size arrays
  
  ! NBO maximum basis angular momentum:
  integer(KINT), parameter :: maxl = 6

  ! maximum linear combination of Cartesians in spherical basis fcts.
  ! for simplicity, we set this to the number of Cartesian polynomials
  ! corresponding to maxl, which is (l+1)(l+2)/2 or 1,3,6,10,15,21,28
  ! for l=0,1,2,3,4,5,6) although the actual number of relevant linear
  ! combos is smaller:
  integer(KINT), parameter :: maxbcomp = (maxl+1)*(maxl+2)/2
  
  ! limit for exponentiation, to avoid floating point exceptions:
  real(KREAL), parameter :: explim=300

  ! file I/O units etc.

  integer(KINT), parameter :: inp=5, out=6, err=0
  
  
end module constants_parameters
  
