
module molecule_basis 

  use definitions
  use constants_parameters
  
  implicit none
  
  integer(KINT), save :: natoms, nspin, nbas, nshell, nexp, lmax

  logical, save :: upper, bodm, form, unrest, bohr

  character(len=LCHARS), save :: jobtitle
  
  real(KREAL), allocatable, save :: coord(:,:), znuc(:), zeta(:), & 
       ccoe(:,:)
  
  integer(KINT), allocatable, save :: center(:), label(:), ncomp(:), nprim(:), &
     nptr(:)

  real(KREAL), allocatable, save :: bas_norm(:), prim_norm(:)

  ! the following arrays are used in bas_powers() to generate the needed
  ! linear combinations of Cartesian basis functions for a spherical
  ! basis (5d,7f,etc)
  ! maxbcomp is defined in module constants_parameters
  integer(KINT) :: nc1, nx1(maxbcomp), ny1(maxbcomp), nz1(maxbcomp)
  integer(KINT) :: nc2, nx2(maxbcomp), ny2(maxbcomp), nz2(maxbcomp)
  real(KREAL) :: cc1(maxbcomp), cc2(maxbcomp)

  ! provide initial values for some of the variables in namelists
  ! &GENNBO and &BASIS. Running the program will fail if those are not
  ! set explicitly in the namelist input. This is on purpose.

  data natoms, nspin, nbas, nshell, nexp, lmax / 0,0,0,0,0,0/

  data &
       upper, bodm, form, unrest, bohr /&
       .false., .false., .false., .false., .false. /
  
end module molecule_basis
  
