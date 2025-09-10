module namelists

  use definitions
  use molecule_basis

  implicit none

  real(KREAL), allocatable :: overlap(:), fock(:,:), density(:,:), &
       dipole(:,:), nuclear(:), lcaomo(:,:,:), kinetic(:)
  
  
  namelist /gennbo/ natoms, nbas, upper, bodm, form, unrest, bohr
  namelist /molecule/ znuc, coord, jobtitle
  namelist /basis/ center, label, nshell, nexp
  namelist /contract/ ncomp, nprim, nptr, zeta, ccoe
  namelist /matrices/ overlap, density, fock, lcaomo, nuclear, kinetic, dipole

end module namelists
