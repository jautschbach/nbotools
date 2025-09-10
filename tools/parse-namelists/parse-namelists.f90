program parse_namelists
  
  ! process a file with data in namelist variables equivalent to
  ! an NBO FILE47. Then test the basis set convention and calculate
  ! the overlap matrix & check if it is the same as the overlap
  ! coming from FILE47.

  ! module constants_parameters contains definitions of some
  ! hardcoded variables.
  ! module factorials has a hardcoded limit of 40! Increase if needed.

  ! (c) 2024 Jochen Autschbach

  use definitions
  use factorials
  use constants_parameters
  use molecule_basis
  use namelists

  implicit none

  integer(KINT) :: i,ii, j,k,l,n,nx,ios, iconv, dbg, found(10), itmp
  
  character(len=LCHARS) :: fmt

  !real(KREAL), external :: dnewt

  real(KREAL), allocatable :: smat(:)


  !integer idigits
  !idigits(i) = nint(log10(float(i)))

  ! ---------------------------------------------------------------------------

  ! debug level
  dbg = 1

  ! --------------
  ! initialization
  ! --------------

  call init_factorials

  !if (dbg>0) write(out,*) fact(0:mxfact)
  !if (dbg>0) write(out,*) dnewt(15,3)
  
  ! -------------------------------------------
  ! provide defaults for basic namelist entries
  ! -------------------------------------------
  
  natoms = 0       ! number of atoms in the molecule
  nbas   = 0       ! number of contracted basis functions
  upper  = .false. ! true if symmetric matrices are in packed storage
  bodm   = .false. ! true if N=tr[P S] with S not being a unit matrix
  form   = .false. ! unused
  unrest = .false. ! true if we have separate alpha and beta spin data

  ! --------------------------------------
  ! read namelist data from STDIN. Arrays
  ! are allocated as needed as we go along
  ! --------------------------------------

  if (dbg>0) write(out,*) 'reading namelist ''&gennbo'''
  rewind(inp)
  read(inp,nml=gennbo,iostat=ios)
  if (dbg>0) write(out,*) 'namelist &gennbo status:',ios
  if (ios.ne.0) stop 'error parsing namelist &gennbo'

  if (natoms.le.0) stop 'natoms not set'
  if (nbas.le.0) stop 'nbas not set'
  
  allocate(znuc(natoms))
  allocate(coord(3,natoms))

  rewind(inp)
  if (dbg>0) write(out,*) 'reading namelist ''&molecule'''
  read(inp,nml=molecule,iostat=ios)
  if (dbg>0) write(out,*) 'namelist &molecule status',ios
  if (ios.ne.0) stop 'error parsing namelist &molecule'

  allocate(center(nbas))
  allocate(label(nbas))

  rewind(inp)
  if (dbg>0) write(out,*) 'reading namelist ''&basis'''
  read(inp,nml=basis,iostat=ios)
  if (dbg>0) write(out,*) 'namelist &basis status',ios
  if (ios.ne.0) stop 'error parsing namelist &basis'

  if (nexp.le.0) stop 'nexp not set'
  if (nshell.le.0) stop 'nshell not set'

  ! determine largest angular momentum in the basis: lmax

  lmax = maxval(label)
  lmax = lmax/100 ! div(lmax,100)
  if (lmax.gt.maxl) stop 'error: max. L of basis too large. aborting'

  allocate(ccoe(nexp,0:lmax))
  ccoe = 0
  allocate(ncomp(nshell))
  allocate(nprim(nshell))
  allocate(nptr(nshell))
  allocate(zeta(nexp))

  rewind(inp)
  if (dbg>0) write(out,*) 'reading namelist ''&contract'''
  read(inp,nml=contract,iostat=ios)
  if (dbg>0) write(out,*) 'namelist &contract status',ios
  if (ios.ne.0) stop 'error parsing namelist &contract'


  ! ---------------------------------------
  ! print the data that we just parsed from
  ! the namelist input
  ! ---------------------------------------

  !write(out,*) 'znuc: ', znuc(1:natoms)
  write(out,'(/1x,a)') 'molecule specification: '
  write(out,'(1x,''NATOMS='',i0)') natoms
  write(out,*) trim(jobtitle)
  do i = 1,natoms
     write(out,'(1x,f3.0,1x,3(f15.6))') znuc(i),coord(1:3,i)
  end do

  write(out,'(/1x,a)') 'basis set specification'
  write(out,'(1x,''NBAS  ='',i0)') nbas
  write(out,'(1x,''NEXP  ='',i0)') nexp
  write(out,'(1x,''NSHELL='',i0)') nshell
  write(out,'(1x,''Largest angular momentum l='',i0)') lmax

  write(out,'(1x,''LABEL = basis function labels:'')')
  write(out,'(1x,10(i3,1x))') label(:)
  
  write(out,'(1x,''Exponents:'')')
  write(out,'(1x,3e20.12)') zeta(:)

  write(out,'(1x,''NPRIM = Primitives per shell:'')')
  write(out,'(1x,10(i3,1x))') nprim(:)

  write(out,'(1x,''NPTR = pointer array:'')')
  write(out,'(1x,10(i3,1x))') nptr(:)

  do l = 0, lmax
     write(out,'(1x,''Contraction coefficients l='',i0)') l
     write(out,'(1x,3e20.12)') ccoe(:,l)
  end do

  ! ---------------------------------------------------
  ! variable 'nspin' determines how many spins we have.
  ! 'unrest' is the same as the OPEN variable in FILE47
  ! ---------------------------------------------------

  nspin = 1
  if (unrest) nspin = 2

  ! ---------------------------------
  ! allocate arrays and read matrices
  ! from &matrices namelist
  ! ---------------------------------

  n = nbas
  if(upper) then
     nx = n*(n+1)/2
  else
     nx = n**2
  end if

  allocate(overlap(nx))
  allocate(kinetic(nx))
  allocate(nuclear(nx))
  allocate(dipole(nx,3))
  allocate(fock(nx,nspin))
  allocate(density(nx,nspin))
  allocate(lcaomo(n,n,nspin))

  rewind(inp)
  if (dbg>0) write(out,*) 'reading namelist ''&matrices'''
  read(inp,nml=matrices,iostat=ios)
  if (dbg>0) write(out,*) 'namelist &matrices status',ios
  if (ios.ne.0) stop 'error parsing namelist &matrices'

  ! convert nuclear coords to Bohr units if in Angstrom
  if(.not.bohr) coord = coord / angstrom

  ! -----------------------------
  ! check and print the basis set
  ! -----------------------------

  found(:) = 0 ! an array of flags for different basis set conventions
  allocate (bas_norm(nbas)) ! defined in molecule_basis_module

  ! test basis set convention:
  ! 1 = Orca, 2 = Gaussian, 3 = Molden2AIM
  ! the last argument .false. in the call means no basis info is printed.
  iconv = 1 
  call check_print_basis(iconv,found(iconv),.false.)
  iconv = 2
  call check_print_basis(iconv,found(iconv),.false.)
  iconv = 3
  call check_print_basis(iconv,found(iconv),.false.)

  !write(out,*) 'found =',found

  itmp = sum(found(:))

  if (itmp.gt.1) then ! we found more than one basis convention
                       ! leading to normalized basis functions
    write(out,*) 'Basis data appears to conform to several conventions ! '
    write(out,*) 'Using Gaussian convention for the following list ...  '
    iconv = 2
  else if (itmp.eq.1) then ! found one matching basis convention
    iconv = maxloc(found,1) ! which one did we find?
    select case (iconv)
    case(1) ! ORCA
      write(out,*) 'Basis seems to be in ORCA 47 convention'
    case(2)
      write(out,*) 'Basis seems to be in Gaussian 47 convention'
    case(3)
      write(out,*) 'Basis seems to be in Molden2AIM convention'
    case default
      stop 'Error in iconv assignment'
    end select
  else if (itmp.lt.1) then ! found no basis convention that works
    write(out,*) '**** WARNING ****'
    write(out,*) &
       'Basis data appears to conform to NONE of the conventions I know ! '
    !write(out,*) 'Using Gaussian convention for the following list ...  '
    iconv = 3
    stop 'error iconv'
  else
    stop 'error iconv: unknown error'
  end if ! itmp
  !write(out,*) iconv

  call check_print_basis(iconv,found(iconv),.true.)

    

  ! ----------------------------------------------------
  ! calculate the overlap matrix of the contracted basis
  ! ----------------------------------------------------
  
  allocate (smat(nx))
  call overlap_matrix(iconv,smat,nx)

  ! compare the overlap matrix elements we just calculated to
  ! the overlap from FILE47, and print WARNING next to the data if they
  ! don't appear to be the same.

  write(out,'(5(/1x,a))') &
     'I have now calculated the overlap matrix S for the AO basis.', &
     'If one or more elements of S do not agree with the overlap from',&
     'FILE47, the elements for which there is a discrepancy (> 1E-5)',&
     'will be printed (order: mine, FILE47) along with the string WARNING.',&
     'If there is no output, then the elements agree to within 1E-5.'
  
  
  ii = 0
  do i = 1,nbas
    do j = 1,i
      ii = ii+1
      if (upper) then
        k = ii
      else
        k = i + nbas*(j-1)
      end if
      if (abs(smat(ii)-overlap(k)).gt.small) then
        write (fmt,*) '(1x,''i,j,Sij :'',i3,1x,i3,f15.8,f15.8,'' WARNING'')'
        write(out,fmt) i, j, smat(ii), overlap(k)
      else
        !write (fmt,*) '(1x,''i,j,Sij :'',i3,i3,f15.8,f15.8)'
        !write(out,fmt) i, j, smat(ii), overlap(k)
        continue
      end if
    end do
  end do

  write(out,'(2(/1x,a)/)') '---','End of overlap matrix comparison.'
  deallocate (smat)

  ! -----------------
  ! deallocate arrays
  ! -----------------
  
  deallocate(znuc,coord)
  deallocate(center,label)
  deallocate(ccoe,ncomp,nprim,nptr,zeta)
  deallocate(overlap,kinetic,nuclear,dipole,fock,density,lcaomo)
  deallocate(bas_norm)
  
  ! ---------------------------------------------------------------------------
  
  stop 'normal termination'

  
end program parse_namelists
