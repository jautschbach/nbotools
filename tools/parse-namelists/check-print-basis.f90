subroutine check_print_basis(iconv,found,lp)

  ! calculate diagonal elements of overlap matrix for the basis
  ! according to different basis conventions (iconv). If according to
  ! a chosen convention the basis is not normalized, then found=0. if
  ! the basis is normalized, found=1. if lp=.true. then print basis set
  ! information to output

  use definitions
  use constants_parameters
  use molecule_basis

  implicit none

  integer(KINT), intent(in) :: iconv ! convention for basis normalization
  
  integer(KINT), intent(out) :: found ! =1 if basis is normalized, 0 otherwise
  
  logical, intent(in) :: lp ! print switch, 'true' if we want to print
                            ! the basis info
  
  logical :: warn

  integer(KINT) :: ishell, ibas, iexp, iprim, icomp, iptr, inuc, ityp, li

  integer(KINT) :: jexp, jprim

  real(KREAL) zeta1, coef1, zeta2, coef2

  integer dbg, i,j, nptot, ii, x1, y1, z1, x2, y2, z2

  real(KREAL) :: sii, rtmp, t1, t2

  character(17) pname


  ! ---------------------------------------------------------------------------

  dbg = 0  ! debug switch
  pname = 'check_print-basis' ! name of subroutine in debug output
  warn = .false.

  bas_norm(:) = zero

  if (.not.lp) write (out,'(1x,a,i0)') 'testing basis convention no. ',iconv
  if (lp) write(out,'(/1x,21(''='')/1x,''Basis Set Information''/1x,21(''=''))')

  select case(iconv)
  case (1) ! ORCA convention
    if (lp) write(out,*) '(using ORCA file-47 convention)'
    continue
  case(2) ! Gaussian convention
    if (lp) write(out,*) '(using Gaussian file-47 convention)'
    continue
  case(3) ! Molden2AIM convention
    if (lp) write(out,*) '(using Molden2AIM file-47 convention)'
    continue
  case default
    stop pname//': basis convention flag out of range'
  end select

  ! sanity check re. basis spec.:
  ibas = 0
  do ishell = 1,nshell
    ibas = ibas+ncomp(ishell)
  end do
  !write(out,*) 'nbas, ibas = ',nbas,ibas
  if (ibas.ne.nbas) stop pname//': ncomp array broken'

  ! start loop over basis functions: counter ibas
  
  ibas = 0
  nptot = 0 ! total number of primitives
  
  do ishell = 1,nshell
    do icomp = 1, ncomp(ishell)
      
      ibas = ibas + 1
      if (ibas.gt.nbas) stop pname//': error ibas'
      
      ityp = label(ibas)
      li = ityp/100 ! div(ityp,100) = angular momentum of AO
      if (li < 0 .or. li>lmax) stop pname//': error li'
      
      inuc = center(ibas)

      ! get linear combinations of Cartesian primitives for ibas
      call bas_powers(iconv,ityp,nc1,nx1,ny1,nz1,cc1)

      if (lp) write(out,'(1x,79(''-'')/1x,a,t9,a,t14,a,t17,a,&
         t24,a,t30,a,t34,a,t38,a,t47,a/1x,79(''-''))') &
         'Fct.#', 'Type', 'L', 'Center','Pow.','x','y','z','Coef.'
      
      if (lp) write(out,'(1x,i0,t9,i0,t14,i0,t17,i0, &
         t30,i0,t34,i0,t38,i0,t42,f15.8)') ibas, ityp, li, inuc, &
         nx1(1), ny1(1), nz1(1), cc1(1)

      if (nc1.gt.1 .and. lp) then
        do i = 2,nc1
          write (out,'(t30,i0,t34,i0,t38,i0,t42,f15.8)') &
             nx1(i), ny1(i), nz1(i), cc1(i)
        end do
      end if

      if (lp) write (out,'(1x,a,t18,a,t43,a)') &
         'Contraction:','Exponent','Coefficient'

      sii = zero ! self-overlap of basis function
      
      do iprim = 1,nprim(ishell)
        iexp = nptr(ishell) + iprim - 1
        if (iexp.le.0 .or. iexp.gt.nexp) stop pname//': error iexp'
        zeta1 = zeta(iexp)
        coef1 = ccoe(iexp,li)
        if (lp) write(out,'(t10,f20.10,t35,f20.10)') zeta1, coef1

        do jprim = 1,nprim(ishell)
          jexp = nptr(ishell) + jprim - 1
          if (jexp.le.0 .or. jexp.gt.nexp) stop pname//': error jexp'
          zeta2 = zeta(jexp)
          coef2 = ccoe(jexp,li)

          do i = 1,nc1
            do j = 1,nc1
              x1 = nx1(i); y1 = ny1(i); z1 = nz1(i)
              x2 = nx1(j); y2 = ny1(j); z2 = nz1(j)
              
              call  ovrlap_1c_prim(zeta1,zeta2, x1,y1,z1, x2,y2,z2,&
                 & rtmp)
              
              if (dbg>1) write(out,*) zeta1,zeta2, x1,y1,z1,&
                 x2,y2,z2,rtmp

              select case(iconv)
              case (1)
                ! ORCA convention: primitive normalization absorbed in
                ! the contraction coefficients
                t1 = one; t2 = one
                
              case(2)
                ! Gaussian convention: contraction coefficients assume
                ! normalized primitives
                call  ovrlap_1c_prim(zeta1,zeta1, x1,y1,z1, x1,y1,z1,&
                   & t1)
                call  ovrlap_1c_prim(zeta2,zeta2, x2,y2,z2, x2,y2,z2,&
                  & t2)

                if (t1.lt.zero .or. t2.lt.zero) then
                  write(err,*) ibas, iprim,jprim,i,j,t1,t2
                  stop pname//': primitive self-overlap negative'
                end if

                t1 = sqrt(t1)
                t2 = sqrt(t2)

                ! the check below was previously made prior to taking
                ! the square roots (they were taken later). However,
                ! for all-electron relativistic basis sets with very
                ! large exponents, this could trigger the error
                ! exit. Therefore, we're now testing the smallness for
                ! the square root of the self-overlap
                
                if (t1.lt.tnsy .or. t2.lt.tnsy) then
                  write(err,*) ibas, iprim,jprim,i,j,t1,t2
                  stop pname//': primitive overlap too small'
                end if
                
                t1 = one/t1
                t2 = one/t2

              case (3)
                ! Molden2AIM convention: primitive normalization absorbed in
                ! the contraction coefficients. Same as ORCA, but some
                ! different conventions for d, f, etc.
                t1 = one; t2 = one
              end select
             
              sii = sii + cc1(i)*cc1(j) * coef1 * coef2 * t1 * t2 * rtmp
            end do ! j 
            
          end do ! i

        end do ! jprim
        
      end do ! iprim

      if (lp) write(out,'(1x,a,f9.5)') 'Normalization: <f|f> = ',sii

      bas_norm(ibas) = sii
      if (abs(sii-one).gt.small) warn = .true. ! basis fct. not normalized!
        
      if (lp) write(out,*)
      
    end do ! icomp
  end do ! ishell
  
  if (lp) write(out,'(/1x,79(''-''))')

  if (.not.warn) then
    found = 1 ! found normalized basis convention according to iconv value
  else
    found = 0 ! did not find normalized basis convention
  end if
        
  return

end subroutine check_print_basis



