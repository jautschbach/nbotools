subroutine overlap_matrix(iconv,smat,nsmat)

  ! calculate overlap matrix of AO basis functions

  ! literature: "Fundamentals of Molecular Integral Evaluation",
  !             JT Ferman, EF Valeev (FV), PDF found at Google.com 08/05

  use definitions
  use constants_parameters
  use molecule_basis

  implicit none

  integer(KINT), intent(in) :: iconv ! basis set convention used
  real(KREAL), intent(out) :: smat(nsmat)
  integer(KINT), intent(in) :: nsmat

  integer(KINT) :: ishell, ibas, iexp, iprim, icomp, iptr, inuc, ityp, li
  integer(KINT) :: jshell, jbas, jexp, jprim, jcomp, jptr, jnuc, jtyp, lj

  real(KREAL) zeta1, zeta2, coef1, coef2

  integer(KINT) :: dbg, i,j, nx, ii, x1, y1, z1, x2, y2, z2

  real(KREAL) :: sij, rtmp, t1, t2
  
  real(KREAL) :: xA, yA, zA, xB, yB, zB, abx, aby, abz, ab2, gamma

  real(KREAL) :: xP, yP, zP, pax, pay, paz, pbx, pby, pbz

  character(7) :: pname

  ! ---------------------------------------------------------------------------

  pname='overlap'
  dbg = 0 ! debug switch
  nx = nbas*(nbas+1)/2 ! number of unique overlap matrix elements

  if (iconv.le.0 .or. iconv.gt.3) &
     stop pname//': basis convention flag out of range'

  if(nx.gt.nsmat) stop 'overlap: nsmat < nx'

  !allocate (smat(nx)) ! output matrix: Overlap of contracted basis fcts.

  smat(:) = zero
  
  ii = 0 ! index for smat(:)

  ! sanity check re. basis spec.:
  ibas = 0
  do ishell = 1,nshell
    ibas = ibas+ncomp(ishell)
  end do
  !write(out,*) 'nbas, ibas = ',nbas,ibas
  if (ibas.ne.nbas) stop 'overlap: ncomp array broken'

  ! start outer loop over basis functions: counter ibas
  
  ibas = 0
  
  do ishell = 1,nshell
    do icomp = 1, ncomp(ishell)
      
      ibas = ibas + 1
      if (ibas.gt.nbas) stop 'overlap: error ibas'
      
      ityp = label(ibas)
      li = ityp/100 ! div(ityp,100) = angular momentum of AO
      if (li < 0 .or. li>lmax) stop 'overlap: error li'
      
      inuc = center(ibas)
      
      if (dbg>0)write(out,*) pname//': basis no., type, angmom, center: ',&
         & ibas, ityp, li, inuc

      ! get linear combinations of Cartesian primitives for ibas
      call bas_powers(iconv,ityp,nc1,nx1,ny1,nz1,cc1)
      

      ! for debug purposes, print primitive info from outer basis loop:
      if(dbg>0) then
        write(out,*) 'basis powers for x, y, z, and coeffs'
        do i = 1,nc1
          write (out,*) nx1(i), ny1(i), nz1(i), cc1(i)
        end do
        write (out,*) 'exponents and contraction coeffs'
        do iprim = 1,nprim(ishell)
          iexp = nptr(ishell) + iprim - 1
          if (iexp.le.0 .or. iexp.gt.nexp) stop 'overlap: error iexp'
          zeta1 = zeta(iexp)
          coef1 = ccoe(iexp,li)
          if (dbg>0) write(out,*) zeta1, coef1
        end do ! iprim
        write(out,*)
      end if ! dbg

      ! start inner loop over basis functions: counter jbas

      jbas = 0
      
      do jshell = 1,nshell
        do jcomp = 1, ncomp(jshell)
          
          jbas = jbas + 1
          if (jbas.gt.nbas) stop 'overlap: error jbas'
          
          jtyp = label(jbas)
          lj = jtyp/100 ! div(jtyp,100) = angular momentum of AO
          if (lj < 0 .or. lj>lmax) stop 'overlap: error lj'
          
          jnuc = center(jbas)
          
          ! skip this basis function if jbas > ibas, otherwise increment ii
          
          if (jbas > ibas) then
            cycle
          else
            ii = ii+1
          end if
          
          ! get linear combinations of Cartesian primitives for jbas
          call bas_powers(iconv,jtyp,nc2,nx2,ny2,nz2,cc2)
          
          ! apply Gaussian product theorem pt. I
          
          xA = coord(1,inuc)
          yA = coord(2,inuc)
          zA = coord(3,inuc)
          
          xB = coord(1,jnuc)
          yB = coord(2,jnuc)
          zB = coord(3,jnuc)
          
          if (dbg>1) write (out,*) pname//': centers',xa,ya,za,xb,yb,zb
          
          abx = xA - xB
          aby = yA - yB
          abz = zA - zB
          ab2 = abx**2 + aby**2 + abz**2
          
          ! calculate overlap in double loop over primitives,
          ! then contract to form the overlap between the pair
          ! ibas and jbas of contracted functions.
          ! skip inner loops if contraction coefficients are effectively 0
          
          do iprim = 1,nprim(ishell)
            iexp = nptr(ishell) + iprim - 1
            if (iexp.le.0 .or. iexp.gt.nexp) stop 'overlap: error iexp'
            zeta1 = zeta(iexp)
            coef1 = ccoe(iexp,li)
            if (abs(coef1).lt.tiny) cycle
            
            do jprim = 1,nprim(jshell)
              jexp = nptr(jshell) + jprim - 1
              if (jexp.le.0 .or. jexp.gt.nexp) stop 'overlap: error jexp'
              zeta2 = zeta(jexp)
              coef2 = ccoe(jexp,lj)
              if (abs(coef2).lt.tiny) cycle
              
              ! Gaussian product theorem pt. II
              
              gamma = (zeta1 + zeta2) ! (FV:2.33)
              if(gamma.le.zero) then
                write(out,'(a,e20.10)') pname//': gamma = ',gamma
                stop 'overlap: gamma < 0'
              endif
              
              xP = (zeta1*xA + zeta2*xB)/gamma ! (FV:2:33 and FV:2.41)
              yP = (zeta1*yA + zeta2*yB)/gamma
              zP = (zeta1*zA + zeta2*zB)/gamma
              
              if (dbg>1) write (out,*) &
                 ' product center:',xP, yP, zP
              
              pax = xP - xA
              pay = yP - yA
              paz = zP - zA
              
              pbx = xP - xB
              pby = yP - yB
              pbz = zP - zB
              
              ! double loop over basis function spherical components
              ! (if applicable)
              
              sij = zero
              
              do i = 1,nc1
                x1 = nx1(i); y1 = ny1(i); z1 = nz1(i)
                
                do j = 1,nc2
                  x2 = nx2(j); y2 = ny2(j); z2 = nz2(j)
                  
                  if (inuc.eq.jnuc) then
                    ! one-center case
                    
                    call  ovrlap_1c_prim(zeta1,zeta2, x1,y1,z1, x2,y2,z2,&
                       & rtmp)
                    if (dbg>1) write(out,*) zeta1,zeta2, x1,y1,z1,&
                       x2,y2,z2,rtmp
                    
                  else
                    ! two-center case:
                                      
                    call ovrlap_2c_prim(ab2,pax,pay,paz,pbx,pby,pbz,gamma, &
                       zeta1,zeta2,x1,y1,z1,x2,y2,z2,rtmp)
                    
                  end if ! one- or two-center case?

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
                    t1 = one/sqrt(t1)
                    t2 = one/sqrt(t2)

                  case (3)
                    ! Molden file convention: primitive normalization absorbed in
                    ! the contraction coefficients. Same as ORCA, but some
                    ! different conventions for d, f, etc.
                    t1 = one; t2 = one
                  end select
            
                  sij = sij + cc1(i)*cc2(j)*t1*t2*rtmp
                  
                end do ! j
              end do ! i
              
              smat(ii) = smat(ii) + coef1*coef2*sij
              
            end do ! jprim          
          end do ! iprim
          
        end do ! jcomp
      end do! jshell
      
    end do ! icomp
  end do ! ishell
        
  return

end subroutine overlap_matrix


! =============================================================================

subroutine ovrlap_1c_prim(zeta1,zeta2,l1,m1,n1,l2,m2,n2,s12)

  use definitions
  use constants_parameters
  use factorials
  implicit none
  
  !  purpose: calculate one-center overlap s12 = <f1|f2> between two 
  !  unnormalized Gaussian primitives, can also be used for normalization
  !  integral
  !
  !  the functions are of the type f = x**l y**m z**n exp(-zeta r**2)
  !
  !  ref:     (FV:2.22))
  !           
  !  output:  s12 - overlap integral
  !  input:   zeta1, zeta2  -  exponents of Gaussians
  !           l1, m1, n1    -  powers for x,y,z of f1
  !           l2, m2, n2    -  powers for x,y,z of f2

  integer(KINT), intent(in) :: l1,m1,n1,l2,m2,n2
  real(KREAL), intent(in) :: zeta1,zeta2
  real(KREAL), intent(out) :: s12
  
  integer(KINT) :: ll, mm, nn, ll2, mm2, nn2, ltot
  real(KREAL) :: rtmp

  real(KREAL), parameter :: r32 = 1.5_KREAL, pi32 = pi**r32
  logical, external :: odd, even

  ! --------------------------------------------------------------------------

  ! if the power of either x, y, or z is odd the integral is zero
  
  s12 = zero
  
  ll = l1 + l2
  mm = m1 + m2
  nn = n1 + n2
  
  if (odd(ll) .or. odd(mm) .or. odd(nn)) then
    return
  end if
  
  ! if all combined powers of x, y, z are even we can apply
  ! (FV:2.22) without the factors of a_i and a_j

  ll2 = ll/2
  mm2 = mm/2
  nn2 = nn/2
  ltot = ll2 + mm2 + nn2

  !      write (6,'(a,4i4,2e15.7)') 
  !     +   'll,mm,nn,ltot,a1,a2',ll,mm,nn,ltot,zeta1,zeta2

  rtmp = pi32 * dblfact(ll-1) * dblfact(mm-1) * dblfact(nn-1)
  !      write(6,*) pi32,rtmp
  rtmp = rtmp / (two**(ltot) * (zeta1+zeta2)**(float(ltot)+r32))

  s12 = rtmp
       
  return

end subroutine ovrlap_1c_prim

! =============================================================================

subroutine ovrlap_2c_prim(ab2,pax,pay,paz,pbx,pby,pbz,gamma,              &
                    zeta1,zeta2,l1,m1,n1,l2,m2,n2,s12)
!
!     purpose: calculate two-center overlap s12 = <f1|f2> between two 
!     Gaussian primitives
!
!     ref:     (FV:3.15)
!              
!     output:  s12 - overlap integral
!     input:   zeta1, zeta2  -  exponents of Gaussians
!              l1, m1, n1      -  powers for x,y,z of f1
!              l2, m2, n2      -  powers for x,y,z of f2
!              ab2             -  (FV:2:37) distance**2 between centers
!              gamma           -  (FV:2.33) sum of exponents
!              pax, pay, ...   -  (FV:2.41) location of product function
!                                 relative to function centers
!
  use definitions
  use constants_parameters
  use factorials
  implicit none

  real(KREAL), intent(in) ::  ab2,pax,pay,paz,pbx,pby,pbz,&
     zeta1,zeta2,gamma

  integer(KINT), intent(in) ::  l1,m1,n1,l2,m2,n2

  real(KREAL), intent(out) :: s12

  real(KREAL) :: prod,fk,exponent

  integer(KINT) :: i,iii,k

  real(KREAL) :: xIx, yIy, zIz, pg12 

  logical, external :: odd, even
  character(LCHARS) :: line

  ! ---------------------------------------------------------------------------

  xIx = zero
  yIy = zero
  zIz = zero
  
  pg12 = sqrt(PI/gamma) ! common factor in all integrals
  
  ! note: if (l1+l2) is odd the loop should go to
  ! (l1 + l2 -1)/2. Needs to be tested here, also for the
  ! other powers of x, y, z.
  
  ! what we also need to test is if the power of x, y, or z of
  ! the product function is zero. In that case we do not apply the
  ! clumsy expansion formula which might then lead to
  ! undetermined expressions involving 0**0
  
  if ((l1+l2) .eq.0) then
    xIx = pg12  
  else
    if (odd(l1+l2)) then
      iii = (l1 + l2 -1)/2
    else
      iii = (l1 + l2)/2
    end if
  
    do i = 0, iii
      
      if ((2*i-1).gt.mxfact) then
        stop 'ovrlap_2c_prim: factorial blown'
      end if
      
      prod = pg12 * dblfact(2*i-1) /( (two * gamma)**(i))
      
      fk = zero
      k = 2 * i
      call FkOverlap(k,l1,l2,pax,pbx,fk)
      
      xIx = xIx + prod * fk
    end do
  end if
  
  if ((m1+m2) .eq.0) then
    yIy = pg12
  else
    if (odd(m1+m2)) then
      iii = (m1 + m2 -1)/2
    else
      iii = (m1 + m2)/2
    end if
    
    do i = 0, iii

      if ((2*i-1).gt.40) then
        stop 'ovrlap_2c_prim: factorial blown'
      end if
      
      prod = pg12 * dblfact(2*i-1) /( (two * gamma)**(i))
    
      fk = zero
      k = 2 * i
      call FkOverlap(k,m1,m2,pay,pby,fk)
    
      yIy = yIy + prod * fk
    end do
  end if

  if ((n1+n2) .eq.0) then
    zIz = pg12
  else
    if (odd(n1+n2)) then
      iii = (n1 + n2 -1)/2
    else
      iii = (n1 + n2)/2
    end if
    
    do i = 0, iii
      
      if ((2*i-1).gt.40) then
        stop 'ovrlap_2c_prim: factorial blown'
      end if
      
      prod = pg12 * dblfact(2*i-1) /( (two * gamma)**(i))
    
      fk = zero
      k = 2 * i
      call FkOverlap(k,n1,n2,paz,pbz,fk)    
      zIz = zIz + prod * fk
    end do
  end if
  
  ! assemble the overlap integral (FV:3.12)
  exponent = zeta1 * zeta2 * ab2 / gamma
  
  if (exponent.lt.explim) then
    s12 = exp(-exponent)* xIx * yIy * zIz
  else
    s12 = zero
  end if
    
  return
    
end subroutine ovrlap_2c_prim

! =============================================================================

subroutine FkOverlap(k,l1,l2,pa,pb,fk)

  ! purpose: Evaluate helper function f_k for calculation of
  !          2-center overlap integrals between primitive Gaussians.
  !          f_k yields the coefficient of the power of x, y, or z
  !          of the product function according to the Gaussian product
  !          theorem (GPT).
  !
  ! ref:     (FV:2.46)
  !
  ! input:   k: index of Gaussian product expansion
  !          l1, l2: powers of x, y, or, z of the two functions
  !          pa, pb: distance x, y, or z of the 2 fcts. relative to GPT center
  ! output:  fk: value of f_k function = expansion coefficient
  
  use definitions
  use constants_parameters
  use factorials
  implicit none

  integer(KINT), intent(in) :: k, l1, l2
  real(KREAL), intent(in) ::  pa, pb
  real(KREAL), intent(out) :: fk

  real(KREAL) :: dnewt, fka, fkb
  integer(KINT) :: i, j  
  integer(KINT) :: q, qlow, qhigh

  ! ---------------------------------------------------------------------------

      qlow = max(-k,(k-2*l2))
      qhigh = min(k,(2*l1 - k))

      fka = zero
      fkb = zero

      ! next : equation (FV:2.46)
      do q = qlow,qhigh,2
        i = (k+q)/2
        j = (k-q)/2
        ! note: function dnewt tests internally if arguments are out of range
        fka = fka + dnewt(l1,i)*dnewt(l2,j)*(pa**(l1-i))*pb**(l2-j)
      end do
      
      fk = fka

      return

!     debug code, do not use for production:
      !     next: equation (FV:2.45)
      
      do i = 0, l1
        do j = 0, l2
          if (i + j .eq. k) then
            fkb = fkb + dnewt(l1,i)*dnewt(l2,j)* &
                (pa**(l1-i)) * (pb**(l2-j))
          end if
        end do
      end do

      write (out,*) 'fk:', fk
      if (fka .ne. fkb) then 
        write (6,*) 'fk not equal', fka, fkb
      end if
      
      return
      end

! =============================================================================

logical function odd(n)
  use definitions
  implicit none
  integer(KINT) ::  n

  if (n .eq. 0) then
    odd = .false.
  else
    if (mod(n,2) .eq.0) then
      odd = .false.
    else
      odd = .true.
    end if
  end if
end function odd

! =============================================================================

logical function even(n)
  use definitions
  implicit none
  integer(KINT) ::  n

  if (n .eq. 0) then
    even = .true.
  else
    if (mod(n,2) .eq.0) then
      even = .true.
    else
      even = .false.
    end if
  end if
end function even

! =============================================================================

