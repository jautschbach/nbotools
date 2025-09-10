module factorials

  use definitions
  
  implicit none

  integer(KINT), parameter :: mxfact = 40
  
  real(KREAL), save :: fact(0:mxfact), dblfact(-1:mxfact)

end module factorials

! ==============================================================================

subroutine init_factorials

  use definitions
  use factorials

  integer(KINT) :: i

  fact(0) = 1.d0
  
  do i=1,mxfact
    fact(i) = fact(i-1)*float(i)
  end do
 
  dblfact(-1) = 1.0d0
  dblfact(0)  = 1.0d0
  dblfact(1)  = 1.0d0
  dblfact(2)  = 2.0d0
  do i=3,mxfact
    dblfact(i) = float(i)*dblfact(i-2)
  end do

  return
end subroutine init_factorials

! ==============================================================================

real(KREAL) function dnewt(l1,l2)
  use definitions
  use factorials
  implicit none
  integer(KINT) :: l1,l2

  if (l1.lt.0 .or. l1.gt.mxfact .or. l2.lt.0 .or. l2.gt.mxfact) &
     stop 'dnewt: l1 or l2 out of range'
  if ((l1-l2).lt.0 .or. (l1-l2).gt.mxfact) &
     stop 'dnewt: l1-l2 out of range'
  
  dnewt = fact(l1)/fact(l2)/fact(l1-l2)
  
  return
end function dnewt
    

