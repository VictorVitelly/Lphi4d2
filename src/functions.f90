module functions
    use iso_fortran_env, only : dp => real64, i4 => int32
    use parameters
    implicit none

contains

  function iv(i)
    integer(i4), intent(in) :: i
    integer(i4) :: iv
    if(i==N+1) then
      iv=1
    else if(i==0) then
      iv=N
    else
      iv=i
    end if
  end function iv

  function lagrangian(m02,lamb0,phi,i1,i2)
    real(dp), intent(in) :: m02,lamb0
    real(dp), dimension(:,:), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp) :: lagrangian
    lagrangian=( (phi(iv(i1+1),i2)-phi(i1,i2) )**2 +(phi(i1,iv(i2+1)) -phi(i1,i2))**2 &
              &+m02*phi(i1,i2)**2 +lamb0*phi(i1,i2)**4 /2._dp )/2._dp
  end function lagrangian

  function S(m02,lamb0,phi)
    real(dp), intent(in) :: m02,lamb0
    real(dp), dimension(:,:), intent(in) :: phi
    real(dp) :: S
    integer(i4) :: i1,i2,Narr
    Narr=size(phi,dim=1)
    S=0._dp
    do i1=1,Narr
      do i2=1,Narr
        S=S+lagrangian(m02,lamb0,phi,i1,i2)
      end do
    end do
  end function S

  function DeltaS(m02,lamb0,phi,i1,i2,phi2)
    real(dp), intent(in) :: m02,lamb0
    real(dp), dimension(:,:), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp), intent(in) :: phi2
    real(dp) :: DeltaS
    real(dp) :: DSa,DSb
    DSa=(2._dp+m02/2._dp )*(phi2**2-phi(i1,i2)**2 )+lamb0*(phi2**4-phi(i1,i2)**4)/4._dp
    DSb=-(phi2-phi(i1,i2))*(phi(iv(i1+1),i2)+phi(i1,iv(i2+1))+phi(iv(i1-1),i2)+phi(i1,iv(i2-1)))
    DeltaS=DSa +DSb
  end function DeltaS

  function DeltaS2(m02,lamb0,phi,i1,i2,phi2)
    real(dp), intent(in) :: m02,lamb0
    real(dp), dimension(:,:), intent(in) :: phi
    real(dp), dimension(size(phi,dim=1),size(phi,dim=1)) :: phi_y
    integer(i4), intent(in) :: i1,i2
    real(dp), intent(in) :: phi2
    real(dp) :: DeltaS2
    phi_y=phi
    phi_y(i1,i2)=phi2
    DeltaS2=S(m02,lamb0,phi_y)-S(m02,lamb0,phi)
  end function DeltaS2

  function mean(phi)
    real(dp), dimension(N,N), intent(in) :: phi
    integer(i4):: i1,i2
    real(dp) :: mean
    mean=0._dp
    do i1=1,N
      do i2=1,N
        mean=mean+phi(i1,i2)
      end do
    end do
  end function mean

  recursive function find(x,parent) result(out)
    integer(i4), intent(in) :: x
    integer(i4), intent(inout) :: parent(:)
    integer(i4) :: out
    if(parent(x) /= x) then
      parent(x)=find(parent(x),parent )
    end if
    out=parent(x)
  end function find

  subroutine union(x,y,parent)
    integer(i4),intent(in) :: x,y
    integer(i4),intent(inout) :: parent(:)
    integer :: root_x, root_y
    root_x=find(x,parent)
    root_y=find(y,parent)
    if (root_x /= root_y) then
      parent(root_y)=root_x
    end if
  end subroutine union

end module functions
