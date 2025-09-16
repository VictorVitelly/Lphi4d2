module statistics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  implicit none

contains

  subroutine random_phi(x,bound)
    implicit none
    real(dp),intent(out) :: x
    real(dp), intent(in) :: bound
    real(dp) :: y
    call random_number(y)
    x = 2._dp*bound*y -bound
  end subroutine random_phi

  subroutine cold_start(phi)
    real(dp), dimension(:,:), intent(out) :: phi
    phi=0.0_dp
  end subroutine cold_start

  subroutine hot_start(phi,hotphi)
    real(dp), dimension(N,N), intent(out) :: phi
    real(dp), intent(in) :: hotphi
    integer(i4) :: i1,i2
    do i1=1,N
      do i2=1,N
        call random_phi(phi(i1,i2),hotphi)
      end do
    end do
  end subroutine hot_start

  subroutine montecarlo(m0,dphi,phi,AR)
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(N,N), intent(inout) :: phi
    real(dp), intent(out) :: AR
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1,i2
    AR=0._dp
    do i1=1,N
      do i2=1,N
        call random_phi(deltaphi,dphi)
        phi2=phi(i1,i2)+deltaphi
        DS=DeltaS(m0,phi,i1,i2,phi2)
        !DS=DeltaS2(m0,phi,i1,i2,phi2)
        if(DS .le. 0._dp) then
          phi(i1,i2)=phi2
          AR=AR+1._dp
        else
          call random_number(r)
          p=Exp(-DS)
          AR=AR+p
          if(r < p ) then
            phi(i1,i2)=phi2
          end if
        end if
      end do
    end do
    AR=AR/real(N**2,dp)
  end subroutine montecarlo

  subroutine metropolis(m0,phi)
    real(dp), intent(in) :: m0
    real(dp), dimension(N,N), intent(inout) :: phi
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1,i2
    do i1=1,N
      do i2=1,N
        call random_phi(deltaphi,dphi_m)
        phi2=phi(i1,i2)+deltaphi
        DS=DeltaS(m0,phi,i1,i2,phi2)
        if(DS .le. 0._dp) then
          phi(i1,i2)=phi2
        else
          call random_number(r)
          p=Exp(-DS)
          if(r < p ) then
            phi(i1,i2)=phi2
          end if
        end if
      end do
    end do
  end subroutine metropolis
  
  subroutine flip_sign(phi,i)
    real(dp), dimension(N,N), intent(inout) :: phi
    integer(i4), intent(in) :: i
    integer(i4) :: i1,i2
    if(i>thermalization .and. mod(i,10)==0) then
      do i1=1,N
        do i2=1,N
          phi(i1,i2)=-phi(i1,i2)
        end do
      end do
    end if
  end subroutine flip_sign
  
  subroutine cluster(phi) 
    real(dp), dimension(N,N),intent(inout) :: phi
    integer(i4), dimension(N,N) :: spin
    logical, dimension(N,N) :: bond_x,bond_y
    integer(i4) :: i,j,label(N,N),parent(N*N),next_label,left_label,up_label
    logical, allocatable :: flip_cluster(:)
    real(dp) :: beta,r,p
    
    spin(:,:)=nint(sign(1._dp,phi(:,:)),i4)
    
    do i=1,N
      do j=1,N
        if(spin(i,j)==spin(mod(i,N)+1,j) ) then
          beta=abs(phi(i,j))*abs(phi(mod(i,N)+1,j))
          p=1._dp-exp(-2._dp*beta )
          call random_number(r)
          bond_x(i,j)=(r<p)
        else
          bond_x(i,j)=.false.
        end if
        if(spin(i,j)==spin(i,mod(j,N)+1) ) then
          beta=abs(phi(i,j))*abs(phi(i,mod(j,N)+1))
          p=1._dp-exp(-2._dp*beta )
          call random_number(r)
          bond_y(i,j)=(r<p)
        else
          bond_y(i,j)=.false.
        end if
      end do
    end do

    label(:,:)=0
    do i=1,N*N
      parent(i)=i
    end do
    next_label=1
    left_label=0
    up_label=0

    do i=1,N
      do j=1,N
        left_label=0
        up_label=0
        if(i>1 .and. bond_x(i-1,j) ) then
          left_label=label(i-1,j)
        end if
        if(j>1 .and. bond_y(i,j-1) ) then
          up_label=label(i,j-1)
        end if
        if(left_label==0 .and. up_label==0) then
          label(i,j)=next_label
          next_label=next_label+1  
        else if(left_label /= 0 .and. up_label==0) then
          label(i,j)=left_label
        else if(left_label== 0 .and. up_label/=0) then
          label(i,j)=up_label
        else
          label(i,j)=min(left_label,up_label)
          call union(left_label,up_label,parent)
        end if
      end do
    end do
    
    do j=1,N
      if(bond_x(N,j) ) then
        call union(label(1,j),label(N,j),parent )
      end if
    end do
    do i=1,N
      if(bond_y(i,N) ) then
        call union(label(i,1),label(i,N),parent )
      end if 
    end do
    do i=1,N
      do j=1,N
        label(i,j)=find(label(i,j),parent)
      end do
    end do

    allocate(flip_cluster(next_label) )
    flip_cluster(:)=.false.
    do i=1,next_label-1
      call random_number(r)
      flip_cluster(i)=(r<0.5_dp)
    end do
    
    do i=1,N
      do j=1,N
        if(flip_cluster(label(i,j))) then
          phi(i,j)=-phi(i,j)
        end if
      end do
    end do
    deallocate(flip_cluster)
    
  end subroutine cluster
  
  subroutine cycles(m0,phi,k)
    real(dp), intent(in) :: m0
    real(dp), dimension(N,N), intent(inout) :: phi
    integer(i4), intent(in) :: k
    integer(i4) :: i
    call cluster(phi)
    do i=1,k
      call metropolis(m0,phi)
    end do
    !call cluster(phi)
  end subroutine cycles


  !Statistics for measurements
  !
  !
  !
  
  subroutine jackknife(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp), allocatable :: xmean(:), delta_y(:)
    integer(i4) :: k,Narr,i,j
      Narr=size(x)
      allocate(delta_y(size(Mbin)))
      do j=1,size(Mbin)
        allocate(xmean(Mbin(j)))
        jackk=0._dp
        xmean=0._dp
        do i=1,Mbin(j)
          do k=1,Narr
            if(k .le. (i-1)*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            else if(k > i*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            end if
          end do
          xmean(i)=xmean(i)/(real(Narr,dp) -real(Narr/Mbin(j),dp))
        end do
        do k=1,Mbin(j)
          jackk=jackk+(xmean(k)-y )**2
        end do
        delta_y(j)=Sqrt(real(Mbin(j)-1,dp)*jackk/real(Mbin(j),dp))
        deallocate(xmean)
      end do
      deltay=maxval(delta_y)
  end subroutine jackknife

  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i4) :: k,Narr
    Narr=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,Narr
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(Narr-1,dp)
    deltay=Sqrt(variance/real(Narr,dp))
  end subroutine standard_error

  subroutine mean_0(x,y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i4) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
  end subroutine mean_0

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    integer(i4) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
    !call standard_error(x,y,deltay)
    call jackknife(x,y,deltay)
  end subroutine mean_scalar

  subroutine mean_vector(x,y,deltay)
    real(dp), dimension(:,:), intent(in) :: x !x(N,Nmsrs or Nmsrs2)
    real(dp), dimension(size(x,dim=1)), intent(out) :: y,deltay
    integer(i4) :: i1,Narr
    y=0._dp
    deltay=0._dp
    Narr=size(x,dim=1)
    do i1=1,Narr
      call mean_scalar(x(i1,:),y(i1),deltay(i1))
    end do
  end subroutine mean_vector

  subroutine mean_matrix(x,y,deltay)
    real(dp), dimension(:,:,:), intent(in) :: x !x(N,N,Nmsrs or Nmsrs2)
    real(dp), dimension(N,N), intent(out) :: y,deltay
    integer(i4) :: i1,i2
    y=0._dp
    deltay=0._dp
    do i1=1,N
      do i2=1,N
      call mean_scalar(x(i1,i2,:),y(i1,i2),deltay(i1,i2))
      end do
    end do
  end subroutine mean_matrix

  subroutine mean_phi(x,y)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i4) :: i1,i2
    y=0._dp
    do i1=1,N
      do i2=1,N
        y=y+x(i1,i2)
      end do
    end do
    y=y/real(size(x,dim=1)*size(x,dim=2),dp)
  end subroutine mean_phi

  subroutine histogram(x,A1,A2)
  real(dp), dimension(:,:), intent(in) :: x
  integer(i4), dimension(bins), intent(inout) :: A1
  real(dp), dimension(bins), intent(in) :: A2
  integer(i4) :: i,j,k
  do i=1,bins
    do j=1,size(x,dim=1)
      do k=1,size(x,dim=2)
        if(x(j,k) .le. real(A2(i),dp)+binwidth/2._dp .and. x(j,k)>real(A2(i),dp)-binwidth/2._dp ) then
          A1(i)=A1(i)+1
          cycle
        end if
      end do
    end do
  end do
  end subroutine histogram

  subroutine histogram2(x,A1,A2)
  real(dp), dimension(:,:), intent(in) :: x
  integer(i4), dimension(bins), intent(inout) :: A1
  real(dp), dimension(bins), intent(in) :: A2
  integer(i4) :: i
  real(dp) :: y
  call mean_phi(x,y)
  do i=1,bins
    if(y .le. real(A2(i),dp)+binwidth/2._dp .and. y>real(A2(i),dp)-binwidth/2._dp ) then
      A1(i)=A1(i)+1
    cycle
    end if
  end do
  end subroutine histogram2

end module statistics
