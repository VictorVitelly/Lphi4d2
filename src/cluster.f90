  subroutine cluster(phi)
    real(dp), dimension(N,N), intent(inout) :: phi
    integer(i4), dimension(N,N) :: spin
    logical, dimension(N,N) :: bond_x,bond_y
    integer(i4) :: i,j,label(N,N),parent(N*N),next_label,left_label,up_label
    logical, allocatable :: flip_cluster(:)
    real(dp) :: beta,r,p

    spin(:,:)=nint(sign(1._dp,phi(:,:)),32)

    do i=1,N
      do j=1,N
        if(spin(i,j)==spin(mod(i,N)+1,j) ) then
          beta=abs(phi(i,j))*abs(phi(mod(i,N)+1,j))
          p=1._dp-exp(-beta*2._dp )
          call random_number(r)
          bond_x(i,j)=(r<p)
        else
          bond_x(i,j)=.false.
        end if
        if(spin(i,j)==spin(i,mod(j,N)+1) ) then
          beta=abs(phi(i,j))*abs(phi(i,mod(j,N)+1))
          p=1._dp-exp(-beta*2._dp )
          call random_number(r)
          bond_y(i,j)=(r<p)
        else
          bond_y(i,j)=.false.
        end if
      end do
    end do

    label(:,:)=0
    do i=1,N
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
        if(j>1 .and. bond_x(i,j-1) ) then
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
    do i=1,N
      do j=1,N
        if(label(i,j)/=0) then
          label(i,j)=find(label(i,j),parent)
        end if
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


  end subroutine cluster
