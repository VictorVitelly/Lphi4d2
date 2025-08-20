module measurements
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine correlation(phi,corr1,corr2)
    real(dp), dimension(N,N), intent(in) :: phi
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp), dimension(N,N), intent(inout) :: corr2
    real(dp), dimension(N) :: varphi
    integer(i4) :: i1,i2
    varphi=0._dp
    do i1=1,N
      do i2=1,N
        varphi(i1)=varphi(i1)+phi(i1,i2)
      end do
    end do
    varphi(:)=varphi(:)/real(N,dp)
    do i1=1,N
      corr1(i1)=corr1(i1)+varphi(i1)
      do i2=1,N
        corr2(i1,i2)=corr2(i1,i2)+varphi(i1)*varphi(i2)
      end do
    end do
  end subroutine correlation

  subroutine correlation_function(corr1,corr2,CF,CFprom)
    real(dp), dimension(N,Nmsrs), intent(in) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(in) :: corr2
    real(dp), dimension(N,N), intent(out) :: CF,CFprom
    real(dp), dimension(N) :: corr1prom,corr1delta
    real(dp), dimension(N,N) :: corr2prom,corr2delta
    integer(i4) :: i1,i2
    corr1prom=0._dp
    corr2prom=0._dp
    corr1delta=0._dp
    corr2delta=0._dp
    call mean_vector(corr1,corr1prom,corr1delta)
    call mean_matrix(corr2,corr2prom,corr2delta)
    do i1=1,N
      do i2=1,N
        CF(i1,i2)=corr2prom(i1,i2)-corr1prom(i1)*corr1prom(i2)
        CFprom(i1,i2)=Sqrt((corr2delta(i1,i2))**2+(corr1prom(i1)*corr1delta(i2))**2 +(corr1prom(i2)*corr1delta(i1) )**2)
      end do
    end do
  end subroutine correlation_function
  
  subroutine autocorrelation(m0,tmax,phi,montecarlos)
    integer(i4), intent(in) :: tmax,montecarlos
    real(dp), intent(in) :: m0
    real(dp), dimension(N,N), intent(inout) :: phi
    real(dp), dimension(tmax+1) :: auto,auto_delta
    real(dp) :: E(Nmsrs+tmax), auto1(Nmsrs)
    real(dp) :: E_ave,auto1_ave,autoj(tmax+1,Nauto)
    integer(i4) :: i,j,tt
    character(len=32) :: filename
    write(filename, '("data/autocorrmcs", I0, ".dat")') montecarlos
    !tt=tmax
    open(10, file = filename, status = 'replace')
    print*, 'File created at:', filename
    do j=1,Nauto
      do i=1,Nmsrs+tmax
        call cycles(m0,phi,montecarlos)
        !E(i)=S(m0,phi)/(real(N,dp)**2)
        E(i)=abs(Magnet(phi))/(real(N,dp)**2)
      end do
      call mean_0(E,E_ave )
      
      do tt=0,tmax
        do i=1,Nmsrs
          auto1(i)=E(i)*E(i+tt)
        end do
        call mean_0(auto1,auto1_ave)
        auto=auto1_ave-(E_ave**2)
        autoj(tt+1,j)=auto1_ave-(E_ave**2)
      end do
    end do
    do tt=0,tmax
      call mean_scalar(autoj(tt+1,:),auto(tt+1),auto_delta(tt+1))
      write(10,*) tt,auto(tt+1),auto_delta(tt+1)
    end do
    close(10)
  end subroutine autocorrelation

!Alternative autocorrelation without overlapping
  subroutine autocorrelation2(m0,tmax,phi)
    integer(i4), intent(in) :: tmax
    real(dp), intent(in) :: m0
    real(dp), dimension(N,N), intent(inout) :: phi
    integer(i4) :: i,j,tt,jcount
    real(dp) :: E(tmax), auto1(tmax,Nmsrs),auto2(Nmsrs)
    real(dp) :: auto1_ave(tmax),auto2_ave, auto(tmax)
    real(dp) :: auto1bin(tmax,Mbins),auto2bin(Mbins)
    real(dp) :: autoc(tmax), jackk(tmax), auto_delta(tmax)
    real(dp) :: Nbins
    open(10, file = 'data/autocorrv2.dat', status = 'replace')
    Nbins=real(Nmsrs,dp)-real(Nmsrs,dp)/real(Mbins,dp)

    jcount=0
    do j=1,sweeps
      call metropolis(m0,phi)
      if(j>2*thermalization .and. mod(j,eachsweep)==0) then
        jcount=jcount+1
        do i=1,tmax
          call metropolis(m0,phi)
          E(i)=S(m0,phi)/(real(N,dp)**2)
        end do
        do tt=0,tmax-1
          auto1(tt+1,jcount)=E(1)*E(1+tt)
        end do
        call mean_0(E(:),auto2(jcount) )
      end if
    end do

    do tt=0,tmax-1
      call mean_0(auto1(tt+1,:),auto1_ave(tt+1) )
    end do

    call mean_0(auto2(:),auto2_ave )
    auto(:)=auto1_ave(:)!-(auto2_ave**2)

    auto1bin=0._dp
    auto2bin=0._dp
    do j=1,Mbins
      do i=1,Nmsrs
        if(i .le. (j-1)*Nmsrs/Mbins) then
          auto1bin(:,j)=auto1bin(:,j)+auto1(:,i)
          auto2bin(j)=auto2bin(j)+auto2(i)
        else if(i > j*Nmsrs/Mbins) then
          auto1bin(:,j)=auto1bin(:,j)+auto1(:,i)
          auto2bin(j)=auto2bin(j)+auto2(i)
        end if
      end do
    end do
    auto1bin(:,:)=auto1bin(:,:)/(Nbins )
    auto2bin(:)=auto2bin(:)/(Nbins )
    jackk=0._dp
    do j=1,Mbins
      autoc(:)=auto1bin(:,j)!-(auto2bin(j)**2)
      jackk(:)=jackk(:)+(autoc(:)-auto(:) )**2
    end do
    auto_delta(:)=Sqrt(real(Mbins-1,dp)*jackk(:)/real(Mbins,dp) )

    do tt=0,tmax-1
      write(10,*) tt,auto(tt+1),auto_delta(tt+1)
    end do
    close(10)
  end subroutine autocorrelation2

end module measurements
