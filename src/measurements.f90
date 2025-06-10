module measurements
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine initialize1(susc1,susc2)
    real(dp), intent(inout),dimension(Nmsrs) :: susc1,susc2
      susc1=0._dp
      susc2=0._dp
  end subroutine initialize1

  subroutine initialize2(corr1,corr2)
    real(dp), dimension(N,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(inout) :: corr2
      corr1=0._dp
      corr2=0._dp
  end subroutine initialize2

  subroutine measure(M,Sprom,susc1,susc2,heat1,heat2)
    real(dp), intent(in) :: M, Sprom
    real(dp), intent(out) :: susc1,susc2,heat1,heat2
    susc1=(M**2)
    susc2=abs(M)
    heat1=(Sprom**2)
    heat2=Sprom !abs(Sprom)
  end subroutine measure

  subroutine measure2(M,susc1,susc2)
    real(dp), intent(in) :: M
    real(dp), intent(out) :: susc1,susc2
    susc1=(M**2)
    susc2=abs(M)
  end subroutine measure2

  subroutine divideN(x1,x2)
    real(dp), intent(inout) :: x1,x2
    x1=x1/real(Nmsrs,dp)
    x2=x2/real(Nmsrs,dp)
  end subroutine divideN

  subroutine correlation(phi,k,corr1,corr2)
    real(dp), dimension(N,N), intent(in) :: phi
    integer(i4), intent(in) :: k
    real(dp), dimension(N,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(inout) :: corr2
    real(dp), dimension(N) :: varphi
    integer(i4) :: i1,i2
    varphi=0._dp
    do i1=1,N
      do i2=1,N
        varphi(i1)=varphi(i1)+phi(i1,i2)
      end do
    end do
    do i1=1,N
      corr1(i1,k)=varphi(i1)
      do i2=1,N
        corr2(i1,i2,k)=varphi(i1)*varphi(i2)
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

   subroutine autocorrelation(m0,tmax,phi)
    integer(i4), intent(in) :: tmax
    real(dp), intent(in) :: m0
    real(dp), dimension(N,N), intent(inout) :: phi
    real(dp) :: auto,auto_delta
    real(dp), dimension(Nmsrs+tmax) :: E
    real(dp), dimension(Nmsrs) :: auto1
    real(dp), allocatable :: auto1j(:),Ej(:)
    real(dp) :: E_ave,auto1_ave,autoj,jackk
    real(dp) :: ARR
    integer(i4) :: i,j,tt
    open(70, file = 'data/autocorr.dat', status = 'replace')

    do i=1,Nmsrs+tmax
      call montecarlo(m0,dphi,phi,ARR)
      E(i)=S(m0,phi)/(real(N,dp)**2)
    end do
    call mean_0(E,E_ave)

    do tt=0,tmax
      auto1=0._dp
      auto=0._dp
      auto_delta=0._dp

      do i=1,Nmsrs
        auto1(i)=E(i)*E(i+tt)
      end do
      call mean_0(auto1,auto1_ave)
      auto=auto1_ave-(E_ave**2)
      allocate(auto1j(Mbins) )
      allocate(Ej(Mbins) )
      auto1j=0._dp
      Ej=0._dp
      do j=1,Mbins
        do i=1,Nmsrs
          if(i .le. (j-1)*Nmsrs/Mbins) then
            auto1j(j)=auto1j(j)+auto1(i)
            Ej(j)=Ej(j)+E(i)
          else if(i > j*Nmsrs/Mbins) then
            auto1j(j)=auto1j(j)+auto1(i)
            Ej(j)=Ej(j)+E(i)
          end if
        end do
      end do
      auto1j=auto1j/(real(Nmsrs,dp)-real(Nmsrs,dp)/real(Mbins,dp) )
      Ej=Ej/(real(Nmsrs,dp)-real(Nmsrs,dp)/real(Mbins,dp) )
      jackk=0._dp
      do j=1,Mbins
        autoj=0._dp
        autoj=auto1j(j)-(Ej(j)**2)
        jackk=jackk+(autoj-auto )**2
      end do
      auto_delta=Sqrt(real(Mbins-1,dp)*jackk/real(Mbins,dp) )

      write(70,*) tt, auto, auto_delta
      deallocate(auto1j,Ej)
    end do
    close(70)
  end subroutine autocorrelation


  subroutine autocorrelation2(m0,tmax,phi)
    integer(i4), intent(in) :: tmax
    real(dp), intent(in) :: m0
    real(dp), dimension(N,N), intent(inout) :: phi
    integer(i4) :: i,j,tt
    real(dp) :: E(tmax), auto1(tmax,Nmsrs),auto2(Nmsrs)
    real(dp) :: auto1_ave(tmax),auto2_ave, auto(tmax)
    real(dp) :: auto1bin(tmax,Mbins),auto2bin(Mbins)
    real(dp) :: autoc(tmax), jackk(tmax), auto_delta(tmax)
    real(dp) :: Nbins
    open(70, file = 'data/autocorr.dat', status = 'replace')
    Nbins=real(Nmsrs,dp)-real(Nmsrs,dp)/real(Mbins,dp)

    do j=1,Nmsrs
      do i=1,tmax
        call metropolis(m0,phi)
        E(i)=S(m0,phi)/(real(N,dp)**2)
      end do
      do tt=0,tmax-1
        auto1(tt+1,j)=E(1)*E(1+tt)
      end do
      call mean_0(E(:),auto2(j) )
    end do
    do tt=0,tmax-1
      call mean_0(auto1(tt+1,:),auto1_ave(tt+1) )
    end do
    call mean_0(auto2(:),auto2_ave )
    auto(:)=auto1_ave(:)-(auto2_ave**2)

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
      autoc(:)=auto1bin(:,j)-(auto2bin(j)**2)
      jackk(:)=jackk(:)+(autoc(:)-auto(:) )**2
    end do
    auto_delta(:)=Sqrt(real(Mbins-1,dp)*jackk(:)/real(Mbins,dp) )

    do tt=0,tmax-1
      write(70,*) tt,auto(tt+1),auto_delta(tt+1)
    end do
    close(70)
  end subroutine autocorrelation2

end module measurements
