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

end module measurements
