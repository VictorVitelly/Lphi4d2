program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none
  !integer(i4) kkkk

  call cpu_time(starting)
  
  !Only measure the acc. rate, must modify dphi
  !call acceptance_rate(-0.1_dp,lambda0)

  !Thermalization history and autocorretion functions
  !write(*,*) 'How many Metropolis per cycle?' 
  !read(*,*) kkkk
  !write(*,*) 'Correctly read k=',kkkk
  !call thermalize(-0.18_dp,lambda0,4)
  !call time_test(-1.25_dp,lambda0,kkkk)

  !Histogram
  !call make_histogram(-1.4_dp,lambda0)

  !Measure action, magnetization, susceptibility and heat cap.
  !call vary_m0(-1.6_dp,-0.8_dp,lambda0,40)
  !call vary_m0(-1.6_dp,-0.6_dp,lambda0,79)
  !call vary_m0(-1.4_dp,-1.15_dp,lambda0,42)
  !call vary_m0(-1.31_dp,-1.21_dp,lambda0,10)
  !call vary_m0(-1.4_dp,-1.0_dp,lambda0,11)  
  
  !Measure correlation function
  call correlate(-1.6_dp,-1.0_dp,lambda0,10)
  !call correlate(-0.5_dp,0.15_dp,0.45_dp,16)
  

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

contains

  subroutine thermalize(m0,lamb0,montecarlos)
  real(dp), intent(in) :: m0,lamb0
  integer(i4),intent(in) :: montecarlos
  real(dp), allocatable :: phi(:,:)
  integer(i4) :: i
  open(10, file = 'data/history.dat', status = 'replace')
  allocate(phi(N,N))
  call hot_start(phi,hotphi)
  !call cold_start(phi)
  do i=1,2*thermalization
    !50 sweeps for L=8, 500 sweeps for L=64
    if(i==1 .or. mod(i,eachsweep)==0) then
      write(10,*) i,",", S(m0,lamb0,phi)/real(N**2,dp)
    end if
    call cycles(m0,lamb0,phi,montecarlos)
  end do
  call autocorrelation(m0,lamb0,201,phi,montecarlos)
  close(10)
  deallocate(phi)
  end subroutine thermalize

  !Measure acceptance rate for different DeltaPhi, recommended start at dphi=0.1
  subroutine acceptance_rate(m0,lamb0)
  real(dp), allocatable :: phi(:,:)
  real(dp), intent(in) :: m0,lamb0
  real(dp) :: ARR,AR_ave,AR_delta
  real(dp), dimension(10000) :: AR
  real(dp),dimension(201,10) :: hist
  integer(i4) :: i,k,k2,j
  open(10, file = 'data/history.dat', status = 'replace')
  dphi=0.1_dp
  hotphi=2._dp*dphi
  do j=1,10
    k2=0
    allocate(phi(N,N))
    call cold_start(phi)
    !call hot_start(phi,hotphi)
    do i=1,20000
      if(i==1 .or. mod(i,500)==0) then
        k2=k2+1
        hist(k2,j)=S(m0,lamb0,phi)/real(N**2,dp)
      end if
      call montecarlo(m0,lamb0,dphi,phi,ARR)
    end do
    k=0
    AR=0._dp
    do i=1,10000000
      call montecarlo(m0,lamb0,dphi,phi,ARR)
      if( mod(i,1000)==0) then
        k=k+1
        AR(k)=ARR
      end if
      call flip_sign(phi,i)
    end do
    write(*,*) k
    call mean_scalar(AR,AR_ave,AR_delta)
    write(*,*) dphi, ",", AR_ave,",", AR_delta
    deallocate(phi)
    dphi=dphi+0.1_dp
    hotphi=2._dp*dphi
  end do
  k2=0
  do i=1,5*thermalization
    if(i==1 .or. mod(i,500)==0) then
      k2=k2+1
      write(10,*) i, hist(k2,:)
    end if
  end do
  close(10)
  end subroutine acceptance_rate

  subroutine make_histogram(m0,lamb0)
  real(dp), intent(in) :: m0,lamb0
  real(dp), allocatable :: phi(:,:)
  real(dp), allocatable :: A2(:)
  integer(i4), allocatable :: A1(:)
  real(dp) :: ARR,norm
  integer(i4) :: i,k
  open(50, file = 'data/histogram.dat', status = 'replace')
  allocate(phi(N,N))
  allocate(A1(bins))
  allocate(A2(bins))
  do i=1,bins
    A2(i)=minn+binwidth/2._dp+real(i-1,dp)*binwidth
  end do
  A1=0
  !call cold_start(phi)
  call hot_start(phi,hotphi)
  k=0

  do i=1,sweeps
    call montecarlo(m0,lamb0,dphi,phi,ARR)
    if(i>thermalization .and. mod(i,eachsweep)==0) then
      k=k+1
      call histogram(phi,A1,A2)
      call flip_sign(phi,i)
    end if
  end do

  norm=0._dp
  do i=1,bins
    norm=norm+A1(i)
  end do
  norm=norm*(real(maxx-minn,dp) )/real(bins,dp)
  do i=1,bins
    write(50,*) A2(i), ",", A1(i)/norm, ",", sqrt( real(A1(i),dp) )/norm
  end do
  deallocate(phi)
  close(50)
  end subroutine make_histogram

  subroutine correlate(mi,mf,lamb0,Nts)
  real(dp), intent(in) :: mi,mf,lamb0
  !subroutine correlate(m0,lambi,lambf,Nts)
  !real(dp), intent(in) :: m0,lambi,lambf
  real(dp) :: x0
  integer(i4) :: Nts
  real(dp), allocatable :: phi(:,:),corr1(:),corr2(:,:),CF(:,:),CF_ave(:,:),CF_delta(:,:)
  integer(i4) :: i,j,k,i2
  open(60, file = 'data/corrfunc.dat', status = 'replace')
  allocate(phi(N,N))
  allocate(corr1(N))
  allocate(corr2(N,N))
  allocate(CF(N,Nmsrs2))
  allocate(CF_ave(N,Nts))
  allocate(CF_delta(N,Nts))

  !call cold_start(phi)
  do k=1,Nts
    CF(:,:)=0._dp
    x0=mi+(mf-mi)*real(k-1,dp)/real(Nts-1,dp)
    !x0=lambi+(lambf-lambi)*real(k-1,dp)/real(Nts-1,dp)
    write(*,*) x0
    call cold_start(phi)
    do j=1,2*thermalization
      call cycles(x0,lamb0,phi,4)
      !call cycles(m0,x0,phi,4)
    end do
    do j=1,Nmsrs2
      corr1(:)=0._dp
      corr2(:,:)=0._dp
      do i=1,Nmsrs
        do i2=1,eachsweep
          call cycles(x0,lamb0,phi,4)
          !call cycles(m0,x0,phi,4)
        end do
        call correlation(phi,corr1,corr2)
      end do
      corr1(:)=corr1(:)/real(Nmsrs,dp)
      corr2(:,:)=corr2(:,:)/real(Nmsrs,dp)
      do i=1,N
        CF(i,j)=corr2(i,1)-(corr1(1)*corr1(1))
      end do
    end do
    do j=1,N
      call mean_scalar(CF(j,:),CF_ave(j,k),CF_delta(j,k))
    end do
  end do
  
  do k=1,N+1
    write(60,*) abs(k-1), CF_ave(iv(k),:), CF_delta(iv(k),:)
  end do

  deallocate(corr1,corr2,CF,CF_ave,CF_delta)
  deallocate(phi)
  close(60)
  end subroutine correlate
  
 subroutine vary_m0(mi,mf,lamb0,Nts)
  real(dp), intent(in) :: mi,mf,lamb0
  integer(i4), intent(in) :: Nts
  real(dp), dimension(N,N) :: phi
  integer(i4) :: i,j,k,ie
  real(dp), dimension(Nmsrs2) :: E,M,suscep,heat,U4
  real(dp) :: m0,vol,norm,EE,MM,E_ave,E_delta,M_ave,M_delta,E2,M2,M4
  real(dp) :: suscep_ave,suscep_delta,heat_ave,heat_delta,U4_ave,U4_delta
  !real(dp) :: ARR,ar(Nmsrs2),ar_ave,ar_delta
  !real(dp) :: csx,csx2,cs(Nmsrs2),cs2(Nmsrs2),cs_ave,cs_delta,cs2_ave,cs2_delta
  !open(1, file = 'data/history.dat', status = 'replace')
  open(10, file = 'data/action.dat', status = 'replace')
  open(20, file = 'data/magnetization.dat', status = 'replace')
  open(30, file = 'data/susceptibility.dat', status = 'replace')
  open(40, file = 'data/heat.dat', status = 'replace')
  open(50, file = 'data/binder.dat', status = 'replace')
  !open(60, file = 'data/rank.dat', status = 'replace')
  !open(70, file = 'data/rank2.dat', status = 'replace')
  !open(80, file = 'data/accrate.dat', status = 'replace')
  norm=real(Nmsrs,dp)
  vol=real(N**2,dp)
  do k=1,Nts
  write(*,*) k
  call cold_start(phi)
    m0=mi+(mf-mi)*real(k-1,dp)/real(Nts-1)
    E(:)=0._dp
    M(:)=0._dp
    !cs(:)=0._dp
    !cs2(:)=0._dp
    !ar(:)=0._dp
    do j=1,thermalization
      call cycles(m0,lamb0,phi,4)
      !write(1,*) m0, j, S(m0,phi)/vol
    end do
    do j=1,Nmsrs2
      E2=0._dp
      M2=0._dp
      M4=0._dp
      do i=1,Nmsrs
        do ie=1,eachsweep
          call cycles(m0,lamb0,phi,4)
        end do
        MM=mean(phi)
        EE=S(m0,lamb0,phi)
        E(j)=E(j)+EE
        M(j)=M(j)+abs(MM)
        E2=E2+EE**2
        M2=M2+MM**2
        M4=M4+MM**4
        !cs(j)=cs(j)+csx
        !cs2(j)=cs2(j)+csx2
        !ar(j)=ar(j)+ARR      
      end do
      E(j)=E(j)/norm
      M(j)=M(j)/norm
      E2=E2/norm
      M2=M2/norm
      M4=M4/norm
      suscep(j)=M2-M(j)**2
      heat(j)=E2-E(j)**2
      U4(j)=1._dp-M4/(3._dp*M2**2)
      !cs(j)=cs(j)/norm
      !cs2(j)=cs2(j)/norm
      !ar(j)=ar(j)/norm
    end do
    call mean_scalar(E,E_ave,E_delta)
    call mean_scalar(M,M_ave,M_delta)
    call mean_scalar(suscep,suscep_ave,suscep_delta)
    call mean_scalar(heat,heat_ave,heat_delta)
    call mean_scalar(U4,U4_ave,U4_delta)
    !call mean_scalar(cs,cs_ave,cs_delta)
    !call mean_scalar(cs2,cs2_ave,cs2_delta)
    !call mean_scalar(ar,ar_ave,ar_delta)
    write(10,*) m0, E_ave/vol, E_delta/vol
    write(20,*) m0, M_ave/vol, M_delta/vol
    write(30,*) m0, suscep_ave/vol, suscep_delta/vol
    write(40,*) m0, heat_ave/vol, heat_delta/vol
    write(50,*) m0, U4_ave, U4_delta
    !write(60,*) m0, cs_ave,cs_delta
    !write(70,*) m0, cs2_ave/(vol**2), cs2_delta/(vol**2)
    !write(80,*) m0, ar_ave, ar_delta
  end do
  
  !close(1)
  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  !close(60)
  !close(70)
  !close(80)
  end subroutine vary_m0
  
  subroutine time_test(m0,lamb0,montecarlos)
  integer(i4), intent(in) :: montecarlos
  real(dp), intent(in) :: m0,lamb0
  integer(i4) :: i,j
  real(dp) :: ti,tf,time_ave,time_delta
  real(dp) :: phi(N,N),time(120)
  call cold_start(phi)
  do j=1,120
    call cpu_time(ti)
    do i=1,100000
      call cycles(m0,lamb0,phi,montecarlos)
    end do
    call cpu_time(tf)
    time(j)=tf-ti
  end do
  call mean_scalar(time,time_ave,time_delta)
  write(*,*) time_ave, time_delta
  end subroutine time_test

end program main
