program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call cpu_time(starting)

  !Only measure the acc. rate, must modify dphi
  !call acceptance_rate(-1.5_dp)

  !Thermalization history and autocorretion functions
  !call thermalize(-1.4_dp)

  !Histogram
  call make_histogram(1.0_dp)

  !Measure action, magnetization, susceptibility and heat cap.
  !call vary_m0(-1.3_dp,-1.1_dp,50)

  !Measure correlation function
  !call correlate(-2.5_dp)

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

contains

  subroutine thermalize(m0)
  real(dp), intent(in) :: m0
  real(dp), allocatable :: phi(:,:)
  integer(i4) :: i
  real(dp) :: ARR
  open(10, file = 'data/history.dat', status = 'replace')
  allocate(phi(N,N))
  !call hot_start(phi,hotphi)
  call cold_start(phi)
  do i=1,5*thermalization
    call montecarlo(m0,dphi,phi,ARR)
    if(i==1 .or. mod(i,100)==0) then
      write(10,*) i,",", S(m0,phi)/real(N**2,dp)
    end if
  end do
  call autocorrelation(m0,200,phi)
  close(10)
  deallocate(phi)
  end subroutine thermalize

  subroutine acceptance_rate(m0)
  real(dp), allocatable :: phi(:,:)
  real(dp), intent(in) :: m0
  real(dp) :: ARR,AR_ave,AR_delta
  real(dp), dimension(Nmsrs) :: AR
  integer(i4) :: i,k,j
  do j=1,19
    k=0
    allocate(phi(N,N))
    call cold_start(phi)
    !call hot_start(phi,hotphi)
    AR=0._dp
    do i=1,sweeps
      call montecarlo(m0,dphi,phi,ARR)
      if(i>thermalization .and. mod(i,eachsweep)==0) then
        k=k+1
        AR(k)=ARR
      end if
      call flip_sign(phi,i)
    end do
    call mean_scalar(AR,AR_ave,AR_delta)
    write(*,*) AR_ave,",", AR_delta
    deallocate(phi)
    dphi=dphi+0.1_dp/2._dp
    hotphi=2._dp*dphi
  end do
  end subroutine acceptance_rate

  subroutine vary_m0(mi,mf,Nms)
  real(dp), intent(in) :: mi,mf
  integer(i4), intent(in) :: Nms
  real(dp), allocatable :: phi(:,:)
  real(dp) :: M,m0,ARR,M_ave,M_delta!,S_ave,S_delta
  real(dp) :: susc_ave,susc_delta!,heat_ave,heat_delta
  real(dp), dimension(Nmsrs) :: susc1,susc2!,heat1,heat2,Sprom
  integer(i4) :: i,j,k
  open(10, file = 'data/action.dat', status = 'replace')
  open(20, file = 'data/magnetization.dat', status = 'replace')
  open(30, file = 'data/susceptibility.dat', status = 'replace')
  !open(40, file = 'data/heat.dat', status = 'replace')

  do j=1,Nms
    m0=mi+(mf-mi)*real(j-1,dp)/real(Nms-1,dp)
    k=0
    allocate(phi(N,N))
    !call cold_start(phi)
    call hot_start(phi,hotphi)
    call initialize1(susc1,susc2)

    do i=1,sweeps
      call montecarlo(m0,dphi,phi,ARR)
      if(i>thermalization .and. mod(i,eachsweep)==0) then
        k=k+1
        !Sprom(k)=S(m0,phi)
        M=mean(phi)
        call measure2(M,susc1(k),susc2(k))
      end if
      call flip_sign(phi,i)
    end do

    !call mean_scalar(Sprom,S_ave,S_delta)
    call mean_scalar(susc2,M_ave,M_delta)
    !call heat_jackk(heat1,heat2,heat_ave,heat_delta)
    call heat_jackk(susc1,susc2,susc_ave,susc_delta)

    !write(10,*) m0, S_ave/real(N**2,dp), S_delta/real(N**2,dp)
    write(20,*) m0, M_ave/real(N**2,dp), M_delta/real(N**2,dp)
    write(30,*) m0,",", susc_ave/real(N**2,dp),",", susc_delta/real(N**2,dp)
   ! write(40,*) m0, heat_ave/real(N**2,dp), heat_delta/real(N**2,dp)
    deallocate(phi)
  end do
  close(10)
  close(20)
  close(30)
  close(40)
  end subroutine vary_m0

  subroutine make_histogram(m0)
  real(dp), intent(in) :: m0
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
    call montecarlo(m0,dphi,phi,ARR)
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

  subroutine correlate(m0)
  real(dp), intent(in) :: m0
  real(dp), allocatable :: phi(:,:), corr1(:,:), corr2(:,:,:), CF(:,:), CFprom(:,:)
  real(dp) :: ARR
  integer(i4) :: i,k
  open(60, file = 'data/corrfunc.dat', status = 'replace')
  allocate(phi(N,N))
  allocate(corr1(N,Nmsrs))
  allocate(corr2(N,N,Nmsrs))
  allocate(CF(N,N))
  allocate(CFprom(N,N))

  !call cold_start(phi)
  call hot_start(phi,hotphi)
  k=0
  do i=1,sweeps
    call montecarlo(m0,dphi,phi,ARR)
    if(i>thermalization .and. mod(i,eachsweep)==0) then
      k=k+1
      call correlation(phi,k,corr1,corr2)
      call flip_sign(phi,i)
    end if
  end do

  call correlation_function(corr1,corr2,CF,CFprom)
  do i=1,N+1
    write(60,*) abs(i-1), CF(iv(i),1), CFprom(iv(i),1)
  end do

  deallocate(corr1,corr2,CF,CFprom)
  deallocate(phi)
  close(60)
  end subroutine correlate

end program main
