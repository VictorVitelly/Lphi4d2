module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    !integer(i4), parameter :: N=64,thermalization=40000,Nmsrs=20000,eachsweep=500
    integer(i4), parameter :: N=16,thermalization=10000,Nmsrs=40000,eachsweep=500,Nmsrs2=1500
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    integer(i4), parameter :: Mbins=10,bins=101,Nauto=15000,Mbin(4)=(/5,10,15,20/)

    real(dp), parameter :: lambda0=1._dp, dphi_m=0.5_dp!, dphi=0.5_dp, hotphi=2._dp*dphi
    real(dp) :: dphi=0.5_dp, hotphi=1._dp
    real(dp), parameter :: maxx=1.5_dp, minn=-1.5_dp
    real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
