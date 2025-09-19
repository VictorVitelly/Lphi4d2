module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    !integer(i4), parameter :: N=8,thermalization=1300,eachsweep=80
    !integer(i4), parameter :: N=8,thermalization=1300,eachsweep=27
    integer(i4), parameter :: N=16,thermalization=1800,eachsweep=36
    !integer(i4), parameter :: N=32,thermalization=2300,eachsweep=45
    !integer(i4), parameter :: N=64,thermalization=2700,eachsweep=55
    !integer(i4), parameter :: N=128,thermalization=3000,eachsweep=60
    integer(i4), parameter :: Nmsrs=1, Nmsrs2=120000
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    integer(i4), parameter :: bins=101,Nauto=150000,Mbin(5)=(/4,5,10,15,20/)
    !integer(i4),parameter :: Mbin(3)=(/5,10,15/)!

    real(dp), parameter :: lambda0=1.0_dp, dphi_m=0.5_dp!, dphi=0.5_dp, hotphi=2._dp*dphi
    real(dp) :: dphi=0.5_dp, hotphi=1._dp
    real(dp), parameter :: maxx=1.5_dp, minn=-1.5_dp
    real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
