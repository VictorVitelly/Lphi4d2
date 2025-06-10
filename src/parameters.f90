module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    !integer(i4), parameter :: N=64,thermalization=6000,Nmsrs=5000,eachsweep=500
    integer(i4), parameter :: N=8,thermalization=20000,Nmsrs=1000000,eachsweep=500
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    integer(i4), parameter :: Mbins=10,bins=101

    real(dp), parameter :: lambda0=1._dp, dphi_m=0.5_dp!, dphi=0.5_dp, hotphi=2._dp*dphi
    real(dp) :: dphi=0.5_dp, hotphi=1._dp
    real(dp), parameter :: maxx=1.5_dp, minn=-1.5_dp
    real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
