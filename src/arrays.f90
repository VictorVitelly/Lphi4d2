module arrays
    use iso_fortran_env, only : dp => real64, i4 => int32
    !use parameters, only : N
    implicit none

    real(dp), allocatable :: phi(:,:), CF(:,:), CFprom(:,:)
    real(dp), allocatable :: corr1(:,:)
    real(dp), allocatable :: corr2(:,:,:)

contains

  function int2str(i)
    integer, intent(in) :: i
    character(:), allocatable :: int2str
    character(20) :: k
    write(k,*) i
    int2str = trim(adjustl(k))
  end function int2str

  subroutine write_array(rhop)
    real(dp), dimension(:,:), intent(in) :: rhop
    integer(i4) :: Narr,n2
    Narr=size(rhop,dim=1)
    write(*,*) "_________________________________________________________"
    do n2=1, Narr
      write(*,'("|"'//int2str(10)//'(f20.10,x)"|")') rhop(n2,:)
    end do
    write(*,*) "_________________________________________________________"
  end subroutine write_array

end module arrays
