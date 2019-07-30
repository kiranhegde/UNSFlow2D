module data_type
implicit none
!integer, parameter :: dp = kind(1.0d0)
integer, parameter :: i1=selected_int_kind(2)
integer, parameter :: i2=selected_int_kind(4)
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(18)
integer, parameter :: sp=selected_real_kind(6,37)
integer, parameter :: dp=selected_real_kind(15,307)
integer, parameter :: qp=selected_real_kind(31,307)
! Small number
!real(kind=dp) :: EPSILON
!parameter(EPSILON=1.0d-16)
real(kind=dp),parameter::eps=epsilon(1.0_dp)
!integer(kind=i4), parameter :: ndim=2
end module data_type

