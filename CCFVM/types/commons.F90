! Common
module commons
use data_type
implicit none
integer(kind=i4), parameter :: no=0
integer(kind=i4), parameter :: yes=1

! viscosity
integer(kind=i4), parameter :: nvar=4

!     GAMMA = ratio of specific heats
!     GAMMA1= GAMMA - 1
!     GAS_CONST = gas constant, this can be set to 1.0
!     M_PI = value of pi
real(kind=dp) :: gamma, gamma1, gas_const, pi
end module commons

