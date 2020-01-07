!---------------------------------------------------------
module output
use data_type
implicit none

real(kind=dp)::Fx=0.0_dp,Fy=0.0_dp,cl=0.0_dp,cd=0.0_dp,cm=0.0_dp
real(kind=dp)::Fx1=0.0_dp,Fy1=0.0_dp


type wallbound
     real(kind=dp)::x,y 
     real(kind=dp)::cp=0.0_dp,cf=0.0_dp
end type wallbound 

type(wallbound),allocatable,dimension(:)::wbc 

end module output

