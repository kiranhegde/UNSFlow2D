subroutine rans 
use commons
use visc
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: ie,c1,c2
!real(kind=dp):: flux



print*,'RANS model not available...'
stop

do ie=1,nof
   c1 = fc(ie)%in
   c2 = fc(ie)%out



!   cell(c1)%res(i)=cell(c1)%res(i)+flux*area
!   cell(c2)%res(i)=cell(c2)%res(i)-flux*area


enddo

end subroutine rans 
