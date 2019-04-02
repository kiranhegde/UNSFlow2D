subroutine solid_flux(ie,c1,c2)
use commons
use pri
use grid
implicit none
integer(kind=i4):: ie,c1,c2
real(kind=dp) :: qc(nvar),p1,p2
real(kind=dp) :: nx, ny,con(nvar),area




!area=0.5d0*fc(ie)%area
!con(:)=pt(cc)%qc(:)
!call con2prim(con)

!nx=fc(ie)%ldx
!ny=fc(ie)%ldy
!if(fc(ie)%ldx==0.d0.and.fc(ie)%ldy==0.d0) then
!nx=fc(ie)%rdx
!ny=fc(ie)%rdy
!endif

!con(:)=pt(cc)%qp(:)!+(pt(cc)%qx(:)*nx+pt(cc)%qy(:)*ny)

p1 = 0.75d0*pt(c1)%qp(4) + 0.25d0*pt(c2)%qp(4)
p2 = 0.25d0*pt(c1)%qp(4) + 0.75d0*pt(c2)%qp(4)
!p1 = 5.d0*pt(c1)%qp(4) +      pt(c2)%qp(4)
!p2 =      pt(c1)%qp(4) + 5.d0*pt(c2)%qp(4)

nx= 0.5d0*(pt(c2)%y-pt(c1)%y)
ny=-0.5d0*(pt(c2)%x-pt(c1)%x)
area=dsqrt(nx*nx+ny*ny)

pt(c1)%la=pt(c1)%la+dsqrt(gamma*pt(c1)%qp(4)/pt(c1)%qp(1))*area
pt(c2)%la=pt(c2)%la+dsqrt(gamma*pt(c2)%qp(4)/pt(c2)%qp(1))*area

pt(c1)%res(3) =pt(c1)%res(3)+p1*nx
pt(c1)%res(4) =pt(c1)%res(4)+p1*ny
pt(c2)%res(3) =pt(c2)%res(3)+p2*nx
pt(c2)%res(4) =pt(c2)%res(4)+p2*ny
!print*,'wall',ie,cc
end
