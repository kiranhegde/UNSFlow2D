subroutine slip_wall(ie,cc,qcl,qcr)
use commons
!use pri
use grid
implicit none
integer(kind=i4):: ie,cc
real(kind=dp) :: ds,nx,ny,con(nvar),un
real(kind=dp) :: qcl(nvar), qcr(nvar)


!call con2prim(con)
qcl(:)=cell(cc)%qp(:)

nx=fc(ie)%sx
ny=fc(ie)%sy
ds=dsqrt(nx*nx+ny*ny)
nx=nx/ds
ny=ny/ds

un=qcl(2)*nx+qcl(3)*ny

qcl(2)=qcl(2)-un*nx
qcl(3)=qcl(3)-un*ny

qcr=qcl
!con(:)=cell(cc)%qp(:)!+(cell(cc)%qx(:)*nx+cell(cc)%qy(:)*ny)

!p=con(4)


!cell(cc)%res(3) =cell(cc)%res(3)+p*nx
!cell(cc)%res(4) =cell(cc)%res(4)+p*ny
!print*,'wall',ie,cc
end
