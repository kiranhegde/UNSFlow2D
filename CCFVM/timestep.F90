!-----------------------------------------------------------------------------
! Time-step from cfl condition
!-----------------------------------------------------------------------------
subroutine time_step1
use param
use pri
use grid
!use commons
implicit none
integer(kind=i4) :: i

call spectral1

do i=1,noc
cell(i)%dt=0.d0
enddo

dtglobal = 1.0d20
do i=1,noc
   cell(i)%dt  = cfl*cell(i)%cv/cell(i)%la ! local  timestep
   dtglobal = dmin1(dtglobal, cell(i)%dt)   ! global timestep
   !if(cell(i)%dt<=0.0) print*,i,cell(i)%dt
enddo
end subroutine time_step1

!======================================================================================
!  time  step
subroutine time_step02
use param
use grid
!use pri
use commons
implicit none
integer(kind=i4) :: i
real(kind=dp)    :: ll

call spectral02
!call spectral1

cell(:)%dt=0.d0
dtglobal = 1.0d20
do i=1,noc
   !call con2prim(cell(i)%qc(:))
   !ll  = dabs(u*cell(i)%dx+v*cell(i)%dy)+a*cell(i)%ds
   ll  = cell(i)%la
   cell(i)%dt  = cfl*cell(i)%cv/ll          ! local  timestep
   dtglobal = dmin1(dtglobal, cell(i)%dt)   ! global timestep
   !if(cell(i)%dt<=0.0) print*,i,cell(i)%dt
enddo

end
!======================================================================================

subroutine spectral02 
use param
!use pri
use visc
use grid
!use commons
implicit none
integer(kind=i4) :: ie,in,out
real(kind=dp)    :: r0, u0, v0, p0, a0, nx, ny,un,area,ll,vll,mu
real(kind=dp)    :: con(nvar),Kv
!real(kind=dp)    :: factor=4.0_dp/3.0_dp,dist(ndim),cdist,ratio

do ie=1,noc
cell(ie)%la=0.d0
enddo

Kv=0.25_dp
vll=0.0_dp
! if viscous ...
if(flow_type /= 'inviscid') vll = 2.0_dp*gamma/(prandtl*Rey*Kv)
!if(flow_type /= 'inviscid') vll = gamma/(prandtl*Kv)
!if(flow_type /= 'inviscid') vll=2.0_dp*dmax1(factor,gamma/prandtl)!*m_inf/Rey/Kv

do ie=1,nof
   in=fc(ie)%in
   out=fc(ie)%out
   nx  = fc(ie)%nx
   ny  = fc(ie)%ny
   area = fc(ie)%area 
   !dist(:)=cell(out)%cen(:)-cell(in)%cen(:)
   !cdist=dsqrt(dist(1)*dist(1)+dist(2)*dist(2)) 
   mu  = fc(ie)%mu
   con(1:4)=fc(ie)%qp(1:4)
   r0   = con(1)
   u0   = con(2)
   v0   = con(3)
   p0   = con(4)
   a0   = dsqrt(gamma*p0/r0)
   un  = u0*nx + v0*ny
   ! inviscid spectral radius
   ll  = area*(dabs(un) + a0)
   ! adding viscous spectral radius
   ll=ll+vll*mu/r0
   !vll=dmax1(factor,gamma/prandtl)/(Rey*r0)
   !ratio = mu/(prandtl*r0*Rey*kv)
   !vll = ratio*dmax1(factor,gamma)
 
   fc(ie)%la=ll
   cell(in)%la=cell(in)%la+ll
   if(out<=noc) cell(out)%la=cell(out)%la+ll
enddo


end subroutine spectral02 

subroutine spectral1 
use param
use pri
use visc
use grid
!use commons
implicit none
integer(kind=i4) :: i
real(kind=dp)    :: r0, u0, v0, p0, a0, sx, sy,area,ll,vll,mu
real(kind=dp)    :: con(nvar),Kv,factor

do i=1,noc
cell(i)%la=0.d0
enddo

factor=4.0_dp/3.0_dp
Kv=2_dp
vll=0.0_dp

do i=1,noc
   sx=0.5_dp*cell(i)%sx
   sy=0.5_dp*cell(i)%sy
   area=cell(i)%cv
   con(1:nvar)=cell(i)%qp(1:nvar)
   mu=cell(i)%mul
   r0   = con(1)
   u0   = con(2)
   v0   = con(3)
   p0   = con(4)
   a0   = dsqrt(gamma*p0/r0)
   ! inviscid spectral radius
   ll  = dabs(u0)*sx+dabs(v0)*sy+a0*(sx+sy)
   ! adding viscous spectral radius
   if(flow_type /= 'inviscid')   vll = dmax1(factor,gamma)*mu*(sx*sx+sy*sy)/(prandtl*r0*area)
   cell(i)%la=ll+kv*vll
enddo


end subroutine spectral1 
