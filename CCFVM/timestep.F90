!-----------------------------------------------------------------------------
! Time-step from cfl condition
!-----------------------------------------------------------------------------
subroutine time_step2
use param
use pri
use grid
!use commons
implicit none
integer(kind=i4) :: i

call spectral

do i=1,noc
cell(i)%dt=0.d0
enddo

dtglobal = 1.0d20
do i=1,noc
   !call con2prim(cell(in)%qc(:))
   !sx = fc(i)%sx
   !sy = fc(i)%sy
   !ds =  dsqrt(sx*sx + sy*sy)
   !lam =  dabs(u*sx1 + v*sy1) + a*ds1
 
   cell(i)%dt  = cfl*cell(i)%cv/cell(i)%la ! local  timestep
   dtglobal = dmin1(dtglobal, cell(i)%dt)   ! global timestep
   !if(cell(i)%dt<=0.0) print*,i,cell(i)%dt
enddo
end subroutine time_step2

!======================================================================================
!  time  step
subroutine time_step02
use param
use grid
use pri
use commons
implicit none
integer(kind=i4) :: i
real(kind=dp)    :: ll

call spectral

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

subroutine spectral 
use param
use pri
use visc
use grid
!use commons
implicit none
integer(kind=i4) :: ie,in,out
real(kind=dp)    :: r0, u0, v0, p0, a0, nx, ny,un,nl,ll,vll,mu
real(kind=dp)    :: con(nvar),Kv

do ie=1,noc
cell(ie)%la=0.d0
enddo

Kv=0.25_dp
vll=0.0_dp
! if viscous ...
if(flow_type /= 'inviscid') vll = 2.0_dp*gamma/(prandtl*Rey*Kv)

do ie=1,nof
   in=fc(ie)%in
   out=fc(ie)%out
   nx  = fc(ie)%sx
   ny  = fc(ie)%sy
   nl  = dsqrt(nx*nx + ny*ny)
   mu  = fc(ie)%mu
   con(1:4)=fc(ie)%qp(1:4)
   r0   = con(1)
   u0   = con(2)
   v0   = con(3)
   p0   = con(4)
   a0   = dsqrt(gamma*p0/r0)
   un  = u0*nx + v0*ny
   ! inviscid spectral radius
   ll  = dabs(un) + a0*nl
   ! adding viscous spectral radius
   ll = ll + vll*mu/r0
   fc(ie)%la=ll
   cell(in)%la=cell(in)%la+ll
   if(out<=noc) cell(out)%la=cell(out)%la+ll
enddo


end subroutine spectral 
