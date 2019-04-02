!-----------------------------------------------------------------------------
! Time-step from cfl condition: refined version
!-----------------------------------------------------------------------------
!======================================================================================
subroutine time_step2
use param
use pri
use grid
!use commons
implicit none
integer(kind=i4) :: i
real(kind=dp)    :: qq,aa 

!call spectral

do i=1,nop
pt(i)%dt=0.d0
enddo

dtglobal = 1.0d20
do i=1,nop
   !call con2prim(pt(i)%qc(:))
   !sx = fc(i)%sx
   !sy = fc(i)%sy
   !ds =  dsqrt(sx*sx + sy*sy)
   !lam =  dabs(u*sx1 + v*sy1) + a*ds1
   qq = dsqrt(q)
   !aa = dsqrt(GAMMA*p/rho)
   !pt(i)%dt = CFL*pt(i)%cv/(qq+a)
 
   pt(i)%dt  = cfl*pt(i)%cv/pt(i)%la ! local  timestep
   dtglobal = dmin1(dtglobal, pt(i)%dt)   ! global timestep
   !if(cell(i)%dt<=0.0) print*,i,cell(i)%dt
enddo
!print*,'Global DT:',dtglobal

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
real(kind=dp)    :: ll,sx,sy,ds

!call spectral

pt(:)%dt=0.d0
dtglobal = 1.0d20
do i=1,nop
   call con2prim(pt(i)%qc(:))
   sx=pt(i)%sx
   sy=pt(i)%sy
   ll  = (dabs(u)+a)*sx+(dabs(v)+a)*sy
   pt(i)%dt  = cfl*pt(i)%cv/ll          ! local  timestep
   !pt(i)%la  = ll 
   dtglobal = dmin1(dtglobal, pt(i)%dt)   ! global timestep
   if(pt(i)%dt<=0.0) print*,i,pt(i)%dt
enddo


end
!======================================================================================

subroutine spectral 
use param
use pri
use grid
!use commons
implicit none
integer(kind=i4) :: i,j,p1,p2,c
real(kind=dp)    :: nx,ny,ds,ll,r,dx,dy,vprod
real(kind=dp)    :: prim1(1:nvar),prim2(1:nvar) 
real(kind=dp)    :: un1,un2,a1,a2,a0,un

do i=1,nop
pt(i)%la=0.d0
enddo
fc(:)%la=0.d0

do i=1,nof

   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   nx=fc(i)%nx
   ny=fc(i)%ny
   ds=fc(i)%area
   
   call con2prim(pt(p1)%qc(:))
   un1  = dabs(u*nx + v*ny)
   a1=a

   call con2prim(pt(p2)%qc(:))
   un2  = dabs(u*nx + v*ny)
   a2=a

   a0=0.5d0*(a1+a2)
   un=0.5d0*(un1+un2)

   ll  = (un + a0)*ds
   pt(p1)%la=pt(p1)%la+ll
   pt(p2)%la=pt(p2)%la+ll
   fc(i)%la=ll
enddo



end subroutine spectral 
