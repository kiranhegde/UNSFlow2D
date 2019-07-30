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
   call con2prim(cell(i)%qc(:))

   ll  = dabs(u*cell(i)%dx+v*cell(i)%dy)+a*cell(i)%ds
   cell(i)%dt  = cfl*cell(i)%cv/ll          ! local  timestep
   dtglobal = dmin1(dtglobal, cell(i)%dt)   ! global timestep
   if(cell(i)%dt<=0.0) print*,i,cell(i)%dt
enddo



end
!======================================================================================

subroutine spectral 
use param
use pri
use grid
!use commons
implicit none
integer(kind=i4) :: i,j,in,out
real(kind=dp)    :: r0, u0, v0, p0, a0, nx, ny,un,nl,ll
real(kind=dp)    :: con(nvar),prim1(npvar),prim2(npvar)

do i=1,noc
cell(i)%la=0.d0
enddo

do i=1,nof

   in=fc(i)%in
   out=fc(i)%out
   nx  = fc(i)%sx
   ny  = fc(i)%sy
   nl  = dsqrt(nx*nx + ny*ny)
    
   if(in>0.and.out>0) then
      call con2prim(cell(in)%qc(:))
      prim1(:)=prim(:)
      call con2prim(cell(out)%qc(:))
      prim2(:)=prim(:)
      do j=1,nvar
         con(j) = 0.5d0*(prim1(j)+prim2(j))
      enddo

      r0   = con(1)
      u0   = con(2)
      v0   = con(3)
      p0   = con(4)
      a0   = dsqrt(gamma*p0/r0)
      un  = u0*nx + v0*ny
      ll  = dabs(un) + a0*nl
      cell(in)%la=cell(in)%la+ll
      cell(out)%la=cell(out)%la+ll
      fc(i)%la=ll
   elseif(in>0.and.out<0) then
      con(:)=cell(in)%qc(:)
      call con2prim(con)
      un  = u*nx + v*ny
      ll  = dabs(un) + a*nl
      cell(in)%la=cell(in)%la+ll
      fc(i)%la=ll
!   elseif(in==0.and.out/=0) then
!      con(:)=cell(out)%qc(:)
!      call con2prim(con)
!      un  = u*nx + v*ny
!      ll  = dabs(un) + a*nl
!      cell(out)%la=cell(out)%la+ll
!      fc(i)%la=ll
   endif

enddo


end subroutine spectral 
