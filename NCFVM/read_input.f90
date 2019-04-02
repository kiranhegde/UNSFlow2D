!-----------------------------------------------------------------------------
!.....Read parameters from an input file and set freestream values
!-----------------------------------------------------------------------------
subroutine read_input
use param
implicit none
integer(kind=i4) :: inp, iargc, n, inpstatus
character sdummy*32

n = iargc()
if(n .eq. 0)then
   print*,'You must specify an input file.'
   stop
endif

call getarg(1,inpfile)

inpfile=trim(inpfile)

inp=11
open(unit=inp, file=inpfile, status='old')
print*,'Reading parameters from ',inpfile
read(inp,*)sdummy, istart
read(inp,*)sdummy, flow_type
read(inp,*)sdummy, m_inf
read(inp,*)sdummy, aoa_deg
read(inp,*)sdummy, Rey
read(inp,*)sdummy, cfl_max
read(inp,*)sdummy, timemode
read(inp,*)sdummy, gmaxiter, prectype, gerrtol
read(inp,*)sdummy, iterlast
read(inp,*)sdummy, maxiter
read(inp,*)sdummy, minres
read(inp,*)sdummy, saveinterval
read(inp,*)sdummy, scrinterval
read(inp,*)sdummy, niso
read(inp,*)sdummy, flux_type
read(inp,*)sdummy, ILIMIT
read(inp,*)sdummy, vortex, xref, yref
read(inp,*)sdummy, gridfile
close(inp)

inpstatus = yes

if(istart .ne. scratch .and. istart .ne. restart)then
   print*,'Unknown start option',istart
   print*,'Possible values: 1=scratch or 2=restart'
   inpstatus = no
endif

if(flow_type /= 'inviscid' .and. flow_type /= 'laminar' .and. &
   flow_type /= 'rans')then
   print*,'Unknown flow type',flow_type
   print*,'Possible values: inviscid, laminar, rans'
   inpstatus = no
endif

if(timemode /= 'rk3' .and. timemode /= 'lusgs' .and. &
   timemode /= 'gmres')then
   print*,'Unknown time integration scheme ', timemode
   print*,'Possible values: rk3, lusgs, gmres'
   inpstatus = no
endif

if(flux_type /= 'roe' .and. flux_type == 'kfvs')then
   print*,'Unknown flux ', flux_type
   print*,'Possible values: roe, kfvs'
   inpstatus = no
endif


if(ilimit .ne. no .and. ilimit .ne. yes)then
   print*,'Unknown limiter option',ilimit
   print*,'Possible values: 0=no, 1=yes'
   inpstatus = no
endif

if(vortex .ne. yes .and. vortex .ne. no)then
   print*,'Unknown vortex option',vortex
   print*,'Possible values: 0=no, 1=yes'
   inpstatus = no
endif

if(inpstatus .eq. no) stop

end

!========================================================
subroutine read_grid
use param
use grid
implicit none

integer(kind=i4) :: i,j,k
integer(kind=i4) :: n1,n2,m1,m2,c1,c2
integer(kind=i4) :: p1,p2,e1,e2,c,f1
real(kind=dp)    :: x1,y1,x2,y2
real(kind=dp)    :: xc,yc,dx,dy,ds
real(kind=dp)    :: check,xfc,yfc,xc1,xc2,yc1,yc2

open(3,file='geometry.inp')
read(3,*)nop,noc,nof
allocate(pt(nop),fc(nof),cell(noc))
do i=1,nop
 read(3,*)j,pt(i)%x,pt(i)%y,pt(i)%bc
enddo
do i=1,noc
 read(3,*)j,cell(i)%cen(1),cell(i)%cen(2),cell(i)%cv
enddo
do i=1,nof
 read(3,*)j,fc(i)%pt(1),fc(i)%pt(2),fc(i)%in,fc(i)%out,fc(i)%nx,fc(i)%ny,fc(i)%bc
enddo
close(3)

write(*,'(" Bounding box:")')
write(*,'(10x, "xmin =", f18.3)') minval(pt(:)%x) 
write(*,'(10x, "xmax =", f18.3)') maxval(pt(:)%x) 
write(*,'(10x, "ymin =", f18.3)') minval(pt(:)%y)
write(*,'(10x, "ymax =", f18.3)') maxval(pt(:)%y) 


end subroutine read_grid
