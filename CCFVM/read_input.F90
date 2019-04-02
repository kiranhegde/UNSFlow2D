!-----------------------------------------------------------------------------
!.....Read parameters from an input file and set freestream values
!-----------------------------------------------------------------------------
subroutine read_input
use data_type,only:i4
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
!read(inp,*)sdummy, gmaxiter, prectype, gerrtol
read(inp,*)sdummy, grad_type
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

if(grad_type /= 'gg'.and. grad_type /= 'lsqr'.and.grad_type /= 'ggfc')then
   print*,'Unknown  gradient type :', grad_type
   print*,'Possible values: gg , lsqr , ggfc '
   inpstatus = no
endif

if(flux_type /= 'roe' .and. flux_type /= 'ausm'.and.flux_type /= 'vanleer')then
   print*,'Unknown flux ', flux_type
   print*,'Possible values: roe, ausm , vanleer'
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
use data_type,only:i4
!use param
use grid,only:cells,fc,pt,cell,nop,noc,nof,noc_bc
implicit none

integer(kind=i4) :: i,j
integer(kind=i4) :: ghostcell
type(cells),allocatable,dimension(:)::elm

open(3,file='geometry.inp')
read(3,*)nop,noc,nof
allocate(pt(nop),fc(nof),elm(noc))
do i=1,nop
 read(3,*)j,pt(i)%x,pt(i)%y,pt(i)%bc
enddo
do i=1,noc
 !read(3,*)j,cell(i)%cen(1),cell(i)%cen(2),cell(i)%cv
 read(3,*)j,elm(i)%cen(1),elm(i)%cen(2),elm(i)%cv
enddo
do i=1,nof
 read(3,*)j,fc(i)%pt(1),fc(i)%pt(2),fc(i)%in,fc(i)%out,fc(i)%sx,fc(i)%sy,fc(i)%bc
enddo
close(3)

ghostcell=0
noc_bc=noc
if(ghostcell==1) then
do i=1,nof
   if(fc(i)%in==0) then 
     noc_bc=noc_bc+1
     fc(i)%in=noc_bc
   elseif(fc(i)%out==0) then 
     noc_bc=noc_bc+1
     fc(i)%out=noc_bc
   endif  
enddo
endif

allocate(cell(noc_bc))
do i=1,noc
   cell(i)%cen(1)=elm(i)%cen(1) 
   cell(i)%cen(2)=elm(i)%cen(2) 
   cell(i)%cv    =elm(i)%cv 
enddo

deallocate(elm)


write(*,'( "Number of Ghost cells =",i10)')noc_bc-noc

write(*,'(" Bounding box:")')
write(*,'(10x, "xmin =", f18.3)') minval(pt(:)%x) 
write(*,'(10x, "xmax =", f18.3)') maxval(pt(:)%x) 
write(*,'(10x, "ymin =", f18.3)') minval(pt(:)%y)
write(*,'(10x, "ymax =", f18.3)') maxval(pt(:)%y) 


end subroutine read_grid
