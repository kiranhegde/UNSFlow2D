!-----------------------------------------------------------------------------
!.....Read parameters from an input file and set freestream values
!-----------------------------------------------------------------------------
subroutine read_input
use data_type,only:i4
use param,only:  istart,flow_type,iwall,t_wall,m_inf,aoa_deg,&  
               & Rey,timemode,grad_type, & 
               & iterlast,maxiter,minres,saveinterval,&
               & scrinterval,niso,flux_type,ILIMIT, irs,& 
               & cfl_min,cfl_max,cfl_type,CFL_ramp_steps, &
               & vortex, xref, yref,gridfile,scratch,inpfile,restart
use commons,only:yes,no

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
read(inp,*)sdummy, flow_type, iwall, t_wall 
read(inp,*)sdummy, m_inf
read(inp,*)sdummy, aoa_deg
read(inp,*)sdummy, Rey
read(inp,*)sdummy, cfl_type,cfl_min,cfl_max,CFL_ramp_steps
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
read(inp,*)sdummy, irs 
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

if(iwall.ne. no .and. iwall.ne. yes)then
   print*,'Unknown wall conditions :',iwall
   print*,'Possible values: 0=Isothermal, 1=Adiabatic'
   inpstatus = no
endif

if(timemode /= 'rk3' .and. timemode /= 'lusgs' .and. &
   timemode /= 'gmres')then
   print*,'Unknown time integration scheme ', timemode
   print*,'Possible values: rk3, lusgs, gmres'
   inpstatus = no
endif

if(cfl_type /= 'tanh' .and. cfl_type /= 'ramp')then
   print*,'Unknown flow type', cfl_type
   print*,'Possible values: tanh, ramp'
   inpstatus = no
endif

if(grad_type /= 'gg'.and. grad_type /= 'lsqr'.and.grad_type /= 'ggfc')then
   print*,'Unknown  gradient type :', grad_type
   print*,'Possible values: gg , lsqr , ggfc '
   inpstatus = no
endif

if(flux_type /= 'roe' .and. flux_type /= 'ausm'.and.flux_type /= &
    & 'vanleer'.and.flux_type /= 'rusanov')then
   print*,'Unknown flux ', flux_type
   print*,'Possible values: roe, ausm, vanleer, rusanov'
   inpstatus = no
endif

if(irs.ne. no .and. irs.ne. yes)then
   print*,'Unknown Implicit Residual Smoothening option',irs
   print*,'Possible values: 0=no, 1=yes'
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
use output,only:wbc
use grid,only:cells,fc,pt,cell,nop,noc,nof,nogc,nbf,nwbc,startBC,endBC,startfc,endfc, &
                & startGC,endGC
use tools
implicit none

integer(kind=i4) :: i,j,p1,p2,in
integer(kind=i4) :: ghostcell
type (vector) :: a
type (coord)  :: p,px

!type(cells),allocatable,dimension(:)::elm

open(3,file='geometry.inp')
read(3,*)nop,noc,nof,nbf
! start index of internal faces
startFC=1
! end index of internal faces
endFC=nof

! start index of boundary faces
startBC=nof+1
!  end  index of boundary faces
endBC=nof+nbf
! total no. of face
nof=nof+nbf

! start index of Ghost cells
startGC=noc+1
!  end  index of Ghost cells
endGC=noc+nbf

!allocate(pt(nop),fc(nof),cell(noc))
allocate(pt(nop),fc(nof),cell(endGC))

do i=1,nop
 read(3,*)j,pt(i)%x,pt(i)%y,pt(i)%bc
enddo
do i=1,noc
 read(3,*)j,cell(i)%cen(1),cell(i)%cen(2),cell(i)%cv
enddo
do i=startFC,endFC
 read(3,*)j,fc(i)%pt(1),fc(i)%pt(2),fc(i)%in,fc(i)%out,fc(i)%sx,fc(i)%sy,fc(i)%bc
enddo
nwbc=0
do i=startBC,endBC
 read(3,*)j,fc(i)%pt(1),fc(i)%pt(2),fc(i)%in,fc(i)%sx,fc(i)%sy,fc(i)%bc
 if(fc(i)%bc==1001) nwbc=nwbc+1
enddo
close(3)

print*,'Wall boundary faces :',nwbc
allocate(wbc(nwbc))

ghostcell=1
if(ghostcell==1) then
  nogc=nbf
  !allocate(gcell(nogc))
  j=startGC
  do i=startBC,endBC
     in=fc(i)%in
     fc(i)%out=j
     p1=fc(i)%pt(1)
     p2=fc(i)%pt(2)
     a%s%x=pt(p1)%x
     a%s%y=pt(p1)%y
     a%e%x=pt(p2)%x
     a%e%y=pt(p2)%y
     p%x=cell(in)%cen(1)
     p%y=cell(in)%cen(2)
     px%x=0.0_dp
     px%y=0.0_dp
     call pt2linemirror(p,a,px)  
     cell(j)%cen(1)=px%x
     cell(j)%cen(2)=px%y
     write(45,*)px%x,px%y
     j=j+1
  enddo
endif

write(*,'( "Number of Ghost cells =",i10)')nogc
write(*,'(" Bounding box:")')
write(*,'(10x, "xmin =", f18.3)') minval(pt(:)%x) 
write(*,'(10x, "xmax =", f18.3)') maxval(pt(:)%x) 
write(*,'(10x, "ymin =", f18.3)') minval(pt(:)%y)
write(*,'(10x, "ymax =", f18.3)') maxval(pt(:)%y) 


end subroutine read_grid
