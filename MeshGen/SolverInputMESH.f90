module  data_type
implicit none
integer, parameter :: i1=selected_int_kind(2)
integer, parameter :: i2=selected_int_kind(4)
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(18)
integer, parameter :: sp=selected_real_kind(6,37)
integer, parameter :: dp=selected_real_kind(15,307)
integer, parameter :: qp=selected_real_kind(31,307)

end module data_type

module grid
use data_type
implicit none

integer(kind=i4):: np,nf,nc 

type points 
     real(kind=dp) :: x,y,z
     integer(kind=i4):: bc,flag,nv2c 
     integer(kind=i4),dimension(:),pointer::v2c
end type points 

type faces  
     integer(kind=i4):: pt(2) 
     integer(kind=i4):: in 
     integer(kind=i4):: out 
     integer(kind=i4):: bc
     integer(kind=i4):: flag
     real(kind=dp)   :: sx,sy 
end type faces  

type cells  
     integer(kind=i4):: nc2v,nc2f,nc2c  
     real(kind=dp)   :: xc,yc,zc,cv
     integer(kind=i4),dimension(:),pointer::c2v
     integer(kind=i4),dimension(:),pointer::c2f
     integer(kind=i4),dimension(:),pointer::c2c
end type cells  

type(points),allocatable,dimension(:)::pt
type(points),allocatable,dimension(:,:)::mesh
type(faces),allocatable,dimension(:)::fc
type(cells),allocatable,dimension(:)::cell

end module grid
!========================================================
program SolverMesh
use grid
implicit none

integer(kind=i4) :: i,j,k,nx,ny,nz,nbk,cc,c1
real(kind=dp) :: dummy

open(3,file='grid.grd')
read(3,*)nbk
if(nbk>1) stop 'Only Single block is allowed'
read(3,*) nx, ny, nz
print*
print*,'nx,ny,nz :',nx, ny, nz
print*
allocate(mesh(nx,ny))
read(3,*) (((mesh(i,j)%x, i=1,nx), j=1,ny), k=1,nz), &
          (((mesh(i,j)%y, i=1,nx), j=1,ny), k=1,nz), &
          (((dummy      , i=1,nx), j=1,ny), k=1,nz)
close(3)

!For O-grid
!nx=nx-1

np=(nx-1)*ny
nf=(nx-1)*(2*ny-1)
nc=(nx-1)*(ny-1)

print*,'No. of vertices =',np
print*,'No. of faces    =',nf
print*,'No. of cells    =',nc

allocate(pt(np),fc(nf),cell(nc))
cc=0
do j=1,ny
   do i=1,nx-1
      cc=cc+1   
      pt(cc)%x=mesh(i,j)%x      
      pt(cc)%y=mesh(i,j)%y      
      pt(cc)%bc=0
      if(j==1) pt(cc)%bc=1001
      if(j==ny) pt(cc)%bc=2001
      pt(cc)%flag=0
   enddo
enddo
print*
print*,'Node count =',cc

cc=0
do j=1,ny
   do i=1,nx-1
      cc=cc+1   
      fc(cc)%pt(1)=cc
      fc(cc)%pt(2)=cc+1
      if(i==nx-1) then 
        fc(cc)%pt(2)=cc-(nx-2)
      endif 
      fc(cc)%bc=0
      fc(cc)%out=cc-nx+1
      if(j==1) then 
        fc(cc)%bc=1001
        fc(cc)%out=0
      endif
      fc(cc)%in=cc
      if(j==ny) then 
        fc(cc)%bc=2001
        fc(cc)%in=0
      endif
      fc(cc)%flag=0
   enddo
enddo

print*,'Zhi-face count=',cc

c1=0
do j=1,ny-1
   do i=1,nx-1
      cc=cc+1   
      c1=c1+1   
      fc(cc)%pt(1)=c1
      fc(cc)%pt(2)=c1+nx-1
      fc(cc)%bc=0
      fc(cc)%out=c1
      fc(cc)%in=c1-1
      if(i==1) then 
!        fc(cc)%bc=2001
        fc(cc)%in=(nx-1)*j
      endif 
      fc(cc)%flag=0
   enddo
enddo

print*
print*,'face count =',cc

call solver_compatible
!call connectivity

end program SolverMesh


!=====================================================
subroutine solver_compatible
!=====================================================
use grid
implicit none
integer(kind=i4) :: i,j,in,out,count(nc),c1,c2
integer(kind=i4) :: n1,n2
real(kind=dp)    :: dx,dy,ds,x1,y1,x2,y2,nx,ny
real(kind=dp)    :: rx,ry,dotprod

do i=1,nc
   count(i) = 0
   cell(i)%xc = 0.0d0
   cell(i)%yc = 0.0d0
   cell(i)%cv = 0.0d0
enddo

do j=1,nf
   n1  = fc(j)%pt(1)
   n2  = fc(j)%pt(2)
   in  = fc(j)%in
   out = fc(j)%out
   !print*, j,in,out
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y

   if(in/=0) then
    cell(in)%xc = cell(in)%xc+x1
    cell(in)%yc = cell(in)%yc+y1
    count(in) = count(in) + 1
    cell(in)%xc = cell(in)%xc+x2
    cell(in)%yc = cell(in)%yc+y2
    count(in) = count(in) + 1
   endif
   if(out/=0) then
    cell(out)%xc = cell(out)%xc+x1
    cell(out)%yc = cell(out)%yc+y1
    count(out) = count(out) + 1
    cell(out)%xc = cell(out)%xc+x2
    cell(out)%yc = cell(out)%yc+y2
    count(out) = count(out) + 1
   endif
enddo

do i = 1,nc
   cell(i)%xc=cell(i)%xc/count(i)
   cell(i)%yc=cell(i)%yc/count(i)
enddo

!  faces forming  a cell
cell(:)%nc2f = 0
do i=1,nf
   in = fc(i)%in
   out = fc(i)%out
   if(in/=0) then
     cell(in)%nc2f = cell(in)%nc2f + 1
   endif
   if(out/=0) then
     cell(out)%nc2f = cell(out)%nc2f + 1
   endif
enddo

do i=1,nc
   allocate(cell(i)%c2f(cell(i)%nc2f))
enddo

cell(:)%nc2f = 0
do i=1,nf
   in = fc(i)%in
   out = fc(i)%out
   if(in/=0) then
     cell(in)%nc2f = cell(in)%nc2f + 1
     cell(in)%c2f(cell(in)%nc2f) = i
   endif
   if(out/=0) then
     cell(out)%nc2f = cell(out)%nc2f + 1
     cell(out)%c2f(cell(out)%nc2f) = i
   endif
enddo



do i=1,nf
   n1 = fc(i)%pt(1)
   n2 = fc(i)%pt(2)
   c1 = fc(i)%in
   c2 = fc(i)%out
 
   if(n1>np.or.n2>np.or.n1<1.or.n2<1) then
      print*,'fc,n1,n2,in,out:',i,n1,n2,c1,c2
      stop  
   endif    

   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds
   x1 = (x2+x1)/2.0d0 ; y1= (y2+y1)/2.0d0
   if( c1/=0 ) then
      x2 = cell(c1)%xc ; y2 = cell(c1)%yc
      dx = x2-x1 ; dy = y2-y1
      ds = dsqrt(dx*dx+dy*dy)
      rx = dx/ds
      ry = dy/ds
      dotprod = rx*nx+ry*ny
      if(dotprod<0.0d0) then
         fc(i)%in = c1
         fc(i)%out = c2
      else
         fc(i)%in = c2
         fc(i)%out = c1
      endif
   elseif( c1==0 ) then
      x2 = cell(c2)%xc ; y2 = cell(c2)%yc
      dx = x2-x1 ; dy = y2-y1
      ds = dsqrt(dx*dx+dy*dy)
      rx = dx/ds
      ry = dy/ds
      dotprod = rx*nx+ry*ny
      if(dotprod<0.0d0) then
       fc(i)%in = c2
       fc(i)%out = c1
      else
       fc(i)%in = c1
       fc(i)%out = c2
      endif

   endif
enddo


do i=1,nf
   n1 = fc(i)%pt(1)
   n2 = fc(i)%pt(2)
   in = fc(i)%in
   out = fc(i)%out
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   dx = 0.5d0*(x2+x1) ; dy = y2-y1

   if(in/=0) then
    cell(in)%cv = cell(in)%cv + dx*dy 
   endif
   if(out/=0) then
    cell(out)%cv = cell(out)%cv - dx*dy 
   endif

   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds
   !fc(i)%sx = nx
   !fc(i)%sy = ny
   fc(i)%sx = dy
   fc(i)%sy = -dx
enddo

ds=0.0
do i=1,nc
   ds=ds+cell(i)%cv
enddo
print*,'Domain volume = ',ds

open(3,file='Vec_plot.dat')
open(4,file='Vec_InCell.dat')
open(5,file='Vec_OutCell.dat')
do i=1,nf
   n1=fc(i)%pt(1)
   n2=fc(i)%pt(2)
   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   write(3,100)x1,y1,x2-x1,y2-y1
   !write(21,*)pt(n2)%x,pt(n2)%y
   !write(21,*)
   x1 = (x2+x1)/2.0d0 ; y1= (y2+y1)/2.0d0
   in=fc(i)%in
   out=fc(i)%out
   if(in>0) write(4,100)x1,y1,cell(in)%xc-x1,cell(in)%yc-y1
   if(out>0) write(5,100)x1,y1,cell(out)%xc-x1,cell(out)%yc-y1
enddo
close(3)
close(4)
close(5)

100 format(1x,4(f15.6,1x))

open(3,file='CC_plot.dat')
do i=1,nc
   write(3,*)cell(i)%xc,cell(i)%yc
   write(3,*)
enddo
close(3)


open(13,file='geometry.inp')
write(13,*)np,nc,nf
do i=1,np
 write(13,*)i,pt(i)%x,pt(i)%y,pt(i)%bc
enddo
do i=1,nc
 write(13,*)i,cell(i)%xc,cell(i)%yc,cell(i)%cv
enddo
do i=1,nf
 write(13,*)i,fc(i)%pt(1),fc(i)%pt(2),fc(i)%in,fc(i)%out,&
                       fc(i)%sx,fc(i)%sy,fc(i)%bc
enddo
close(13)

end subroutine solver_compatible
