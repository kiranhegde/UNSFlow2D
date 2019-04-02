subroutine geometric  
use grid
implicit none
!integer(kind=i4) :: i

!    Read grid from file
call read_grid

! Mesh connectivity generation

call cell_2_cell

call renumber

call cell_2_face

!call face_co_volume
!call cell_face_norm_c2f_dist

call vertex_2_cell
call cell_2_vertex
call avg_cell_face_length

call face_co_volume
call cell_face_norm_c2f_dist

!    writes cell face  normal vector's in gnuplot format to files
call write_2_file
call LSQR_Coeff

!========================================================
contains
!========================================================

subroutine cell_2_cell
use param
use grid
implicit none

integer(kind=i4) :: i,j,c
integer(kind=i4) :: in,out

! cells surrounding a cell
print*
print*,"==> Cells surrounding a cell "

cell(:)%nc2c=0
do i=1,nof
   in=fc(i)%in
   out=fc(i)%out
 
   if(in/=0.and.out/=0) then

     c=0
     do j=1,cell(in)%nc2c
        if(cell(in)%c2c(j)==out) c=c+1  
     enddo

     if(c==0) then 
       cell(in)%nc2c=cell(in)%nc2c+1
       call alloc_int_ptr(cell(in)%c2c,cell(in)%nc2c)      
       cell(in)%c2c(cell(in)%nc2c)=out      
     endif


     c=0
     do j=1,cell(out)%nc2c
        if(cell(out)%c2c(j)==in) c=c+1  
     enddo

     if(c==0) then 
       cell(out)%nc2c=cell(out)%nc2c+1
       call alloc_int_ptr(cell(out)%c2c,cell(out)%nc2c)      
       cell(out)%c2c(cell(out)%nc2c)=in      
     endif


   endif
   
enddo

print*, 'Max. cell nbr  =',maxval(cell(1:noc)%nc2c),'at node',maxloc(cell(1:noc)%nc2c)
print*, 'Min. cell nbr  =',minval(cell(1:noc)%nc2c),'at node',minloc(cell(1:noc)%nc2c)

!print*,'cell'
!do i=1,noc
!   print*,i,(cell(i)%c2c(j), j=1,cell(i)%nc2c) 
!enddo
!stop

end  subroutine cell_2_cell
!===========================================================
!Cells to face   connectivity 
subroutine cell_2_face
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: in,out

print*
print*,"==> Faces surrounding a cell "

cell(:)%nc2f=0
do i=1,nof
   in=fc(i)%in
   out=fc(i)%out

   if(in/=0) then
   cell(in)%nc2f=cell(in)%nc2f+1
   call alloc_int_ptr(cell(in)%c2f,cell(in)%nc2f)      
   cell(in)%c2f(cell(in)%nc2f)=i      
   endif

   if(out/=0) then
   cell(out)%nc2f=cell(out)%nc2f+1
   call alloc_int_ptr(cell(out)%c2f,cell(out)%nc2f)      
   cell(out)%c2f(cell(out)%nc2f)=i      
   endif
enddo

print*, 'Max. cell faces =',maxval(cell(1:noc)%nc2f),'at node',maxloc(cell(1:noc)%nc2f)
print*, 'Min. cell faces =',minval(cell(1:noc)%nc2f),'at node',minloc(cell(1:noc)%nc2f)


!print*,'cell'
!do i=1,noc
!   write(80,*),i,(cell(i)%c2f(j), j=1,cell(i)%nc2f) 
!enddo
!stop

end  subroutine cell_2_face
!===========================================================
! Finding face co-volume for face & Cell gradient
subroutine face_co_volume
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: p1,p2
integer(kind=i4) :: in,out
real(kind=dp)    :: x1,y1,x2,y2
real(kind=dp)    :: dx,dy


do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   in=fc(i)%in
   out=fc(i)%out
   fc(i)%cov=0.d0
   
   if(in/=0.and.out/=0) then
      x1 = pt(p1)%x    ; y1 = pt(p1)%y
      x2 = cell(out)%cen(1) ; y2 = cell(out)%cen(2)
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   
      x1 = cell(out)%cen(1) ; y1 = cell(out)%cen(2)
      x2 = pt(p2)%x    ; y2 = pt(p2)%y
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   
      x1 = pt(p2)%x    ; y1 = pt(p2)%y
      x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   
      x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
      x2 = pt(p1)%x    ; y2 = pt(p1)%y
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   endif

   if(in==0.and.out/=0) then
      x1 = pt(p1)%x    ; y1 = pt(p1)%y
      x2 = cell(out)%cen(1) ; y2 = cell(out)%cen(2)
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   
      x1 = cell(out)%cen(1) ; y1 = cell(out)%cen(2)
      x2 = pt(p2)%x    ; y2 = pt(p2)%y
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   
      x1 = pt(p2)%x ; y1 = pt(p2)%y
      x2 = pt(p1)%x ; y2 = pt(p1)%y
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   endif

   if(out==0.and.in/=0) then
      x1 = pt(p1)%x ; y1 = pt(p1)%y
      x2 = pt(p2)%x ; y2 = pt(p2)%y
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy

      x1 = pt(p2)%x    ; y1 = pt(p2)%y
      x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   
      x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
      x2 = pt(p1)%x    ; y2 = pt(p1)%y
      dx = 0.5d0*(x2+x1) ; dy = y2-y1
      fc(i)%cov=fc(i)%cov+dx*dy
   endif
enddo

do i=1,noc
   cell(i)%cov=0.d0
enddo

do i=1,nof
   in=fc(i)%in
   out=fc(i)%out
   if(fc(i)%cov <= 0.d0) then  
      print*,i,'-ve face co-volume'
      write(75,*)i,in,out
      stop
   endif

   if(in/=0) cell(in)%cov=cell(in)%cov+fc(i)%cov
   if(out/=0) cell(out)%cov=cell(out)%cov+fc(i)%cov
enddo

do i=1,noc
   if(cell(i)%cov<=0.d0) then
      print*,i,'-ve cell co-volume'
      write(75,*)i,in,out
      stop
   endif
enddo

end subroutine face_co_volume
!===========================================================
! Finding face co-volume for face & Cell gradient
subroutine cell_face_norm_c2f_dist
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: p1,p2,in,out
real(kind=dp)    :: x1,y1,x2,y2
real(kind=dp)    :: xc,yc,dx,dy
real(kind=dp)    :: check,dotprod

! face normal components
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   dx = x2-x1 ; dy = y2-y1
   fc(i)%sx = dy
   fc(i)%sy = -dx
enddo


! making boundary face normal outward
do i=1,nof
   if(fc(i)%bc==2001) then
   in=fc(i)%in
   out=fc(i)%out
   if(in==0) then
      xc=cell(out)%cen(1)
      yc=cell(out)%cen(2)
      p1=fc(i)%pt(1)
      p2=fc(i)%pt(2)
      check=dotprod(p1,p2,xc,yc)
      if(check>0.0d0) then
        fc(i)%pt(1)=p2
        fc(i)%pt(2)=p1
        x1 = pt(p2)%x ; y1 = pt(p2)%y
        x2 = pt(p1)%x ; y2 = pt(p1)%y
        dx = x2-x1 ; dy = y2-y1
        fc(i)%sx = dy
        fc(i)%sy = -dx
      endif  
   endif  

   if(out==0) then
      xc=cell(in)%cen(1)
      yc=cell(in)%cen(2)
      p1=fc(i)%pt(1)
      p2=fc(i)%pt(2)
      check=dotprod(p1,p2,xc,yc)
      if(check>0.0d0) then
        fc(i)%pt(1)=p2
        fc(i)%pt(2)=p1
        x1 = pt(p2)%x ; y1 = pt(p2)%y
        x2 = pt(p1)%x ; y2 = pt(p1)%y
        dx = x2-x1 ; dy = y2-y1
        fc(i)%sx = dy
        fc(i)%sy = -dx
      endif  
   endif  
   endif  
enddo

! cell center to face distance
do i=1,nof
   !in=fc(i)%in
   !out=fc(i)%out
   
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   
   fc(i)%cen(1)=0.5d0*(pt(p1)%x+pt(p2)%x)
   fc(i)%cen(2)=0.5d0*(pt(p1)%y+pt(p2)%y)
   
!   ! Left cell distance 
!   if(in/=0) then
!   xc1=cell(in)%cen(1)
!   yc1=cell(in)%cen(2)
!   fc(i)%ldx = xfc-xc1 
!   fc(i)%ldy = yfc-yc1 
!   endif
!
!   ! Right cell distance 
!   if(out/=0) then
!   xc2=cell(out)%cen(1)
!   yc2=cell(out)%cen(2)
!   fc(i)%rdx = xfc-xc2 
!   fc(i)%rdy = yfc-yc2 
!   endif

enddo

!do i=1,noc
!   print*,i,(cell(i)%c2c(j), j=1,cell(i)%nc2c) 
!enddo

!print*, 'Max. cell nbr=',maxval(cell(:)%nc2c),'at node',maxloc(cell(:)%nc2c)
!print*, 'Min. cell nbr=',minval(cell(:)%nc2c),'at node',minloc(cell(:)%nc2c)
end subroutine cell_face_norm_c2f_dist

!==============================================================================
subroutine vertex_2_cell
use param
use grid
implicit none

integer(kind=i4) :: i,j
integer(kind=i4) :: p1,p2,in,out,c
real(kind=dp)    :: xc,yc,dx,dy



! cells surrounding a node
print*
print*,"==> Cells surrounding a node"

pt(:)%nv2c=0
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   in=fc(i)%in
   out=fc(i)%out

   if(in/=0) then
         c=0
         do j=1,pt(p1)%nv2c
            if(pt(p1)%v2c(j)==in) c=c+1  
         enddo
         if(c==0) then 
           pt(p1)%nv2c=pt(p1)%nv2c+1
           call alloc_int_ptr(pt(p1)%v2c,pt(p1)%nv2c)      
           pt(p1)%v2c(pt(p1)%nv2c)=in      
         endif
      
         c=0
         do j=1,pt(p2)%nv2c
            if(pt(p2)%v2c(j)==in) c=c+1  
         enddo
         if(c==0) then 
           pt(p2)%nv2c=pt(p2)%nv2c+1
           call alloc_int_ptr(pt(p2)%v2c,pt(p2)%nv2c)      
           pt(p2)%v2c(pt(p2)%nv2c)=in      
         endif
   endif

   if(out/=0) then
         c=0
         do j=1,pt(p1)%nv2c
            if(pt(p1)%v2c(j)==out) c=c+1  
         enddo
         if(c==0) then 
           pt(p1)%nv2c=pt(p1)%nv2c+1
           call alloc_int_ptr(pt(p1)%v2c,pt(p1)%nv2c)      
           pt(p1)%v2c(pt(p1)%nv2c)=out      
         endif
      
         c=0
         do j=1,pt(p2)%nv2c
            if(pt(p2)%v2c(j)==out) c=c+1  
         enddo
         if(c==0) then 
           pt(p2)%nv2c=pt(p2)%nv2c+1
           call alloc_int_ptr(pt(p2)%v2c,pt(p2)%nv2c)      
           pt(p2)%v2c(pt(p2)%nv2c)=out      
         endif
   endif

enddo   

!print*,'kkr'

do i=1,nop
   allocate(pt(i)%wt(pt(i)%nv2c))
   do j=1,pt(i)%nv2c
      c=pt(i)%v2c(j)
      xc=cell(c)%cen(1) 
      yc=cell(c)%cen(2) 
      dx=xc-pt(i)%x
      dy=yc-pt(i)%y
      pt(i)%wt(j)=1.d0/dsqrt(dx*dx+dy*dy)
      if(pt(i)%wt(j)<0.d0) then
        print*,"'Distance negative..."
        stop
      endif    
   enddo 
   !print*,i,(pt(i)%v2c(j), j=1,pt(i)%nv2c) 
enddo



print*, 'Max. nbr=',maxval(pt(1:nop)%nv2c),'at node',maxloc(pt(1:nop)%nv2c)
print*, 'Min. nbr=',minval(pt(1:nop)%nv2c),'at node',minloc(pt(1:nop)%nv2c)
end subroutine vertex_2_cell
!====================================================================
subroutine avg_cell_face_length
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: in,out
real(kind=dp)    :: ds,sx,sy


print*
print*,"==> Cell face length  "

! average cell face length
!do i=1,noc
!   cell(i)%ds=0.d0
!   cell(i)%dx=0.d0
!   cell(i)%dy=0.d0
!   do j=1,cell(i)%nc2f
!      c=cell(i)%c2f(j) 
!      dx=dabs(fc(c)%sx)
!      dy=dabs(fc(c)%sy)
!      ds=dsqrt(dx*dx+dy*dy)
!      cell(i)%ds=cell(i)%ds+ds
!      cell(i)%dx=cell(i)%dx+dx
!      cell(i)%dy=cell(i)%dy+dy
!   enddo 
!   cell(i)%ds=cell(i)%ds/cell(i)%nc2f
!   cell(i)%dx=0.5d0*cell(i)%dx
!   cell(i)%dy=0.5d0*cell(i)%dy
!   if(cell(i)%dx<=0.d0.or.cell(i)%dy<=0.d0) then 
!     print*, 'Time step length scale  zero'
!     stop 
!   endif
!enddo

cell(:)%ds=0.d0
cell(:)%dx=0.d0
cell(:)%dy=0.d0
do i=1,nof
   in=fc(i)%in
   out=fc(i)%out
   sx=0.5d0*fc(i)%sx
   sy=0.5d0*fc(i)%sy
   ds=dsqrt(sx*sx+sy*sy)
   if(in/=0) then 
     cell(in)%ds=cell(in)%ds+ds
     cell(in)%dx=cell(in)%dx+sx
     cell(in)%dy=cell(in)%dy+sy
   endif 

   if(out/=0) then 
     cell(out)%ds=cell(out)%ds+ds
     cell(out)%dx=cell(out)%dx+sx
     cell(out)%dy=cell(out)%dy+sy
   endif 
enddo


!do i=1,noc
!   !cell(i)%ds=cell(i)%ds/cell(i)%nc2f
!   cell(i)%ds=0.5d0*cell(i)%ds
!   cell(i)%dx=0.5d0*cell(i)%dx
!   cell(i)%dy=0.5d0*cell(i)%dy
   !cell(i)%dx=cell(i)%dx/cell(i)%nc2f
   !cell(i)%dy=cell(i)%dy/cell(i)%nc2f
!enddo
print*, 'Max. cell face length =',maxval(cell(1:noc)%ds),'at node',maxloc(cell(1:noc)%ds)
print*, 'Min. cell face length =',minval(cell(1:noc)%ds),'at node',minloc(cell(1:noc)%ds)

end subroutine avg_cell_face_length
!====================================================================
subroutine cell_2_vertex
use param
use grid
implicit none

integer(kind=i4) :: i,j,k,m1,m2
integer(kind=i4) :: p1,p2,f1
real(kind=dp)    :: xc,yc
real(kind=dp)    :: check,dotprod


!Cells to vertex connectivity 
print*
print*,"==> Vertex surrounding a cell "

!cell(:)%nc2v=0
!do i=1,nop
!   do j=1,pt(i)%nv2c
!      in=pt(i)%v2c(j)
!      cell(in)%nc2v=cell(in)%nc2v+1
!      call alloc_int_ptr(cell(in)%c2v,cell(in)%nc2v)      
!      cell(in)%c2v(cell(in)%nc2v)=i      
!   enddo    
!enddo

cell(:)%nc2v=0
do i=1,noc
   xc=cell(i)%cen(1)
   yc=cell(i)%cen(2)
      p1=fc(cell(i)%c2f(1))%pt(1)
      p2=fc(cell(i)%c2f(1))%pt(2)
      check=dotprod(p1,p2,xc,yc)
      !print*,'check',check

      if(check<0.0d0) then
         cell(i)%nc2v=cell(i)%nc2v+1
         call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)      
         cell(i)%c2v(cell(i)%nc2v)=p1      

         cell(i)%nc2v=cell(i)%nc2v+1
         call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)      
         cell(i)%c2v(cell(i)%nc2v)=p2      
      else
         cell(i)%nc2v=cell(i)%nc2v+1
         call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)      
         cell(i)%c2v(cell(i)%nc2v)=p2      

         cell(i)%nc2v=cell(i)%nc2v+1
         call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)      
         cell(i)%c2v(cell(i)%nc2v)=p1      
      endif
      
   do k=2,cell(i)%nc2f-1
      p2=cell(i)%c2v(cell(i)%nc2v)
      !print*,(cell(i)%c2v(j), j=1,cell(i)%nc2v) 

      do j=k,cell(i)%nc2f
         m1=fc(cell(i)%c2f(j))%pt(1)
         m2=fc(cell(i)%c2f(j))%pt(2)
         if(m1==p2.or.m2==p2) then 
            check=dotprod(m1,m2,xc,yc)
            !print*,i,k,j,c,check
            !print*,p1,p2,m1,m2 
            if(check<0.0d0) then
               cell(i)%nc2v=cell(i)%nc2v+1
               call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)      
               cell(i)%c2v(cell(i)%nc2v)=m2      
            else
               cell(i)%nc2v=cell(i)%nc2v+1
               call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)      
               cell(i)%c2v(cell(i)%nc2v)=m1      
            endif
              
            f1=cell(i)%c2f(j)      
            cell(i)%c2f(j)=cell(i)%c2f(k)      
            cell(i)%c2f(k)=f1      
            !p2=cell(i)%c2v(c)
            cycle 
           ! exit
         endif
      enddo 
   enddo 
   cell(i)%nc2v=cell(i)%nc2v+1
   call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)
   cell(i)%c2v(cell(i)%nc2v)= cell(i)%c2v(1)

   !print*,i,(cell(i)%c2v(j), j=1,cell(i)%nc2v) 
   !stop
enddo

print*, 'Max. cell vertex=',maxval(cell(1:noc)%nc2v)-1,'at node',maxloc(cell(1:noc)%nc2v)
print*, 'Min. cell vertex=',minval(cell(1:noc)%nc2v)-1,'at node',minloc(cell(1:noc)%nc2v)

end subroutine cell_2_vertex
!==========================================================
subroutine write_2_file
use param
use grid
implicit none

integer(kind=i4) :: i,j
integer(kind=i4) :: p1,p2
real(kind=dp)    :: x1,y1,x2,y2

x1=0.0_dp
x2=0.0_dp
y1=0.0_dp
y2=0.0_dp

do i=1,noc
!   print*,i
   do j=1,cell(i)%nc2v-1
      p1=cell(i)%c2v(j)
      p2=cell(i)%c2v(j+1)
      x1 = pt(p1)%x ; y1 = pt(p1)%y
      x2 = pt(p2)%x ; y2 = pt(p2)%y
      !write(25,101)x1,y1,x2-x1,y2-y1
      !write(25,*)pt(p1)%x,pt(p1)%y
      !write(25,*)pt(p2)%x,pt(p2)%y
      !write(25,*)
   enddo
   !print*,i,(cell(i)%c2v(j), j=1,cell(i)%nc2v) 
   !exit
enddo

open(3,file='Vec_plot.dat')
open(4,file='Vec_InCell.dat')
open(5,file='Vec_OutCell.dat')
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   write(3,100)x1,y1,x2-x1,y2-y1
   !write(21,*)pt(p2)%x,pt(p2)%y
   !write(21,*)
   x1 = (x2+x1)/2.0d0 ; y1= (y2+y1)/2.0d0
   !in=fc(i)%in
   !out=fc(i)%out
   !if(fc(i)%bc==1001.and.in>0) write(4,100)x1,y1,cell(in)%cen(1)-x1,cell(in)%cen(2)-y1
   if(fc(i)%bc==1001) write(4,100)x1,y1,fc(i)%sx,fc(i)%sy
   !if(fc(i)%bc==1001.and.out>0) write(4,100)x1,y1,cell(out)%cen(1)-x1,cell(out)%cen(2)-y1
   if(fc(i)%bc==2001) write(5,100)x1,y1,fc(i)%sx,fc(i)%sy
enddo
close(3)
close(4)
close(5)

do i=1,nof
   if(fc(i)%bc==1001) then
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   !write(36,*)x1,y1
   !write(36,*)x2,y2
   !write(36,*)
   endif
   if(fc(i)%bc==2001) then
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   !write(37,*)x1,y1
   !write(37,*)x2,y2
   !write(37,*)
   endif
enddo

open(3,file="Normals.dat")  
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   x1 = (x2+x1)/2.0d0 ; y1= (y2+y1)/2.0d0
   write(3,100)x1,y1,fc(i)%sx,fc(i)%sy
enddo
close(3)

100 format(1x,4(f15.8,1x))

open(3,file='CC_plot.dat')
do i=1,noc
   write(3,*)cell(i)%cen(1),cell(i)%cen(2)
   write(3,*)
enddo
close(3)
!101 format(1x,4(f15.6,1x))
end subroutine write_2_file
!=================================================================================

!===========================================================================
!                      Allocate/extend array
!===========================================================================
subroutine alloc_int_ptr1(x,n)
use data_type
implicit none
integer(kind=i4),intent(in) ::n 
integer(kind=i4)::i,j
integer(kind=i4),dimension(:),pointer::temp
integer(kind=i4),dimension(:),pointer::x

if (n <= 0) then
 write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
 stop
endif

! If not allocated, allocate and return
if (.not.(associated(x))) then
!if (.not.(allocated(x))) then
 allocate(x(n))
 return
endif

! If reallocation, create a pointer with a target of new dimension.
allocate(temp(n))
temp=0

! (1) Expand the array dimension
if ( n > size(x) ) then
   do i = 1, size(x)
      temp(i) = x(i)
   end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
else
   j=ubound(x,1,i4)
   if(j<=0) then
      print*, "allocation issue in alloc_int_ptr..."
      print*,j,size(x)
      stop
   endif
   do i = 1, n
      temp(i) = x(i)
   end do
endif

! Destroy the target of x
  deallocate(x)

! Re-assign the pointer
 x => temp

!deallocate(temp)

return

end subroutine alloc_int_ptr1

!pure subroutine add_to(vec,val,n,chunk_size,finished)
!http://degenerateconic.com/dynamically-sizing-arrays/
pure subroutine alloc_int_ptr(vec,n)
use data_type
  
 implicit none
  
 integer(kind=i4),dimension(:),allocatable,intent(inout) :: vec
    !! the vector to add to
! integer(kind=i4) ,intent(in) :: val  
    !! the value to add
 integer(kind=i4) ,intent(inout) :: n  
    !! counter for last element added to vec.
    !! must be initialized to size(vec)
    !! (or 0 if not allocated) before first call
 !integer(kind=i4) ,intent(in) :: chunk_size  
 integer(kind=i4) :: chunk_size  
    !! allocate vec in blocks of this size (>0)
! logical,intent(in) :: finished 
    !! set to true to return vec
    !! as its correct size (n)
  
 integer(kind=i4) ,dimension(:),allocatable :: tmp

chunk_size=1 
n=n-chunk_size
 
 if (allocated(vec)) then
     if (n==size(vec)) then
         ! have to add another chunk:
         !allocate(tmp(size(vec)+chunk_size))
         allocate(tmp(size(vec)+chunk_size))
         tmp(1:size(vec)) = vec
         call move_alloc(tmp,vec)
     end if
     n = n + 1
 else
     ! the first element:
     allocate(vec(chunk_size))
     n = 1
 end if
  
! vec(n) = val
  
! if (finished) then
     ! set vec to actual size (n):
     if (allocated(tmp)) deallocate(tmp)
     allocate(tmp(n))
     tmp = vec(1:n)
     call move_alloc(tmp,vec)
! end if
  
 end subroutine alloc_int_ptr



integer(kind=i4) function append(n, array)
  integer(kind=i4),pointer,dimension(:) :: array
  integer(kind=i4),pointer,dimension(:) :: tmp_arr
  integer(kind=i4) ::n

  if (size(array) .eq. n) then
     allocate(tmp_arr(2*size(array)))
     tmp_arr(1:size(array)) = array
     deallocate(array)
     array => tmp_arr
  end if
  n = n + 1
  append = n
end function



!=====================================================================

subroutine  renumber
use grid
use param
implicit none
integer(kind=i4) :: i,j, in,out,c,it,t1
integer(kind=i4),allocatable:: oldnum(:), newnum(:)

type(cells),allocatable,dimension(:)::elmn
type(faces),allocatable,dimension(:)::fac

!do i=1,noc_bc
!   write(46,*)'old',i,(cell(i)%c2c(j),j=1,cell(i)%nc2c)
!enddo

!allocate(oldnum(noc),newnum(noc))
!allocate(elmn(noc))


allocate(oldnum(noc_bc),newnum(noc_bc))
allocate(elmn(noc_bc))

do i=1,noc_bc
   oldnum(i) = 0
   newnum(i) = 0
enddo

!c=1
!cc=minloc(cell(1:noc)%cen(1),dim=1)
!oldnum(1)=cc
!newnum(cc)=1
!is=c
!ie=c
!do !while (c<=noc)
!   mem=i2
!   do i=is,ie
!      it=oldnum(i)
!      if(it .eq. 0)then
!         print*,'renumber: Fatal error. it is zero for i=',i
!         stop
!      endif
!      do j=1,cell(it)%nc2c
!         t1=cell(it)%c2c(j)
!         if(t1.gt.0.and.newnum(t1).eq.0) then
!            c=c+1
!            oldnum(c)=t1   
!            newnum(t1)=c   
!         endif 
!      enddo
!    enddo
!    is=mem+1
!    ie=c
!    !print*,c,noc
!if (c==noc)  exit
!enddo


!do i=1,noc
!   oldnum(i) = 0
!   newnum(i) = 0
!enddo
!
c=1
!c=minloc(cell(1:noc)%cen(1),dim=1)
oldnum(1)=c
newnum(1)=c
!print*,'minloc x=',c

do i=1,noc_bc
   it=oldnum(i)
   if(it .eq. 0)then
      print*,'renumber: Fatal error. it is zero for i=',i
      stop
   endif

   do j=1,cell(it)%nc2c
      t1=cell(it)%c2c(j)
      if(t1.gt.0)  then 
      !if(t1.le.noc)  then 
      if(newnum(t1).eq.0) then
        c=c+1
        oldnum(c)=t1   
        newnum(t1)=c   
      endif 
      endif 
   enddo 

enddo

!print*,'here'
if(noc_bc.ne. c)then
  print*,'renumber: count does not match no. of cells.'
  print*,noc_bc,noc,noc_bc-noc
  print*,c,noc_bc-c
  print*,'          Possible bug'
  stop
endif


do i=1,noc
   elmn(i)%nc2c=0
   elmn(i)%cen(1)=cell(i)%cen(1)
   elmn(i)%cen(2)=cell(i)%cen(2)
   elmn(i)%cv=cell(i)%cv
   do j=1,cell(i)%nc2c
      elmn(i)%nc2c=elmn(i)%nc2c+1
      call alloc_int_ptr(elmn(i)%c2c,elmn(i)%nc2c)      
      elmn(i)%c2c(elmn(i)%nc2c)=cell(i)%c2c(j)      
   enddo 
enddo


do i=1,nof
   in=fc(i)%in
   out=fc(i)%out

   if(in/=0) then
     fc(i)%in=newnum(in)
   endif

   if(out/=0) then
     fc(i)%out=newnum(out)
   endif
enddo    


do i=1,noc
   it=oldnum(i)
   cell(i)%nc2c=0
   cell(i)%cen(1)=elmn(it)%cen(1)
   cell(i)%cen(2)=elmn(it)%cen(2)
   cell(i)%cv=elmn(it)%cv
   do j=1,elmn(it)%nc2c
      t1=elmn(it)%c2c(j)
      if(t1.gt.0) then  
         cell(i)%nc2c=cell(i)%nc2c+1
         call alloc_int_ptr(cell(i)%c2c,cell(i)%nc2c)      
         cell(i)%c2c(cell(i)%nc2c)=newnum(t1)
      else
         cell(i)%nc2c=cell(i)%nc2c+1
         call alloc_int_ptr(cell(i)%c2c,cell(i)%nc2c)      
         cell(i)%c2c(cell(i)%nc2c)=t1
      endif
   enddo 
enddo



!do i=1,noc
!   write(47,*)'new',i,(cell(i)%c2c(j),j=1,cell(i)%nc2c)
!enddo


deallocate(oldnum,newnum,elmn)

!===============================================================
return
!===============================================================


call cell_2_face

! Edge renumbering

allocate(oldnum(nof),newnum(nof))
allocate(fac(nof))

do i=1,nof
   oldnum(i) = 0
   newnum(i) = 0
enddo

c=0
!oldnum(1)=c
!newnum(1)=c

do i=1,noc
   do j=1,cell(i)%nc2f
      t1=cell(i)%c2f(j)
      !print*,i,j,t1
      !pause
      if(t1>0.and.newnum(t1)==0) then 
        c=c+1
        oldnum(c)=t1
        newnum(t1)=c
      endif 
   enddo
enddo

if(nof .ne. c)then
  print*,'renumber: count does not match nof.'
  print*,'          Possible bug'
  print*,nof,c
  stop
endif


do i=1,nof
   fac(i)%pt(1)=fc(i)%pt(1)
   fac(i)%pt(2)=fc(i)%pt(2)
   fac(i)%in   =fc(i)%in   
   fac(i)%out  =fc(i)%out  
   fac(i)%bc   =fc(i)%bc   
   fac(i)%sx   =fc(i)%sx   
   fac(i)%sy   =fc(i)%sy   
enddo

do i=1,nof
   it=oldnum(i) 
   fc(i)%pt(1)=fac(it)%pt(1)
   fc(i)%pt(2)=fac(it)%pt(2)
   fc(i)%in   =fac(it)%in   
   fc(i)%out  =fac(it)%out  
   fc(i)%bc   =fac(it)%bc   
   fc(i)%sx   =fac(it)%sx   
   fc(i)%sy   =fac(it)%sy   
enddo

deallocate(oldnum,newnum,fac)


end subroutine  renumber

!========================================================

subroutine LSQR_Coeff
use param
use grid
implicit none

integer(kind=i4) :: i,j,c
real(kind=dp)    :: dx,dy,xc,yc
real(kind=dp)    :: r11,r12,r22

do i=1,noc
   xc=cell(i)%cen(1) 
   yc=cell(i)%cen(2) 
   dx=0.d0
   dy=0.d0
   r11=0.d0
   r12=0.d0
   r22=0.d0

   do j=1,cell(i)%nc2c
      c=cell(i)%c2c(j)
      dx=cell(c)%cen(1)-xc
      dy=cell(c)%cen(2)-yc
      r11=r11+dx*dx
      r12=r12+dx*dy
      r22=r22+dy*dy
   enddo

   r11=dsqrt(r11)
   r12=r12/r11
   r22=dsqrt(r22-r12*r12 )

   cell(i)%r11=r11
   cell(i)%r12=r12
   cell(i)%r22=r22
enddo

end subroutine LSQR_Coeff


!========================================================


end subroutine geometric  
