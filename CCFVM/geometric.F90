subroutine geometric  
use grid
implicit none
integer(kind=i4) :: i,in,out

!    Read grid from file
call read_grid


do i=startFC,endFC
   in=fc(i)%in
   out=fc(i)%out
 
   if(out==0.or.in ==0)  then 
      print*,in,out 
      Stop 'Inside faces in/out==0' 
   endif
enddo

do i=startBC,endBC
   in=fc(i)%in
   out=fc(i)%out
 
   if(in <=0)  then 
      print*,in,out 
      Stop 'BC faces in/out==0' 
   endif
enddo
! Mesh connectivity generation

call cell_face_norm_c2f_dist

call cell_2_cell

call renumber

call cell_2_face

!call cell_face_norm_c2f_dist

call vertex_2_cell

call cell_2_vertex

call avg_cell_face_length

call face_co_volume

! compute  geometrical co-efficient for least square 
call LSQR_Coeff

!    writes cell face  normal vector's in gnuplot format to files
call write_2_file

!========================================================
contains
!========================================================
! 
subroutine cell_face_norm_c2f_dist
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: p1,p2,in,out
real(kind=dp)    :: x1,y1,x2,y2
real(kind=dp)    :: xc,yc,dx,dy
real(kind=dp)    :: check,dotprod

! face normal components & face midpoint
do i=startFC,endFC
   in=fc(i)%in
   out=fc(i)%out
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   xc=cell(out)%cen(1)
   yc=cell(out)%cen(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   dx = x2-x1 ; dy = y2-y1
   fc(i)%sx = dy
   fc(i)%sy = -dx
   check=dotprod(p1,p2,xc,yc)

   if(check<0.0_dp) then
!     print*,'Domain face normal in->out:',i,in,out 
      print*,in,out,check,fc(i)%bc
      fc(i)%in=out
      fc(i)%out=in
      fc(i)%pt(1)=p2
      fc(i)%pt(2)=p1
      x1 = pt(p2)%x ; y1 = pt(p2)%y
      x2 = pt(p1)%x ; y2 = pt(p1)%y
      dx = x2-x1 ; dy = y2-y1
      fc(i)%sx = dy
      fc(i)%sy = -dx
   endif 
   
   fc(i)%cen(1)=0.5_dp*(pt(p1)%x+pt(p2)%x)
   fc(i)%cen(2)=0.5_dp*(pt(p1)%y+pt(p2)%y)
enddo

! making boundary face normal outward
do i=startBC,endBC

   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   in=fc(i)%in
   out=fc(i)%out

   if(fc(i)%bc==2001) then
      xc=cell(in)%cen(1)
      yc=cell(in)%cen(2)
      p1=fc(i)%pt(1)
      p2=fc(i)%pt(2)
      check=dotprod(p1,p2,xc,yc)
      if(check>0.0_dp) then
        !print*,i,'in->out',2001 
        fc(i)%pt(1)=p2
        fc(i)%pt(2)=p1
        x1 = pt(p2)%x ; y1 = pt(p2)%y
        x2 = pt(p1)%x ; y2 = pt(p1)%y
        dx = x2-x1 ; dy = y2-y1
        fc(i)%sx = dy
        fc(i)%sy = -dx
        !stop
      endif  
   endif  

   if(fc(i)%bc==1001) then
      xc=cell(in)%cen(1)
      yc=cell(in)%cen(2)
      p1=fc(i)%pt(1)
      p2=fc(i)%pt(2)
      check=dotprod(p1,p2,xc,yc)
      if(check>0.0_dp) then
        !print*,i,'in->out',1001 
        fc(i)%pt(1)=p2
        fc(i)%pt(2)=p1
        x1 = pt(p2)%x ; y1 = pt(p2)%y
        x2 = pt(p1)%x ; y2 = pt(p1)%y
        dx = x2-x1 ; dy = y2-y1
        fc(i)%sx = dy
        fc(i)%sy = -dx
        !stop
      endif  
   endif  
   fc(i)%cen(1)=0.5_dp*(pt(p1)%x+pt(p2)%x)
   fc(i)%cen(2)=0.5_dp*(pt(p1)%y+pt(p2)%y)

enddo


do i=1,nof
   in=fc(i)%in
   out=fc(i)%out

   if(in<0) then 
     fc(i)%pt(1)=p2
     fc(i)%pt(2)=p1
     x1 = pt(p2)%x ; y1 = pt(p2)%y
     x2 = pt(p1)%x ; y2 = pt(p1)%y
     dx = x2-x1 ; dy = y2-y1
     fc(i)%sx = dy
     fc(i)%sy = -dx
     fc(i)%in=out
     fc(i)%out=in
      write(33,*)x1,y1 
      write(33,*)x2,y2 
      write(33,*)
      print*,'====>',i,in,out,fc(i)%bc 
      !Stop 'in or out=0....-ve'
   endif

 
!   if(out>=0.or.in<=0)  then 
!      p1=fc(i)%pt(1)
!      p2=fc(i)%pt(2)
!      x1 = pt(p1)%x ; y1 = pt(p1)%y
!      x2 = pt(p2)%x ; y2 = pt(p2)%y
!      write(33,*)x1,y1 
!      write(33,*)x2,y2 
!      write(33,*)
!      print*,'====>',i,in,out,fc(i)%bc 
!      !Stop 'in or out=0....' 
!   endif

enddo



!print*, 'Max. cell nbr=',maxval(cell(:)%nc2c),'at node',maxloc(cell(:)%nc2c)
!print*, 'Min. cell nbr=',minval(cell(:)%nc2c),'at node',minloc(cell(:)%nc2c)
end subroutine cell_face_norm_c2f_dist
!===========================================================
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
do i=startFC,endFC
   in=fc(i)%in
   out=fc(i)%out
 
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

   if(in>0) then
   cell(in)%nc2f=cell(in)%nc2f+1
   call alloc_int_ptr(cell(in)%c2f,cell(in)%nc2f)      
   cell(in)%c2f(cell(in)%nc2f)=i      
   endif

   if(out>0) then
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
integer(kind=i4),parameter :: nn=3
real(kind=dp)    :: xy(ndim,nn)
real(kind=dp)    :: area

fc(:)%cov=0.0_dp
cell(:)%cov=0.0_dp
do i=startFC,endFC
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   in=fc(i)%in
   out=fc(i)%out

   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   !xy(1,2)=cell(out)%cen(1) ; xy(2,2)=cell(out)%cen(2)
   xy(1,2)=pt(p2)%x         ; xy(2,2)=pt(p2)%y
   xy(1,3)=cell(in)%cen(1)  ; xy(2,3)=cell(in)%cen(2)
   area=0.0_dp
   call triarea(xy,area,nn)
   fc(i)%cov=area

   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   xy(1,2)=cell(out)%cen(1) ; xy(2,2)=cell(out)%cen(2)
   xy(1,3)=pt(p2)%x         ; xy(2,3)=pt(p2)%y
   area=0.0_dp
   call triarea(xy,area,nn)
   fc(i)%cov=fc(i)%cov+area
 
!   x1 = pt(p1)%x    ; y1 = pt(p1)%y
!   x2 = cell(out)%cen(1) ; y2 = cell(out)%cen(2)
!   dx = 0.5_dp*(x2+x1) ; dy = y2-y1
!   fc(i)%cov=fc(i)%cov+dx*dy
!
!   x1 = cell(out)%cen(1) ; y1 = cell(out)%cen(2)
!   x2 = pt(p2)%x    ; y2 = pt(p2)%y
!   dx = 0.5_dp*(x2+x1) ; dy = y2-y1
!   fc(i)%cov=fc(i)%cov+dx*dy
!
!   x1 = pt(p2)%x    ; y1 = pt(p2)%y
!   x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
!   dx = 0.5_dp*(x2+x1) ; dy = y2-y1
!   fc(i)%cov=fc(i)%cov+dx*dy
!
!   x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
!   x2 = pt(p1)%x    ; y2 = pt(p1)%y
!   dx = 0.5_dp*(x2+x1) ; dy = y2-y1
!   fc(i)%cov=fc(i)%cov+dx*dy

   if(fc(i)%cov <= 0.0_dp) then  
      print*,i,'0 ? -ve face co-volume @ inner domain'
      write(75,*)i,in,out
      stop
   endif

   cell(in)%cov=cell(in)%cov+fc(i)%cov
   cell(out)%cov=cell(out)%cov+fc(i)%cov
enddo

do i=startBC,endBC
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   in=fc(i)%in
   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   xy(1,2)=pt(p2)%x         ; xy(2,2)=pt(p2)%y
   xy(1,3)=cell(in)%cen(1)  ; xy(2,3)=cell(in)%cen(2)
   area=0.0_dp
   call triarea(xy,area,nn)
   fc(i)%cov=area

!   x1 = pt(p1)%x ; y1 = pt(p1)%y
!   x2 = pt(p2)%x ; y2 = pt(p2)%y
!   dx = 0.5_dp*(x2+x1) ; dy = y2-y1
!   fc(i)%cov=fc(i)%cov+dx*dy
!
!   x1 = pt(p2)%x    ; y1 = pt(p2)%y
!   x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
!   dx = 0.5_dp*(x2+x1) ; dy = y2-y1
!   fc(i)%cov=fc(i)%cov+dx*dy
!
!   x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
!   x2 = pt(p1)%x    ; y2 = pt(p1)%y
!   dx = 0.5_dp*(x2+x1) ; dy = y2-y1
!   fc(i)%cov=fc(i)%cov+dx*dy
   if(fc(i)%cov <= 0.0_dp) then  
      print*,i,'? 1 -ve face co-volume @ BC'
      write(75,*)i,in,out
      stop
   endif
   cell(in)%cov=cell(in)%cov+fc(i)%cov
enddo

!do i=1,noc
!   cell(i)%cov=0.0_dp
!enddo
!
!do i=1,nof
!   in=fc(i)%in
!   out=fc(i)%out
!   if(fc(i)%cov <= 0.0_dp) then  
!      print*,i,'-ve face co-volume @ BC or inner domain'
!      write(75,*)i,in,out
!      stop
!   endif
!
!   if(in>0) cell(in)%cov=cell(in)%cov+fc(i)%cov
!   if(out>0) cell(out)%cov=cell(out)%cov+fc(i)%cov
!enddo

do i=1,noc
   if(cell(i)%cov<=0.0_dp) then
      print*,i,'-ve cell co-volume'
      write(75,*)i,in,out
      stop
   endif
enddo

end subroutine face_co_volume

subroutine triarea(xy,area,nn)
implicit none
integer(kind=i4) :: i,nn,ip1
real(kind=dp)    :: xy(ndim,nn)
real(kind=dp)    :: area,dx,dy

! area of a polygon
area=0.0_dp
do i=1,nn
   ip1=i+1
   if(i==nn) ip1=1
   dx=0.5_dp*(xy(1,ip1)+xy(1,i))
   dy=        xy(2,ip1)-xy(2,i)
   area=area+dx*dy
enddo

if (area <= 0.0_dp) stop 'geom,triarea,  -ve area'
end subroutine triarea



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

   if(in>0) then
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

   if(out>0) then
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
      pt(i)%wt(j)=1._dp/dsqrt(dx*dx+dy*dy)
      if(pt(i)%wt(j)<0.0_dp) then
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
!   cell(i)%ds=0.0_dp
!   cell(i)%dx=0.0_dp
!   cell(i)%dy=0.0_dp
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
!   cell(i)%dx=0.5_dp*cell(i)%dx
!   cell(i)%dy=0.5_dp*cell(i)%dy
!   if(cell(i)%dx<=0.0_dp.or.cell(i)%dy<=0.0_dp) then 
!     print*, 'Time step length scale  zero'
!     stop 
!   endif
!enddo

cell(:)%ds=0.0_dp
cell(:)%dx=0.0_dp
cell(:)%dy=0.0_dp
do i=1,nof
   in=fc(i)%in
   out=fc(i)%out
   sx=0.5_dp*fc(i)%sx
   sy=0.5_dp*fc(i)%sy
   ds=dsqrt(sx*sx+sy*sy)
   if(in>0) then 
     cell(in)%ds=cell(in)%ds+ds
     cell(in)%dx=cell(in)%dx+sx
     cell(in)%dy=cell(in)%dy+sy
   endif 

   if(out>0) then 
     cell(out)%ds=cell(out)%ds+ds
     cell(out)%dx=cell(out)%dx+sx
     cell(out)%dy=cell(out)%dy+sy
   endif 
enddo


!do i=1,noc
!   !cell(i)%ds=cell(i)%ds/cell(i)%nc2f
!   cell(i)%ds=0.5_dp*cell(i)%ds
!   cell(i)%dx=0.5_dp*cell(i)%dx
!   cell(i)%dy=0.5_dp*cell(i)%dy
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

      if(check<0.0_dp) then
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
            if(check<0.0_dp) then
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

integer(kind=i4) :: i
integer(kind=i4) :: p1,p2
real(kind=dp)    :: x1,y1,x2,y2

x1=0.0_dp
x2=0.0_dp
y1=0.0_dp
y2=0.0_dp

open(3,file='BC-Vec_plot.dat')
open(4,file='BC_plot.dat')
do i=startBC,endBC
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   write(4,*)x1,y1
   write(4,*)x2,y2
   write(4,*)
   x1 = (x2+x1)/2.0_dp ; y1= (y2+y1)/2.0_dp
   write(3,100)x1,y1,fc(i)%sx,fc(i)%sy
enddo
close(3)
close(4)


open(3,file='Mesh_plot.dat')
open(4,file='MeshVec_plot.dat')
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   write(3,101)x1,y1
   write(3,101)x2,y2
   write(3,*)
   write(4,100)x1,y1,x2-x1,y2-y1
enddo
close(4)
close(3)


open(4,file='InCell.dat')
open(5,file='OutCell.dat')
do i=startFC,endFC
   in=fc(i)%in
   out=fc(i)%out
   write(4,101)cell(in)%cen(1),cell(in)%cen(2)
   write(5,101)cell(out)%cen(1),cell(out)%cen(2)
enddo
close(4)
close(5)

open(3,file="NormalsOut.dat")  
open(4,file="NormalsIn.dat")  
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   x1 = pt(p1)%x ; y1 = pt(p1)%y
   x2 = pt(p2)%x ; y2 = pt(p2)%y
   x1 = (x2+x1)/2.0_dp ; y1= (y2+y1)/2.0_dp
   write(3,100)x1,y1,fc(i)%sx*0.5,fc(i)%sy*0.5
   write(4,100)x1,y1,-fc(i)%sx*0.5,-fc(i)%sy*0.5
enddo
close(4)
close(3)

100 format(1x,4(f15.8,1x))
101 format(1x,2(f15.8,1x))

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
!subroutine alloc_int_ptr1(x,n)
!use data_type
!implicit none
!integer(kind=i4),intent(in) ::n 
!integer(kind=i4)::i,j
!integer(kind=i4),dimension(:),pointer::temp
!integer(kind=i4),dimension(:),pointer::x
!
!if (n <= 0) then
! write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
! stop
!endif
!
!! If not allocated, allocate and return
!if (.not.(associated(x))) then
!!if (.not.(allocated(x))) then
! allocate(x(n))
! return
!endif
!
!! If reallocation, create a pointer with a target of new dimension.
!allocate(temp(n))
!temp=0
!
!! (1) Expand the array dimension
!if ( n > size(x) ) then
!   do i = 1, size(x)
!      temp(i) = x(i)
!   end do
!
!! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
!else
!   j=ubound(x,1,i4)
!   if(j<=0) then
!      print*, "allocation issue in alloc_int_ptr..."
!      print*,j,size(x)
!      stop
!   endif
!   do i = 1, n
!      temp(i) = x(i)
!   end do
!endif
!
!! Destroy the target of x
!  deallocate(x)
!
!! Re-assign the pointer
! x => temp
!
!!deallocate(temp)
!
!return
!
!end subroutine alloc_int_ptr1

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

!do i=1,noc
!   write(46,*)'old',i,(cell(i)%c2c(j),j=1,cell(i)%nc2c)
!enddo

!allocate(oldnum(noc),newnum(noc))
!allocate(elmn(noc))


allocate(oldnum(noc),newnum(noc))
allocate(elmn(noc))

do i=1,noc
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

do i=1,noc
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
if(noc.ne. c)then
  print*,'renumber: count does not match no. of cells.'
  print*,noc,noc,noc-noc
  print*,c,noc-c
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

   if(in>0) then
     fc(i)%in=newnum(in)
   endif

   if(out>0) then
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

integer(kind=i4) :: i
real(kind=dp)    :: dx,dy,xo,yo,xi,yi,ds
real(kind=dp)    :: r11,r12,r22

cell(:)%r11=0.0_dp
cell(:)%r12=0.0_dp
cell(:)%r22=0.0_dp
ds=1.0_dp

do i=startFC,endFC
   in  =fc(i)%in   
   out =fc(i)%out  
   xi = cell(in)%cen(1)
   yi = cell(in)%cen(2)
   xo = cell(out)%cen(1)
   yo = cell(out)%cen(2)
   dx = xo-xi 
   dy = yo-yi 
   !ds=1.0_dp/dsqrt(dx*dx+dy*dy)
   cell(in )%r11=cell(in )%r11+ds*dx*dx
   cell(in )%r12=cell(in )%r12+ds*dx*dy
   cell(in )%r22=cell(in )%r22+ds*dy*dy

   dx = xi-xo 
   dy = yi-yo 
   cell(out)%r11=cell(out)%r11+ds*dx*dx
   cell(out)%r12=cell(out)%r12+ds*dx*dy
   cell(out)%r22=cell(out)%r22+ds*dy*dy
enddo

do i=startBC,endBC
   in =fc(i)%in   
   xi = cell(in)%cen(1)
   yi = cell(in)%cen(2)
   xo = fc(i)%cen(1)
   yo = fc(i)%cen(2)
   dx = xo-xi 
   dy = yo-yi 
   !ds=1.0_dp/dsqrt(dx*dx+dy*dy)
   cell(in)%r11=cell(in)%r11+ds*dx*dx
   cell(in)%r12=cell(in)%r12+ds*dx*dy
   cell(in)%r22=cell(in)%r22+ds*dy*dy 
enddo

do i=1,noc
   r11=dsqrt(cell(i)%r11)
   cell(i)%r11=r11
   r12=cell(i)%r12/r11 
   cell(i)%r12=r12
   r22=cell(i)%r22
   cell(i)%r22=dsqrt(r22-r12*r12)
enddo


end subroutine LSQR_Coeff


!subroutine LSQR_Coeff
!use param
!use grid
!implicit none
!
!integer(kind=i4) :: i,j,c,in,out
!real(kind=dp)    :: dx,dy,xc,yc
!real(kind=dp)    :: r11,r12,r22
!real(kind=dp),dimension(2*ndim,noc):: dmat
!
!cell(:)%r11=0.0_dp
!cell(:)%r12=0.0_dp
!cell(:)%r22=0.0_dp
!cell(:)%det=0.0_dp
!
!do i=startFC,endFC
!   in   =fc(i)%in   
!   out  =fc(i)%out  
!   dx=cell(out)%cen(1)-cell(in)%cen(1)
!   dy=cell(out)%cen(2)-cell(in)%cen(2)
!   cell(in )%r11=cell(in )%r11+dx*dx
!   cell(in )%r12=cell(in )%r12+dx*dy
!   cell(in )%r22=cell(in )%r22+dy*dy
!
!   dx=cell(in)%cen(1)-cell(out)%cen(1)
!   dy=cell(in)%cen(2)-cell(out)%cen(2)
!   cell(out)%r11=cell(out)%r11+dx*dx
!   cell(out)%r12=cell(out)%r12+dx*dy
!   cell(out)%r22=cell(out)%r22+dy*dy
!enddo
!
!do i=startBC,endBC
!   in   =fc(i)%in   
!   !out  =fc(i)%out  
!   !if(out<0) print*,i,in,out 
!     
!   dx=fc(i)%cen(1)-cell(in)%cen(1)
!   dy=fc(i)%cen(2)-cell(in)%cen(2)
!
!   cell(in)%r11=cell(in)%r11+dx*dx
!   cell(in)%r12=cell(in)%r12+dx*dy
!   cell(in)%r22=cell(in)%r22+dy*dy
!enddo
!
!do i=1,noc
!   r11=cell(i)%r11
!   r12=cell(i)%r12
!   r22=cell(i)%r22
!   cell(i)%det=r11*r22-r12*r12
!   if(dabs(cell(i)%det)<=eps) then
!     print*,i," det -->  0",r11,r12,r22
!     stop
!   endif  
!enddo 
!
!
!end subroutine LSQR_Coeff

end subroutine geometric  
