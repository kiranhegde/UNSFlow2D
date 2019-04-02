subroutine geometric  
use grid
implicit none
integer(kind=i4) :: i,p1,p2,in,out


!do i=1,nof
!   p1=fc(i)%pt(1)
!   p2=fc(i)%pt(2)
!   in=fc(i)%in
!   out=fc(i)%out
!
!   if(p1>p2) then
!      !write(83,*)i,p1,p2
!      fc(i)%pt(1)=p2
!      fc(i)%pt(2)=p1
!      fc(i)%in   =out
!      fc(i)%out  =in
!   endif 
!enddo

! Mesh connectivity generation

call vertex_2_vertex

call renumber

call cell_face_norm_c2f_dist

call vertex_2_face 

call cell_2_face

call vertex_2_cell

call cell_2_vertex

call avg_cell_face_length

call face_co_volume

call LSQR_Coeff

!    writes cell face  normal vector's in gnuplot format to files
call write_2_file
!stop

!========================================================
contains
!========================================================

subroutine vertex_2_vertex
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: p1,p2

! nodes surrounding a node
print*
print*,"==> nodes surrounding a node"

pt(:)%nv2v=0
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)

   pt(p1)%nv2v=pt(p1)%nv2v+1
   call alloc_int_ptr(pt(p1)%v2v,pt(p1)%nv2v)      
   pt(p1)%v2v(pt(p1)%nv2v)=p2      

   pt(p2)%nv2v=pt(p2)%nv2v+1
   call alloc_int_ptr(pt(p2)%v2v,pt(p2)%nv2v)      
   pt(p2)%v2v(pt(p2)%nv2v)=p1      
enddo

print*, 'Max. nbr=',maxval(pt(:)%nv2v),'at node',maxloc(pt(:)%nv2v)
print*, 'Min. nbr=',minval(pt(:)%nv2v),'at node',minloc(pt(:)%nv2v)
 
end subroutine vertex_2_vertex

!====================================================================

subroutine  renumber
use grid
use param
implicit none
integer(kind=i4) :: i,j, in,out,c,cc,it,t1,max_bw,min_bw,is,ie,mem
integer(kind=i4),allocatable:: oldnum(:), newnum(:)

type(points),allocatable,dimension(:)::ptn
type(faces),allocatable,dimension(:)::fac

!do i=1,nop
!   write(46,*)'old',i,(pt(i)%v2v(j),j=1,pt(i)%nv2v)
!enddo

allocate(oldnum(nop),newnum(nop))
allocate(ptn(nop))

c=1
oldnum(1)=c
newnum(1)=c

do i=1,nop
   it=oldnum(i)
   if(it .eq. 0)then
      print*,'renumber: Fatal error. it is zero for i=',i
      stop
   endif

   do j=1,pt(it)%nv2v
      t1=pt(it)%v2v(j)
      if(t1.gt.0)  then 
      if(newnum(t1).eq.0) then
        c=c+1
        oldnum(c)=t1   
        newnum(t1)=c   
      endif 
      endif 
   enddo 

enddo

if(nop .ne. c)then
  print*,'renumber: count does not match nt.'
  print*,'          Possible bug'
  stop
endif

ptn(:)%nv2v=0
do i=1,nop
   ptn(i)%nv2v=0
   ptn(i)%x=pt(i)%x
   ptn(i)%y=pt(i)%y
   ptn(i)%bc=pt(i)%bc
   do j=1,pt(i)%nv2v
      ptn(i)%nv2v=ptn(i)%nv2v+1
      call alloc_int_ptr(ptn(i)%v2v,ptn(i)%nv2v)      
      ptn(i)%v2v(ptn(i)%nv2v)=pt(i)%v2v(j)      
   enddo 
enddo


do i=1,nop
   it=oldnum(i)
   pt(i)%nv2v=0
   pt(i)%x=ptn(it)%x
   pt(i)%y=ptn(it)%y
   pt(i)%bc=ptn(it)%bc
   do j=1,ptn(it)%nv2v
      t1=ptn(it)%v2v(j)
      if(t1.gt.0) then  
         pt(i)%nv2v=pt(i)%nv2v+1
         call alloc_int_ptr(pt(i)%v2v,pt(i)%nv2v)      
         pt(i)%v2v(pt(i)%nv2v)=newnum(t1)
      else
         pt(i)%nv2v=pt(i)%nv2v+1
         call alloc_int_ptr(pt(i)%v2v,pt(i)%nv2v)      
         pt(i)%v2v(pt(i)%nv2v)=t1
      endif
   enddo 
enddo

do i=1,nof
   p1 = fc(i)%pt(1)
   p2 = fc(i)%pt(2)
   fc(i)%pt(1)=newnum(p1)
   fc(i)%pt(2)=newnum(p2)
enddo

!do i=1,nop
!   write(47,*)'new',i,(pt(i)%v2v(j),j=1,pt(i)%nv2v)
!enddo

deallocate(oldnum,newnum,ptn)


end subroutine  renumber

!====================================================================
! Finding face normals 
subroutine cell_face_norm_c2f_dist
use param
use grid
implicit none

integer(kind=i4) :: i,j
integer(kind=i4) :: p1,p2,in,out
real(kind=dp)    :: x1,y1,x2,y2
real(kind=dp)    :: xc,yc,dx,dy,ds,xm,ym,sx,sy,vprod
real(kind=dp)    :: check,xfc,yfc,xc1,xc2,yc1,yc2,dotprod

fc(:)%nx=0.d0
fc(:)%ny=0.d0

! face normal components
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   xm=0.5d0*(pt(p1)%x+pt(p2)%x)
   ym=0.5d0*(pt(p1)%y+pt(p2)%y)
   dx=(pt(p2)%x-pt(p1)%x)
   dy=(pt(p2)%y-pt(p1)%y)
   fc(i)%cen(1)=xm
   fc(i)%cen(2)=ym

   in=fc(i)%in
   out=fc(i)%out

   if(in/=0) then
      xc=cell(in)%cen(1)
      yc=cell(in)%cen(2)
      sx =  yc-ym 
      sy = -(xc-xm)
      vprod = sx*dx + sy*dy
      if(vprod<0.d0) then 
        sx = -sx 
        sy = -sy 
        !fc(i)%pt(1)=p2  
        !fc(i)%pt(2)=p1  
        !fc(i)%in=out
        !fc(i)%out=in
      endif
      fc(i)%nx = sx 
      fc(i)%ny = sy 
   endif    

   if(out/=0) then
      xc=cell(out)%cen(1)
      yc=cell(out)%cen(2)

      sx = (ym-yc) 
      sy =-(xm-xc)
      vprod = sx*dx + sy*dy
      if(vprod<0.d0) then 
        sx = -sx 
        sy = -sy 
        !fc(i)%pt(1)=p2  
        !fc(i)%pt(2)=p1  
        !fc(i)%in=out
        !fc(i)%out=in
      endif

      fc(i)%nx = fc(i)%nx+sx
      fc(i)%ny = fc(i)%ny+sy
   !elseif(out==0) then
      
   endif    

enddo

do i=1,nof
   sx = fc(i)%nx
   sy = fc(i)%ny
   ds = dsqrt(sx*sx+sy*sy)
   fc(i)%nx = fc(i)%nx/ds
   fc(i)%ny = fc(i)%ny/ds
   fc(i)%area = ds
enddo 

pt(:)%cv=0.d0
! Dual-volume for node centered fvm
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   in=fc(i)%in
   out=fc(i)%out
   xm=fc(i)%cen(1)
   ym=fc(i)%cen(2)
   if(p2==0.or.p1==0) then
     print*, 'p1 or p2 is zero'
   endif

   dx=0.d0
   dy=0.d0

   if(in/=0) then
      xc=cell(in)%cen(1)
      yc=cell(in)%cen(2)
      dy = yc-ym ; dx = 0.5d0*(xc+xm)
      pt(p1)%cv=pt(p1)%cv+dx*dy 
      pt(p2)%cv=pt(p2)%cv-dx*dy 
   elseif(in==0) then
      xc=pt(p1)%x
      yc=pt(p1)%y
      dy = yc-ym ; dx = 0.5d0*(xc+xm)
      pt(p1)%cv=pt(p1)%cv+dx*dy 

      xc=pt(p2)%x
      yc=pt(p2)%y
      dy = yc-ym ; dx = 0.5d0*(xc+xm)
      pt(p2)%cv=pt(p2)%cv-dx*dy 
   endif    

   if(out/=0) then
      xc=cell(out)%cen(1)
      yc=cell(out)%cen(2)
      dy = ym-yc ; dx = 0.5d0*(xc+xm)
      pt(p1)%cv=pt(p1)%cv+dx*dy 
      pt(p2)%cv=pt(p2)%cv-dx*dy 
   elseif(out==0) then
      xc=pt(p1)%x
      yc=pt(p1)%y
      dy = ym-yc ; dx = 0.5d0*(xc+xm)
      pt(p1)%cv=pt(p1)%cv+dx*dy 

      xc=pt(p2)%x
      yc=pt(p2)%y
      dy = ym-yc ; dx = 0.5d0*(xc+xm)
      pt(p2)%cv=pt(p2)%cv-dx*dy 
   endif    

enddo

x1=0.d0
do i=1,nop
   if(pt(i)%cv<=0.d0) then
      write(84,*)i,pt(i)%x,pt(i)%y,pt(i)%cv
   endif
   x1=x1+pt(i)%cv
enddo 

x2=0.d0
do i=1,noc
   if(cell(i)%cv<=0.d0) then
      write(85,*)i,pt(i)%x,pt(i)%y,pt(i)%cv
   endif
   x2=x2+cell(i)%cv
enddo 

print*,x1,x2

!stop



!do i=1,noc
!   print*,i,(cell(i)%c2c(j), j=1,cell(i)%nc2c) 
!enddo

!print*, 'Max. cell nbr=',maxval(cell(:)%nc2c),'at node',maxloc(cell(:)%nc2c)
!print*, 'Min. cell nbr=',minval(cell(:)%nc2c),'at node',minloc(cell(:)%nc2c)
end subroutine cell_face_norm_c2f_dist
!===========================================================
subroutine vertex_2_face 
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: p1,p2



! nodes surrounding a node
print*
print*,"==> faces connecting a node"

pt(:)%nv2f=0
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)

   pt(p1)%nv2f=pt(p1)%nv2f+1
   call alloc_int_ptr(pt(p1)%v2f,pt(p1)%nv2f)      
   pt(p1)%v2f(pt(p1)%nv2f)=i       

   pt(p2)%nv2f=pt(p2)%nv2f+1
   call alloc_int_ptr(pt(p2)%v2f,pt(p2)%nv2f)      
   pt(p2)%v2f(pt(p2)%nv2f)=i      
enddo


print*, 'Max. nbr=',maxval(pt(:)%nv2f),'at node',maxloc(pt(:)%nv2f)
print*, 'Min. nbr=',minval(pt(:)%nv2f),'at node',minloc(pt(:)%nv2f)
!stop 
end subroutine vertex_2_face
!====================================================================
!Cells to face   connectivity 
subroutine cell_2_face
use param
use grid
implicit none

integer(kind=i4) :: i,j
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

print*, 'Max. cell faces =',maxval(cell(:)%nc2f),'at node',maxloc(cell(:)%nc2f)
print*, 'Min. cell faces =',minval(cell(:)%nc2f),'at node',minloc(cell(:)%nc2f)


!print*,'cell'
!do i=1,noc
!   write(80,*),i,(cell(i)%c2f(j), j=1,cell(i)%nc2f) 
!enddo
!stop

end  subroutine cell_2_face
!===========================================================

!==============================================================================
subroutine vertex_2_cell
use param
use grid
implicit none

integer(kind=i4) :: i,j
integer(kind=i4) :: p1,p2,in,out,c
real(kind=dp)    :: xc,yc,dx,dy,ds



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



print*, 'Max. nbr=',maxval(pt(:)%nv2c),'at node',maxloc(pt(:)%nv2c)
print*, 'Min. nbr=',minval(pt(:)%nv2c),'at node',minloc(pt(:)%nv2c)
end subroutine vertex_2_cell
!====================================================================
subroutine avg_cell_face_length
use param
use grid
implicit none

integer(kind=i4) :: i
integer(kind=i4) :: p1,p2
real(kind=dp)    :: sx,sy,ds



! average cell face length

pt(:)%ds=0.d0
pt(:)%sx=0.d0
pt(:)%sy=0.d0
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)

   ds=fc(i)%area
   sx=dabs(fc(i)%nx)*ds
   sy=dabs(fc(i)%ny)*ds

   pt(p1)%ds=pt(p1)%ds+ds
   pt(p1)%sx=pt(p1)%sx+sx
   pt(p1)%sy=pt(p1)%sy+sy

   pt(p2)%ds=pt(p2)%ds+ds
   pt(p2)%sx=pt(p2)%sx+sx
   pt(p2)%sy=pt(p2)%sy+sy
enddo


do i=1,nop
   !pt(i)%ds=pt(i)%ds/pt(i)%nc2f
   pt(i)%ds=0.5d0*pt(i)%ds
   pt(i)%sx=0.5d0*pt(i)%sx
   pt(i)%sy=0.5d0*pt(i)%sy
   !pt(i)%nx=pt(i)%nx/pt(i)%nc2f
   !pt(i)%ny=pt(i)%ny/pt(i)%nc2f
enddo

end subroutine avg_cell_face_length
!====================================================================
subroutine cell_2_vertex
use param
use grid
implicit none

integer(kind=i4) :: i,j,k,m1,m2
integer(kind=i4) :: p1,p2,in,out,c,f1
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
   !cell(i)%nc2v=cell(i)%nc2v+1
   !call alloc_int_ptr(cell(i)%c2v,cell(i)%nc2v)
   !cell(i)%c2v(cell(i)%nc2v)= cell(i)%c2v(1)

!   print*,i,(cell(i)%c2v(j), j=1,cell(i)%nc2v) 
enddo

print*, 'Max. cell vertex=',maxval(cell(:)%nc2v),'at node',maxloc(cell(:)%nc2v)
print*, 'Min. cell vertex=',minval(cell(:)%nc2v),'at node',minloc(cell(:)%nc2v)

end subroutine cell_2_vertex
!==========================================================
! Finding face co-volume for face & Cell gradient
subroutine face_co_volume
use param
use grid
implicit none

integer(kind=i4) :: i,j,k
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

!do i=1,noc
!   cell(i)%cov=0.d0
!enddo

do i=1,nof
   in=fc(i)%in
   out=fc(i)%out
   if(fc(i)%cov <= 0.d0) then  
      print*,i,'-ve face co-volume'
      write(75,*)i,in,out
      stop
   endif

!   if(in/=0) cell(in)%cov=cell(in)%cov+fc(i)%cov
!   if(out/=0) cell(out)%cov=cell(out)%cov+fc(i)%cov
enddo

!do i=1,noc
!   if(cell(i)%cov<=0.d0) then
!      print*,i,'-ve cell co-volume'
!      write(75,*)i,in,out
!      stop
!   endif
!enddo

end subroutine face_co_volume
!===========================================================
subroutine LSQR_Coeff
use param
use grid
implicit none

integer(kind=i4) :: i,j,c
real(kind=dp)    :: dx,dy,xc,yc
real(kind=dp)    :: r11,r12,r22,alfa1,alfa2

do i=1,nop
   xc=pt(i)%x 
   yc=pt(i)%y 
   dx=0.d0
   dy=0.d0
   r11=0.d0
   r12=0.d0
   r22=0.d0

   do j=1,pt(i)%nv2v
      c=pt(i)%v2v(j)
      dx=pt(c)%x-xc
      dy=pt(c)%y-yc
      r11=r11+dx*dx
      r12=r12+dx*dy
      r22=r22+dy*dy
   enddo

   r11=dsqrt(r11)
   r12=r12/r11
   r22=dsqrt(r22-r12*r12 )

   allocate(pt(i)%wx(pt(i)%nv2v),pt(i)%wy(pt(i)%nv2v))

   do j=1,pt(i)%nv2v
      c=pt(i)%v2v(j)
      dx=pt(c)%x-xc
      dy=pt(c)%y-yc

      alfa1=dx/r11/r11
      alfa2=(dy-dx*r12/r11)/r22/r22
      pt(i)%wx(j)=alfa1-alfa2*r12/r11
      pt(i)%wy(j)=alfa2
   enddo

enddo

end subroutine LSQR_Coeff

!===========================================================
subroutine write_2_file
use param
use grid
implicit none

integer(kind=i4) :: i,j
integer(kind=i4) :: p1,p2
real(kind=dp)    :: x1,y1,x2,y2



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
   if(fc(i)%bc==1001) write(4,100)x1,y1,fc(i)%nx,fc(i)%ny
   !if(fc(i)%bc==1001.and.out>0) write(4,100)x1,y1,cell(out)%cen(1)-x1,cell(out)%cen(2)-y1
   if(fc(i)%bc==2001) write(5,100)x1,y1,fc(i)%nx,fc(i)%ny
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
   write(3,100)x1,y1,fc(i)%nx,fc(i)%ny
enddo
close(3)

100 format(1x,4(f15.8,1x))

open(3,file='CC_plot.dat')
do i=1,noc
   write(3,*)cell(i)%cen(1),cell(i)%cen(2)
   write(3,*)
enddo
close(3)
101 format(1x,4(f15.6,1x))
end subroutine write_2_file
!=================================================================================

!===========================================================================
!                      Allocate/extend array
!===========================================================================
subroutine alloc_int_ptr(x,n)
use data_type
implicit none
integer(kind=i4),intent(in) ::n 
integer(kind=i4)::i
integer(kind=i4),dimension(:),pointer::temp
integer(kind=i4),dimension(:),pointer::x

if (n <= 0) then
 write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
 stop
endif

! If not allocated, allocate and return
if (.not.(associated(x))) then
 allocate(x(n))
 return
endif

! If reallocation, create a pointer with a target of new dimension.
allocate(temp(n))

! (1) Expand the array dimension
if ( n > size(x) ) then
   do i = 1, size(x)
      temp(i) = x(i)
   end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
else
   do i = 1, n
      temp(i) = x(i)
   end do
endif

! Destroy the target of x
!  deallocate(x)

! Re-assign the pointer
 x => temp

!deallocate(temp)

return

end subroutine alloc_int_ptr

!=====================================================================

end subroutine geometric  
