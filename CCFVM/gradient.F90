!======================================================================================
subroutine Gradient_LSQR
use grid 
use commons
use param,only:ilimit
implicit none
integer(kind=i4) :: i,k,in,out,p1,p2
real(kind=dp)    :: xn,yn,xc,yc,dx,dy,ds,var
real(kind=dp)    :: wx0,wy0,wx1,wy1     
real(kind=dp)    :: r11,r12,r22,alfa1,alfa2


do i=1,noc
   cell(i)%grad=0.0_dp
enddo

do i=startFC,endFC
   in =fc(i)%in
   out=fc(i)%out
   xc=cell(in)%cen(1)
   yc=cell(in)%cen(2)
   xn=cell(out)%cen(1)
   yn=cell(out)%cen(2)
   dx=xn-xc
   dy=yn-yc
   ds=1.0_dp/dsqrt(dx*dx+dy*dy)
   r11=cell(in)%r11
   r12=cell(in)%r12
   r22=cell(in)%r22

   alfa1=dx/r11/r11
   alfa2=(dy-dx*r12/r11)/(r22*r22)
   wx0=alfa1-alfa2*r12/r11
   wy0=alfa2

   dx=xc-xn
   dy=yc-yn
   r11=cell(out)%r11
   r12=cell(out)%r12
   r22=cell(out)%r22

   alfa1=dx/r11/r11
   alfa2=(dy-dx*r12/r11)/(r22*r22)
   wx1=alfa1-alfa2*r12/r11
   wy1=alfa2

   do k=1,ngrad
      var=cell(out)%qp(k)-cell(in)%qp(k) 
      cell(in)%grad(1,k)=cell(in)%grad(1,k)+var*wx0
      cell(in)%grad(2,k)=cell(in)%grad(2,k)+var*wy0
      cell(out)%grad(1,k)=cell(out)%grad(1,k)-var*wx1
      cell(out)%grad(2,k)=cell(out)%grad(2,k)-var*wy1
   enddo

enddo


do i=startBC,endBC
   in =fc(i)%in
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   dx=fc(i)%cen(1)-cell(in)%cen(1)
   dy=fc(i)%cen(2)-cell(in)%cen(2)
   ds=1.0_dp/dsqrt(dx*dx+dy*dy)
   r11=cell(in)%r11
   r12=cell(in)%r12
   r22=cell(in)%r22
   alfa1=dx/r11/r11
   alfa2=(dy-dx*r12/r11)/r22/r22
   wx0=alfa1-alfa2*r12/r11
   wy0=alfa2

   do k=1,ngrad
      var=0.5_dp*(pt(p1)%prim(k)+pt(p2)%prim(k))-cell(in)%qp(k)
      cell(in)%grad(1,k)=cell(in)%grad(1,k)+var*wx0
      cell(in)%grad(2,k)=cell(in)%grad(2,k)+var*wy0
   enddo
enddo




if(ILIMIT==1) call limit

end subroutine Gradient_LSQR
!======================================================================================
subroutine Gradient_GG
use grid 
use commons
use param,only:ilimit

implicit none
integer(kind=i4) :: i,j,in,out,p1,p2
real(kind=dp) :: var,dx,dy

do i=1,noc
   cell(i)%grad(:,:)=0.0_dp
enddo

do i=startFC,endFC
   in = fc(i)%in
   out = fc(i)%out
   p1=fc(i)%pt(1)       
   p2=fc(i)%pt(2)       
   !dx= pt(p2)%y-pt(p1)%y    
   !dy=-(pt(p2)%x-pt(p1)%x)    
   dx=fc(i)%sx
   dy=fc(i)%sy

   do j=1,ngrad
      var=fc(i)%qp(j)
      cell(in)%grad(1,j)=cell(in)%grad(1,j)+var*dx
      cell(in)%grad(2,j)=cell(in)%grad(2,j)+var*dy
      cell(out)%grad(1,j)=cell(out)%grad(1,j)-var*dx
      cell(out)%grad(2,j)=cell(out)%grad(2,j)-var*dy
   enddo
enddo

do i=startBC,endBC
   in = fc(i)%in
   p1=fc(i)%pt(1)       
   p2=fc(i)%pt(2)       
   !dx= pt(p2)%y-pt(p1)%y    
   !dy=-(pt(p2)%x-pt(p1)%x)    
   dx=-fc(i)%sx
   dy=-fc(i)%sy

   do j=1,ngrad
      var=fc(i)%qp(j)
      cell(in)%grad(1,j)=cell(in)%grad(1,j)+var*dx
      cell(in)%grad(2,j)=cell(in)%grad(2,j)+var*dy
   enddo
enddo


do i=1,noc
   var=cell(i)%cv
   do j=1,ngrad
      cell(i)%grad(1,j)=cell(i)%grad(1,j)/var
      cell(i)%grad(2,j)=cell(i)%grad(2,j)/var
   enddo
enddo

if(ILIMIT==1) call limit

end subroutine Gradient_GG
!======================================================================================
! GG gradient using diamond path around an edge 
subroutine Gradient_GG_FC
use grid 
use commons
use param,only:ilimit
implicit none
integer(kind=i4) :: i,j,in,out,p1,p2
integer(kind=i4),parameter :: nn=3
real(kind=dp) :: dx,dy

real(kind=dp) :: a1,a2
real(kind=dp) :: grad1(ndim,ngrad)
real(kind=dp) :: grad2(ndim,ngrad)
real(kind=dp) :: prim(npvar,nn),xy(ndim,nn)


do i=1,noc
   cell(i)%grad(:,:)=0.0_dp
enddo

do i=1,nof
   fc(i)%grad(:,:) =0.0_dp
enddo

! Face gradient and cell gradient is calculated 
do i=startFC,endFC
   a1=0.0_dp
   a2=0.0_dp

   in = fc(i)%in
   out = fc(i)%out
   p1=fc(i)%pt(1)       
   p2=fc(i)%pt(2)       

   ! triangular cell to the right of the edge
   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   xy(1,2)=cell(out)%cen(1) ; xy(2,2)=cell(out)%cen(2)
   xy(1,3)=pt(p2)%x         ; xy(2,3)=pt(p2)%y
   prim(:,1)=pt(p1)%prim(:)       
   prim(:,2)=cell(out)%qp(:)
   prim(:,3)=pt(p2)%prim(:)       
   grad1=0.0_dp 
   ! gradient to the right of cell
   call gradtrixy(prim,xy,grad1,a1,nn)

   ! triangular cell to the left  of the edge
   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   xy(1,2)=pt(p2)%x         ; xy(2,2)=pt(p2)%y
   xy(1,3)=cell(in)%cen(1)  ; xy(2,3)=cell(in)%cen(2)
   prim(:,1)=pt(p1)%prim(:)       
   prim(:,2)=pt(p2)%prim(:)       
   prim(:,3)=cell(in)%qp(:)

   ! gradient to the left  of cell
   grad2=0.0_dp
   call gradtrixy(prim,xy,grad2,a2,nn)

   ! face gradient
   fc(i)%grad(1,:)=(a1*grad1(1,:)+a2*grad2(1,:))/(a1+a2) 
   fc(i)%grad(2,:)=(a1*grad1(2,:)+a2*grad2(2,:))/(a1+a2) 
   
   ! adding contribution to the cells sharing the face 
   cell(in)%grad(1,:)=cell(in)%grad(1,:)+ fc(i)%grad(1,:)*(a1+a2)
   cell(in)%grad(2,:)=cell(in)%grad(2,:)+ fc(i)%grad(2,:)*(a1+a2)
   cell(out)%grad(1,:)=cell(out)%grad(1,:)+ fc(i)%grad(1,:)*(a1+a2)
   cell(out)%grad(2,:)=cell(out)%grad(2,:)+ fc(i)%grad(2,:)*(a1+a2)
enddo
 
! gradient at the boundary faces and to the cells   
do i=startBC,endBC
   in = fc(i)%in
   p1=fc(i)%pt(1)       
   p2=fc(i)%pt(2)       

   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   xy(1,2)=pt(p2)%x         ; xy(2,2)=pt(p2)%y
   xy(1,3)=cell(in)%cen(1)  ; xy(2,3)=cell(in)%cen(2)

   prim(:,1)=pt(p1)%prim(:)       
   prim(:,2)=pt(p2)%prim(:)       
   prim(:,3)=cell(in)%qp(:)
   grad1=0.0_dp
   call gradtrixy(prim,xy,grad1,a1,nn)

   fc(i)%grad(1,:)=grad1(1,:)
   fc(i)%grad(2,:)=grad1(2,:)  
   cell(in)%grad(1,:)=cell(in)%grad(1,:)+fc(i)%grad(1,:)*a1
   cell(in)%grad(2,:)=cell(in)%grad(2,:)+fc(i)%grad(2,:)*a1
enddo

! Finally cell gradient, using cell co-volume
do i=1,noc
   do j=1,npvar
      cell(i)%grad(1,j)= cell(i)%grad(1,j)/cell(i)%cov
      cell(i)%grad(2,j)= cell(i)%grad(2,j)/cell(i)%cov
   enddo
enddo

if(ILIMIT==1) call limit


!100 format(1x,4(f15.8,1x))
contains

subroutine gradtrixy(prim,xy,grad,tarea,nn)
implicit none
integer(kind=i4) :: i,nn,ip1
real(kind=dp)    :: prim(npvar,nn),grad(ndim,ngrad),xy(ndim,nn)
real(kind=dp)    :: tarea,var(npvar) 

! area of a polygon
tarea=0.0_dp
do i=1,nn
   ip1=i+1
   if(i==nn) ip1=1
   dx=0.5_dp*(xy(1,ip1)+xy(1,i)) 
   dy=        xy(2,ip1)-xy(2,i)
   tarea=tarea+dx*dy
enddo

if (tarea <= 0.0_dp) stop 'ggfc,gradtrixy,  -ve area'

! gradient in a triangular cell
!grad(1,:)= 0.5_dp*(prim(:,1)*(xy(2,2)-xy(2,3))+prim(:,2)*(xy(2,3)-xy(2,1)) & 
!          &      +prim(:,3)*(xy(2,1)-xy(2,2)))
!grad(2,:)=-0.5_dp*(prim(:,1)*(xy(1,2)-xy(1,3))+prim(:,2)*(xy(1,3)-xy(1,1)) &
!          &      +prim(:,3)*(xy(1,1)-xy(1,2)))

do i=1,nn
   ip1=i+1
   if(i==nn) ip1=1
   var(:)=0.5_dp*(prim(:,i)+prim(:,ip1)) 
   dx=  xy(2,ip1)-xy(2,i) 
   dy=-(xy(1,ip1)-xy(1,i))
!   dx=  xy(1,ip1)-xy(1,i) 
!   dy=  xy(2,ip1)-xy(2,i) 
   grad(1,:)=grad(1,:)+var(:)*dx
   grad(2,:)=grad(2,:)+var(:)*dy
enddo

grad(1,:)=grad(1,:)/tarea
grad(2,:)=grad(2,:)/tarea

 
end subroutine gradtrixy

end subroutine Gradient_GG_FC

!======================================================================================
!subroutine Gradient_LSQR
!use grid 
!use commons
!use param,only:ilimit
!implicit none
!integer(kind=i4) :: i,j,k,c,in,out,p1,p2
!real(kind=dp)    :: xc,yc,xn,yn,dx,dy,wt,var,det
!real(kind=dp)    :: wx,wy     
!real(kind=dp)    :: r11,r12,r22,alfa1,alfa2
!real(kind=dp)    :: qxy(2,nvar,noc)
!
!do i=1,noc
!   cell(i)%grad(:,:)=0.0_dp
!enddo
!
!qxy=0.0_dp
!
!do i=startFC,endFC
!   in =fc(i)%in
!   out=fc(i)%out
!   xc=cell(in)%cen(1)
!   yc=cell(in)%cen(2)
!   xn=cell(out)%cen(1)
!   yn=cell(out)%cen(2)
!   dx=xn-xc
!   dy=yn-yc
!
!   do k=1,nvar
!      var=cell(out)%qp(k)-cell(in)%qp(k)
!      qxy(1,k,in)=qxy(1,k,in)+var*dx    
!      qxy(2,k,in)=qxy(2,k,in)+var*dy    
!      var=-var
!      dx=-dx
!      dy=-dy
!      qxy(1,k,out)=qxy(1,k,out)+var*dx    
!      qxy(2,k,out)=qxy(2,k,out)+var*dy    
!   enddo
!enddo
!
!
!do i=startBC,endBC
!   in =fc(i)%in
!   p1=fc(i)%pt(1)
!   p2=fc(i)%pt(2)
!   xc=cell(in)%cen(1)
!   yc=cell(in)%cen(2)
!   xn=fc(i)%cen(1)
!   yn=fc(i)%cen(2)
!   dx=xn-xc
!   dy=yn-yc
!
!   do k=1,nvar
!      var=0.5_dp*(pt(p1)%prim(k)+pt(p2)%prim(k))
!      var=var-cell(in)%qp(k)
!      qxy(1,k,in)=qxy(1,k,in)+var*dx    
!      qxy(2,k,in)=qxy(2,k,in)+var*dy    
!   enddo
!enddo
!
!do i=1,noc
!   det=cell(i)%det
!   r11=cell(i)%r11/det
!   r12=cell(i)%r12/det
!   r22=cell(i)%r22/det
!   do k=1,nvar
!      cell(i)%grad(1,k)= r22*qxy(1,k,i)-r12*qxy(2,k,i)
!      cell(i)%grad(2,k)=-r12*qxy(1,k,i)+r11*qxy(2,k,i)
!   enddo
!enddo
!
!if(ILIMIT==1) call limit
!
!end subroutine Gradient_LSQR
!
!!======================================================================================
!subroutine Gradient_GG_FC1
!use grid 
!use commons
!use param,only:ilimit
!implicit none
!integer(kind=i4) :: i,j,in,out,p1,p2
!real(kind=dp)    :: pv(nvar)
!real(kind=dp)    :: x1,y1,x2,y2,dx,dy
!
!
!do i=1,noc
!   cell(i)%grad(:,:)=0.d0
!enddo
!
!do i=startFC,endFC
!      in = fc(i)%in
!      out = fc(i)%out
!      p1=fc(i)%pt(1)       
!      p2=fc(i)%pt(2)       
!      fc(i)%grad(1,:)=0.d0
!      fc(i)%grad(2,:)=0.d0
!
!      do j=1,nvar    
!                  pv(j)=0.5d0*(pt(p1)%prim(j)+cell(out)%qp(j)) 
!                  x1 = pt(p1)%x    ; y1 = pt(p1)%y
!                  x2 = cell(out)%cen(1) ; y2 = cell(out)%cen(2)
!                  dx = y2-y1       ;  dy = -(x2-x1) 
!                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
!                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
!         
!                  pv(j)=0.5d0*(cell(out)%qp(j)+pt(p2)%prim(j)) 
!                  x1 = cell(out)%cen(1) ; y1 = cell(out)%cen(2)
!                  x2 = pt(p2)%x    ; y2 = pt(p2)%y
!                  dx = y2-y1       ;  dy = -(x2-x1) 
!                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
!                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
!              
!                  pv(j)=0.5d0*(pt(p2)%prim(j)+cell(in)%qp(j)) 
!                  x1 = pt(p2)%x    ; y1 = pt(p2)%y
!                  x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
!                  dx = y2-y1       ;  dy = -(x2-x1) 
!                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
!                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
!              
!                  pv(j)=0.5d0*(cell(in)%qp(j)+pt(p1)%prim(j)) 
!                  x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
!                  x2 = pt(p1)%x    ; y2 = pt(p1)%y
!                  dx = y2-y1       ;  dy = -(x2-x1) 
!                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
!                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
!     enddo
!     fc(i)%grad(1,:)=fc(i)%grad(1,:)/fc(i)%cov
!     fc(i)%grad(2,:)=fc(i)%grad(2,:)/fc(i)%cov
!
!     cell(in)%grad(1,:)=cell(in)%grad(1,:)+fc(i)%grad(1,:)*fc(i)%cov  
!     cell(in)%grad(2,:)=cell(in)%grad(2,:)+fc(i)%grad(2,:)*fc(i)%cov  
!     cell(out)%grad(1,:)=cell(out)%grad(1,:)+fc(i)%grad(1,:)*fc(i)%cov  
!     cell(out)%grad(2,:)=cell(out)%grad(2,:)+fc(i)%grad(2,:)*fc(i)%cov  
!enddo
!    
!do i=startBC,endBC
!      in = fc(i)%in
!      p1=fc(i)%pt(1)       
!      p2=fc(i)%pt(2)       
!
!      fc(i)%grad(1,:)=0.d0
!      fc(i)%grad(2,:)=0.d0
!
!      do j=1,nvar    
!                 pv(j)=0.5d0*(pt(p1)%prim(j)+pt(p2)%prim(j)) 
!                 x1 = pt(p1)%x ; y1 = pt(p1)%y
!                 x2 = pt(p2)%x ; y2 = pt(p2)%y
!                 dx = y2-y1    ;  dy = -(x2-x1) 
!                 fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
!                 fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
!        
!                 pv(j)=0.5d0*(pt(p2)%prim(j)+cell(in)%qp(j)) 
!                 x1 = pt(p2)%x    ; y1 = pt(p2)%y
!                 x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
!                 dx = y2-y1       ;  dy = -(x2-x1) 
!                 fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
!                 fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
!             
!                 pv(j)=0.5d0*(cell(in)%qp(j)+pt(p1)%prim(j)) 
!                 x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
!                 x2 = pt(p1)%x    ; y2 = pt(p1)%y
!                 dx = y2-y1       ;  dy = -(x2-x1) 
!                 fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
!                 fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
!     enddo
!     fc(i)%grad(1,:)=fc(i)%grad(1,:)/fc(i)%cov
!     fc(i)%grad(2,:)=fc(i)%grad(2,:)/fc(i)%cov
!    
!     cell(in)%grad(1,:)=cell(in)%grad(1,:)+fc(i)%grad(1,:)*fc(i)%cov  
!     cell(in)%grad(2,:)=cell(in)%grad(2,:)+fc(i)%grad(2,:)*fc(i)%cov  
!
!enddo
!
!
!do i=1,noc
!   do j=1,nvar
!   cell(i)%grad(1,j)= cell(i)%grad(1,j)/cell(i)%cov
!   cell(i)%grad(2,j)= cell(i)%grad(2,j)/cell(i)%cov
!   enddo
!enddo
!
!if(ILIMIT==1) call limit
!
!
!!100 format(1x,4(f15.8,1x))
!end subroutine Gradient_GG_FC1


!======================================================================================

!subroutine Gradient_GG_FC
!use grid 
!use commons
!use param,only:ilimit
!implicit none
!integer(kind=i4) :: i,j,k,ie,in,out,p1,p2,c
!integer(kind=i4),parameter :: nn=4
!real(kind=dp) :: q2, pv(nvar)
!real(kind=dp)    :: x1,y1,x2,y2,dx,dy
!
!real(kind=dp) :: ql(nvar),qr(nvar)!,gradx(nvar),grady(nvar)
!real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
!real(kind=dp) :: TOL,alfa,D_L,D_L0
!real(kind=dp) :: prim(nvar,nn),grad(ndim,nvar),xy(ndim,nn)
!
!
!do i=1,noc
!   cell(i)%grad(:,:)=0.0_dp
!enddo
!
!do i=startFC,endFC
!   in = fc(i)%in
!   out = fc(i)%out
!   p1=fc(i)%pt(1)       
!   p2=fc(i)%pt(2)       
!
!   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
!   xy(1,2)=cell(out)%cen(1) ; xy(2,2)=cell(out)%cen(2)
!   xy(1,3)=pt(p2)%x         ; xy(2,3)=pt(p2)%y
!   xy(1,4)=cell(in)%cen(1)  ; xy(2,4)=cell(in)%cen(2)
!
!   prim(:,1)=pt(p1)%prim(:)       
!   prim(:,2)=cell(out)%qp(:)
!   prim(:,3)=pt(p2)%prim(:)       
!   prim(:,4)=cell(in )%qp(:)
!
!   call gradxy(prim,xy,grad,nn)
!
!   cell(in)%grad(1,:)=cell(in)%grad(1,:)+grad(1,:)
!   cell(in)%grad(2,:)=cell(in)%grad(2,:)+grad(2,:)
!   cell(out)%grad(1,:)=cell(out)%grad(1,:)+grad(1,:)
!   cell(out)%grad(2,:)=cell(out)%grad(2,:)+grad(2,:)
!   fc(i)%grad(1,:)=grad(1,:)/fc(i)%cov 
!   fc(i)%grad(2,:)=grad(2,:)/fc(i)%cov     
!enddo
!    
!do i=startBC,endBC
!   in = fc(i)%in
!   p1=fc(i)%pt(1)       
!   p2=fc(i)%pt(2)       
!
!   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
!   xy(1,2)=pt(p2)%x         ; xy(2,2)=pt(p2)%y
!   xy(1,3)=cell(in)%cen(1)  ; xy(2,3)=cell(in)%cen(2)
!
!   prim(:,1)=pt(p1)%prim(:)       
!   prim(:,2)=pt(p2)%prim(:)       
!   prim(:,3)=cell(in)%qp(:)
!
!   call gradxy(prim,xy,grad,nn-1)
!
!   cell(in)%grad(1,:)=cell(in)%grad(1,:)+grad(1,:)
!   cell(in)%grad(2,:)=cell(in)%grad(2,:)+grad(2,:)
!   fc(i)%grad(1,:)=grad(1,:)/fc(i)%cov 
!   fc(i)%grad(2,:)=grad(2,:)/fc(i)%cov     
!enddo
!
!
!do i=1,noc
!   do j=1,nvar
!      cell(i)%grad(1,j)= cell(i)%grad(1,j)/cell(i)%cov
!      cell(i)%grad(2,j)= cell(i)%grad(2,j)/cell(i)%cov
!   enddo
!enddo
!
!if(ILIMIT==1) call limit
!
!
!!100 format(1x,4(f15.8,1x))
!contains
!
!subroutine gradxy(prim,xy,grad,nn)
!implicit none
!integer(kind=i4) :: i,j,nn,ip1
!real(kind=dp)    :: prim(nvar,nn),grad(ndim,nvar),xy(ndim,nn),var(nvar)
!
!grad=0.0_dp
!do i=1,nn
!   ip1=i+1
!   if(i==nn) ip1=1
!   var(:)=0.5_dp*(prim(:,i)+prim(:,ip1)) 
!   dx=  xy(2,ip1)-xy(2,i) 
!   dy=-(xy(1,ip1)-xy(1,i))
!   grad(1,:)=grad(1,:)+var(:)*dx
!   grad(2,:)=grad(2,:)+var(:)*dy
!enddo
!
!end subroutine gradxy
!
!
!end subroutine Gradient_GG_FC
!
