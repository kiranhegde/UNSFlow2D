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
   cell(i)%grad=0.d0
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

   do k=1,nvar
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

   do k=1,nvar
      var=0.5_dp*(pt(p1)%prim(k)+pt(p2)%prim(k))-cell(in)%qp(k)
      cell(in)%grad(1,k)=cell(in)%grad(1,k)+var*wx0
      cell(in)%grad(2,k)=cell(in)%grad(2,k)+var*wy0
   enddo
enddo




if(ILIMIT==1) call limit

end subroutine Gradient_LSQR

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
!   cell(i)%grad(:,:)=0.d0
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

!======================================================================================
subroutine Gradient_GG
use grid 
use commons
use param,only:ilimit

implicit none
integer(kind=i4) :: i,j,in,out
real(kind=dp) :: var,dx,dy

do i=1,noc
   cell(i)%grad(:,:)=0.d0
enddo

do i=startFC,endFC
   in = fc(i)%in
   out = fc(i)%out
   dx=fc(i)%sx
   dy=fc(i)%sy

   do j=1,nvar
      var=fc(i)%qp(j)
      cell(in)%grad(1,j)=cell(in)%grad(1,j)+var*dx
      cell(in)%grad(2,j)=cell(in)%grad(2,j)+var*dy
      cell(out)%grad(1,j)=cell(out)%grad(1,j)-var*dx
      cell(out)%grad(2,j)=cell(out)%grad(2,j)-var*dy
   enddo
enddo

do i=startBC,endBC
   in = fc(i)%in
   dx=fc(i)%sx
   dy=fc(i)%sy

   do j=1,nvar
      var=fc(i)%qp(j)
      cell(in)%grad(1,j)=cell(in)%grad(1,j)+var*dx
      cell(in)%grad(2,j)=cell(in)%grad(2,j)+var*dy
   enddo
enddo


do i=1,noc
   var=cell(i)%cv
   do j=1,nvar
   cell(i)%grad(1,j)=cell(i)%grad(1,j)/var
   cell(i)%grad(2,j)=cell(i)%grad(2,j)/var
   enddo
enddo

if(ILIMIT==1) call limit

end subroutine Gradient_GG

!======================================================================================

subroutine Gradient_GG_FC
use grid 
use commons
use param,only:ilimit
implicit none
integer(kind=i4) :: i,j,in,out,p1,p2
integer(kind=i4),parameter :: nn=3
!real(kind=dp) :: q2, pv(nvar)
real(kind=dp) :: dx,dy

real(kind=dp) :: a1,a2
real(kind=dp) :: gradxy(ndim,nvar)
real(kind=dp) :: prim(nvar,nn),grad(ndim,nvar),xy(ndim,nn)


do i=1,noc
   cell(i)%grad(:,:)=0.d0
enddo

! GG gradient using diamond path around an edge 
! Face gradient and cell gradient is calculated 
do i=startFC,endFC
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
   ! gradient to the right of cell
   grad=0.0_dp
   call gradtrixy(prim,xy,grad,a1,nn)
   gradxy=grad

   ! triangular cell to the left  of the edge
   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   xy(1,2)=pt(p2)%x         ; xy(2,2)=pt(p2)%y
   xy(1,3)=cell(in)%cen(1)  ; xy(2,3)=cell(in)%cen(2)
   prim(:,1)=pt(p1)%prim(:)       
   prim(:,2)=pt(p2)%prim(:)       
   prim(:,3)=cell(in)%qp(:)

   ! gradient to the left  of cell
   call gradtrixy(prim,xy,grad,a2,nn)
   gradxy=gradxy+grad
   
   ! adding contribution to the cells sharing the face 
   cell(in)%grad(1,:)=cell(in)%grad(1,:)+gradxy(1,:)
   cell(in)%grad(2,:)=cell(in)%grad(2,:)+gradxy(2,:)
   cell(out)%grad(1,:)=cell(out)%grad(1,:)+gradxy(1,:)
   cell(out)%grad(2,:)=cell(out)%grad(2,:)+gradxy(2,:)
   ! face gradient
   fc(i)%grad(1,:)=gradxy(1,:)/fc(i)%cov 
   fc(i)%grad(2,:)=gradxy(2,:)/fc(i)%cov     
enddo
 
! gradient at the boundary faces and to the cells   
do i=startBC,endBC
   in = fc(i)%in
   p1=fc(i)%pt(1)       
   p2=fc(i)%pt(2)       

   grad=0.0_dp
   xy(1,1)=pt(p1)%x         ; xy(2,1)=pt(p1)%y
   xy(1,2)=pt(p2)%x         ; xy(2,2)=pt(p2)%y
   xy(1,3)=cell(in)%cen(1)  ; xy(2,3)=cell(in)%cen(2)

   prim(:,1)=pt(p1)%prim(:)       
   prim(:,2)=pt(p2)%prim(:)       
   prim(:,3)=cell(in)%qp(:)
   call gradtrixy(prim,xy,grad,a1,nn)
   gradxy=grad

   cell(in)%grad(1,:)=cell(in)%grad(1,:)+gradxy(1,:)
   cell(in)%grad(2,:)=cell(in)%grad(2,:)+gradxy(2,:)
   fc(i)%grad(1,:)=gradxy(1,:)/fc(i)%cov 
   fc(i)%grad(2,:)=gradxy(2,:)/fc(i)%cov     
enddo

! Finally cell gradient, using cell co-volume
do i=1,noc
   do j=1,nvar
      cell(i)%grad(1,j)= cell(i)%grad(1,j)/cell(i)%cov
      cell(i)%grad(2,j)= cell(i)%grad(2,j)/cell(i)%cov
   enddo
enddo

if(ILIMIT==1) call limit


!100 format(1x,4(f15.8,1x))
contains

subroutine gradtrixy(prim,xy,grad,area,nn)
implicit none
integer(kind=i4) :: i,nn,ip1
real(kind=dp)    :: prim(nvar,nn),grad(ndim,nvar),xy(ndim,nn)
real(kind=dp)    :: area 

! area of a polygon
area=0.0_dp
do i=1,nn
   ip1=i+1
   if(i==nn) ip1=1
   dx=0.5_dp*(xy(1,ip1)+xy(1,i)) 
   dy=        xy(2,ip1)-xy(2,i)
   area=area+dx*dy
enddo

! gradient in a triangular cell
grad(1,:)= 0.5_dp*(prim(:,1)*(xy(2,2)-xy(2,3))+prim(:,2)*(xy(2,3)-xy(2,1)) & 
          &      +prim(:,3)*(xy(2,1)-xy(2,2)))
grad(2,:)=-0.5_dp*(prim(:,1)*(xy(1,2)-xy(1,3))+prim(:,2)*(xy(1,3)-xy(1,1)) &
          &      +prim(:,3)*(xy(1,1)-xy(1,2)))

end subroutine gradtrixy

end subroutine Gradient_GG_FC


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
!   cell(i)%grad(:,:)=0.d0
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
!   var(:)=0.5d0*(prim(:,i)+prim(:,ip1)) 
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
