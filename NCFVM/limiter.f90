
!======================================================================================
subroutine LinearReconLimiter 
implicit none

call limiter_VKN
!call limiter_vanAlbada
!call limiter_min_max


end subroutine LinearReconLimiter

!======================================================================================

subroutine limiter_VKN 
use grid 
use commons
implicit none
integer(kind=i4):: i,j,k,c,p1,p2
real(kind=dp)   :: nr,dr,phi,q_min,q_max,kappa
real(kind=dp)   :: TOL,alfa,D_L,var,D_L0,du
real(kind=dp)   :: pv(nvar),dist(ndim)

do i=1,nop
   pt(i)%DUmin(:)=eps
   pt(i)%DUmax(:)=-eps 
enddo

do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   do j=1,nvar
      du=pt(p2)%qp(j)-pt(p1)%qp(j)  
      pt(p1)%dumax(j)=dmax1(pt(p1)%dumax(j), du)
      pt(p1)%dumin(j)=dmin1(pt(p1)%dumin(j), du)
      pt(p2)%dumax(j)=dmax1(pt(p2)%dumax(j),-du)
      pt(p2)%dumin(j)=dmin1(pt(p2)%dumin(j),-du)
   enddo  
enddo

kappa=0.01
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   dist(1)=pt(p2)%x-pt(p1)%x
   dist(2)=pt(p2)%y-pt(p1)%y 
   
   do j=1,nvar 

         !---------------i
         D_L=0.5d0*(sum(pt(p1)%grad(:,j)*dist(:)))

         !TOL=(kappa*(pt(p1)%dumax(j)-pt(p1)%dumin(j)))**2
         TOL=(kappa*dsqrt(pt(p1)%cv))**3
 
         if(D_L>0.d0) then 
            alfa=pt(p1)%dumax(j)
         else
            alfa=pt(p1)%dumin(j)
         endif

         if(dabs(alfa)<eps) alfa=0.d0

         nr=alfa*alfa+2.d0*D_L*alfa+TOL
         dr=alfa*alfa+2.d0*D_L*D_L+D_L*alfa+TOL    
         !phi=nr/dr
         phi=dmin1(1.d0,nr/dr)
         pt(p1)%phi(j)=dmin1(pt(p1)%phi(j),phi)

         !---------------j
         D_L=-0.5d0*(sum(pt(p2)%grad(:,j)*dist(:)))

         !TOL=(kappa*(pt(p2)%dumax(j)-pt(p2)%dumin(j)))**2
         TOL=(kappa*dsqrt(pt(p2)%cv))**3
 
         if(D_L>0.d0) then 
            alfa=pt(p2)%dumax(j)
         else
            alfa=pt(p2)%dumin(j)
         endif

         if(dabs(alfa)<eps) alfa=0.d0

         nr=alfa*alfa+2.d0*D_L*alfa+TOL
         dr=alfa*alfa+2.d0*D_L*D_L+D_L*alfa+TOL    
         phi=dmin1(1.d0,nr/dr)
         pt(p2)%phi(j)=dmin1(pt(p2)%phi(j),phi)

   enddo  
enddo  

!do i=1,nop
!   do j=1,nvar 
!      pt(i)%grad(1:ndim,j)=pt(i)%grad(1:ndim,j)*pt(i)%phi(j)
!   enddo
!enddo

end subroutine limiter_VKN 

!======================================================================================

subroutine limiter_vanAlbada
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,ie,p1,p2,c
real(kind=dp) :: q2, pv(nvar)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar),dist(ndim)
real(kind=dp) :: nr, dr, phi, q_min, q_max
real(kind=dp) :: alfa,D_L,D_L0,rk,tol,kappa,du

do i=1,nop
   pt(i)%DUmin(:)=eps
   pt(i)%DUmax(:)=-eps 
enddo

do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   do j=1,nvar
      du=pt(p2)%qp(j)-pt(p1)%qp(j)  
      pt(p1)%dumax(j)=dmax1(pt(p1)%dumax(j), du)
      pt(p1)%dumin(j)=dmin1(pt(p1)%dumin(j), du)
      pt(p2)%dumax(j)=dmax1(pt(p2)%dumax(j),-du)
      pt(p2)%dumin(j)=dmin1(pt(p2)%dumin(j),-du)
   enddo  
enddo

kappa=0.3
do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   dist(1)=pt(p2)%x-pt(p1)%x
   dist(2)=pt(p2)%y-pt(p1)%y 
   
   do j=1,nvar 

         !---------------i
         D_L=0.5d0*(sum(pt(p1)%grad(:,j)*dist(:)))

         !TOL=(kappa*(pt(p1)%dumax(j)-pt(p1)%dumin(j)))**2
         TOL=(kappa*dsqrt(pt(p1)%cv))**3

         if(D_L>0.d0) then 
            phi=Albada(pt(p1)%dumax(j),D_L,tol)
         else
            phi=Albada(pt(p2)%dumin(j),D_L,tol)
         endif
      
         pt(p1)%phi(j)=dmin1(pt(p1)%phi(j),phi)

         !---------------j
         D_L=-0.5d0*(sum(pt(p2)%grad(:,j)*dist(:)))

         !TOL=(kappa*(pt(p2)%dumax(j)-pt(p2)%dumin(j)))**2
         TOL=(kappa*dsqrt(pt(p2)%cv))**3

         if(D_L>0.d0) then 
            phi=Albada(pt(p2)%dumax(j),D_L,tol)
         else
            phi=Albada(pt(p2)%dumin(j),D_L,tol)
         endif
      
         pt(p2)%phi(j)=dmin1(pt(p2)%phi(j),phi)
   enddo  
enddo  

!do i=1,nop
!   do j=1,nvar 
!      pt(i)%grad(1:ndim,j)=pt(i)%grad(1:ndim,j)*pt(i)%phi(j)
!   enddo
!enddo

contains
!-----------------------------------------------------------------------------
real(kind=dp) function Albada(a,b,tol)
implicit none 
real(kind=dp) :: rk,a,b,tol


Albada=dmax1(0.d0,(2.d0*a*b+eps)/(a*a+b*b+eps))


end function Albada

end subroutine limiter_vanAlbada

!-----------------------------------------------------------------------------

subroutine limiter_min_max
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,p1,p2

real(kind=dp) :: var,dist(ndim)
real(kind=dp) :: phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,du

do i=1,nop
   pt(i)%DUmin(:)=eps
   pt(i)%DUmax(:)=-eps 
enddo

do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   do j=1,nvar
      du=pt(p2)%qp(j)-pt(p1)%qp(j)  
      pt(p1)%dumax(j)=dmax1(pt(p1)%dumax(j), du)
      pt(p1)%dumin(j)=dmin1(pt(p1)%dumin(j), du)
      pt(p2)%dumax(j)=dmax1(pt(p2)%dumax(j),-du)
      pt(p2)%dumin(j)=dmin1(pt(p2)%dumin(j),-du)
   enddo  
enddo


do i=1,nof
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   dist(1)=pt(p2)%x-pt(p1)%x
   dist(2)=pt(p2)%y-pt(p1)%y 
   
   do j=1,nvar 
         !---------------i
         D_L=0.5d0*(sum(pt(p1)%grad(:,j)*dist(:)))

         if(D_L>0.d0 ) then
            phi=dmin1(1.d0,pt(p1)%dumax(j)/D_L)
         elseif(D_L< 0.d0) then
            phi=dmin1(1.d0,pt(p1)%dumin(j)/D_L)
         else
            phi=1.d0
         endif 
      
         pt(p1)%phi(j)=dmin1(pt(p1)%phi(j),phi)

         !---------------i
         D_L=0.5d0*(sum(pt(p2)%grad(:,j)*dist(:)))

         if(D_L>0.d0 ) then
            phi=dmin1(1.d0,pt(p2)%dumax(j)/D_L)
         elseif(D_L< 0.d0) then
            phi=dmin1(1.d0,pt(p2)%dumin(j)/D_L)
         else
            phi=1.d0
         endif 
      
         pt(p2)%phi(j)=dmin1(pt(p2)%phi(j),phi)
   enddo  
enddo  

!do i=1,nop
!   do j=1,nvar 
!      pt(i)%grad(1:ndim,j)=pt(i)%grad(1:ndim,j)*pt(i)%phi(j)
!   enddo
!enddo
end subroutine limiter_min_max


!-----------------------------------------------------------------------------

subroutine muscle_va(ie,left,right)
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,k,c,p1,p2,ie

real(kind=dp) :: var,dist(ndim),left(nvar),right(nvar)
real(kind=dp) :: phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_p,D_m,du,lim

left=0.d0
right=0.d0

kappa=1.d0/3.d0
p1=fc(ie)%pt(1)
p2=fc(ie)%pt(2)
dist(1)=pt(p2)%x-pt(p1)%x
dist(2)=pt(p2)%y-pt(p1)%y 
   
do j=1,nvar 

   du=pt(p2)%qp(j)-pt(p1)%qp(j) 
   D_m=2.0d0*(sum(pt(p1)%grad(:,j)*dist(:)))-du
   D_p=2.0d0*(sum(pt(p2)%grad(:,j)*dist(:)))-du
   !---------------i
   lim=va(D_m,du,pt(p1)%cv)
   left(j) =pt(p1)%qp(j)+0.25d0*((1.d0-kappa*lim)*D_m+(1.d0+kappa*lim)*du)*lim
   !---------------j
   lim=va(D_p,du,pt(p2)%cv)
   right(j)=pt(p2)%qp(j)-0.25d0*((1.d0-kappa*lim)*D_p+(1.d0+kappa*lim)*du)*lim

enddo

contains
real(kind=dp) function va(a,b,cv)
implicit none
real(kind=dp) :: a,b
real(kind=dp) :: nr,dr,tol,ka,cv

!va=dmax1(0.d0,(2.d0*a*b+eps)/(a*a+b*b+eps))
! VenkataKrishnan
ka=0.3
TOL=(ka*dsqrt(cv))**3
nr=a*a+2.d0*b*a+tol
dr=a*a+2.d0*b*b+a*b+tol    
va=dmin1(1.d0,nr/dr)

end function va

end subroutine muscle_va

