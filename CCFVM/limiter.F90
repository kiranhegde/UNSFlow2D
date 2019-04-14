
!======================================================================================
subroutine limit
implicit none

!call limiter_VKN
call limiter_VKN1
!call limiter_vanAlbada
!call limiter_min_max


end subroutine limit

!======================================================================================
subroutine limiter_VKN 
use grid 
use commons
implicit none
integer(kind=i4):: i,j,c1,c2,cc
real(kind=dp)   :: nr,dr,phi0,kappa,onethird
real(kind=dp)   :: TOL,alfa,D_L,dist(ndim),ll,ds,qmax(nvar),qmin(nvar)
real(kind=dp)   :: checkmax(nvar),checkmin(nvar)
real(kind=dp),allocatable,save :: umin(:), umax(:), phi(:)

!allocate(dumin(nvar,noc), dumax(nvar,noc), phi(nvar,noc))
if(.not. allocated(umax))then
    allocate(umax(nvar))
    allocate(umin(nvar))
    allocate(phi(nvar))
endif

onethird=1.0_dp/3.0_dp
!kappa=0.8_dp
kappa=0.05_dp
checkmax=-1e30
checkmin=1e30

do j=1,nvar
   qmax(j)=maxval(cell(:)%qp(j))
   qmin(j)=minval(cell(:)%qp(j))
enddo

do i=1,noc
     ll=0.0_dp
     umax=cell(i)%qp(:) 
     umin=cell(i)%qp(:) 

     do c1=1,cell(i)%nc2c
        cc=cell(i)%c2c(c1) 
        do j=1,nvar
           ! Max/Min b/n cells sharing a face
           umax(j)=dmax1(umax(j),cell(cc)%qp(j))
           umin(j)=dmin1(umin(j),cell(cc)%qp(j))
        enddo 
        dist(:)=cell(cc)%cen(:)-cell(i)%cen(:)
        ds=dsqrt(dist(1)**2+dist(2)**2) 
        !ll=dmax1(ll,ds)
        ll=ll+ds 
     enddo 
        ll=ll/cell(i)%nc2c

     !TOL=(kappa*dsqrt(cell(i)%cv))**onethird
     !TOL=kappa*ll**3
     phi=cell(i)%phi
     !phi=1.0d0
     do c1=1,cell(i)%nc2f

        cc=cell(i)%c2f(c1)
        do j=1,nvar
              TOL=kappa*(qmax(j)-qmin(j))
              dist(:)=fc(cc)%cen(:)-cell(i)%cen(:)
              D_L=sum(cell(i)%grad(:,j)*dist(:))
    
              phi0=1.0_dp 
              if(D_L>0.d0) then
                 alfa=umax(j)-cell(i)%qp(j)
                 phi0=vk_limit(alfa,D_L,tol) 
              elseif(D_L<0.d0) then
                 alfa=umin(j)-cell(i)%qp(j)
                 phi0=vk_limit(alfa,D_L,tol) 
              endif
     
              phi(j)=dmin1(phi(j),phi0)
        enddo

     enddo

     do j=1,nvar 
        cell(i)%phi(j)=phi(j)
        cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*phi(j)
        checkmax(j)=dmax1(checkmax(j),phi(j))
        checkmin(j)=dmin1(checkmin(j),phi(j))
     enddo
enddo

!print*,checkmax(1),checkmin(1)
!print*,checkmax(2),checkmin(2)
!print*,checkmax(3),checkmin(3)
!print*,checkmax(4),checkmin(4)


!deallocate(dumin,dumax,phi)
contains

real(kind=dp) function vk_limit(a,b,tol)
implicit none
real(kind=dp)   :: nr,dr,tol,a,b

  nr=a*a+2.0_dp*a*b+tol*tol
  dr=a*a+2.0_dp*b*b+a*b+tol*tol
  vk_limit=nr/dr       

end function vk_limit


end subroutine limiter_VKN 
!===============================================
subroutine limiter_VKN1 
use grid 
use commons
implicit none
integer(kind=i4):: i,j,in,out,p1,p2
real(kind=dp)   :: nr,dr,phi0,kappa,du,onethird,var
real(kind=dp)   :: TOL,alfa,D_L,dist(ndim),qmax(nvar),qmin(nvar)
real(kind=dp),allocatable,save :: umin(:,:), umax(:,:), phi(:,:)
!real(kind=dp),allocatable :: umin(:,:), umax(:,:), phi(:,:)

onethird=1.0_dp/3.0_dp
!kappa=3.0_dp
kappa=0.05_dp
!kappa=0.8_dp

!allocate(dumin(nvar,noc), dumax(nvar,noc), phi(nvar,noc))
if(.not. allocated(umax))then
    allocate(umax(nvar,noc))
    allocate(umin(nvar,noc))
    allocate(phi(nvar,noc))
endif


do j=1,nvar
   qmax(j)=maxval(cell(:)%qp(j))
   qmin(j)=minval(cell(:)%qp(j))
enddo


do i=1,noc
   phi(:,i)=cell(i)%phi(:)
   umax(:,i)=cell(i)%qp(:)
   umin(:,i)=cell(i)%qp(:)
enddo

!phi=1.0_dp
!do j=1,nvar
!print*,"=====>",maxval(cell(:)%phi(j),mask=(cell(:)%phi(j)<1.0_dp)  ),minval(cell(:)%phi(j),mask=(cell(:)%phi(j)>0.0_dp))
!print*,"=====>",maxval(phi(j,:))
!enddo
!stop

do i=startFC,endFC
   in = fc(i)%in
   out = fc(i)%out

   do j=1,nvar
      ! Max/Min b/n cells sharing a face
      umax(j,in)=dmax1(umax(j,in),cell(out)%qp(j))
      umin(j,in)=dmin1(umin(j,in),cell(out)%qp(j))
      umax(j,out)=dmax1(umax(j,out),cell(in)%qp(j))
      umin(j,out)=dmin1(umin(j,out),cell(in)%qp(j))
   enddo
enddo


do i=startBC,endBC
   in = fc(i)%in
   p1=fc(i)%pt(1)
   p2=fc(i)%pt(2)
   do j=1,nvar
      var=0.5_dp*(pt(p1)%prim(j)+pt(p2)%prim(j))
      umax(j,in)=dmax1(umax(j,in),var)
      umin(j,in)=dmin1(umin(j,in),var)
   enddo
enddo

!print*,maxval(umax(1,:)),minval(umax(1,:))
!print*,maxval(umax(2,:)),minval(umax(2,:))
!print*,maxval(umax(3,:)),minval(umax(3,:))
!print*,maxval(umax(4,:)),minval(umax(4,:))

!print*,maxval(umin(1,:)),minval(umin(1,:))
!print*,maxval(umin(2,:)),minval(umin(2,:))
!print*,maxval(umin(3,:)),minval(umin(3,:))
!print*,maxval(umin(4,:)),minval(umin(4,:))
!stop

do i=startFC,endFC
   in = fc(i)%in
   out = fc(i)%out

   if(in<0.or.out<0) then
      print*,'-ve value in face neigbhor'
      print*,i,in,out  
   endif

      do j=1,nvar
              TOL=kappa*(qmax(j)-qmin(j))
              !TOL=(kappa*dsqrt(cell(in)%cv))**onethird
              dist(:)=fc(i)%cen(:)-cell(in)%cen(:)
              D_L=sum(cell(in)%grad(:,j)*dist(:))
    
              phi0=1.0_dp 
              if(D_L>0.0_dp) then
                 alfa=umax(j,in)-cell(in)%qp(j)
                 phi0=vk_limit(alfa,D_L,tol) 
              elseif(D_L<0.0_dp) then
                 alfa=umin(j,in)-cell(in)%qp(j)
                 phi0=vk_limit(alfa,D_L,tol) 
              endif
     
              phi(j,in)=dmin1(phi(j,in),phi0)

              !TOL=(kappa*dsqrt(cell(out)%cv))**onethird
              dist(:)=fc(i)%cen(:)-cell(out)%cen(:)
              D_L=sum(cell(out)%grad(:,j)*dist(:))
    
              phi0=1.0_dp 
              if(D_L>0.0_dp) then
                 alfa=umax(j,out)-cell(out)%qp(j)
                 phi0=vk_limit(alfa,D_L,tol) 
              elseif(D_L<0.0_dp) then
                 alfa=umin(j,out)-cell(out)%qp(j)
                 phi0=vk_limit(alfa,D_L,tol) 
              endif
     
              phi(j,out)=dmin1(phi(j,out),phi0)
        enddo
enddo

do i=startBC,endBC
   in = fc(i)%in

      dist(:)=fc(i)%cen(:)-cell(in)%cen(:)
      !TOL=(kappa*dsqrt(cell(in)%cv))**onethird

      do j=1,nvar
           D_L=sum(cell(in)%grad(:,j)*dist(:))
           TOL=kappa*(qmax(j)-qmin(j))
             
           phi0=1.0_dp    
           if(D_L>0.0_dp) then
              alfa=umax(j,in)-cell(in)%qp(j)
              phi0=vk_limit(alfa,D_L,tol) 
           elseif(D_L<0.0_dp) then
              alfa=umin(j,in)-cell(in)%qp(j)
              phi0=vk_limit(alfa,D_L,tol) 
           endif
  
           !cell(in)%phi(j)=dmin1(cell(in)%phi(j),phi0)
           phi(j,in)=dmin1(phi(j,in),phi0)
      enddo    
enddo

!print*,maxval(phi(1,:)),minval(phi(1,:))
!print*,maxval(phi(2,:)),minval(phi(2,:))
!print*,maxval(phi(3,:)),minval(phi(3,:))
!print*,maxval(phi(4,:)),minval(phi(4,:))

do i=1,noc
   do j=1,nvar 
      cell(i)%phi(j)=phi(j,i)
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
      !cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*phi(j,i)
   enddo
enddo


!deallocate(umin,umax,phi)

contains

real(kind=dp) function vk_limit(a,b,tol)
implicit none
real(kind=dp)   :: nr,dr,tol,a,b

  nr=a*a+2.0_dp*a*b+tol*tol
  dr=a*a+2.0_dp*b*b+a*b+tol*tol
  vk_limit=nr/dr       

end function vk_limit

end subroutine limiter_VKN1 

!======================================================================================

subroutine limiter_vanAlbada
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,c1,c2
real(kind=dp) :: dist(ndim)
real(kind=dp) :: phi,du
real(kind=dp) :: D_L
real(kind=dp),allocatable :: dumin(:,:), dumax(:,:), phi0(:,:)

!allocate(dumin(nvar,noc), dumax(nvar,noc), phi(nvar,noc))


!DUmin(:,:)=eps
!DUmax(:,:)=-eps

do i=startFC,endFC
   c1 = fc(i)%in
   c2 = fc(i)%out
   if(c1/=0.and.c2/=0) then 
   do j=1,nvar
      du=cell(c2)%qp(j)-cell(c1)%qp(j)
      cell(c1)%dumax(j)=dmax1(cell(c1)%dumax(j), du)
      cell(c1)%dumin(j)=dmin1(cell(c1)%dumin(j), du)
      cell(c2)%dumax(j)=dmax1(cell(c2)%dumax(j),-du)
      cell(c2)%dumin(j)=dmin1(cell(c2)%dumin(j),-du)
   enddo
   endif
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out

   if(c1>0.and.c2>0) then 
   do j=1,nvar

         !---------------i
         dist(:)=fc(i)%cen(:)-cell(c1)%cen(:)
         D_L=sum(cell(c1)%grad(:,j)*dist(:))
         if(D_L>0.d0) then
            phi=Albada(cell(c1)%dumax(j),D_L)
         else
            phi=Albada(cell(c1)%dumin(j),D_L)
         endif

         cell(c1)%phi(j)=dmin1(cell(c1)%phi(j),phi)

         !---------------j
         dist(:)=fc(i)%cen(:)-cell(c2)%cen(:)
         D_L=sum(cell(c2)%grad(:,j)*dist(:))
         !D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
         if(D_L>0.d0) then
            phi=Albada(cell(c2)%dumax(j),D_L)
         else
            phi=Albada(cell(c2)%dumin(j),D_L)
         endif

         cell(c2)%phi(j)=dmin1(cell(c2)%phi(j),phi)

   enddo
   endif
enddo

do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo

contains
!-----------------------------------------------------------------------------
!real(kind=dp) function VanAlbada(rk)
real(kind=dp) function Albada(a,b)
implicit none 
real(kind=dp) :: a,b


!VanAlbada=(rk*rk+rk)/(1.d0+rk*rk)
!VanAlbada=(rk*rk+2.d0*rk)/(rk*rk+rk+2.d0)
Albada=dmax1(0.d0,(2.d0*a*b+eps*eps)/(a*a+b*b+eps*eps) )


end function Albada

end subroutine limiter_vanAlbada

!-----------------------------------------------------------------------------

subroutine limiter_min_max
use grid 
use commons
implicit none
integer(kind=i4) :: i,j,c1,c2

real(kind=dp) :: dist(ndim)
real(kind=dp) :: phi
real(kind=dp) :: D_L,du,D_L0

do i=1,noc
   cell(i)%DUmin(:)=eps
   cell(i)%DUmax(:)=-eps
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out
   if(c1/=0.and.c2/=0) then 
   do j=1,nvar
      du=cell(c2)%qp(j)-cell(c1)%qp(j)
      cell(c1)%dumax(j)=dmax1(cell(c1)%dumax(j), du)
      cell(c1)%dumin(j)=dmin1(cell(c1)%dumin(j), du)
      cell(c2)%dumax(j)=dmax1(cell(c2)%dumax(j),-du)
      cell(c2)%dumin(j)=dmin1(cell(c2)%dumin(j),-du)
   enddo
   endif
enddo

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out

   if(c1/=0.and.c2/=0) then 
   do j=1,nvar

         !---------------i
         dist(:)=fc(i)%cen(:)-cell(c1)%cen(:)
         D_L=sum(cell(c1)%grad(:,j)*dist(:))
         D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
         if(D_L>0.d0) then
            phi=dmin1(1.d0,cell(c1)%dumax(j) /D_L0)
         else
            phi=dmin1(1.d0,cell(c1)%dumin(j) /D_L0)
         endif

         cell(c1)%phi(j)=dmin1(cell(c1)%phi(j),phi)

         !---------------j
         dist(:)=fc(i)%cen(:)-cell(c2)%cen(:)
         D_L=sum(cell(c2)%grad(:,j)*dist(:))
         D_L0=dsign(D_L,1.d0)*(dabs(D_L)+eps) 
         if(D_L>0.d0) then
            phi=dmin1(1.d0,cell(c2)%dumax(j) /D_L0)
         else
            phi=dmin1(1.d0,cell(c2)%dumin(j) /D_L0)
         endif

         cell(c2)%phi(j)=dmin1(cell(c2)%phi(j),phi)

   enddo
   endif
enddo


do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo
end subroutine limiter_min_max
