
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
real(kind=dp)   :: TOL,alfa,D_L,dist(ndim),ll,ds
real(kind=dp),allocatable,save :: umin(:), umax(:), phi(:)

!allocate(dumin(nvar,noc), dumax(nvar,noc), phi(nvar,noc))
if(.not. allocated(umax))then
    allocate(umax(nvar))
    allocate(umin(nvar))
    allocate(phi(nvar))
endif

onethird=1.0_dp/3.0_dp
kappa=0.8_dp

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

     !TOL=(kappa*(cell(i)%cv)**(onethird))**3
     TOL=(kappa*dsqrt(cell(i)%cv))**0.5_dp
     !TOL=kappa*ll**3
     do c1=1,cell(i)%nc2f

        cc=cell(i)%c2f(c1)
        do j=1,nvar
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
        cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*phi(j)
     enddo
enddo


!deallocate(dumin,dumax,phi)
contains

real(kind=dp) function vk_limit(a,b,tol)
implicit none
real(kind=dp)   :: nr,dr,tol,a,b

  nr=a*a+2.0_dp*a*b+tol
  dr=a*a+2.0_dp*b*b+a*b+tol
  vk_limit=nr/dr       

end function vk_limit


end subroutine limiter_VKN 

subroutine limiter_VKN1 
use grid 
use commons
implicit none
integer(kind=i4):: i,j,c1,c2
real(kind=dp)   :: nr,dr,phi0,kappa,du,onethird
real(kind=dp)   :: TOL,alfa,D_L,dist(ndim),qmax(nvar),qmin(nvar)
real(kind=dp),allocatable,save :: umin(:,:), umax(:,:)

onethird=1.0_dp/3.0_dp
kappa=3.0_dp

!allocate(dumin(nvar,noc), dumax(nvar,noc), phi(nvar,noc))
if(.not. allocated(umax))then
    allocate(umax(nvar,noc))
    allocate(umin(nvar,noc))
endif

umax(:,:)=-1e30
umin(:,:)= 1e30
qmax(:)=-1e30
qmin(:)= 1e30

do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out

   if(c1/=0.and.c2/=0) then 
      do j=1,nvar
         ! Max/Min b/n cells sharing a face
         umax(j,c1)=dmax1(umax(j,c1),cell(c1)%qp(j),cell(c2)%qp(j))
         umin(j,c1)=dmin1(umin(j,c1),cell(c1)%qp(j),cell(c2)%qp(j))
         umax(j,c2)=dmax1(umax(j,c2),cell(c1)%qp(j),cell(c2)%qp(j))
         umin(j,c2)=dmin1(umin(j,c2),cell(c1)%qp(j),cell(c2)%qp(j))
         qmax(j)=dmax1(qmax(j),umax(j,c1),umax(j,c2)) 
         qmin(j)=dmin1(qmin(j),umin(j,c1),umin(j,c2)) 
      enddo
   elseif(c1/=0.and.c2==0) then 
      do j=1,nvar
         umax(j,c1)=dmax1(umax(j,c1),cell(c1)%qp(j))
         umin(j,c1)=dmin1(umin(j,c1),cell(c1)%qp(j))
         qmax(j)=dmax1(qmax(j),umax(j,c1)) 
         qmin(j)=dmin1(qmin(j),umin(j,c1)) 
      enddo
   elseif(c1==0.and.c2/=0) then 
      do j=1,nvar
         umax(j,c2)=dmax1(umax(j,c2),cell(c2)%qp(j))
         umin(j,c2)=dmin1(umin(j,c2),cell(c2)%qp(j))
         qmax(j)=dmax1(qmax(j),umax(j,c2)) 
         qmin(j)=dmin1(qmin(j),umin(j,c2)) 
      enddo
   endif
       
enddo


kappa=0.15_dp
do i=1,nof
   c1 = fc(i)%in
   c2 = fc(i)%out


   if(c1/=0) then 
      do j=1,nvar
           dist(:)=fc(i)%cen(:)-cell(c1)%cen(:)
           D_L=sum(cell(c1)%grad(:,j)*dist(:))
  
           TOL=kappa*(qmax(j)-qmin(j))
        
           phi0=1.0_dp    
           if(D_L>0.d0) then
              alfa=umax(j,c1)-cell(c1)%qp(j)
              phi0=vk_limit(alfa,D_L,tol) 
           elseif(D_L<0.d0) then
              alfa=umin(j,c1)-cell(c1)%qp(j)
              phi0=vk_limit(alfa,D_L,tol) 
           endif
  
           cell(c1)%phi(j)=dmin1(cell(c1)%phi(j),phi0)
      enddo    
   endif

   if(c2/=0) then 
      do j=1,nvar
           dist(:)=fc(i)%cen(:)-cell(c2)%cen(:)
           D_L=sum(cell(c2)%grad(:,j)*dist(:))
  
           TOL=kappa*(qmax(j)-qmin(j))
  
           phi0=1.0_dp    
           if(D_L>0.d0) then
              alfa=umax(j,c2)-cell(c2)%qp(j)
              phi0=vk_limit(alfa,D_L,tol) 
           elseif(D_L<0.d0) then
              alfa=umin(j,c2)-cell(c2)%qp(j)
              phi0=vk_limit(alfa,D_L,tol) 
           endif
  
           cell(c2)%phi(j)=dmin1(cell(c2)%phi(j),phi0)
      enddo
   endif
enddo

do i=1,noc
   do j=1,nvar 
      cell(i)%grad(1:ndim,j)=cell(i)%grad(1:ndim,j)*cell(i)%phi(j)
   enddo
enddo

!deallocate(dumin,dumax,phi)
contains

real(kind=dp) function vk_limit(a,b,tol)
implicit none
real(kind=dp)   :: nr,dr,tol,a,b

  nr=a*a+2.0_dp*a*b+tol
  dr=a*a+2.0_dp*b*b+a*b+tol
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
