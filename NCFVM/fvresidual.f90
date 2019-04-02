!------------------------------------------------------------------------------
! Computes finite volume residual
!------------------------------------------------------------------------------
subroutine fvresidual
use grid 
implicit none
integer(kind=i4)  :: i, j, ie, p1, p2,recon
real(kind=dp) :: nx, ny,x,y,left(nvar),right(nvar)
real(kind=dp) :: dist(ndim),grad(nvar)

<<<<<<< HEAD
call Gradient_LSQR
!call Gradient_GG
=======
!call Gradient_LSQR
call Gradient_GG
>>>>>>> a688afb70276d1c80f99d1548c4045d94d229af4

recon=1

do i=1,nop
   do j=1,nvar
      pt(i)%res(j) = 0.0d0
   enddo
enddo
pt(:)%la=0.d0
fc(:)%la=0.d0

!     Compute flux for interior edges

if(Recon==1) then
call LinearReconLimiter

do ie=1,nof
   p1 = fc(ie)%pt(1)
   p2 = fc(ie)%pt(2)
   dist(1)=pt(p2)%x-pt(p1)%x
   dist(2)=pt(p2)%y-pt(p1)%y

!  Left state
   grad(:)=0.d0
   do i=1,ndim
      grad(:)=grad(:)+pt(p1)%grad(i,:)*dist(i)
   enddo

   left(:)=pt(p1)%qp(:)+0.5d0*grad(:)*pt(p1)%phi(:)

!  Right state
   grad(:)=0.d0
   do i=1,ndim
      grad(:)=grad(:)+pt(p2)%grad(i,:)*dist(i)
   enddo

   right(:)=pt(p2)%qp(:)-0.5d0*grad(:)*pt(p2)%phi(:)


   call vanleer_flux(ie,p1,p2,left,right)
   !if(fc(ie)%bc==0) call vanleer_flux(ie,p1,p2)
   !if(p1 /=0 .and.p2 /=0) call vanleer_flux(ie,p2,p1)
   !if(p1 /=0 .and.p2 /=0) call roe_flux(ie,p1,p2)
   !if(p1 /=0 .and.p2 /=0) call ausmPlus_flux(ie,p1,p2)
   if(fc(ie)%bc==1001)  call solid_flux(ie,p1,p2)
   if(fc(ie)%bc==2001)  call farfield_flux1(ie,p1,p2)
enddo

else
do ie=1,nof
   p1 = fc(ie)%pt(1)
   p2 = fc(ie)%pt(2)

   call muscle_va(ie,left,right)
   call vanleer_flux(ie,p1,p2,left,right)
   !if(fc(ie)%bc==0) call vanleer_flux(ie,p1,p2)
   !if(p1 /=0 .and.p2 /=0) call vanleer_flux(ie,p2,p1)
   !if(p1 /=0 .and.p2 /=0) call roe_flux(ie,p1,p2)
   !if(p1 /=0 .and.p2 /=0) call ausmPlus_flux(ie,p1,p2)
   if(fc(ie)%bc==1001)  call solid_flux(ie,p1,p2)
   if(fc(ie)%bc==2001)  call farfield_flux1(ie,p1,p2)
enddo
endif

call time_step2
!call time_step02

end
