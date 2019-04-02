!------------------------------------------------------------------------------
! Computes finite volume residual
!------------------------------------------------------------------------------
subroutine fvresidual
use commons
use grid 
use param
implicit none
integer(kind=i4) :: i, j, ie, c1, c2
real(kind=dp)    :: qcl(nvar), qcr(nvar),flux(nvar)

do i=1,noc
   do j=1,nvar
      cell(i)%res(j) = 0.0d0
   enddo
enddo

!     Compute flux for interior edges
if (flux_type == 'vanleer') then
   do ie=1,nof
      c1 = fc(ie)%in
      c2 = fc(ie)%out
      if(c1/=0.and.c2/=0) then 
         call recon(ie,c1,c2,qcl,qcr)
         call vanleer_flux(ie,c1,c2,qcl,qcr)
      if(fc(ie)%bc==1001.and.c1/=0) then 
         call slip_wall(ie,c1,qcl,qcr)
         call roe_flux(ie,qcl,qcr,flux)
         cell(c1)%res(:)=cell(c1)%res(:)+flux(:)
      endif
      if(fc(ie)%bc==1001.and.c2/=0) then 
         call slip_wall(ie,c2,qcl,qcr)
         call roe_flux(ie,qcl,qcr,flux)
         cell(c2)%res(:)=cell(c2)%res(:)+flux(:)
      endif
      endif
   enddo
elseif (flux_type == 'roe') then
   do ie=1,nof
      c1 = fc(ie)%in
      c2 = fc(ie)%out
      if(c1/=0.and.c2/=0) then 
         call recon(ie,c1,c2,qcl,qcr)
         call roe_flux(ie,qcl,qcr,flux)
         cell(c1)%res(:)=cell(c1)%res(:)+flux(:)
         cell(c2)%res(:)=cell(c2)%res(:)-flux(:)
      endif
      if(fc(ie)%bc==1001.and.c1/=0) then 
         call slip_wall(ie,c1,qcl,qcr)
         call roe_flux(ie,qcl,qcr,flux)
         cell(c1)%res(:)=cell(c1)%res(:)+flux(:)
      endif
      if(fc(ie)%bc==1001.and.c2/=0) then 
         call slip_wall(ie,c2,qcl,qcr)
         call roe_flux(ie,qcl,qcr,flux)
         cell(c2)%res(:)=cell(c2)%res(:)+flux(:)
      endif

!   if(fc(ie)%bc==1001.and.c1/=0)  call solid_flux(ie,c1)
!   if(fc(ie)%bc==1001.and.c2/=0)  call solid_flux(ie,c2)
   enddo
elseif (flux_type == 'ausm') then
   do ie=1,nof
      c1 = fc(ie)%in
      c2 = fc(ie)%out
      if(c1/=0.and.c2/=0) then 
         call recon(ie,c1,c2,qcl,qcr)
         call ausmPlus_flux(ie,c1,c2,qcl,qcr)
      endif 
      if(fc(ie)%bc==1001.and.c1/=0) then 
         call slip_wall(ie,c1,qcl,qcr)
         call roe_flux(ie,qcl,qcr,flux)
         cell(c1)%res(:)=cell(c1)%res(:)+flux(:)
      endif
      if(fc(ie)%bc==1001.and.c2/=0) then 
         call slip_wall(ie,c2,qcl,qcr)
         call roe_flux(ie,qcl,qcr,flux)
         cell(c2)%res(:)=cell(c2)%res(:)+flux(:)
      endif
   enddo
endif

!     Compute flux for solid wall edges

!do ie=1,nof
!   c1 = fc(ie)%in
!   c2 = fc(ie)%out
!   if(fc(ie)%bc==1001.and.c1/=0)  call solid_flux(ie,c1)
!   if(fc(ie)%bc==1001.and.c2/=0)  call solid_flux(ie,c2)
!enddo

!     Flux for far-field points

do ie=1,nof
   c1 = fc(ie)%in
   c2 = fc(ie)%out
   if(fc(ie)%bc==2001.and.c1/=0)  call farfield_flux(ie,c1)
   if(fc(ie)%bc==2001.and.c2/=0)  call farfield_flux(ie,c2)
enddo



return
if (flow_type=="laminar") then
    call viscous_flux
elseif (flow_type=="turbulent") then
    call viscous_flux
    call rans
else
   print*,'Unknown flow type:',flow_type
   stop
endif

end


subroutine recon(ie,c1,c2,qcl,qcr)
use commons
use pri
use grid
implicit none
integer(kind=i4) :: i,ie,c1,c2
real(kind=dp)    :: qcl(nvar), qcr(nvar)
real(kind=dp)    :: distl(ndim),distr(ndim) 


qcl(:)=0.d0
qcr(:)=0.d0
distl=0.0_dp
distr=0.0_dp

distl(:)=fc(ie)%cen(:)-cell(c1)%cen(:)
distr(:)=fc(ie)%cen(:)-cell(c2)%cen(:)

do i=1,nvar
   !Left state
   qcl(i)=cell(c1)%qp(i)+sum(cell(c1)%grad(:,i)*distl(:))
   !Right state
   qcr(i)=cell(c2)%qp(i)+sum(cell(c2)%grad(:,i)*distr(:))
enddo

end subroutine recon 
