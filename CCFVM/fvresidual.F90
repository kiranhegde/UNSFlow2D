!------------------------------------------------------------------------------
! Computes finite volume residual
!------------------------------------------------------------------------------
subroutine fvresidual
use commons
use grid 
use param
use inf
use pri
use output
implicit none
integer(kind=i4) :: i, j, ie, in, out,key
real(kind=dp)    :: qcl(nvar), qcr(nvar),flux(nvar),phi
real(kind=dp)    :: distl(ndim)!,distr(ndim)

phi=1.0_dp
do i=1,noc
   do j=1,nvar
      cell(i)%res(j) = 0.0d0
   enddo
enddo

key=0
Fx1=0.0_dp
Fy1=0.0_dp
! Compute flux for boundary edges
do ie=startBC,endBC
   qcl=0.d0
   qcr=0.d0
   distl=0.0_dp
   in = fc(ie)%in
   out = fc(ie)%out

   ! Farfield
   !if(fc(ie)%bc==2001.and.in>0) then  
   if(fc(ie)%bc==2001) then  
      if(key==0) then
         distl(:)=fc(ie)%cen(:)-cell(in)%cen(:)
         !distr(:)=fc(ie)%cen(:)-cell(out)%cen(:)
         do i=1,nvar
            phi=cell(in)%phi(i)
            qcl(i)=cell(in)%qp(i)+phi*sum(cell(in)%grad(:,i)*distl(:))
         enddo
         call farfield_flux(ie,qcl,qcr)
         if (flux_type == 'ausm')    call ausmPlus_flux(ie,qcl,qcr,flux)
         if (flux_type == 'roe')     call roe_flux(ie,qcl,qcr,flux)
         if (flux_type == 'vanleer') call vanleer_flux(ie,qcl,qcr,flux)
         if (flux_type == 'rusanov') call rusanov_flux(ie,qcl,qcr,flux)
         cell(in)%res(:)=cell(in)%res(:)+flux(:)
         call prim2con(qcr)  
         cell(out)%qc(:)=conv(:)
         call con2prim(conv)
         cell(out)%qp(:)=prim(:)
      else
         call prim2con(fs_inf)  
         cell(out)%qc(1:nvar)=conv(1:nvar) 
         call con2prim(conv)
         cell(out)%qp(:)=prim(:)
      endif
   endif

   ! Slip BC on wall (Euler)
   if(fc(ie)%bc==1001) then 
      if(key==1) then
         distl(:)=fc(ie)%cen(:)-cell(in)%cen(:)
         !distr(:)=fc(ie)%cen(:)-cell(out)%cen(:)
         do i=1,nvar
            phi=cell(in)%phi(i)
            qcl(i)=cell(in)%qp(i)+phi*sum(cell(in)%grad(:,i)*distl(:))
         enddo
         call slip_wall1(ie,qcl,qcr)
         if (flux_type == 'ausm')    call ausmPlus_flux(ie,qcl,qcr,flux)
         if (flux_type == 'roe')     call roe_flux(ie,qcl,qcr,flux)
         if (flux_type == 'vanleer') call vanleer_flux(ie,qcl,qcr,flux)
         if (flux_type == 'rusanov') call rusanov_flux(ie,qcl,qcr,flux)
         cell(in)%res(:)=cell(in)%res(:)+flux(:)
      else 
         distl(:)=fc(ie)%cen(:)-cell(in)%cen(:)
         !distr(:)=fc(ie)%cen(:)-cell(out)%cen(:)
         do i=1,nvar
            phi=cell(in)%phi(i)
            qcl(i)=cell(in)%qp(i)+phi*sum(cell(in)%grad(:,i)*distl(:))
         enddo

         ! Defining Ghost cell state 
         if (flow_type/="inviscid") then
             call NoSlip_wall(ie)   
         else
             call slip_wall(ie)
         endif
         qcr(1:nvar)=cell(out)%qp(1:nvar)   

         if (flux_type == 'ausm')    call ausmPlus_flux(ie,qcl,qcr,flux)
         if (flux_type == 'roe')     call roe_flux(ie,qcl,qcr,flux)
         if (flux_type == 'vanleer') call vanleer_flux(ie,qcl,qcr,flux)
         if (flux_type == 'rusanov') call rusanov_flux(ie,qcl,qcr,flux)
         cell(in)%res(:)=cell(in)%res(:)+flux(:)
      endif  
         Fx1=Fx1+flux(3) 
         Fy1=Fy1+flux(4) 
   endif


   !distl(:)=fc(ie)%cen(:)-cell(in)%cen(:)
   !distr(:)=fc(ie)%cen(:)-cell(out)%cen(:)

!   do i=1,nvar
!      phi=cell(in)%phi(i)
!      qcl(i)=cell(in)%qp(i)+phi*sum(cell(in)%grad(:,i)*distl(:))
!   enddo  
!
!   if(fc(ie)%bc==2001.and.in>0) then  
!      call farfield_flux(ie,qcl,qcr)
!      !qcr=fs_inf
!      if (flux_type == 'ausm')    call ausmPlus_flux(ie,qcl,qcr,flux)
!      if (flux_type == 'roe')     call roe_flux(ie,qcl,qcr,flux)
!      if (flux_type == 'vanleer') call vanleer_flux(ie,qcl,qcr,flux)
!      if (flux_type == 'rusanov') call rusanov_flux(ie,qcl,qcr,flux)
!      cell(in)%res(:)=cell(in)%res(:)+flux(:)
!   endif
!   ! Slip BC on wall (Euler)
!   if(fc(ie)%bc==1001.and.in>0) then 
!      call slip_wall(ie,qcl,qcr)
!      if (flux_type == 'ausm')    call ausmPlus_flux(ie,qcl,qcr,flux)
!      if (flux_type == 'roe')     call roe_flux(ie,qcl,qcr,flux)
!      if (flux_type == 'vanleer') call vanleer_flux(ie,qcl,qcr,flux)
!      if (flux_type == 'rusanov') call rusanov_flux(ie,qcl,qcr,flux)
!      cell(in)%res(:)=cell(in)%res(:)+flux(:)
!   endif
enddo


! Compute flux for interior edges
if (flux_type == 'vanleer') then
   do ie=startFC,endFC
      in = fc(ie)%in
      out = fc(ie)%out
      call recon(ie,in,out,qcl,qcr)
      call vanleer_flux(ie,qcl,qcr,flux)
      cell(in)%res(:)=cell(in)%res(:)+flux(:)
      cell(out)%res(:)=cell(out)%res(:)-flux(:)
   enddo
elseif (flux_type == 'roe') then
   do ie=startFC,endFC
   !do ie=1,nof
      in = fc(ie)%in
      out = fc(ie)%out
      call recon(ie,in,out,qcl,qcr)
      call roe_flux(ie,qcl,qcr,flux)
      cell(in)%res(:)=cell(in)%res(:)+flux(:)
      cell(out)%res(:)=cell(out)%res(:)-flux(:)
   enddo
elseif (flux_type == 'ausm') then
   do ie=startFC,endFC
      in = fc(ie)%in
      out = fc(ie)%out
      call recon(ie,in,out,qcl,qcr)
      call ausmPlus_flux(ie,qcl,qcr,flux)
      cell(in)%res(:)=cell(in)%res(:)+flux(:)
      cell(out)%res(:)=cell(out)%res(:)-flux(:)
   enddo
elseif (flux_type == 'rusanov') then
   do ie=startFC,endFC
      in = fc(ie)%in
      out = fc(ie)%out
      call recon(ie,in,out,qcl,qcr)
      call rusanov_flux(ie,qcl,qcr,flux)
      cell(in)%res(:)=cell(in)%res(:)+flux(:)
      cell(out)%res(:)=cell(out)%res(:)-flux(:)
   enddo
endif


if (flow_type=="laminar") then
   !do ie=startBC,endBC
   !   call NoSlip_wall(ie)   
   !enddo 
   call viscous_flux
elseif (flow_type=="turbulent") then
    call viscous_flux
    call rans
endif

! central implicit residual smoothening iterations
if (irs==yes) call cirs

end


subroutine recon(ie,in,out,qcl,qcr)
use commons
use pri
use grid
implicit none
integer(kind=i4) :: i,ie,in,out
real(kind=dp)    :: qcl(nvar), qcr(nvar)
real(kind=dp)    :: distl(ndim),distr(ndim) 
real(kind=dp)    :: zhi,beta,alfa,phi 

!zhi=0.0_dp       ! 
!zhi=-1.0_dp      ! 2nd order fully upwind
zhi=1.0_dp/3.0_dp ! 3rd order fully upwind
!zhi=2.0_dp/3.0_dp ! 3rd order fully upwind
phi=1.0_dp

qcl(:)=0.d0
qcr(:)=0.d0
distl=0.0_dp
distr=0.0_dp

distl(:)=fc(ie)%cen(:)-cell(in)%cen(:)
distr(:)=fc(ie)%cen(:)-cell(out)%cen(:)

!do i=1,nvar
   !Left state
!   qcl(i)=cell(in)%qp(i)+sum(cell(in)%grad(:,i)*distl(:))
   !Right state
!   qcr(i)=cell(out)%qp(i)+sum(cell(out)%grad(:,i)*distr(:))
!enddo

!return

!do i=1,nvar
!   !Left state
!   alfa=sum(cell(in)%grad(:,i)*distl(:))
!   beta=0.5_dp*(cell(out)%qp(i)-cell(in)%qp(i))
!   qcl(i)=cell(in)%qp(i)+zhi*beta+(1.0_dp-zhi)*alfa
!   !Right state
!   alfa=sum(cell(out)%grad(:,i)*distr(:))
!   beta=0.5_dp*(cell(in)%qp(i)-cell(out)%qp(i))
!   qcr(i)=cell(out)%qp(i)+zhi*beta+(1.0_dp-zhi)*alfa
!enddo


do i=1,nvar
   !Left state
   phi=cell(in)%phi(i)
   alfa=sum(cell(in)%grad(:,i)*distl(:))
   beta=cell(out)%qp(i)-cell(in)%qp(i)
   qcl(i)=cell(in)%qp(i)+0.5_dp*phi*(zhi*beta+(1.0_dp-zhi)*alfa)

   !Right state
   phi=cell(out)%phi(i)
   alfa=sum(cell(out)%grad(:,i)*distr(:))
   beta=cell(in)%qp(i)-cell(out)%qp(i)
   qcr(i)=cell(out)%qp(i)+0.5_dp*phi*(zhi*beta+(1.0_dp-zhi)*alfa)
enddo




end subroutine recon 
