!-----------------------------------------------------------------------------
! LUSGS implicit scheme - matrix free
!-----------------------------------------------------------------------------
subroutine lusgs
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: i,j,it,iv,p1,p2,c,itr,inner_itr
real(kind=dp):: omega,sign 
real(kind=dp):: D(nop),cres(nvar), flux1(nvar), flux2(nvar) 
real(kind=dp):: nx,ny,ds
type fvms
       real(kind=dp),dimension(nvar):: dqc
       real(kind=dp),dimension(nvar):: qflux,qold
end type fvms
type(fvms),allocatable,dimension(:)::fvm



allocate(fvm(nop))
! Over-relaxation factor: higher value improves stability but retards
! convergence. Needs to be tuned.
inner_itr=1
omega = 1.5d0

! Compute residual vector
call fvresidual




! Compute diagonal term
do it=1,nop
   D(it)=pt(it)%cv/pt(it)%dt+0.5d0*omega*pt(it)%la
   fvm(it)%dqc(:)=0.d0
   fvm(it)%qold(:)=pt(it)%qold(:)
enddo

!lower(:,:)=0.d0
!upper(:,:)=0.d0

do itr=1,inner_itr



! Forward loop
do it=1,nop
   cres(:) = 0.0d0

   do j=1,pt(it)%nv2f

      c=pt(it)%v2f(j) 
      nx = fc(c)%nx
      ny = fc(c)%ny
      ds = fc(c)%area
      p1 = fc(c)%pt(1)
      p2 = fc(c)%pt(2)
      sign=1.d0
      if(p2==it) then  
         p2=p1  
         sign=-1.d0 
      endif

      if(p2 .lt. it .and. p2 .gt. 0)then
         flux1(:)=pt(p2)%qold(:)
         call normalflux(flux1)

         flux2(:)=pt(p2)%qc(:)
         call normalflux(flux2)

         do iv=1,nvar
            cres(iv)=cres(iv)+ 0.5d0*( sign*(flux2(iv)-flux1(iv))-omega*fc(c)%la*fvm(p2)%dqc(iv))
         enddo
      endif

   enddo

   do iv=1,nvar
      fvm(it)%dqc(iv)=(-pt(it)%res(iv)-cres(iv))/D(it)
      pt(it)%qc(iv) = pt(it)%qold(iv) + fvm(it)%dqc(iv)
   enddo

enddo

! Reverse loop
do it=nop,1,-1
   cres(:) = 0.0d0

   do j=1,pt(it)%nv2f

      c=pt(it)%v2f(j) 
      nx = fc(c)%nx
      ny = fc(c)%ny
      ds = fc(c)%area
      p1 = fc(c)%pt(1)
      p2 = fc(c)%pt(2)
      sign=1.d0
      if(p2==it) then  
         p2=p1  
         sign=-1.d0
      endif

      if(p2 .gt. it .and. p2 .gt. 0)then
         flux1(:)=pt(p2)%qold(:)
         call normalflux(flux1)

         flux2(:)=pt(p2)%qc(:)
         call normalflux(flux2)

         do iv=1,nvar
            cres(iv)=cres(iv)+ 0.5d0*( sign*(flux2(iv)-flux1(iv))-omega*fc(c)%la*fvm(p2)%dqc(iv))
         enddo
      endif

   enddo

   do iv=1,nvar
      fvm(it)%dqc(iv)=fvm(it)%dqc(iv)-cres(iv)/D(it)
      pt(it)%qc(iv) = pt(it)%qold(iv) + fvm(it)%dqc(iv)
   enddo

enddo


enddo


deallocate(fvm)

contains

!-----------------------------------------------------------------------------
! Computes flux along (nx,ny)
!-----------------------------------------------------------------------------
subroutine normalflux(flux)
use param
use pri
implicit none
real(kind=dp) :: flux(nvar)
real(kind=dp) :: un,et

call con2prim(flux)
!e       = p/GAMMA1 + 0.5d0*rho*(u*u + v*v)

et=e*rho
un      = u*nx + v*ny
flux(1) = (Et+ p)*un
flux(2) = rho*un
flux(3) = flux(2)*u + p*nx
flux(4) = flux(2)*v + p*ny

flux(:)=flux(:)*ds
end subroutine normalflux

end
