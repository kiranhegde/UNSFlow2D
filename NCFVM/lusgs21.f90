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
real(kind=dp):: omega,lam,qq,aa 
real(kind=dp):: cres(nvar)
real(kind=dp):: nx,ny,area,dflux(nvar)
type fvms
       real(kind=dp)                :: diag
       real(kind=dp),dimension(nvar):: lower,upper
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
   lam=pt(it)%la
   fvm(it)%diag=pt(it)%cv/pt(it)%dt+0.5*omega*lam
   fvm(it)%lower(:)=0.d0
   fvm(it)%upper(:)=0.d0
enddo


do itr=1,inner_itr

! Forward loop
do it=1,nop
   cres(:) = 0.0d0
   do j=1,pt(it)%nv2f
      c=pt(it)%v2f(j) 
      nx = fc(c)%nx
      ny = fc(c)%ny
      p1 = fc(c)%pt(1)
      p2 = fc(c)%pt(2)
      if(p2==it) then  
         p2=p1  
         nx =-nx
         ny =-ny
      endif
      area = fc(c)%area
    
      if(p2 .lt. it)then
         call normal_flux(dflux,fvm(p2)%lower(1:nvar),p2)
         !call maxeig(pt(p2)%qc(1:nvar), lam)
         lam=fc(p2)%la
         do iv=1,nvar
            cres(iv)=cres(iv)+(dflux(iv) &
                             -omega*lam*fvm(p2)%lower(iv))
         enddo
      endif
   enddo

   do iv=1,nvar
      fvm(it)%lower(iv)=(-pt(it)%res(iv)-0.5d0*cres(iv))/fvm(it)%diag
   enddo
enddo


! Reverse loop
do it=nop,1,-1
   cres(:) = 0.0d0
   do j=1,pt(it)%nv2f
      c=pt(it)%v2f(j) 
      nx = fc(c)%nx
      ny = fc(c)%ny
      p1 = fc(c)%pt(1)
      p2 = fc(c)%pt(2)
      if(p2==it) then  
         p2=p1  
         nx =-nx
         ny =-ny
      endif
      area = fc(c)%area

      if(p2 .gt. it)then
         call normal_flux(dflux,fvm(p2)%upper(1:nvar),p2)
         !call maxeig(pt(p2)%qc(1:nvar), lam)
         lam=fc(p2)%la
         do iv=1,nvar
            cres(iv)=cres(iv)+(dflux(iv) &
                             -omega*lam*fvm(p2)%upper(iv))
         enddo
      endif
   enddo

   do iv=1,nvar
      fvm(it)%upper(iv)=fvm(it)%lower(iv)-0.5d0*cres(iv)/fvm(it)%diag
   enddo
enddo

enddo

do it=1,nop
   do iv=1,nvar
      pt(it)%qc(iv) = pt(it)%qold(iv) + fvm(it)%upper(iv) 
   enddo
enddo


deallocate(fvm)

contains

!-----------------------------------------------------------------------------
! Computes flux along (nx,ny)
!-----------------------------------------------------------------------------
subroutine normal_flux(flux,dq,p2)
use param
use pri
implicit none
integer(kind=i4):: p2
real(kind=dp) :: flux(nvar),dq(nvar)
real(kind=dp) :: flux1(nvar)
real(kind=dp) :: flux2(nvar)
real(kind=dp) :: un,et,nx,ny

flux(:)=pt(p2)%qold(:)
call con2prim(flux)
!e       = p/GAMMA1 + 0.5d0*rho*(u*u + v*v)

et=e*rho
un       = u*nx + v*ny
flux1(1) = (Et+ p)*un
flux1(2) = rho*un
flux1(3) = flux1(2)*u + p*nx
flux1(4) = flux1(2)*v + p*ny


flux(:)=pt(p2)%qold(:)+dq(:)
call con2prim(flux)

et=e*rho
un       = u*nx + v*ny
flux2(1) = (Et+ p)*un
flux2(2) = rho*un
flux2(3) = flux2(2)*u + p*nx
flux2(4) = flux2(2)*v + p*ny

flux(:)=(flux2(:)-flux1(:))*area

end subroutine normal_flux

!-----------------------------------------------------------------------------
! Computes maximum eigenvalue normal to a face with normal (nx, ny)
! and ds = dsqrt(nx*nx + ny*ny) is face length
!-----------------------------------------------------------------------------
subroutine maxeig(qc,lam)
!use commons
!use pri
!use grid
implicit none
real(dp) :: qc(nvar), lam

call con2prim(qc)
!a   = dsqrt(GAMMA*p/rho)
lam = (dabs(u*nx + v*ny) + a)*area

end subroutine maxeig 
              

end subroutine lusgs
