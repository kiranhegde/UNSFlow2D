!-----------------------------------------------------------------------------
! LUSGS implicit scheme - matrix free
!-----------------------------------------------------------------------------
subroutine lusgs
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: j,it,iv,in,out,c,itr,inner_itr
real(kind=dp):: omega,lam 
real(kind=dp):: cres(nvar)
real(kind=dp):: sx,sy,ds,dflux(nvar)
type fvms
       real(kind=dp)                :: diag
       real(kind=dp),dimension(nvar):: lower,upper
end type fvms
type(fvms),allocatable,dimension(:)::fvm


allocate(fvm(noc))
! Over-relaxation factor: higher value improves stability but retards
! convergence. Needs to be tuned.
inner_itr=1
omega = 1.5d0

! Compute residual vector
call fvresidual

! Compute diagonal term
do it=1,noc
   fvm(it)%diag=cell(it)%cv/cell(it)%dt+0.5*omega*cell(it)%la
   fvm(it)%lower(:)=0.d0
   fvm(it)%upper(:)=0.d0
enddo


do itr=1,inner_itr

! Forward loop
do it=1,noc
   fvm(it)%lower(:)=0.d0
   cres(:) = 0.0d0
   do j=1,cell(it)%nc2f
      c=cell(it)%c2f(j) 
      sx = fc(c)%sx
      sy = fc(c)%sy
      ds = dsqrt(sx*sx + sy*sy)
      in = fc(c)%in
      out = fc(c)%out
      if(out==it) then  
         out=in  
         sx =-fc(c)%sx
         sy =-fc(c)%sy
      endif
      sx=sx/ds
      sy=sy/ds

      if(out .lt. it .and. out .gt. 0.and.out<=noc)then
         call normal_flux(dflux,fvm(out)%lower(1:nvar),out)
         !call maxeig(cell(out)%qc(1:nvar), lam)
         lam=fc(c)%la     
         !lam=cell(out)%la
         do iv=1,nvar
            cres(iv)=cres(iv)+(dflux(iv) &
                             -omega*lam*fvm(out)%lower(iv))
         enddo
      endif
   enddo

   do iv=1,nvar
      fvm(it)%lower(iv)=(-cell(it)%res(iv)-0.5d0*cres(iv))/fvm(it)%diag
   enddo
enddo

! Reverse loop
do it=noc,1,-1
   cres(:) = 0.0d0
   fvm(it)%upper(:)=0.d0
   do j=1,cell(it)%nc2f
      c=cell(it)%c2f(j) 
      sx = fc(c)%sx
      sy = fc(c)%sy
      ds = dsqrt(sx*sx + sy*sy)
      in = fc(c)%in
      out = fc(c)%out
      if(out==it) then  
         out=in  
         sx =-fc(c)%sx
         sy =-fc(c)%sy
      endif
      sx=sx/ds
      sy=sy/ds

      if(out .gt. it .and. out .gt. 0.and.out<=noc)then
         call normal_flux(dflux,fvm(out)%upper(1:nvar),out)
         !call maxeig(cell(out)%qc(1:nvar), lam)
         lam=fc(c)%la     
         !lam=cell(out)%la
         do iv=1,nvar
            cres(iv)=cres(iv)+(dflux(iv) &
                             -omega*lam*fvm(out)%upper(iv))
         enddo
      endif
   enddo

   do iv=1,nvar
      fvm(it)%upper(iv)=fvm(it)%lower(iv)-0.5d0*cres(iv)/fvm(it)%diag
   enddo
enddo

enddo

do it=1,noc
   do iv=1,nvar
      cell(it)%qc(iv) = cell(it)%qold(iv) + fvm(it)%upper(iv) 
   enddo
enddo


deallocate(fvm)

contains

!-----------------------------------------------------------------------------
! Computes flux along (sx,sy)
!-----------------------------------------------------------------------------
subroutine normal_flux(flux,dq,out)
use param
use pri
implicit none
integer(kind=i4):: out
real(kind=dp) :: flux(nvar),dq(nvar)
real(kind=dp) :: flux1(nvar)
real(kind=dp) :: flux2(nvar)
real(kind=dp) :: un,et

flux(:)=cell(out)%qold(:)
call con2prim(flux)
!e       = p/GAMMA1 + 0.5d0*rho*(u*u + v*v)

et=e*rho
un       = u*sx + v*sy
flux1(1) = (Et+ p)*un
flux1(2) = rho*un
flux1(3) = flux1(2)*u + p*sx
flux1(4) = flux1(2)*v + p*sy


flux(:)=cell(out)%qold(:)+dq(:)
call con2prim(flux)

et=e*rho
un       = u*sx + v*sy
flux2(1) = (Et+ p)*un
flux2(2) = rho*un
flux2(3) = flux2(2)*u + p*sx
flux2(4) = flux2(2)*v + p*sy

flux(:)=ds*(flux2(:)-flux1(:))

end subroutine normal_flux

!-----------------------------------------------------------------------------
! Computes maximum eigenvalue normal to a face with normal (sx, sy)
! and ds = dsqrt(sx*sx + sy*sy) is face length
!-----------------------------------------------------------------------------
!subroutine maxeig(qc,lam)
!use commons
!use pri
!use grid
!implicit none
!real(dp) :: qc(nvar), lam

!call con2prim(qc)
!a   = dsqrt(GAMMA*p/rho)
!lam = dabs(u*sx + v*sy) + a*ds

!end
              

end
