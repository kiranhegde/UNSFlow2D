!-----------------------------------------------------------------------------
! LUSGS implicit scheme - matrix free
!-----------------------------------------------------------------------------

subroutine lusgs

!call lusgs1
call lusgs2

end subroutine lusgs
!-----------------------------------------------------------------------------
subroutine lusgs1
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: j,it,iv,in,out,c,itr,inner_itr
real(kind=dp):: omega,rA,lamA,viscA,rAdiag 
real(kind=dp):: cres(nvar),diag
real(kind=dp):: nx,ny,area,dflux(nvar)
real(kind=dp)   :: dist(ndim),cdist
type fvms
!       real(kind=dp)                :: diag
       real(kind=dp),dimension(nvar):: lower,upper
end type fvms
type(fvms),allocatable,dimension(:)::fvm


allocate(fvm(noc))
! Over-relaxation factor: higher value improves stability but retards
! convergence. Needs to be tuned.
! omega=1.5
inner_itr=1
omega = 1.5d0

! Compute residual vector
call fvresidual

! Compute diagonal term
!do it=1,noc
!   fvm(it)%diag=cell(it)%cv/cell(it)%dt+0.5*omega*cell(it)%la
!   fvm(it)%lower(:)=0.d0
!   fvm(it)%upper(:)=0.d0
!enddo


do itr=1,inner_itr

! Forward loop
do it=1,noc
   fvm(it)%lower(:)=0.d0
   cres(:) = 0.0d0
   rAdiag = 0.0_dp

   do j=1,cell(it)%nc2f
      c=cell(it)%c2f(j) 
      nx = fc(c)%nx
      ny = fc(c)%ny
      area=fc(c)%area
      in = fc(c)%in
      out = fc(c)%out
      dist(:)=cell(out)%cen(:)-cell(in)%cen(:)
      cdist=dsqrt(dist(1)*dist(1)+dist(2)*dist(2)) 
      if(out==it) then  
         out=in  
         nx =-fc(c)%nx
         ny =-fc(c)%ny
      endif

      if(out .lt. it .and. out .gt. 0.and.out<=noc)then
         call normal_flux(dflux,fvm(out)%lower(1:nvar),out)
         call maxeig(cell(out)%qc(1:nvar),lamA,viscA,cdist,out)
         rA = omega*lamA+viscA
         rAdiag = rAdiag + 0.5_dp*omega*lamA+viscA
         do iv=1,nvar
            cres(iv)=cres(iv)+(dflux(iv)-rA*fvm(out)%lower(iv))
         enddo
      endif
   enddo
   
   diag=cell(it)%cv/cell(it)%dt+rAdiag
   do iv=1,nvar
      fvm(it)%lower(iv)=(-cell(it)%res(iv)-0.5d0*cres(iv))/diag
   enddo
enddo


! Reverse loop
do it=noc,1,-1
   cres(:) = 0.0_dp
   rAdiag = 0.0_dp
   fvm(it)%upper(:)=0.d0

   do j=1,cell(it)%nc2f
      c=cell(it)%c2f(j) 
      nx = fc(c)%nx
      ny = fc(c)%ny
      area=fc(c)%area
      in = fc(c)%in
      out = fc(c)%out
      dist(:)=cell(out)%cen(:)-cell(in)%cen(:)
      cdist=dsqrt(dist(1)*dist(1)+dist(2)*dist(2)) 
      if(out==it) then  
         out=in  
         nx =-fc(c)%nx
         ny =-fc(c)%ny
      endif

      if(out .gt. it .and. out .gt. 0.and.out<=noc)then
         call normal_flux(dflux,fvm(out)%upper(1:nvar),out)
         call maxeig(cell(out)%qc(1:nvar),lamA,viscA,cdist,out)
         rA = omega*lamA+viscA
         rAdiag = rAdiag + 0.5_dp*omega*lamA+viscA
         do iv=1,nvar
            cres(iv)=cres(iv)+(dflux(iv)-rA*fvm(out)%upper(iv))
         enddo
      endif
   enddo
   diag=cell(it)%cv/cell(it)%dt+rAdiag
   do iv=1,nvar
      fvm(it)%upper(iv)=fvm(it)%lower(iv)-0.5d0*cres(iv)/diag
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
! Computes flux along (nx,ny)
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
et=e*rho
un       = u*nx + v*ny
flux1(1) = (Et+ p)*un
flux1(2) = rho*un
flux1(3) = flux1(2)*u + p*nx
flux1(4) = flux1(2)*v + p*ny


flux(:)=cell(out)%qold(:)+dq(:)
call con2prim(flux)
et=e*rho
un       = u*nx + v*ny
flux2(1) = (Et+ p)*un
flux2(2) = rho*un
flux2(3) = flux2(2)*u + p*nx
flux2(4) = flux2(2)*v + p*ny

flux(:)=area*(flux2(:)-flux1(:))

end subroutine normal_flux

!-----------------------------------------------------------------------------
! Computes maximum eigenvalue normal to a face with normal (nx, ny)
!-----------------------------------------------------------------------------
subroutine maxeig(qc,lamA,viscA,cdist,out)
use commons
use pri
use grid
use visc
implicit none
integer(kind=i4):: out
real(kind=dp)   :: qc(nvar),lamA,viscA,cdist
real(kind=dp)   :: factor,ratio,mu,kv

kv=0.25_dp
factor=4.0_dp/3.0_dp
mu  = cell(out)%mul

call con2prim(qc)
a = dsqrt(GAMMA*p/rho)
! inviscid spectral radii
lamA  =(dabs(u*nx + v*ny) + a)*area

! viscous spectral radii
ratio = area*mu/(cdist*prandtl*rho*Rey) 
viscA = ratio*dmax1(factor,gamma) 

end
              

end subroutine lusgs1

!------------------------------------------------------------------------------
subroutine lusgs2
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: j,it,iv,in,out,c,itr,inner_itr
real(kind=dp):: omega,lam 
real(kind=dp):: cres(nvar)
real(kind=dp):: nx,ny,area,dflux(nvar)
type fvms
       real(kind=dp)                :: diag
       real(kind=dp),dimension(nvar):: lower,upper
end type fvms
type(fvms),allocatable,dimension(:)::fvm


allocate(fvm(noc))
! Over-relaxation factor: higher value improves stability but retards
! convergence. Needs to be tuned.
! omega=1.5
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
      nx = fc(c)%nx
      ny = fc(c)%ny
      area=fc(c)%area
      in = fc(c)%in
      out = fc(c)%out
      if(out==it) then  
         out=in  
         nx =-fc(c)%nx
         ny =-fc(c)%ny
      endif

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
      nx = fc(c)%nx
      ny = fc(c)%ny
      area=fc(c)%area
      in = fc(c)%in
      out = fc(c)%out
      if(out==it) then  
         out=in  
         nx =-fc(c)%nx
         ny =-fc(c)%ny
      endif

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
! Computes flux along (nx,ny)
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
et=e*rho
un       = u*nx + v*ny
flux1(1) = (Et+ p)*un
flux1(2) = rho*un
flux1(3) = flux1(2)*u + p*nx
flux1(4) = flux1(2)*v + p*ny


flux(:)=cell(out)%qold(:)+dq(:)
call con2prim(flux)
et=e*rho
un       = u*nx + v*ny
flux2(1) = (Et+ p)*un
flux2(2) = rho*un
flux2(3) = flux2(2)*u + p*nx
flux2(4) = flux2(2)*v + p*ny

flux(:)=area*(flux2(:)-flux1(:))

end subroutine normal_flux

!-----------------------------------------------------------------------------
! Computes maximum eigenvalue normal to a face with normal (nx, ny)
! and ds = dsqrt(nx*nx + ny*ny) is face length
!-----------------------------------------------------------------------------
!subroutine maxeig(qc,lam)
!use commons
!use pri
!use grid
!implicit none
!real(dp) :: qc(nvar), lam

!call con2prim(qc)
!a   = dsqrt(GAMMA*p/rho)
!lam = dabs(u*nx + v*ny) + a*area

!end
              

end subroutine lusgs2


