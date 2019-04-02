!-----------------------------------------------------------------------------
! Flux for a farfield edge: Steger-Wardmin1g flux splitting
!-----------------------------------------------------------------------------
subroutine farfield_flux1(ie,cc)
use grid
use pri
use inf
implicit none
integer(kind=i4) :: i,ie,cc
real(kind=dp) :: con(nvar)
real(kind=dp) :: dr, nx, ny, q2,un, l2, l3
real(kind=dp) :: uinf, vinf, pinf, rinf, ainf, &
            q2inf, l1p, l2p, l3p, l1n, l2n, l3n, fp(4),  &
            fn(4), f1, f2, a1, a2

uinf  = u_inf
vinf  = v_inf
q2inf = uinf**2 + vinf**2
pinf  = p_inf
rinf  = r_inf
ainf  = a_inf

nx = fc(ie)%sx
ny = fc(ie)%sy
dr =  dsqrt(nx*nx + ny*ny)
nx = nx/dr
ny = ny/dr

! Positive flux
!con(:)=cell(cc)%qc(:)
!call con2prim(con)

!dx=fc(ie)%ldx
!dy=fc(ie)%ldy
!if(fc(ie)%ldx==0.d0.and.  fc(ie)%ldy==0.d0) then
!dx=fc(ie)%rdx
!dy=fc(ie)%rdy
!endif
con(:)=cell(cc)%qp(:)!+(cell(cc)%qx(:)*dx+cell(cc)%qy(:)*dy)


rho=con(1)
u  =con(2)
v  =con(3)
p  =con(4)

q2= u*u + v*v 
a= dsqrt(GAMMA*p/rho)

un = u*nx + v*ny
l2 = un + a
l3 = un - a
l1p= dmax1(un, 0.0d0)
l2p= dmax1(l2, 0.0d0)
l3p= dmax1(l3, 0.0d0)
a1 = 2.0d0*(GAMMA-1.0d0)*l1p + l2p + l3p
f1 = 0.5d0*rho/GAMMA

fp(2) = f1*a1
fp(3) = f1*( a1*u + a*(l2p - l3p)*nx )
fp(4) = f1*( a1*v + a*(l2p - l3p)*ny )
fp(1) = f1*( 0.5d0*a1*q2 + a*un*(l2p - l3p) + &
             a**2*(l2p+l3p)/GAMMA1 )

! Negative flux
un = uinf*nx + vinf*ny
l2 = un + ainf
l3 = un - ainf
l1n= dmin1(un, 0.0d0)
l2n= dmin1(l2, 0.0d0)
l3n= dmin1(l3, 0.0d0)
a2 = 2.0d0*(GAMMA-1.0d0)*l1n + l2n + l3n
f2 = 0.5d0*rinf/GAMMA

fn(2) = f2*a2
fn(3) = f2*( a2*uinf + ainf*(l2n - l3n)*nx )
fn(4) = f2*( a2*vinf + ainf*(l2n - l3n)*ny )
fn(1) = f2*( 0.5d0*a2*q2inf + ainf*un*(l2n - l3n) + &
             ainf**2*(l2n+l3n)/GAMMA1 )

do i=1,nvar
cell(cc)%res(i)=cell(cc)%res(i)+dr*( fp(i) + fn(i) )
enddo


end

subroutine farfield_flux(ie,cc)
use grid
use pri
use inf
implicit none
integer(kind=i4)  :: i,ie,cc
real(kind=dp) :: dr, nx, ny, q2,un,con(nvar)
real(kind=dp) :: uinf, vinf, pinf, rinf, ainf 
real(kind=dp) :: Rp,Rm,rhof,uf,vf,pf,af,Unf,Sf,Ef 
real(kind=dp) :: q2inf,un_inf,flux(nvar) 
            

nx = fc(ie)%sx
ny = fc(ie)%sy
dr =  dsqrt(nx*nx + ny*ny)
nx = nx/dr
ny = ny/dr

uinf  = u_inf
vinf  = v_inf
q2inf = uinf**2 + vinf**2
pinf  = p_inf
rinf  = r_inf
ainf  = a_inf
un_inf = uinf*nx + vinf*ny

! Positive flux
!con(:)=cell(cc)%qc(:)
!call con2prim(con)

!dx=fc(ie)%ldx
!dy=fc(ie)%ldy
!if(fc(ie)%ldx==0.d0.and.  fc(ie)%ldy==0.d0) then
!dx=fc(ie)%rdx
!dy=fc(ie)%rdy
!endif
con(:)=cell(cc)%qp(:)
!do i=1,ndim
!dist=fc(ie)%cen(i)-cell(cc)%cen(i)
!con(:)=con(:)+cell(cc)%grad(i,:)*dist
!enddo

!con(:)=cell(cc)%qp(:)+(cell(cc)%qx(:)*dx+cell(cc)%qy(:)*dy)


rho=con(1)
u  =con(2)
v  =con(3)
p  =con(4)

q2= u*u + v*v 
a= dsqrt(GAMMA*p/rho)

un = u*nx + v*ny

Rp=un+2.d0*a/gamma1
Rm=un_inf-2.d0*ainf/gamma1


Unf=0.5d0*(Rp+Rm)
af=0.25d0*(Rp-Rm)*gamma1

if(Unf>0.d0) then
  uf=u+nx*(Unf-un)
  vf=v+ny*(Unf-un)
  Sf=p/(rho**gamma)
else
  uf=uinf+nx*(Unf-un_inf)
  vf=vinf+ny*(Unf-un_inf)
  Sf=pinf/(rinf**gamma)
endif

rhof=(af*af/gamma/Sf)**(1.d0/gamma1)
pf=rhof*af*af/gamma
ef=pf/gamma1 + 0.5d0*rhof*(uf*uf+vf*vf)

flux(1) = (ef + pf)*unf
flux(2) = rhof*unf
flux(3) = flux(2)*uf + pf*nx
flux(4) = flux(2)*vf + pf*ny


do i=1,nvar
cell(cc)%res(i)=cell(cc)%res(i)+dr*flux(i)
enddo


end


subroutine farfield_flux0(ie,cc)
use grid
use pri
use inf
implicit none
integer(kind=i4)  :: i,ie,cc
real(kind=dp) :: dr, nx, ny, q2,un,con(nvar)
real(kind=dp) :: uinf, vinf, pinf, rinf, ainf 
real(kind=dp) :: Rp,Rm,rhof,uf,vf,pf,af,Unf,Sf,Ef 
real(kind=dp) :: q2inf,un_inf,flux(nvar),la0(nvar),la(nvar) 
            

nx = fc(ie)%sx
ny = fc(ie)%sy
dr =  dsqrt(nx*nx + ny*ny)
nx = nx/dr
ny = ny/dr

uinf  = u_inf
vinf  = v_inf
q2inf = uinf**2 + vinf**2
pinf  = p_inf
rinf  = r_inf
ainf  = a_inf
un_inf = uinf*nx + vinf*ny

la(1)=un_inf-a_inf
la(2)=un_inf
la(3)=un_inf
la(4)=un_inf+a_inf


! Positive flux
!con(:)=cell(cc)%qc(:)
!call con2prim(con)

!dx=fc(ie)%ldx
!dy=fc(ie)%ldy
!if(fc(ie)%ldx==0.d0.and.  fc(ie)%ldy==0.d0) then
!dx=fc(ie)%rdx
!dy=fc(ie)%rdy
!endif
con(:)=cell(cc)%qp(:)!+(cell(cc)%qx(:)*dx+cell(cc)%qy(:)*dy)


rho=con(1)
u  =con(2)
v  =con(3)
p  =con(4)

q2= u*u + v*v 
a= dsqrt(GAMMA*p/rho)

un = u*nx + v*ny


la0(1)=un-a
la0(2)=un
la0(3)=un
la0(4)=un+a



Rp=un+2.d0*a/gamma1
Rm=un_inf-2.d0*ainf/gamma1


Unf=0.5d0*(Rp+Rm)
af=0.25d0*(Rp-Rm)*gamma1

if(Unf>0.d0) then
  uf=u+nx*(Unf-un)
  vf=v+ny*(Unf-un)
  Sf=p/(rho**gamma)
else
  uf=uinf+nx*(Unf-un_inf)
  vf=vinf+ny*(Unf-un_inf)
  Sf=pinf/(rinf**gamma)
endif

rhof=(af*af/gamma/Sf)**(1.d0/gamma1)
pf=rhof*af*af/gamma
ef=pf/gamma1 + 0.5d0*rhof*(uf*uf+vf*vf)

flux(1) = (ef + pf)*unf
flux(2) = rhof*unf
flux(3) = flux(2)*uf + pf*nx
flux(4) = flux(2)*vf + pf*ny


do i=1,nvar
cell(cc)%res(i)=cell(cc)%res(i)+dr*flux(i)
enddo


end




