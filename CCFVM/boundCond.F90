subroutine slip_wall(ie,qcl,qcr)
use commons
!use pri
use grid
implicit none
integer(kind=i4):: ie,cc
real(kind=dp) :: ds,nx,ny,con(nvar),un
real(kind=dp) :: qcl(nvar), qcr(nvar)

nx=fc(ie)%sx
ny=fc(ie)%sy
ds=dsqrt(nx*nx+ny*ny)
nx=nx/ds
ny=ny/ds

un=qcl(2)*nx+qcl(3)*ny

! Weak-Riemann : By setting ghost normal velocity to zero
!qcl(2)=qcl(2)-un*nx
!qcl(3)=qcl(3)-un*ny
qcr=qcl

! Weak-Riemann : By setting same tangential component in ghost as in the  interior
! V_ = Vin-2(Vin.n)*n
qcr(2)=(qcl(2)*qcl(1)-2.0_dp*un*nx)/qcl(1)
qcr(3)=(qcl(3)*qcl(1)-2.0_dp*un*ny)/qcl(1)

end

subroutine farfield_flux(ie,qcl,qcr)
use grid
use pri
use inf
implicit none
integer(kind=i4)  :: i,ie,cc
real(kind=dp) :: qcl(nvar), qcr(nvar)
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

rho=qcl(1)
u  =qcl(2)
v  =qcl(3)
p  =qcl(4)

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

qcr(1)=rhof
qcr(2)=uf
qcr(3)=vf
qcr(4)=pf

end
