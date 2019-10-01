! No Slip Wall BC
subroutine NoSlip_wall(ie)
use commons
use pri
use grid
implicit none
integer(kind=i4):: ie,in,out
!real(kind=dp) :: ds,nx,ny,uni,tg,pvar(nvar),rhog
real(kind=dp) :: tg,rhog
!real(kind=dp) :: qcl(nvar), qcr(nvar)

in  = fc(ie)%in
out = fc(ie)%out
if(out<noc) print*,out,'s'

call con2prim(cell(in)%qc(1:nvar))

! setting up Ghost state
if(iwall==0) then 
  ! Isothermal wall BC
  tg=2.0_dp*t_wall-t 
  if(tg<=0.0_dp) tg=0.01_dp*t_wall
  rhog=p/(rho*znd) 
  cell(out)%qp(1)=rhog
  cell(out)%qp(2)=-u
  cell(out)%qp(3)=-v
  cell(out)%qp(4)=p
  cell(out)%qp(5)=tg
  call prim2con(cell(out)%qp(1:nvar)) 
  cell(out)%qc(:)=conv(:)
elseif(iwall==1) then 
  ! Adiabatic wall  BC
  cell(out)%qc(1:2)=cell(in)%qc(1:2)
  cell(out)%qc(3)=-cell(in)%qc(3)
  cell(out)%qc(4)=-cell(in)%qc(4)
  call con2prim(cell(out)%qc(1:nvar))
  cell(out)%qp(:)=prim(:)
else
    Stop 'Unknown no-slip wall BC'
endif
   

end  subroutine   NoSlip_wall
!============================================================================
! Slip Wall BC
subroutine slip_wall(ie)
use commons
use pri
use grid
implicit none
integer(kind=i4):: ie,in,out
real(kind=dp) :: area,nx,ny,un
!real(kind=dp) :: qcl(nvar), qcr(nvar)

in  = fc(ie)%in
out = fc(ie)%out
if(out<noc) print*,out,'s'
nx=fc(ie)%sx
ny=fc(ie)%sy
area=fc(ie)%area
nx=nx/Area
ny=ny/Area

call con2prim(cell(in)%qc(1:nvar))
un=u*nx+v*ny

cell(out)%qc(1:2)=cell(in)%qc(1:2)
! Weak-Riemann : By setting ghost normal velocity to zero
!qcl(2)=qcl(2)-un*nx
!qcl(3)=qcl(3)-un*ny
!qcr=qcl

! Weak-Riemann : By setting same tangential component in ghost as in the  interior
! V_ = Vin-2(Vin.n)*n
cell(out)%qc(3)=rho*u-2.0_dp*un*nx
cell(out)%qc(4)=rho*v-2.0_dp*un*ny

call con2prim(cell(out)%qc(1:nvar))
cell(out)%qp(:)=prim(:)

end  subroutine  slip_wall
!============================================================================
subroutine slip_wall1(ie,qcl,qcr)
use commons
!use pri
use grid
implicit none
integer(kind=i4):: ie
real(kind=dp) :: area,nx,ny,un
real(kind=dp) :: qcl(nvar), qcr(nvar)

nx=fc(ie)%sx
ny=fc(ie)%sy
area=fc(ie)%area
nx=nx/Area
ny=ny/Area

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
! ============================================================================
subroutine farfield_flux(ie,qcl,qcr)
use grid
use pri
use inf
implicit none
integer(kind=i4)  :: ie
real(kind=dp) :: qcl(nvar), qcr(nvar)
real(kind=dp) :: area, nx, ny, q2,un,tx,ty,tn
real(kind=dp) :: uinf, vinf, pinf, rinf, ainf
real(kind=dp) :: Rp,Rm,rhob,ub,vb,pb,ab,Sb
real(kind=dp) :: q2inf,un_inf,Vnrm,Vtan,Vninf,Vtinf
real(kind=dp) :: lam1,lam2,lam3,mach


nx = fc(ie)%sx
ny = fc(ie)%sy
area=fc(ie)%area
nx=nx/Area
ny=ny/Area
tx=-ny
ty= nx

rinf  = r_inf
uinf  = u_inf
vinf  = v_inf
pinf  = p_inf
q2inf = uinf**2 + vinf**2
!ainf  = a_inf
ainf= dsqrt(GAMMA*pinf/rinf)
un_inf = uinf*nx + vinf*ny

rho=qcl(1)
u  =qcl(2)
v  =qcl(3)
p  =qcl(4)
q2= u*u + v*v
a= dsqrt(GAMMA*p/rho)
un = u*nx + v*ny
tn = u*tx + v*ty

rhob=rinf
ub=uinf
vb=vinf
pb=pinf

lam1=un+a
lam2=un-a
lam3=un
mach=un/a


if(mach>=1.0_dp) then
   Vnrm=un
   Vtan=tn 
   ab=a
   rhob=rho
   ub=Vnrm*nx+Vtan*tx 
   vb=Vnrm*ny+Vtan*ty 
   pb=p
elseif(0.0_dp<=mach.and.mach<1.0_dp) then
   Vninf=uinf*nx+vinf*ny
   Vtinf=uinf*tx+vinf*ty
   Rp=un+2.0_dp*a/gamma1
   Rm=Vninf-2.0_dp*ainf/gamma1
   Sb=p/(rho**gamma)
   Vnrm=0.5_dp*(Rp+Rm)
   Vtan=tn 
   ub=Vnrm*nx+Vtan*tx 
   vb=Vnrm*ny+Vtan*ty 
   ab=0.25_dp*(Rp-Rm)*gamma1
   rhob=(ab*ab/gamma/Sb)**(1.0_dp/gamma1)
   pb=rhob*ab*ab/gamma
elseif(-1.0_dp<mach.and.mach<0.0_dp) then
   Vninf=uinf*nx+vinf*ny
   Vtinf=uinf*tx+vinf*ty
   Rp=un+2.0_dp*a/gamma1
   Rm=Vninf-2.0_dp*ainf/gamma1
   Sb=pinf/(rinf**gamma)
   Vnrm=0.5_dp*(Rp+Rm)
   Vtan=Vtinf 
   ub=Vnrm*nx+Vtan*tx 
   vb=Vnrm*ny+Vtan*ty 
   ab=0.25_dp*(Rp-Rm)*gamma1
   rhob=(ab*ab/gamma/Sb)**(1.0_dp/gamma1)
   pb=rhob*ab*ab/gamma
elseif(-1.0_dp<=mach) then
   rhob=rinf
   ub=uinf   
   vb=vinf   
   pb=pinf
endif

qcr(1)=rhob
qcr(2)=ub
qcr(3)=vb
qcr(4)=pb

end

!subroutine farfield_flux(ie,qcl,qcr)
!use grid
!use pri
!use inf
!implicit none
!integer(kind=i4)  :: ie
!real(kind=dp) :: qcl(nvar), qcr(nvar)
!real(kind=dp) :: dr, nx, ny, q2,un
!real(kind=dp) :: uinf, vinf, pinf, rinf, ainf
!real(kind=dp) :: Rp,Rm,rhof,uf,vf,pf,af,Unf,Sf,Ef
!real(kind=dp) :: q2inf,un_inf
!real(kind=dp) :: lamda1,lamda3
!
!
!nx = fc(ie)%sx
!ny = fc(ie)%sy
!dr =  dsqrt(nx*nx + ny*ny)
!nx = nx/dr
!ny = ny/dr
!
!rinf  = r_inf
!uinf  = u_inf
!vinf  = v_inf
!pinf  = p_inf
!q2inf = uinf**2 + vinf**2
!ainf  = a_inf
!un_inf = uinf*nx + vinf*ny
!
!rho=qcl(1)
!u  =qcl(2)
!v  =qcl(3)
!p  =qcl(4)
!q2= u*u + v*v
!a= dsqrt(GAMMA*p/rho)
!un = u*nx + v*ny
!
!lamda1=un-a
!lamda3=un+a
!if(lamda1<0.0_dp) then
!   Rm=un_inf-2.0_dp*ainf/gamma1
!else   
!   Rm=un-2.0_dp*a/gamma1
!endif
!
!if(lamda3<0.0_dp) then
!  Rp=un_inf+2.0_dp*ainf/gamma1
!else   
!  Rp=un+2.0_dp*a/gamma1
!endif
!
!!Rp=un+2.0_dp*a/gamma1
!!Rm=un_inf-2.0_dp*ainf/gamma1
!
!Unf=0.5_dp*(Rp+Rm)
!af=0.25_dp*(Rp-Rm)*gamma1
!
!if(Unf<=0.0_dp) then
!  !uf=u+nx*(Unf-un)
!  !vf=v+ny*(Unf-un)
!  uf=u-nx*un
!  vf=v-ny*un
!  Sf=p/(rho**gamma)
!else
!  !uf=uinf+nx*(Unf-un_inf)
!  !vf=vinf+ny*(Unf-un_inf)
!  uf=uinf-nx*un_inf
!  vf=vinf-ny*un_inf
!  Sf=pinf/(rinf**gamma)
!endif
!
!rhof=(af*af/gamma/Sf)**(1.0_dp/gamma1)
!pf=rhof*af*af/gamma
!ef=pf/gamma1 + 0.5_dp*rhof*(uf*uf+vf*vf)
!
!qcr(1)=rhof
!qcr(2)=uf
!qcr(3)=vf
!qcr(4)=pf
!
!end

