subroutine rusanov_flux(ie,qcl,qcr,flux)
! convective flux across a face using Rusanov's scheme
! ie   - interface edge number
! qcl  - left state, primitive variables(rho,u,v,p)
! qcr  - right state, primitive variables(rho,u,v,p)
! flux - interface flux 
use data_type,only:dp,i4
use commons
use pri
use grid
implicit none
integer(kind=i4) :: ie
real(kind=dp):: flux(nvar)
real(kind=dp):: fl(nvar),fr(nvar)
real(kind=dp):: qcl(nvar), qcr(nvar)
real(kind=dp):: area,nx,ny
real(kind=dp):: rol,ul,vl,pl,cl,unl, hl,ql2,al2
real(kind=dp):: ror,ur,vr,pr,cr,unr, hr,qr2,ar2
real(kind=dp):: eigen,spec_rad_r,spec_rad_l

nx = fc(ie)%nx
ny = fc(ie)%ny
area = fc(ie)%area
!nx = nx/area
!ny = ny/area

! Left state
rol = qcl(1)
ul = qcl(2)
vl = qcl(3)
pl = qcl(4)
ql2= ul*ul + vl*vl
al2= GAMMA*pl/rol
hl = al2/GAMMA1 + 0.5d0*ql2
cl  = dsqrt(gamma*pl/rol)

! Right state
ror = qcr(1)
ur = qcr(2)
vr = qcr(3)
pr = qcr(4)
qr2= ur*ur + vr*vr
ar2= GAMMA*pr/ror
hr = ar2/GAMMA1 + 0.5d0*qr2
cr  =dsqrt(gamma*pr/ror)

! face normal velocity 
unl = ul*nx + vl*ny
unr = ur*nx + vr*ny

spec_rad_r=dabs(unl)+cl
spec_rad_l=dabs(unl)+cr

fl(1) = rol*hl*unl!+spec_rad_l*rol*hl
fl(2) = rol*unl!+spec_rad_l*rol
fl(3) = rol*unl*ul+nx*pl
fl(4) = rol*unl*vl+ny*pl

fr(1) = ror*hr*unr!-spec_rad_r*ror*hr
fr(2) = ror*unr!-spec_rad_r*ror
fr(3) = ror*unr*ur+nx*pr
fr(4) = ror*unr*vr+ny*pr

call prim2con(qcl)
qcl=conv 
call prim2con(qcr)
qcr=conv 

!  eigenvalue
eigen = dmax1(dabs(unl)+cl,dabs(unr)+cr)
! Flux
flux(:) = fl(:)+fr(:)
flux(:) = flux(:)+eigen*(qcl(:)-qcr(:))
flux(:) = 0.5*area*flux(:)

!------------------------------------------------------------------------------
end subroutine rusanov_flux
