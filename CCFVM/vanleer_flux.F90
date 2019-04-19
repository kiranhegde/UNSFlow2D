subroutine vanleer_flux(ie,qcl,qcr,flux)
! convective flux across a face using VanLeer's scheme
! ie   - interface edge number
! qcl  - left state, primitive variables(rho,u,v,p)
! qcr  - right state, primitive variables(rho,u,v,p)
! flux - interface flux 
!------------------------------------------------------------------------------
use data_type,only:dp,i4
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: i,ie
real(kind=dp):: rl, ul, vl, pl, cl, unl, ml, hl,ql2,al2
real(kind=dp):: rr, ur, vr, pr, cr, unr, mr, hr,qr2,ar2
real(kind=dp):: mlp, mrm, mrp, mlm,p5p,p5m,m4p,m4m,alpha,beta
real(kind=dp):: p12, plp, prm, fmass12, delta
real(kind=dp):: fmass_p,fmass_m,dissi
real(kind=dp):: fluxL(nvar),fluxR(nvar),fluxP(nvar),flux(nvar)
real(kind=dp):: qcl(nvar), qcr(nvar)
real(kind=dp):: nx,ny,area

!------------------------------------------------------------------------------
nx = fc(ie)%sx 
ny = fc(ie)%sy
area = dsqrt(nx*nx + ny*ny)
nx = nx/area
ny = ny/area


rl = qcl(1)
ul = qcl(2)
vl = qcl(3)
pl = qcl(4)

unl = ul*nx + vl*ny
ql2= ul*ul + vl*vl
al2= GAMMA*pl/rl
hl = al2/GAMMA1 + 0.5d0*ql2
cl =dsqrt(al2)
ml = unl/cl


rr = qcr(1)
ur = qcr(2)
vr = qcr(3)
pr = qcr(4)

unr = ur*nx + vr*ny
qr2= ur*ur + vr*vr
ar2= GAMMA*pr/rr
hr = ar2/GAMMA1 + 0.5d0*qr2
cr =dsqrt(ar2)
mr = unr/cr

beta=0.0d0
alpha=0.0d0
!beta=1.0d0/8.0d0
!alpha=3.0d0/16.0d0
plp=p5p(ml,alpha)
prm=p5m(mr,alpha)
mlp=m4p(ml,beta)
mlm=m4m(ml,beta)
mrp=m4p(mr,beta)
mrm=m4m(mr,beta)

!------>  total  flux
p12 = plp*pl+prm*pr 
fmass_p = rl*cl*mlp
fmass_m = rr*cr*mrm 
fmass12 = fmass_p+fmass_m
delta = -gamma*(cr*mrp*mrm*pr-cl*mlp*mlm*pl)/(gamma-1.0d0)
dissi = fmass_p-fmass_m

fluxL(:) = 0.0d0
fluxR(:) = 0.0d0
flux  = 0.0d0
fluxP(:) = 0.0d0

fluxL(1) = hl
fluxL(2) = 1.0d0
fluxL(3) = ul
fluxL(4) = vl

fluxR(1) = hr
fluxR(2) = 1.0d0
fluxR(3) = ur
fluxR(4) = vr

!------>  Pressure Flux
fluxP(3)=nx*P12
fluxP(4)=ny*P12
fluxP(1)=delta
do i=1,nvar
flux(i) =0.5d0*fmass12*(fluxL(i)+fluxR(i))-0.5d0*dissi*(fluxR(i)-fluxL(i))+fluxP(i)
enddo

flux=flux*area

end subroutine vanleer_flux
!==============================================================================
