!==============================================================================
subroutine ausmPlus_flux(ie,qcl,qcr,flux)
!==============================================================================
! convective flux across a face using AUSM+ -up scheme
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
!------------------------------------------------------------------------------
integer(kind=i4):: ie
real(kind=dp):: flux(nvar),fluxP(nvar)
real(kind=dp):: fluxl(nvar),fluxr(nvar)
real(kind=dp):: rl, ul, vl, pl, cl, unl, ml, hl,al2,ql2
real(kind=dp):: rr, ur, vr, pr, cr, unr, mr, hr,ar2,qr2

real(kind=dp):: m12, mlp, mrm, p5p,p5m,m4p,m4m,alpha,beta,m1m,m1p
real(kind=dp):: p12, plp, prm, ro12, mp
real(kind=dp):: mref, mo,  mbar2,  mstar2, fa, c12
real(kind=dp):: sigma, kp, ku, pu
real(kind=dp):: aL_star,aR_star,aL_hat,aR_hat
real(kind=dp):: qcl(nvar), qcr(nvar)
real(kind=dp):: nx,ny,area
!------------------------------------------------------------------------------

sigma=1.0_dp
ku=0.75_dp
kp=0.25_dp
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
hl = al2/GAMMA1 + 0.5_dp*ql2
cl =dsqrt(al2)

rr = qcr(1)
ur = qcr(2)
vr = qcr(3)
pr = qcr(4)

unr = ur*nx + vr*ny
qr2= ur*ur + vr*vr
ar2= GAMMA*pr/rr
hr = ar2/GAMMA1 + 0.5_dp*qr2
cr =dsqrt(ar2)

aL_star=dsqrt(2.0_dp*(gamma-1.0_dp)*hl/(gamma+1.0_dp))
aR_star=dsqrt(2.0_dp*(gamma-1.0_dp)*hr/(gamma+1.0_dp))
aL_hat = aL_star*aL_star/dmax1(aL_star,dabs(unl))
aR_hat = aR_star*aR_star/dmax1(aR_star,-1.0_dp*dabs(unr))
c12=dmin1(aL_hat,aR_hat)
!c12=0.5_dp*(cl+cr)
ml = unl/c12
mr = unr/c12

mbar2  = 0.5_dp*(unl*unl+unr*unr)/(c12*c12)

!mref = dmax1(m_inf,dsqrt(mbar2))
!mref = dmin1(mref,1.0)
!mref = m_inf
!mref = 1e-2 
mref = dmax1(0.32,0.5*m_inf)

if(mbar2 >=1.0_dp) then
fa = 1.0_dp
else
mstar2 = dmin1(1.0_dp,dmax1(mbar2,mref*mref))
Mo = dsqrt(mstar2)
Mo = mstar2
fa = Mo*(2.0_dp-Mo)
endif


beta=1.0_dp/8.0_dp
alpha=3.0_dp/16.0_dp*(-4.0_dp+5.0_dp*fa*fa)
!alpha=3.0_dp/16.0_dp
plp=p5p(ml,alpha)
prm=p5m(mr,alpha)
mlp=m4p(ml,beta)
mrm=m4m(mr,beta)

ro12= 0.5_dp*(rl+rr)

mp = -(kp/fa)*dmax1((1.0_dp-sigma*mbar2),0.0_dp)*(pr-pl)/(ro12*c12*c12) 
m12 = mlp + mrm + mp

!------->  Mass Flux
!mass12 =  rr*c12*m12
!if(m12>eps) mass12 =  rl*c12*m12

pu = -ku*plp*prm*2.0_dp*ro12*fa*c12*(unr-unl)
p12 = plp*pl+prm*pr+pu

!------>  total  flux
flux(:)  = 0.0_dp
fluxl(:) = 0.0_dp
fluxr(:) = 0.0_dp
fluxP(:) = 0.0_dp

fluxl(1) = hl
fluxl(2) = 1.0_dp
fluxl(3) = ul
fluxl(4) = vl

fluxr(1) = hr
fluxr(2) = 1.0_dp
fluxr(3) = ur
fluxr(4) = vr

fluxP(3)=nx*P12
fluxP(4)=ny*P12

flux=c12*(m1p(m12)*rl*fluxl+m1m(m12)*rr*fluxr)

flux(3)=flux(3)+fluxP(3)
flux(4)=flux(4)+fluxP(4)
flux=flux*area

!return
!
!if(m12>0._dp) then
!flux(1) = hl
!flux(2) = 1.0_dp
!flux(3) = ul
!flux(4) = vl
!else
!flux(1) = hr
!flux(2) = 1.0_dp
!flux(3) = ur
!flux(4) = vr
!endif
!!------>  Pressure Flux
!fluxP(3)=nx*P12
!fluxP(4)=ny*P12
!!----------
!
!
!do i=1,nvar
!flux(i) =mass12*flux(i)+fluxP(i)
!!cell(c1)%res(i)=cell(c1)%res(i)+area*temp
!!cell(c2)%res(i)=cell(c2)%res(i)-area*temp
!enddo
!flux=flux*area

end subroutine ausmPlus_flux 
!==============================================================================
