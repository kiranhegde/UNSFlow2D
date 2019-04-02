!==============================================================================
subroutine ausmPlus_flux(ie,c1,c2,qcl,qcr)
!==============================================================================
!#     computes total convective flux across a face using van leer
!#  flux vector splitting method given left and right conserved states
!#
!#     ql,qr - left & right vector of conserved variables.
!------------------------------------------------------------------------------
use commons
use pri
use grid
implicit none
!------------------------------------------------------------------------------
integer(kind=i4):: i,ie,c1,c2
real(kind=dp):: flux(nvar)
real(kind=dp):: rl, ul, vl, pl, cl, unl, ml, hl,al2,dl,ql2
real(kind=dp):: rr, ur, vr, pr, cr, unr, mr, hr,ar2,dr,qr2
real(kind=dp):: LIMIT 

real(kind=dp):: m12, mlp, mrm, aml, amr, mrp, mlm,p5p,p5m,m4p,m4m,alpha,beta
real(kind=dp):: p12, plp, prm,dmm, dpp, ro12, mp, m1m, m1p, m2p, m2m
real(kind=dp):: fluxL(nvar),fluxR(nvar),fluxP(nvar)
real(kind=dp):: mco, mref, mo, mbar, mbar2, mstar, mstar2, fa, c12, mass_p, mass_m, mass12, dissi
real(kind=dp):: sigma, kp, ku, pu
real(kind=dp):: aL_star,aR_star,aL_hat,aR_hat,atilR,atilL,astarL,astarR
real(kind=dp):: mass,fluxN,fluxT 
real(kind=dp) :: x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar), &
            qvr(nvar), resl(nvar), resr(nvar)
real(kind=dp):: li(nvar),con(nvar)
real(kind=dp):: nx,ny,area,temp,dist
!------------------------------------------------------------------------------

sigma=1.0d0
ku=0.75d0
kp=0.25d0
!------------------------------------------------------------------------------

nx = fc(ie)%sx
ny = fc(ie)%sy
area = dsqrt(nx*nx + ny*ny)
nx = nx/area
ny = ny/area

qcl(:)=0.d0
qcr(:)=0.d0


!     Left state
qcl(:)=cell(c1)%qp(:)
do i=1,ndim
dist=fc(ie)%cen(i)-cell(c1)%cen(i)
qcl(:)=qcl(:)+cell(c1)%grad(i,:)*dist
enddo

rl = qcl(1)
ul = qcl(2)
vl = qcl(3)
pl = qcl(4)

unl = ul*nx + vl*ny
ql2= ul*ul + vl*vl
al2= GAMMA*pl/rl
hl = al2/GAMMA1 + 0.5d0*ql2
cl =dsqrt(al2)

!     Right state
qcr(:)=cell(c2)%qp(:)
do i=1,ndim
dist=fc(ie)%cen(i)-cell(c2)%cen(i)
qcr(:)=qcr(:)+cell(c2)%grad(i,:)*dist
enddo

rr = qcr(1)
ur = qcr(2)
vr = qcr(3)
pr = qcr(4)

unr = ur*nx + vr*ny
qr2= ur*ur + vr*vr
ar2= GAMMA*pr/rr
hr = ar2/GAMMA1 + 0.5d0*qr2
cr =dsqrt(ar2)

aL_star=dsqrt(2.0*(gamma-1.0)*hl/(gamma+1.0))
aR_star=dsqrt(2.0*(gamma-1.0)*hr/(gamma+1.0))
aL_hat = aL_star*aL_star/dmax1(aL_star,unl)
aR_hat = aR_star*aR_star/dmax1(aR_star,-1.0*unr)
c12=dmin1(aL_hat,aR_hat)
!c12=0.5d0*(cl+cr)
ml = unl/c12
mr = unr/c12

mbar2  = 0.5d0*(ml*ml+mr*mr)

!mref = dmax1(m_inf,dsqrt(mbar2))
!mref = dmin1(mref,1.0)
mref = m_inf
!mref = 1e-2 
!mref = dmax1(0.32,0.5*m_inf)

if(mbar2 >=1.d0) then
fa = 1.0
else
mstar2 = dmin1(1.0,dmax1(mbar2,mref*mref))
Mo = dsqrt(mstar2)
Mo = mstar2
fa = Mo*(2.0-Mo)
endif


beta=1.0/8.0d0
alpha=3.0/16.0d0*(-4.0+5.0*fa*fa)
!alpha=3.0d0/16.0d0
plp=p5p(ml,alpha)
prm=p5m(mr,alpha)
mlp=m4p(ml,beta)
mrm=m4m(mr,beta)

ro12= 0.5d0*(rl+rr)

mp = (kp/fa)*max((1-sigma*mbar2),0.0)*(pl-pr)/(ro12*c12*c12) 
m12 = mlp + mrm + mp

!------->  Mass Flux
mass12 =  rr*c12*m12
if(m12>0.d0) mass12 =  rl*c12*m12

pu = ku*plp*prm*2.0*ro12*fa*c12*(unr-unl)
p12 = plp*pl+prm*pr+ pu

!------>  total  flux

flux(:)  = 0.0d0
fluxP(:) = 0.0d0
if(m12>0.d0) then
flux(1) = hl
flux(2) = 1.0d0
flux(3) = ul
flux(4) = vl
else
flux(1) = hr
flux(2) = 1.0d0
flux(3) = ur
flux(4) = vr
endif
!------>  Pressure Flux
fluxP(3)=nx*P12
fluxP(4)=ny*P12
!----------


do i=1,nvar
temp  =mass12*flux(i)+fluxP(i)
cell(c1)%res(i)=cell(c1)%res(i)+area*temp
cell(c2)%res(i)=cell(c2)%res(i)-area*temp
enddo


end subroutine ausmPlus_flux 
!==============================================================================
