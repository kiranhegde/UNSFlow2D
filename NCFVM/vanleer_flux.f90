!==============================================================================
subroutine vanleer_flux(ie,c1,c2,qcl,qcr)
!==============================================================================
!#  computes total convective flux across a face using van leer
!#  flux vector splitting method given left and right conserved states
!#
!#     qcl,qcr - left & right vector of conserved variables.
!------------------------------------------------------------------------------
use commons
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: i,ie,c1,c2
real(kind=dp):: flux
real(kind=dp):: rl, ul, vl, pl, cl, unl, ml, hl,dl,ql(nvar),ql2,al2
real(kind=dp):: rr, ur, vr, pr, cr, unr, mr, hr,dr,qr(nvar),qr2,ar2
real(kind=dp):: m12, mlp, mrm, aml, amr, mrp, mlm,p5p,p5m,m4p,m4m,alpha,beta
real(kind=dp):: p12, plp, prm, fa_lp, fa_rm, fb_lp, fb_rm, fmass12, delta
real(kind=dp):: fmass_p,fmass_m,fener_p,fener_m,dissi,ll,un
real(kind=dp):: fluxL(nvar),fluxR(nvar),fluxP(nvar),LIMIT,limit1
real(kind=dp) :: x1(2), x2(2), qcl(nvar), qcr(nvar), qvl(nvar), &
            qvr(nvar), resl(nvar), resr(nvar)
real(kind=dp):: li(nvar),nx,ny,area,con(nvar),dist(ndim),grad(nvar),rt

!------------------------------------------------------------------------------
nx = fc(ie)%nx 
ny = fc(ie)%ny
area = fc(ie)%area  

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
flux =0.5d0*fmass12*(fluxL(i)+fluxR(i))-0.5d0*dissi*(fluxR(i)-fluxL(i))+fluxP(i)
pt(c1)%res(i)=pt(c1)%res(i)+flux*area
pt(c2)%res(i)=pt(c2)%res(i)-flux*area

enddo

 RT = dsqrt(rR/rL)
rho = RT*rL
  u = (uL+RT*uR)/(1.d0+RT)
  v = (vL+RT*vR)/(1.d0+RT)
  H = (HL+RT* HR)/(1.d0+RT)
  a = dsqrt( (gamma1)*(H-0.5d0*(u*u+v*v)) )
 un = dabs(u*nx+v*ny)

!un=0.5d0*(dabs(unl)+dabs(unr))
!a=0.5d0*(cl+cr)
ll=(un+a)*area
pt(c1)%la=pt(c1)%la+ll
pt(c2)%la=pt(c2)%la+ll
fc(ie)%la=ll


end subroutine vanleer_flux
!==============================================================================
