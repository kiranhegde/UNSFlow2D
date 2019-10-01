subroutine roe_flux(ie,qcl,qcr,flux)
! --------------------------------------------------
! convective flux across a face using Roe's scheme
! ie   - interface edge number
! qcl  - left state, primitive variables(rho,u,v,p)
! qcr  - right state, primitive variables(rho,u,v,p)
! flux - interface flux 
! --------------------------------------------------
use data_type,only:dp,i4
use commons ,only : gamma,gamma1,nvar 
use grid,only : fc
implicit none
integer(kind=i4) :: i,ie
real(kind=dp)    :: qcl(nvar), qcr(nvar),flux(nvar)
real(kind=dp)    :: nx,ny,area
real(kind=dp)    :: rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr, &
                    ua, va, qa2, aa2, aa, ha, &
                    ql2, qr2, rl12, rr12, rd, &
                    unl, unr, una, vna, F_c(4), Fd(4), &
                    m1, m2, a1, a2, a3, a4, l1, l2, l3, l4, &
                    a1l1, a2l2, a3l3, a4l4, aact, aast, &
                    du1, du2, du3, du4, &
                    e1,e2, e4, del
real(kind=dp),parameter :: ETOL=0.1_dp

nx = fc(ie)%sx
ny = fc(ie)%sy
area = fc(ie)%area
nx = nx/area
ny = ny/area

rl = qcl(1)
ul = qcl(2)
vl = qcl(3)
pl = qcl(4)

rr = qcr(1)
ur = qcr(2)
vr = qcr(3)
pr = qcr(4)


ql2= ul*ul + vl*vl
al2= GAMMA*pl/rl
hl = al2/GAMMA1 + 0.5_dp*ql2

qr2= ur*ur + vr*vr
ar2= GAMMA*pr/rr
hr = ar2/GAMMA1 + 0.5_dp*qr2

!     Rotated velocity
unl = ul*nx + vl*ny
unr = ur*nx + vr*ny

!     Centered flux
f_c(1) = rl*hl*unl         + rr*hr*unr
f_c(2) = rl*unl            + rr*unr
f_c(3) = pl*nx + rl*ul*unl + pr*nx + rr*ur*unr
f_c(4) = pl*ny + rl*vl*unl + pr*ny + rr*vr*unr

!     Roe average
rl12 = dsqrt(rl)
rr12 = dsqrt(rr)
rd   = 1.0d0/(rl12 + rr12)

ua   = (ul*rl12 + ur*rr12)*rd
va   = (vl*rl12 + vr*rr12)*rd
ha   = (hl*rl12 + hr*rr12)*rd
qa2  = ua**2 + va**2
aa2  = GAMMA1*(ha - 0.5_dp*qa2)

!#ifdef DEBUG
if(aa2 .le. 0.0d0)then
   print*,'Roe scheme ....'
   print*,'Sonic speed is negative'
   print*,'Left/right cell values'
   print*,qcl(1),qcl(2),qcl(3),qcl(4)
   print*,qcr(1),qcr(2),qcr(3),qcr(4)
   print*,'Left/right vertex values'
!   print*,qvl(1),qvl(2),qvl(3),qvl(4)
!   print*,qvr(1),qvr(2),qvr(3),qvr(4)
!   print*
!   print*,rl,ul,vl,pl
!   print*,rr,ur,vr,pr
!   print*,li
   stop
endif
!#endif
aa  = dsqrt(aa2)
una = ua*nx + va*ny
vna =-ua*ny + va*nx

!     Eigenvalues with entropy fix
e1 = dabs(una - aa)
e2 = dabs(una)
e4 = dabs(una + aa)

del= ETOL*aa
if(e1 .lt. del)then
   l1 = 0.5_dp*(del + e1**2/del)
else
   l1 = e1
endif

if(e2 .lt. del)then
   l2 = 0.5_dp*(del + e2**2/del)
   l3=l2
else
   l2 = e2
   l3 = l2
endif

if(e4 .lt. del)then
   l4 = 0.5_dp*(del + e4**2/del)
else
   l4 = e4
endif

!     Difference of conserved variables
du1 = rr           - rl
du2 = rr*ur        - rl*ul
du3 = rr*vr        - rl*vl
du4 = (rr*hr - pr) - (rl*hl - pl)

!     Amplitudes
m1 = (nx*du2 + ny*du3 - una*du1)/aa
m2 = GAMMA1*(du4 - ua*du2 - va*du3 + qa2*du1)/aa**2

a1 = 0.5_dp*(m2 - m1)
a2 = ( ny*du2 - nx*du3 + vna*du1 )/aa
a4 = 0.5_dp*(m1 + m2)
a3 = du1 - a1 - a4

!     Diffusive flux
a1l1  = a1*l1
a2l2  = a2*l2
a3l3  = a3*l3
a4l4  = a4*l4
aact  = aa*nx
aast  = aa*ny

Fd(1) = a1l1*(ha - una*aa) + a2l2*aa*vna + a3l3*0.5_dp*qa2 + &
        a4l4*(ha + una*aa)
Fd(2) = a1l1               +               a3l3           + a4l4
Fd(3) = a1l1*(ua - aact)   + a2l2*aa*ny  + a3l3*ua        + &
        a4l4*(ua + aact)
Fd(4) = a1l1*(va - aast)   - a2l2*aa*nx  + a3l3*va        + &
        a4l4*(va + aast)

!     Total flux
do i=1,4
   flux(i) = 0.5_dp*area*( f_c(i) - Fd(i) )
enddo

end
