subroutine farfield_flux1(ie,c1,c2)
use grid
use pri
use inf
implicit none
integer(kind=i4)  :: i,ie,c1,c2
real(kind=dp) :: dr, nx, ny, q2,un,dx,dy,con(nvar)
real(kind=dp) :: uinf, vinf, pinf, rinf, ainf 
real(kind=dp) :: Rp,Rm,rhof,uf,vf,pf,af,Unf,Sf,Ef 
real(kind=dp) :: q2inf,un_inf,flux(nvar),dist,area 
            

nx=-0.5*(pt(c2)%y-pt(c1)%y)
ny= 0.5*(pt(c2)%x-pt(c1)%x)

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
!con(:)=pt(cc)%qc(:)
!call con2prim(con)

!dx=fc(ie)%ldx
!dy=fc(ie)%ldy
!if(fc(ie)%ldx==0.d0.and.  fc(ie)%ldy==0.d0) then
!dx=fc(ie)%rdx
!dy=fc(ie)%rdy
!endif
con(:)=0.5d0*(pt(c1)%qp(:)+pt(c2)%qp(:))
!do i=1,ndim
!dist=fc(ie)%cen(i)-pt(cc)%cen(i)
!con(:)=con(:)+pt(cc)%grad(i,:)*dist
!enddo

!con(:)=pt(cc)%qp(:)+(pt(cc)%qx(:)*dx+pt(cc)%qy(:)*dy)


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

area=dsqrt(nx*nx+ny*ny)

pt(c1)%la=pt(c1)%la+(dabs(unf)+dsqrt(gamma*pt(c1)%qp(4)/pt(c1)%qp(1)))*area
pt(c2)%la=pt(c2)%la+(dabs(unf)+dsqrt(gamma*pt(c2)%qp(4)/pt(c2)%qp(1)))*area


do i=1,nvar
pt(c1)%res(i)=pt(c1)%res(i)+dr*flux(i)
pt(c2)%res(i)=pt(c2)%res(i)+dr*flux(i)
enddo


end
!=================================================================
subroutine farfield_flux2(ie,c1,c2)
use grid
use pri
use inf
implicit none
integer(kind=i4)  :: ie,c1,c2
integer(kind=i4) ::  i, j, k
real(kind=dp) :: dx, dy, dref, theta, circ, fact1, fact2, fact3,  &
                fact4, fact, uinf, vinf, pinf, rinf, ainf, &
                q2inf, flux1(nvar), flux2(nvar), &
                ua, va, pa, qa2, aa2, aa, ha, ra, &
                una, vna, ct, st, lent, &
                m2, t1, t2, t3, t4, t5, l1, l2, l3, l4, &
                lp1, lp2, lp3, lp4, lm1, lm2, lm3, lm4, &
                S(nvar,nvar), TT(nvar,nvar), Tp(nvar,nvar), &
                Tm(nvar,nvar), jacp(nvar,nvar), &
                jacm(nvar,nvar), Uin1(nvar), Uin2(nvar), &
                Uout(nvar),prim1(nvar),prim2(nvar)


prim1(:)=pt(c1)%qp(:)
prim2(:)=pt(c2)%qp(:)

   uinf  = u_inf
   vinf  = v_inf
   q2inf = uinf**2 + vinf**2
   pinf  = p_inf
   rinf  = r_inf
   ainf  = a_inf

!    Conserved state for infinity
Uout(1) = rinf
Uout(2) = rinf*uinf
Uout(3) = rinf*vinf
Uout(4) = pinf/gamma1 + 0.5d0*rinf*q2inf


!Average state on edge
ra = 0.5d0*(prim1(1) + prim2(1))
ua = 0.5d0*(prim1(2) + prim2(2))
va = 0.5d0*(prim1(3) + prim2(3))
pa = 0.5d0*(prim1(4) + prim2(4))
qa2= ua**2 + va**2
aa2= GAMMA*pa/ra
aa = dsqrt(gamma*pa/ra)
ha = aa2/gamma1 + 0.5d0*qa2

!    Conserved state for prim1
Uin1(1) = prim1(1)
Uin1(2) = prim1(1)*prim1(2)
Uin1(3) = prim1(1)*prim1(3)
Uin1(4) = prim1(4)/gamma1 + 0.5d0*prim1(1)*(prim1(2)**2 + &
                                                 prim1(3)**2)

!    Conserved state for prim2
Uin2(1) = prim2(1)
Uin2(2) = prim2(1)*prim2(2)
Uin2(3) = prim2(1)*prim2(3)
Uin2(4) = prim2(4)/gamma1 + 0.5d0*prim2(1)*(prim2(2)**2 + &
                                                 prim2(3)**2)

ct=-(pt(c2)%y-pt(c1)%y)
st= (pt(c2)%x-pt(c1)%x)

!ct   = fc(ie)%sx 
!st   = fc(ie)%sy 
lent = dsqrt(ct**2 + st**2)
ct   = ct/lent
st   = st/lent

!   Rotated velocity
una = ua*ct + va*st
vna =-ua*st + va*ct

!   Eigenvalues
l1 = una
l2 = una
l3 = una + aa
l4 = una - aa

!   Positive Eigenvalues
lp1 = dmax1(l1, 0.0d0)
lp2 = lp1
lp3 = dmax1(l3, 0.0d0)
lp4 = dmax1(l4, 0.0d0)

!   Negative Eigenvalues
lm1 = l1 - lp1
lm2 = l2 - lp2
lm3 = l3 - lp3
lm4 = l4 - lp4

!    Right eigenvector matrix
t1 = 0.5d0*ra/aa
m2 = qa2/aa2

TT(1,1) = 1.0d0
TT(2,1) = ua
TT(3,1) = va
TT(4,1) = 0.5d0*qa2

TT(1,2) = 0.0d0
TT(2,2) = ra*st
TT(3,2) = -ra*ct
TT(4,2) = -ra*vna

TT(1,3) = t1
TT(2,3) = t1*(ua + aa*ct)
TT(3,3) = t1*(va + aa*st)
TT(4,3) = t1*(ha + aa*una)

TT(1,4) = t1
TT(2,4) = t1*(ua - aa*ct)
TT(3,4) = t1*(va - aa*st)
TT(4,4) = t1*(ha - aa*una)

!    Inverse of right eigenvector matrix
t1     = 0.5d0*gamma1*m2*aa/ra
t2     = una/ra
t3     = gamma1*ua/aa
t4     = gamma1*va/aa
t5     = gamma1/(ra*aa)

S(1,1) = 1.0d0 - 0.5d0*gamma1*m2
S(2,1) = vna/ra
S(3,1) = t1 - t2
S(4,1) = t1 + t2

S(1,2) = t3/aa
S(2,2) = st/ra
S(3,2) = (ct - t3)/ra
S(4,2) =-(ct + t3)/ra

S(1,3) = t4/aa
S(2,3) = -ct/ra
S(3,3) = (st - t4)/ra
S(4,3) =-(st + t4)/ra

S(1,4) = -gamma1/aa2
S(2,4) = 0.0d0
S(3,4) = t5
S(4,4) = t5

!    Multiply T * lambda
do i=1,nvar
   Tp(i,1) = lp1*TT(i,1)
   Tp(i,2) = lp2*TT(i,2)
   Tp(i,3) = lp3*TT(i,3)
   Tp(i,4) = lp4*TT(i,4)
   Tm(i,1) = lm1*TT(i,1)
   Tm(i,2) = lm2*TT(i,2)
   Tm(i,3) = lm3*TT(i,3)
   Tm(i,4) = lm4*TT(i,4)
enddo

!    Now multiply with S to get pos/neg jacobians
do i=1,nvar
   do j=1,nvar
      jacp(i,j) = 0.0d0
      jacm(i,j) = 0.0d0
      do k=1,nvar
         jacp(i,j) = jacp(i,j) + Tp(i,k)*S(k,j)
         jacm(i,j) = jacm(i,j) + Tm(i,k)*S(k,j)
      enddo
   enddo
enddo


!    Finally the flux jacp*Uin + jacm*Uout
do i=1,nvar
   flux1(i) = 0.0d0
   flux2(i) = 0.0d0
   do j=1,nvar
      flux1(i) = flux1(i) + jacp(i,j)*Uin1(j) + jacm(i,j)*Uout(j)
      flux2(i) = flux2(i) + jacp(i,j)*Uin2(j) + jacm(i,j)*Uout(j)
   enddo
   flux1(i) = 0.5d0*lent*flux1(i)
   flux2(i) = 0.5d0*lent*flux2(i)
   !pt(c1)%res(i) = pt(c1)%res(i) - flux1(i)
   !pt(c2)%res(i) = pt(c2)%res(i) - flux2(i)

enddo

   pt(c1)%res(1) = pt(c1)%res(1) + flux1(4)
   pt(c1)%res(2) = pt(c1)%res(2) + flux1(1)
   pt(c1)%res(3) = pt(c1)%res(3) + flux1(2)
   pt(c1)%res(4) = pt(c1)%res(4) + flux1(3)

   pt(c2)%res(1) = pt(c2)%res(1) + flux2(4)
   pt(c2)%res(2) = pt(c2)%res(2) + flux2(1)
   pt(c2)%res(3) = pt(c2)%res(3) + flux2(2)
   pt(c2)%res(4) = pt(c2)%res(4) + flux2(3)

return
end
