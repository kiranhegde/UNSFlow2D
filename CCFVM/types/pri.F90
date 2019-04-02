!---------------------------------------------------------
module pri
use data_type
use inf
implicit none

real(kind=dp):: rho,u,v,e,q,p,h,a,t,znd
real(kind=dp):: prim(nvar)

contains

subroutine con2prim(qtemp)
!converts conserved variable to primitive variables

implicit none
real(kind=dp):: qtemp(nvar),Et


rho = qtemp(2)
u   = qtemp(3)/qtemp(2)
v   = qtemp(4)/qtemp(2)
Et  = qtemp(1)
e   = Et/rho
q   = u*u + v*v
p   = (gamma-1.d0)*(Et - 0.5d0*rho*q)
h   = (Et+p)/rho
a   = dsqrt(gamma*p/rho)
znd = 1.d0/(gamma*m_inf*m_inf)
t   = p/(rho*znd)

prim(1)=rho
prim(2)=u
prim(3)=v
prim(4)=p
end subroutine con2prim

end module pri

