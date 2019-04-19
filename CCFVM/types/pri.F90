!---------------------------------------------------------
module pri
use data_type
use inf
implicit none

real(kind=dp):: rho,u,v,e,q,p,h,a,t,znd
real(kind=dp):: prim(nvar)
real(kind=dp):: conv(nvar)

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

subroutine prim2con(prim)
!converts primitive variables 2 conserved variable 
use data_type
use inf
implicit none
real(kind=dp) :: prim(nvar)
   
   rho=prim(1)
     u=prim(2)
     v=prim(3)
     p=prim(4)

   ! Et(rho*e)
   conv(1) = p/(gamma-1.0_dp) + 0.5_dp*rho*(u**2 + v**2)
   conv(2)=rho
   conv(3)=rho*u
   conv(4)=rho*v

end subroutine prim2con

end module pri

