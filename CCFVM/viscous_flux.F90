subroutine viscous_flux
use commons
use visc
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: ie,c1,c2
real(kind=dp):: two_third,two
real(kind=dp):: Tauxx,Tauxy,tauyy,mu

return

call Gradient_GG_FC
call sutherland

two_third=2.0_dp/3.0_dp
two=2.0_dp

do ie=1,nof
   c1 = fc(ie)%in
   c2 = fc(ie)%out
   mu=0.5_dp*(cell(c1)%mul+cell(c2)%mul)

   Tauxx=mu*two_third*(two*fc(ie)%grad(1,2)-fc(ie)%grad(2,3))/Rey
   Tauxy=mu*(fc(ie)%grad(2,2)+fc(ie)%grad(1,3))/Rey
   Tauyy=mu*two_third*(two*fc(ie)%grad(2,3)-fc(ie)%grad(1,2))/Rey

   !qx=-mu*       fc(ie)%grad(1,2)


!   cell(c1)%res(i)=cell(c1)%res(i)+flux*area
!   cell(c2)%res(i)=cell(c2)%res(i)-flux*area


enddo

contains
!
!!....Calculate total Reynolds number
!subroutine viscosity(prim, nut, mul, mu)
!!use param
!implicit none
!real(kind=dp) :: prim(nvar,npmax), nut(npmax), mul(npmax), &
!                      mu(npmax)
!
!integer(kind=i4) ::          i
!real(kind=dp) :: mut, mu_turb
!
!if(flow_type .eq. 'laminar')then
!   do i=1,noc
!      mu(i) = mul(i)
!   enddo
!elseif(flow_type .eq. 'turbulent')then
!   do i=1,noc
!      mut   = mu_turb(mul(i), prim(1,i), nut(i))
!      mu(i) = mul(i) + mut
!   enddo
!endif
!
!return
!end
!
!
!!....Calculate Reynolds number based on Sutherland law
subroutine sutherland
use param
use pri
implicit none
integer(kind=i4) ::i 
!real(kind=dp) :: sutherland_viscosity
real(kind=dp) :: flux(nvar) 

do i=1,noc
   flux(:)=cell(i)%qc(:)
   call con2prim(flux)
   cell(i)%mul= sutherland_viscosity(rho,p)
enddo

!return
end
!
!Sutherland viscosity formula
real(kind=dp) function sutherland_viscosity(density, pressure)
use param
implicit none
real(kind=dp) :: density, pressure
real(kind=dp) :: temp, num, den

temp = pressure/density/GAS_CONST
num  = T_inf + SCONST
den  = temp  + SCONST
sutherland_viscosity = (temp/T_inf)**1.5d0*(num/den)/Rey

!return
end

end subroutine viscous_flux
