subroutine viscous_flux
use commons
use visc
use pri
use grid
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: ie,in,out
real(kind=dp):: two_third,two,MubyRey
real(kind=dp):: uvel,vvel,temprature 
real(kind=dp):: ux,uy,vx,vy,tx,ty,qx,qy 
real(kind=dp):: Tauxx,Tauxy,tauyy,mu
real(kind=dp):: flux(nvar),Pr
real(kind=dp):: nx,ny,area,factor


!return

if(trim(grad_type)/='ggfc') call Gradient_GG_FC
call sutherland

two_third=2.0_dp/3.0_dp
two=2.0_dp
flux=0.0_dp

do ie=startFC,endFC
!do ie=1,nof
   nx = fc(ie)%sx
   ny = fc(ie)%sy
   area = dsqrt(nx*nx + ny*ny)
   nx = nx/area
   ny = ny/area
   in = fc(ie)%in
   out= fc(ie)%out
   if(out>0) then 
      mu=0.5_dp*(cell(in)%mul+cell(out)%mul)
   else
      mu=cell(in)%mul
   endif
 
   uvel=fc(ie)%qp(2) 
   vvel=fc(ie)%qp(3) 
   temprature=fc(ie)%qp(5) 
   ux=fc(ie)%grad(1,2)
   uy=fc(ie)%grad(2,2)
   vx=fc(ie)%grad(1,3)
   vy=fc(ie)%grad(2,3)
   tx=fc(ie)%grad(1,5)
   ty=fc(ie)%grad(2,5)
   MubyRey=mu/Rey

   Tauxx=MubyRey*two_third*(two*ux-vy)
   Tauxy=MubyRey*(uy+vx)
   Tauyy=MubyRey*two_third*(two*vy-ux)
    
   Pr=prandtl
   factor=(gamma-1.0_dp)*m_inf*m_inf*Pr 
   qx=-tx/factor
   qy=-ty/factor

   flux(1)=  (uvel*tauxx+vvel*tauxy-qx)*nx &
          & -(uvel*tauxy+vvel*tauyy-qy)*ny 
   flux(3)=tauxx*nx+tauxy*ny  
   flux(4)=tauxy*nx+tauyy*ny  

   cell(in)%res(:)=cell(in)%res(:)+flux(:)*area
   if(out>0) cell(out)%res(:)=cell(out)%res(:)-flux(:)*area
enddo


!do i=startBC,endBC
!   c1 = fc(i)%in
!enddo


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
