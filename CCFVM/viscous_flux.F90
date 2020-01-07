subroutine viscous_flux
use commons
use visc
use pri
use grid
use output
implicit none
!!------------------------------------------------------------------------------
integer(kind=i4):: i,ie,in,out,bc
real(kind=dp):: two_third,two,MubyRey
real(kind=dp):: uvel,vvel
real(kind=dp):: ux,uy,vx,vy,tx,ty,qx,qy 
real(kind=dp):: Tauxx,Tauxy,tauyy,mu
real(kind=dp):: flux(nvar),Pr
real(kind=dp):: nx,ny,area,factor

if(trim(grad_type)/='ggfc') call Gradient_GG_FC
!call sutherland

do i=1,noc
   cell(i)%mul= sutherland_viscosity(i)
enddo

two=2.0_dp
two_third=two/3.0_dp
qx=0.0_dp
qy=0.0_dp

!do ie=startFC,endFC
do ie=1,nof
   flux=0.0_dp
   nx = fc(ie)%nx
   ny = fc(ie)%ny
   area = fc(ie)%area 
   in = fc(ie)%in
   out= fc(ie)%out
   bc= fc(ie)%bc 
   if(out<=noc) then 
      mu=0.5_dp*(cell(in)%mul+cell(out)%mul)
   else
      mu=cell(in)%mul
   endif
   fc(ie)%mu=mu
 
   uvel=fc(ie)%qp(2) 
   vvel=fc(ie)%qp(3) 
   ux=fc(ie)%grad(1,2)
   uy=fc(ie)%grad(2,2)
   vx=fc(ie)%grad(1,3)
   vy=fc(ie)%grad(2,3)
   tx=fc(ie)%grad(1,5)
   ty=fc(ie)%grad(2,5)

   Tauxx=two_third*(two*ux-vy)
   Tauxy=(uy+vx)
   Tauyy=two_third*(two*vy-ux)
    
   Pr=prandtl
   factor=(gamma-1.0_dp)*m_inf*m_inf*Pr 
   qx=-tx/factor
   qy=-ty/factor

   MubyRey=mu/Rey

   flux(1)=MubyRey*( (uvel*tauxx+vvel*tauxy-qx)*nx &
                  & +(uvel*tauxy+vvel*tauyy-qy)*ny )
   flux(3)=MubyRey*(tauxx*nx+tauxy*ny ) 
   flux(4)=MubyRey*(tauxy*nx+tauyy*ny ) 

   cell(in)%res(:)=cell(in)%res(:)-flux(:)*area
   if(out<=noc) cell(out)%res(:)=cell(out)%res(:)+flux(:)*area
   if(bc==1001) then
      Fx1=Fx1-flux(3)*area
      Fy1=Fy1-flux(4)*area
   endif    
enddo


!do ie=startBC,endBC
!   in = fc(ie)%in
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
!subroutine sutherland
!use param
!use pri
!implicit none
!integer(kind=i4) ::i 
!real(kind=dp) :: sutherland_viscosity

!do i=1,noc
!   cell(i)%mul= sutherland_viscosity(i)
!enddo

!return
!end
!
!Sutherland viscosity formula
real(kind=dp) function sutherland_viscosity(i)
use param
implicit none
integer(kind=i4) ::i 
real(kind=dp) :: num, den

call con2prim(cell(i)%qc(:))
num  = 1.0_dp+T_infd/SCONST
den  =   t   +T_infd/SCONST
sutherland_viscosity = (num/den)*t**1.5d0

!return
end

end subroutine viscous_flux
