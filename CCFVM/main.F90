!------------------------------------------------------------------------------
! Main program of Cell-centroid finite volume flow solver 
!------------------------------------------------------------------------------
program main
use param
use pri
use grid
implicit none
integer(kind=i4)  :: i, j, irk
real(kind=dp) :: c1(nvar), c2(nvar), c3(nvar)
real(kind=dp) :: dtbyarea,alfa
real(kind=sp) :: etime, elapsed(2), totaltime,ts,te
!real(kind=dp),allocatable:: ggqx(:,:),ggqy(:,:) 

elapsed=0.0_sp
ts=0.0_sp
te=0.0_sp

!    Set some constants
call math

!    Read input parameters from file, also sets np, nt
call read_input

!    Generating mesh connectivity 
call geometric

! Set initial condition
call initialize

!     Runge-Kutta time stepping
if(timemode=="rk3".and.cfl_max > 1.0d0) then
   print*,'CFL should be less than 1 for rk3 !!!'
   stop
elseif(timemode=="rk3".and.cfl_max < 1.0d0) then
print*,'RK3'
NIRK    = 3
allocate(airk(nirk),birk(nirk))
airk(1) = 0.0d0
airk(2) = 3.0d0/4.0d0
airk(3) = 1.0d0/3.0d0
birk(1) = 1.0d0
birk(2) = 1.0d0/4.0d0
birk(3) = 2.0d0/3.0d0
endif
!print*,'eps=',eps
!stop


iter = 0
fres = 1.0d0
call system('rm -f FLO.RES')

print*,'Beginning of iterations ...'
open(unit=99, file='FLO.RES')

do while(iter .lt. MAXITER .and. fres .gt. MINRES)

   iter = iter + 1
   call cfl_number(CFL)

   do i=1,noc
      call con2prim(cell(i)%qc(1:nvar))
      do j=1,npvar
         cell(i)%qp(j)=prim(j)
      enddo
   enddo

   call avg_c2v
   !call time_step2
   call time_step02
   call save_old
   
   !Gradient calculation using Green-Gauss
   if(trim(grad_type)=='gg') then
      call Gradient_GG
   !Gradient calculation using GG diamond path reconstruction
   elseif(trim(grad_type)=='ggfc') then 
      call Gradient_GG_FC
   !Least square based Gradient calculation 
   elseif(trim(grad_type)=='lsqr') then
      call Gradient_LSQR
   else
      print*,' Unknown  Gradient calculation method:',grad_type
      stop
   endif
     
   if(timemode == 'rk3')then
      alfa=0.d0
      do irk=1,nirk
         call fvresidual
         do i=1,noc
            c1(:)=cell(i)%qold(:)
            c2(:)=cell(i)%qc(:)
            dtbyarea=cell(i)%dt/cell(i)%cv
            do j=1,nvar
               c3(j) = airk(irk)*c1(j) + &
               birk(irk)*(c2(j)-dtbyarea*cell(i)%res(j))
            enddo
            cell(i)%qc(:)=c3(:)
         enddo
!       alfa=1.d0/dble(nirk-irk+1)     
!       do i=1,noc
!          dtbyarea=cell(i)%dt/cell(i)%cv
!          do j=1,nvar
!             cell(i)%qc(j)=cell(i)%qold(j)-alfa*dtbyarea*cell(i)%res(j)
!          enddo
!       enddo
      enddo
   elseif(timemode == 'lusgs')then
      call lusgs
   endif

   call check_positivity
   call residue
   write(99,'(i6,2e16.6)') iter, dlog10(fres),cfl

   if(mod(iter,scrinterval) .eq. 0)then
      call tecplt 
      !call tecplt_del 
      call screen
      te = etime(elapsed)
      te = te/60.0
      print *, 'Time taken :',te,te-ts
      ts=te   
   endif
   flush(6)
enddo
close(99)

call tecplt 
!call tecplt_del 

call check_positivity

totaltime = etime(elapsed)
totaltime = totaltime/60.0
elapsed(1)= elapsed(1)/60.0
elapsed(2)= elapsed(1)/60.0
print *, 'Time: total=', totaltime, ' user=', elapsed(1), &
         ' system=', elapsed(2)

!100 format(1x,i6,1x,4(f15.8,1x))


contains

subroutine cfl_number(cfl_no)
use param
real(kind=dp)   :: cfl1,cfl2,exp_factor,s,cfl_no

   ! Hyperblocally increasing the cfl to cfl_max
   if(cfl_type=='tanh') then 
       !cfl_no = 0.5_dp*cfl_max*(1.0_dp+dtanh(iter*7.0_dp/500.0_dp - 3.5_dp))
       s = real(iter,dp)/real(CFL_ramp_steps-1,dp)
       cfl_no = 0.5_dp*cfl_max*(1.0_dp+dtanh(s*7.0_dp - 3.5_dp))
       !cfl_no = 0.5_dp*cfl_max*(1+dtanh(s*10.0_dp - 5.0_dp))
   elseif(cfl_type=='ramp'.and.iter<=CFL_ramp_steps) then 
   ! Ramping the cfl number to cfl_max in steps
       CFL1 = CFL_min  ! Initial CFL 
       CFL2 = CFL_max  ! Final CFL 
       exp_factor = 0.5_dp
       !s = real(iter,dp)/real(CFL_ramp_steps-1,dp)
       s = real(iter,dp)/real(CFL_ramp_steps,dp)
       CFL_no = CFL1 + (CFL2-CFL1)*(1.0_dp-dexp(-s*exp_factor))/(1.0_dp-dexp(-exp_factor))
   endif

end subroutine cfl_number

end  program


subroutine screen
use param
use pri
use grid
use output
implicit none
real(kind=dp) :: con(nvar)
integer(kind=i4)  :: is, i
real(kind=dp) :: q2, mach, ent 


call find_cl_cd
 
!call  check_positivity

pmin =  1.0d10
pmax = -1.0d10
mmin =  1.0d10
mmax = -1.0d10
rmin =  1.0d10
rmax = -1.0d10
umin =  1.0d10
umax = -1.0d10
vmin =  1.0d10
vmax = -1.0d10
emin =  1.0d10
emax = -1.0d10
tmin =  1.0d10
tmax = -1.0d10
do is=1,noc
    con(:)=cell(is)%qc(:)
    call con2prim(con)
    q2   = u*u  + v*v 
    rmin =dmin1(rmin, rho) 
    rmax =dmax1(rmax, rho) 
    umin =dmin1(umin, u) 
    umax =dmax1(umax, u) 
    vmin =dmin1(vmin, v) 
    vmax =dmax1(vmax, v) 
    pmin =dmin1(pmin, p) 
    pmax =dmax1(pmax, p) 
    mach =dsqrt(q2*rho/(GAMMA*p))
    mmin =dmin1(mmin, mach)
    mmax =dmax1(mmax, mach)
    ent  =dlog10(p/rho**GAMMA/ent_inf)
    emin =dmin1(emin, ent)
    emax =dmax1(emax, ent)
    tmin =dmin1(tmin, t)
    tmax =dmax1(tmax, t)
enddo

!100 format(1x,i6,1x,4(f15.6,1x))
write(*,9)('-',i=1,75)
if (flow_type=="inviscid") then
    write(*,10)m_inf,aoa_deg,CFL
else
   write(*,13)m_inf,aoa_deg,CFL,Rey
endif
write(*,11)irs,ilimit, trim(grad_type), trim(gridfile)
write(*,12)trim(timemode),trim(flow_type),trim(flux_type)
write(*,9)('-',i=1,75)
write(*,'(" Iterations        =",i12)')iter
write(*,'(" Global dt         =",e16.6)') dtglobal
write(*,'(" L2 residue        =",e16.6)')fres
write(*,'(" Linf residue      =",e16.6)')fresi
write(*,'(" Linf cell         =",i12)') iresi
!write(*,*)
write(*,'(" Cl, Cd, Cm            =",3f12.6)')cl,cd,cm
write(*,*)
write(*,'(27x,"Min",8x,"Max")')          
write(*,'(" Density           =",2f12.6)')rmin, rmax
write(*,'(" x velocity        =",2f12.6)')umin, umax
write(*,'(" y velocity        =",2f12.6)')vmin, vmax
write(*,'(" Pressure          =",2f12.6)')pmin, pmax
write(*,'(" Temparature       =",2f12.6)')tmin, tmax
write(*,'(" Mach number       =",2f12.6)')mmin, mmax
write(*,'(" Entropy           =",2f12.6)')emin, emax
write(*,*)
!write(*,*)
!write(*,*)
!write(*,*)
9     format(75a)
10    format(' Mach =',f6.3,'    AOA =',f6.2, '   CFL = ',f8.2)
11    format(' CIRS = ',i2,2x,'  Limiter = ',i2,2x,' Gradient :',a5,2x,'Grid :',a15)
12    format(' TimeMode :',a5,2x,'FlowType :',a12,2x,' FluxScheme :',a7)
13    format(' Mach =',f6.3,'    AOA =',f6.2,'   CFL = ',f8.2,'  Rey = ',ES12.3)
flush(6)

contains

subroutine  find_cl_cd
use output
use inf
use commons
use visc

integer(kind=i4) :: ie,p1,p2
real(kind=dp)    :: pwall,xi,yi,nx,ny,cp,cf,dcx,dcy,xa,ya
real(kind=dp)    :: ux,uy,vx,vy,tx,ty,tw,tauxx,tauxy,tauyy,mu
real(kind=dp)    :: two,two_third,dyna,area

ux=0.0_dp     
uy=0.0_dp     
vx=0.0_dp     
vy=0.0_dp     
tauxx=0.0_dp     
tauxy=0.0_dp     
tauyy=0.0_dp     
!dyna=0.5_dp*r_inf*q_inf*q_inf
dyna=0.5_dp*gamma*p_inf*m_inf*m_inf
 
two=2.0_dp
two_third=two/3.0_dp
fx=0.0_dp
fy=0.0_dp
cm=0.0_dp
i=0
open(3,file='XYCpCf.dat') 
do ie=startBC,endBC
   if(fc(ie)%bc==1001) then
      i=i+1

      pwall=cell(fc(ie)%in)%qp(4)
      !pwall=fc(ie)%qp(4)

      p1=fc(ie)%pt(1) 
      p2=fc(ie)%pt(2) 
      xi=0.5_dp*(pt(p1)%x+pt(p2)%x)
      yi=0.5_dp*(pt(p1)%y+pt(p2)%y)
      xa=xi-xref
      ya=yi-yref
      nx=fc(ie)%sx
      ny=fc(ie)%sy
      area =fc(ie)%area
      nx=nx/area 
      ny=ny/area 
      cp=(pwall-p_inf)/dyna
      !cp=pwall
      dcx=cp*nx
      dcy=cp*ny
      fx=fx+dcx
      fy=fy+dcy

      cf=0.0_dp  
      if(flow_type/="inviscid") then 
         mu=fc(ie)%mu/Rey
         ux=fc(ie)%grad(1,2)
         uy=fc(ie)%grad(2,2)
         vx=fc(ie)%grad(1,3)
         vy=fc(ie)%grad(2,3)
         Tauxx=two_third*(two*ux-vy)
         Tauxy=(uy+vx)
         Tauyy=two_third*(two*vy-ux)
         fx=fx-mu*(Tauxx*nx+Tauxy*ny)
         fy=fy-mu*(Tauxy*nx+Tauyy*ny)

         !ds=dsqrt(nx*nx+ny*ny) 
         !nx=nx/ds
         !ny=ny/ds
         tx=ny
         ty=-nx
         tw = mu*((ux*tx + vx*ty)*nx + (uy*tx + vy*ty)*ny)
         cf = tw/dyna
!        us = dsqrt( dabs(tw)/rho )
!        ys = rho*us*wd1(i)/mu
      endif  

      cm=cm+dcx*ya-dcy*xa
      wbc(i)%cp=-cp
      wbc(i)%cf=-cf
      wbc(i)%x=xi
      wbc(i)%y=yi
      write(3,*)xi,yi,-cp,cf
   endif
enddo
close(3) 
cl=0.0_dp
cd=0.0_dp
fx=fx*area
fy=fy*area
fx1=fx1*area
fy1=fy1*area
print*,'Fx Fx1:',fx,fx1
print*,'Fy Fy1:',fy,fy1
cl = fy*dCos(aoa) - fx*dSin(aoa)
cd = fy*dSin(aoa) + fx*dCos(aoa)
print*,'cl,cd',cl/dyna,cd/dyna
cl = fy1*dCos(aoa) - fx1*dSin(aoa)
cd = fy1*dSin(aoa) + fx1*dCos(aoa)
print*,'cl1,cd1',cl/dyna,cd/dyna

end subroutine  find_cl_cd
end

SUBROUTINE check_positivity
use param
use pri
use grid
implicit none
integer(kind=i4) :: i,j,c
real(kind=dp) :: con(nvar),x1,y1

!i=1
!      print*
!      print*,'Density = ',pt(i)%prim(1)
!      print*,'u vel   = ',pt(i)%prim(2)
!      print*,'v vel   = ',pt(i)%prim(3)
!      print*,'Pressure= ',pt(i)%prim(4)
!      print*
!      print 200
!      do j=1,pt(i)%nv2c
!           c=pt(i)%v2c(j)
!           !write(3,*)cell(c)%cen(1),cell(c)%cen(2)
!           con(:)=cell(c)%qc(:)
!           call con2prim(con)
!           print 201,c,u,v,p,cell(c)%la
!      enddo



do i=1,noc
   if(cell(i)%qp(1).le. 0.0_dp .or.cell(i)%qp(4).le. 0.0_dp) then
   print*
   print*,'Density/pressure is negative at point ', i
   open(3,file='PositivtyCheck.dat') 
      do j=1,cell(i)%nc2v
           c=cell(i)%c2v(j)
           x1 = pt(c)%x ; y1 = pt(c)%y
           write(3,*) x1,y1
      enddo 
      c=cell(i)%c2v(1)
      x1 = pt(c)%x ; y1 = pt(c)%y
      write(3,*) x1,y1
      con(:)=cell(i)%qc(:)
      call con2prim(con)
      write(3,*)
      print*,'Density = ',rho
      print*,'u-vel   = ',u
      print*,'v-vel   = ',v
      print*,'Pressure= ',p
      print*,'La = ',cell(i)%la
   print*
   close(3)


   stop
   endif
enddo
!200 format(1x,'#',20x,'u',15x,'v',15x,'p',14x,'Spec')
!201 format(1x,i8,4(f15.4,1x))
end SUBROUTINE check_positivity


SUBROUTINE tecplt
use param
use grid
implicit none
integer(kind=i4) :: i,j,funit
real(kind=dp) :: q2, mach,entropy
real(kind=dp) :: ro,uo,vo,po,to,Vmag
character(len=15) :: ctype

if(maxval(cell(:)%nc2f) == 3) ctype="TRIANGLE"
if(maxval(cell(:)%nc2f) == 4) ctype="quadrilateral"

funit = 3
call avg_c2v

OPEN(funit,file='tecplt.dat')
WRITE(funit,*) 'TITLE = "flo2d output" '

!WRITE(funit,*) 'VARIABLES="X","Y","rho","u","v","p" '
WRITE(funit,*) 'VARIABLES="X","Y","density","u-velocity","v-velocity","Mach","Pressure","Temprature","Entropy","Velocity" ,"mu"'
WRITE(funit,*) 'ZONE F=FEPOINT,ET=quadrilateral'
!WRITE(funit,*) 'ZONE F=FEPOINT,ET=',trim(ctype)
!WRITE(funit,*) 'ZONE F=FEPOINT,ET=',
WRITE(funit,*) 'N=',nop,',E=',noc
do i=1,nop

    ro=pt(i)%prim(1)
    uo=pt(i)%prim(2)
    vo=pt(i)%prim(3)
    po=pt(i)%prim(4)
    to=pt(i)%prim(5)
    Vmag=dsqrt(uo*uo+vo*vo)

    q2   = uo**2 + vo**2
    mach = dsqrt(q2*ro/(GAMMA*po))
    entropy=log10(po/ro**GAMMA/ent_inf)

    !WRITE(funit,*)pt(i)%x,pt(i)%y,r,u,v,p,mach
    WRITE(funit,100)pt(i)%x,pt(i)%y,ro,uo,vo,mach,po,to,entropy,Vmag,pt(i)%mu
      
END DO
DO i=1,noc
   !WRITE(funit,*)elem(1,i), elem(2,i), elem(3,i) 
   if(cell(i)%nc2f==3) WRITE(funit,*)(cell(i)%c2v(j),j=1,cell(i)%nc2v-1),cell(i)%c2v(3) 
   if(cell(i)%nc2f==4) WRITE(funit,*)(cell(i)%c2v(j),j=1,cell(i)%nc2v-1)
END DO

CLOSE(funit)

!DO i=1,noc
!   do j=1,cell(i)%nc2v-1
!      c=cell(i)%c2v(j)
!      x2 = pt(c)%x ; y2 = pt(c)%y
!      write(34,*)x2,y2
!   END DO
!      x1 = pt(cell(i)%c2v(1))%x ; y1 = pt(cell(i)%c2v(1))%y
!      write(34,*)x1,y1
!      write(34,*)
!END DO

100 format(1x,7(f15.6,2x))


END SUBROUTINE


SUBROUTINE tecplt_del
use param
use grid
implicit none
integer(kind=i4) :: i,j,funit
character(len=15) :: ctype

funit = 4
call avg_Grad_c2v


if(maxval(cell(:)%nc2f) == 3) ctype="TRIANGLE"
if(maxval(cell(:)%nc2f) == 4) ctype="quadrilateral"

OPEN(funit,file='tecplt_del.dat')
WRITE(funit,*) 'TITLE = "flo2d derivatives" '
WRITE(funit,*)'VARIABLES="X","Y","rhodx","rhody","udx","udy","vdx","vdy","pdx","pdy"'
WRITE(funit,*) 'ZONE F=FEPOINT,ET=',trim(ctype)
WRITE(funit,*) 'N=',nop,',E=',noc
do i=1,nop
    WRITE(funit,100)pt(i)%x,pt(i)%y,pt(i)%grad(1,1),pt(i)%grad(2,1),pt(i)%grad(1,2),pt(i)%grad(2,2), & 
    &          pt(i)%grad(1,3),pt(i)%grad(2,3),pt(i)%grad(1,4),pt(i)%grad(2,4)
    WRITE(funit,100)pt(i)%x,pt(i)%y,pt(i)%grad(1,1),pt(i)%grad(2,1),pt(i)%grad(1,2), & 
         & pt(i)%grad(2,2),pt(i)%grad(1,3),pt(i)%grad(2,3),pt(i)%grad(1,4),pt(i)%grad(2,4)
END DO
DO i=1,noc
   !WRITE(funit,*)elem(1,i), elem(2,i), elem(3,i) 
   WRITE(funit,*)(cell(i)%c2v(j),j=1,cell(i)%nc2v-1) 
END DO

CLOSE(funit)


100 format(1x,10(f15.6,2x))


END SUBROUTINE
