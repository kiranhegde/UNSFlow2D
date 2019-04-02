!-----------------------------------------------------------------------------
!.....Definition of some constants
!-----------------------------------------------------------------------------
subroutine math
use param
implicit none


GAMMA        = 1.4d0
GAMMA1       = GAMMA-1.0d0
GAS_CONST    = 1.0d0
PI         = 4.0d0*datan(1.0d0)

end
!-----------------------------------------------------------------------------
!.....Variables stored are primitive - density, u, v, pressure
!.....Initialize primitive variables to free stream values
!-----------------------------------------------------------------------------
subroutine initialize
use param
use grid
implicit none
integer(kind=i4)  :: j

q_inf = 1.d0
aoa     = aoa_deg*PI/180.0d0
qinf(2) = 1.d0
qinf(3) = dcos(aoa)*q_inf
qinf(4) = dsin(aoa)*q_inf
p_inf    = 1.d0/(gamma*m_inf*m_inf)
qinf(1) = p_inf/(gamma-1.d0) + 0.5d0
r_inf = qinf(2)
ent_inf = p_inf/(r_inf**gamma)
a_inf   = dsqrt(gamma*p_inf/r_inf)

u_inf=qinf(3)/qinf(2)
v_inf=qinf(4)/qinf(2)

!     Print some useful info
if(flow_type == 'inviscid') print*,'Euler computation'
if(flow_type == 'laminar')  print*,'Laminar Navier-Stokes computation'
if(flow_type == 'rans')print*,'Turbulent Navier-Stokes computation'
print*,'Free-stream values:'
write(*,'(5x, " Mach number =", f8.4)')m_inf
write(*,'(5x, " AOA         =", f8.4)')aoa_deg
write(*,'(5x, " u velocity  =", f8.4)')u_inf
write(*,'(5x, " v velocity  =", f8.4)')v_inf
write(*,'(5x, " Pressure    =", f15.6)')p_inf

do j=1,nop
   pt(j)%qc(:) = qinf(:)
enddo

end
!-----------------------------------------------------------------------------
! Save flow solution into a file. Conserved variables are saved.
!-----------------------------------------------------------------------------
subroutine save_flow
use grid
implicit none
integer(kind=i4)  :: i, j

open(unit=50, file='FLO.DAT')
do i=1,nop
   write(50,'(4e20.10)') (pt(i)%qc(j), j=1,nvar)
enddo
close(50)

end
!-----------------------------------------------------------------------------
! Read flow solution from file
!-----------------------------------------------------------------------------
subroutine read_flow
use grid
implicit none
integer(kind=i4)  :: i, j

open(unit=50, file='FLO.DAT', status='OLD')
do i=1,nop
   read(50,*)(pt(i)%qc(j), j=1,nvar)
enddo
close(50)

end
!-----------------------------------------------------------------------------
! Save old solution
!-----------------------------------------------------------------------------
subroutine save_old
use param
use grid
use pri
implicit none
integer(kind=i4)  :: i, j

do i=1,nop
      call con2prim(pt(i)%qc(:))
      pt(i)%qp(:)=prim(:)
do j=1,nvar
      pt(i)%qold(j) = pt(i)%qc(j)
enddo
enddo


end
!-----------------------------------------------------------------------------
! L2 and Linf norm of the finite volume residual
!-----------------------------------------------------------------------------
subroutine residue
use param
use grid
implicit none
integer(kind=i4)  :: i, j
real(kind=dp) :: fr,dt

! Store current residual norm into fres_old for gmres
fres_old = fres

fres  = 0.0d0
fresi = 0.0d0
iresi = 0

!do i=1,nop
!      fresi = pt(i)%qc(2)/pt(i)%qold(2) - 1.d0
!      dt=pt(i)%dt
!      fres = fres + fresi*fresi/(dt*dt)
!enddo
!     fres=dsqrt(fres)/nop
!return

do i=1,nop
   fr = 0.0d0
   do j=1,nvar
      fr = fr + pt(i)%res(j)**2
   enddo
   fres = fres + fr
   fr   = dsqrt(fr)
   if(fr .gt. fresi)then
      fresi = fr
      iresi = i
   endif
enddo

fres = dsqrt(fres/nop)

if(iter .eq. 1)then
   fres1 = fres
   print*,'Residue in first iteration =',fres1
endif

if(fres1 .ne. 0.0d0) fres = fres/fres1



end

!-----------------------------------------------------------------------------

!!.....Error function, from Abromovitz and Stegun
!!-----------------------------------------------------------------------------
!double precision function ERRF(X)
!double precision X,ARG,E,VB,T,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
!
!ARG = X*X
!if(ARG .lt. 20.0d0)then
!      E = exp(-ARG)
!else
!      E = 0.0d0
!endif
!VB = abs(X)
!T = 1.0d0/(1.0d0 + 0.3275911d0*VB)
!tmp1 = 1.061405429d0*T
!tmp2 = (tmp1 - 1.453152027d0)*T
!tmp3 = (tmp2 + 1.421413741d0)*T
!tmp4 = (tmp3 - 0.284496736d0)*T
!tmp5 = (tmp4 + 0.254829592d0)*T
!tmp6 = 1.0d0 - tmp5*E
!if(X .lt. 0.0d0)then
!      ERRF = -tmp6
!else
!      ERRF =  tmp6
!endif
!end
!!-----------------------------------------------------------------------------
!!.....Prints data into a file for visualizing contours using gnuplot
!!.....Taken from NSC2KE of Bijan Mohammadi
!!-----------------------------------------------------------------------------
SUBROUTINE isocont(ifile,F,COOR,NVAL,VAL)
implicit none
integer          ifile, nval
double precision F(3),COOR(2,3),VAL(100)

integer          IP1(3), ival, itr, k
double precision epsi, ff1, ff2, ff3, ffma, d12, d23, val1, fk, &
                 fk1, fmi, fma, dif, eps, hh, x, y
!
epsi   = 1.0d-5
IP1(1) = 2
IP1(2) = 3
IP1(3) = 1
FF1    = F(1)
FF2    = F(2)
FF3    = F(3)
FFMA   = DMAX1(DABS(FF1),DABS(FF2))
FFMA   = DMAX1(ffma,DABS(FF3))
D12    = DABS(FF1-FF2)
D23    = DABS(FF2-FF3)
IF(D12+D23.LT.DMAX1(epsi,epsi*FFMA)) GOTO 1000
!  PAS DE RESTRICTION
!  ******************    
DO 100 IVAL=1,NVAL
      VAL1 = VAL(IVAL)
      ITR  = 0
      DO 110 K=1,3
            FK  = F(K)
            FK1 = F(IP1(K))
            FMI = DMIN1(FK,FK1)
            FMA = DMAX1(FK,FK1)
            DIF = FMA-FMI
            IF(DIF.LT.epsi) GOTO 110
            EPS = epsi*DIF
            IF(VAL1.LT.FMI-EPS .OR. VAL1.GT.FMA+EPS) GOTO 110
            HH  = DABS(FK-VAL1)/DIF
            X   = COOR(1,K) + HH*(COOR(1,IP1(K))-COOR(1,K))
            Y   = COOR(2,K) + HH*(COOR(2,IP1(K))-COOR(2,K))
            IF(ITR.EQ.0) GOTO 115
            write(ifile,*) x,y
            write(ifile,*) 
            GOTO 100
115               ITR = 1
            write(ifile,*) x,y
110         CONTINUE
100   CONTINUE
1000  return      
END



!==============================================================================
real(kind=dp) function m1p(mach)
use data_type
implicit none

real(kind=dp):: mach

m1p=0.50*(mach+dabs(mach))

end function m1p
!====================================
real(kind=dp) function m1m(mach)
use data_type
implicit none

real(kind=dp):: mach

m1m=0.5d0*(mach-dabs(mach))

end function m1m
!====================================
real(kind=dp) function m2p(mach)
use data_type
implicit none

real(kind=dp):: mach

m2p=0.25d0*(mach+1.d0)*(mach+1.d0)

end function m2p
!====================================
real(kind=dp) function m2m(mach)
use data_type
implicit none

real(kind=dp):: mach

m2m=-0.25d0*(mach-1.d0)*(mach-1.d0)

end function m2m
!====================================
real(kind=dp) function m4p(mach,beta)
use data_type
implicit none

real(kind=dp):: mach,beta,m2m,m2p

if(dabs(mach)>=1.d0) then
   m4p=0.5d0*(mach+dabs(mach))
else
   m4p=m2p(mach)*(1.d0-16.d0*beta*m2m(mach))
endif

end function m4p
!====================================
real(kind=dp) function m4m(mach,beta)
use data_type
implicit none

real(kind=dp):: mach,beta,m2m,m2p

if(dabs(mach)>=1.d0) then
   m4m=0.5d0*(mach-dabs(mach))
else
   m4m=m2m(mach)*(1.d0+16.d0*beta*m2p(mach))
endif

end function m4m
!====================================
real(kind=dp) function p5p(mach,alpha)
use data_type
implicit none

real(kind=dp):: mach,alpha,m2m,m2p

if(dabs(mach)>=1.d0) then
   p5p=0.5d0*(mach+dabs(mach))/mach
else
   p5p=m2p(mach)*( (2.d0-mach)-16.d0*alpha*mach*m2m(mach))
endif

end function p5p
!====================================
real(kind=dp) function p5m(mach,alpha)
use data_type
implicit none

real(kind=dp):: mach,alpha,m2m,m2p

if(dabs(mach)>=1.d0) then
   p5m=0.5d0*(mach-dabs(mach))/mach
else
   p5m=m2m(mach)*( (-2.d0-mach)+16.d0*alpha*mach*m2p(mach))
endif

end function p5m

!====================================

function dotprod(n1,n2,xc,yc)
   use grid
   implicit none
   integer(kind=i4):: n1,n2
   real(kind=dp)   :: x1,y1,x2,y2
   real(kind=dp)   :: dx,dy,ds,xc,yc
   real(kind=dp)   :: dotprod,nx,ny,rx,ry


   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y
   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds

   dx = xc-x1 ; dy = yc-y1
   ds = dsqrt(dx*dx+dy*dy)
   rx = dx/ds
   ry = dy/ds
   dotprod = rx*nx+ry*ny

end function dotprod


function crossprod(n1,n2,m1,m2)
   use grid
   implicit none
   integer(kind=i4):: n1,n2,m1,m2
   real(kind=dp)   :: x1,y1,x2,y2
   real(kind=dp)   :: x3,y3,x4,y4
   real(kind=dp)   :: dx,dy,ds
   real(kind=dp)   :: crossprod,nx,ny,rx,ry


   x1 = pt(n1)%x ; y1 = pt(n1)%y
   x2 = pt(n2)%x ; y2 = pt(n2)%y

   x3 = pt(m1)%x ; y3 = pt(m1)%y
   x4 = pt(m2)%x ; y4 = pt(m2)%y

   dx = x2-x1 ; dy = y2-y1
   ds = dsqrt(dx*dx+dy*dy)
   nx = dy/ds ; ny = -dx/ds

   dx = x4-x3 ; dy = y4-y3
   ds = dsqrt(dx*dx+dy*dy)
   rx = dx/ds
   ry = dy/ds

   crossprod = nx*ry-ny*rx

end function crossprod

!------------------------------------------------------------------------------
