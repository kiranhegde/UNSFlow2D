! Data size
module data_type
implicit none
!integer, parameter :: dp = kind(1.0d0)
integer, parameter :: i1=selected_int_kind(2)
integer, parameter :: i2=selected_int_kind(4)
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(18)
integer, parameter :: sp=selected_real_kind(6,37)
integer, parameter :: dp=selected_real_kind(15,307)
integer, parameter :: qp=selected_real_kind(31,307)
! Small number
!real(kind=dp) :: EPSILON
!parameter(EPSILON=1.0d-16)
real(kind=sp),parameter::eps=epsilon(1.0)
integer(kind=i4), parameter :: ndim=2
end module data_type

! Common
module commons
use data_type
implicit none
integer(kind=i4), parameter :: no=0
integer(kind=i4), parameter :: yes=1

! viscosity
integer(kind=i4), parameter :: nvar=4

!     GAMMA = ratio of specific heats
!     GAMMA1= GAMMA - 1
!     GAS_CONST = gas constant, this can be set to 1.0
!     M_PI = value of pi
real(kind=dp) :: GAMMA, GAMMA1, GAS_CONST, PI
end module commons

module inf
use data_type
use commons
implicit none
! Freestream values
real(kind=dp) :: q_inf,m_inf, aoa, aoa_deg, u_inf, v_inf,  &
            r_inf, p_inf, T_inf, T_infd, ent_inf,  &
            conv_inf(nvar), qinf(nvar),a_inf, H_inf

! Vortex correction for farfield BC
integer(kind=i4)  :: vortex
real(kind=dp) :: xref, yref
end module inf

module visc
use data_type
use commons
implicit none
! laminar paramters
! Rey = Reynolds number
! SCONST = Constant in Sutherland Law
character(len=24)  :: flow_type
real(kind=dp) :: Rey, prandtl, prandtl_turb, SCONST

! Parameters in Spallart-Allmaras model
real(kind=dp) :: Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1,  &
            Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2

end module visc


module param
use data_type
use commons
use visc
use inf
implicit none

integer(kind=i4) istart, scratch, restart
parameter(scratch=1, restart=2)

! Range of bounding box
real(kind=dp) :: xmin, xmax, ymin, ymax

! Type of grid
! any other value implies hybrid grid
character gridfile*64, inpfile*32

real(kind=dp) :: CFL,cfl_max,MINRES, dtglobal, gerrtol
character(len=24) :: timemode
integer(kind=i4)  :: iter, ITERLAST, MAXITER, saveinterval, &
            gmaxiter, prectype, scrinterval

real(kind=dp),allocatable :: airk(:), birk(:)
integer(kind=i4)  :: NIRK

!     Number of contours
integer(kind=i4)  :: niso

! Range of some variables
! rmin,rmax = density
! pmin,pmax = pressure
! mmin,mmax = mach number
real(kind=dp) :: rmin, rmax, umin, umax, vmin, vmax, pmin, pmax, &
            mmin, mmax, emin, emax, nmin, nmax

real(kind=dp) :: fres, fres_old, fres1, fresi
integer(kind=i4)  :: iresi

!real(kind=dp) :: wd1(nspmax), wdmin, wdmax

!     Define tags for point types
integer(kind=i4)  :: interior, solid, farfield, outflow, boundary
parameter(interior=0)
parameter(solid=3)
parameter(outflow=4)
parameter(farfield=5)
parameter(boundary=-1)


! Limiter factor for MUSCL
integer(kind=i4)  :: GRADTYPE, ILIMIT
real(kind=dp) :: LFACT, ALBADA11, ALBADA12, ALBADA21, ALBADA22

! Size of connectivity list; required by mayavi
integer(kind=i4)  :: lsize

character(len=24) :: flux_type

integer(kind=i4)  :: inviscid, laminar, turbulent
parameter(inviscid=1, laminar=2, turbulent=3)

integer(kind=i4)  :: xvstatus, display

real(kind=dp) :: minelarea, maxelarea, mincvarea, maxcvarea
real(kind=dp) :: minflen, maxflen

end module param

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
!---------------------------------------------------------



!==============================================
module grid
use data_type
use commons
implicit none

integer(kind=i4):: nop,nof,noc,noc_bc=0

type points
     integer(kind=i4):: nv2v,nv2c,nv2f
     real(kind=dp)   :: x,y,z
     real(kind=dp)   :: sx,sy,dt,la,ls,cv,ds
     real(kind=dp)   :: grad(1:ndim,1:nvar)=0.d0 
     real(kind=dp)   :: qp(1:nvar)=0.d0,qc(1:nvar)=0.d0,qold(1:nvar)=0.d0,res(1:nvar)=0.d0
     real(kind=dp)   :: DUmax(1:nvar)=0.d0,DUmin(1:nvar)=0.d0,phi(1:nvar)=1e20
     integer(kind=i4):: bc,flag
     integer(kind=i4),dimension(:),pointer::v2c
     integer(kind=i4),dimension(:),pointer::v2f
     integer(kind=i4),dimension(:),pointer::v2v
     real(kind=dp),dimension(:),pointer::wt
     real(kind=dp),dimension(:),pointer::wx,wy
end type points

type faces
     integer(kind=i4):: pt(ndim)
     integer(kind=i4):: in
     integer(kind=i4):: out
     integer(kind=i4):: bc
     integer(kind=i4):: flag
     real(kind=dp)   :: grad(1:ndim,1:nvar)=0.d0
     real(kind=dp)   :: nx,ny,area,cov,la
     real(kind=dp)   :: cen(1:ndim)=0.d0 ! face center
end type faces

type cells
     integer(kind=i4):: nc2v,nc2f,nc2c
     real(kind=dp)   :: cen(1:ndim)=0.d0,cv
     !real(kind=dp)   :: cov,phi(1:nvar)=1e20,ds
     !real(kind=dp)   :: dx,dy,dt,la,ls
     !real(kind=dp)   :: r11,r12,r22 
     !real(kind=dp)   ::  cv,cov 
     !real(kind=dp)   :: qp(1:nvar)=0.d0,qc(1:nvar)=0.d0,qold(1:nvar)=0.d0,res(1:nvar)=0.d0
     !real(kind=dp)   :: DUmax(1:nvar)=0.d0,DUmin(1:nvar)=0.d0
     !real(kind=dp)   :: grad(1:ndim,1:nvar)=0.d0
     integer(kind=i4),dimension(:),pointer::c2v
     integer(kind=i4),dimension(:),pointer::c2f
     integer(kind=i4),dimension(:),pointer::c2c
end type cells

type(points),allocatable,dimension(:)::pt
type(points),allocatable,dimension(:,:)::mesh
type(faces),allocatable,dimension(:)::fc
type(cells),allocatable,dimension(:)::cell

end module grid


module tools
use data_type
!use commons
implicit none

type coord
     real(kind=dp)   :: x,y
end type coord

type vector
     type(coord)::s,e
end type vector



contains

real(kind=dp) function norm_dist(a,p)

real(kind=dp)   :: NR,DR
type (vector) :: a
type (coord)  :: p

NR=(a%s%y-a%e%y)*p%x+(a%e%x-a%s%x)*p%y+(a%s%x*a%e%y-a%e%x*a%s%y)
DR=dsqrt( (a%e%x-a%s%x)**2+(a%e%y-a%s%y)**2 )
norm_dist=dabs(NR)/DR
end function norm_dist

end module tools

