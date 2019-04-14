!==============================================
module grid
use data_type
use commons
implicit none

integer(kind=i4):: nop,nof,noc,noc_bc=0,startBC=1,endBC=1,startfc=1,endfc=1,nbf 

type points
     real(kind=dp) :: x,y,z
     real(kind=dp) :: prim(1:nvar)=0.d0,grad(1:ndim,1:nvar)=0.d0
     integer(kind=i4):: bc,flag,nv2c
     !integer(kind=i4),dimension(:),pointer::v2c
     !real(kind=dp),dimension(:),pointer::wt
     integer(kind=i4),dimension(:),allocatable::v2c
     real(kind=dp),dimension(:),allocatable::wt
end type points

type faces
     integer(kind=i4):: pt(2)
     integer(kind=i4):: in=-1
     integer(kind=i4):: out=-1
     integer(kind=i4):: bc
     integer(kind=i4):: flag
     real(kind=dp)   :: grad(1:ndim,1:nvar)=0.d0
     real(kind=dp)   :: tgrad(1:ndim)=0.d0
     real(kind=dp)   :: sx,sy,cov,la
     real(kind=dp)   :: cen(1:ndim)=0.d0 ! face center
     real(kind=dp)   :: qp(1:nvar)=0.d0
end type faces

type cells
     integer(kind=i4):: nc2v,nc2f,nc2c
     real(kind=dp)   :: cen(1:ndim)=0.d0,cv,cov,phi(1:nvar),ds
     real(kind=dp)   :: dx,dy,dt,la,ls
     real(kind=dp)   :: r11,r12,r22,det
     real(kind=dp)   :: qp(1:nvar)=0.d0,qc(1:nvar)=0.d0,qold(1:nvar)=0.d0,res(1:nvar)=0.d0
     real(kind=dp)   :: DUmax(1:nvar)=0.d0,DUmin(1:nvar)=0.d0
     real(kind=dp)   :: grad(1:ndim,1:nvar)=0.d0
     real(kind=dp)   :: mul,mut
     !integer(kind=i4),dimension(:),pointer::c2v
     !integer(kind=i4),dimension(:),pointer::c2f
     !integer(kind=i4),dimension(:),pointer::c2c
     integer(kind=i4),dimension(:),allocatable::c2v
     integer(kind=i4),dimension(:),allocatable::c2f
     integer(kind=i4),dimension(:),allocatable::c2c
end type cells

type solution
     ! conserved variables
     real(kind=dp),dimension(:,:),allocatable   :: qc,qold 
     ! primitive variables
     real(kind=dp),dimension(:,:),allocatable   :: qp
     ! gradient of primitive variables
     ! rho,u,v,p,T, and turbulance
     real(kind=dp),dimension(:,:,:),allocatable :: dp
     ! residual            
     real(kind=dp),dimension(:,:),allocatable   :: res
     ! time step           
     real(kind=dp),dimension(:),allocatable     :: dt
     ! limiter            
     real(kind=dp),dimension(:,:),allocatable   :: phi
end type solution



type(points),allocatable,dimension(:)::pt
type(points),allocatable,dimension(:,:)::mesh
type(faces),allocatable,dimension(:)::fc
type(faces),allocatable,dimension(:)::bfc
type(cells),allocatable,dimension(:)::cell

end module grid

