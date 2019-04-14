module param
use data_type,only:i4,dp
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
parameter(boundary=-1)
parameter(interior=0)
parameter(solid=3)
parameter(outflow=4)
parameter(farfield=5)


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

