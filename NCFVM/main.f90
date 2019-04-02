!------------------------------------------------------------------------------
! Main program of Node-centered finite volume flow solver 
!------------------------------------------------------------------------------
program main
use param
use grid
implicit none
integer(kind=i4)  :: i, j, irk,rk
real(kind=dp) :: c1(nvar), c2(nvar), c3(nvar),con(nvar)
real(kind=dp) :: cl, cd,dtbyarea,alfa
real(kind=sp) :: etime, elapsed(2), totaltime
real(kind=dp),allocatable:: ggqx(:,:),ggqy(:,:) 


!    Set some constants
call math

!    Read input parameters from file, also sets np, nt
call read_input

!    Read grid from file
call read_grid

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

iter = 0
fres = 1.0d0
call system('rm -f FLO.RES')
print*,'Beginning of iterations ...'
open(unit=99, file='FLO.RES')
do while(iter .lt. MAXITER .and. fres .gt. MINRES)

   iter = iter + 1
   cfl = 0.5d0*cfl_max*(1+dtanh(iter*7.d0/500.0 - 3.5))
   !cfl = 0.5d0*cfl_max*(1+dtanh(iter*10.d0/1000.0 - 5.0))
   !cfl=cfl_max
   call save_old

   if(timemode == 'rk3')then
      alfa=0.d0
      do irk=1,nirk
         call fvresidual
         do i=1,nop
            c1(:)=pt(i)%qold(:)
            c2(:)=pt(i)%qc(:)
            dtbyarea=pt(i)%dt/pt(i)%cv
            do j=1,nvar
               c3(j) = airk(irk)*c1(j) + &
               birk(irk)*(c2(j)-dtbyarea*pt(i)%res(j))
            enddo
            pt(i)%qc(:)=c3(:)
         enddo
!       alfa=1.d0/dble(nirk-irk+1)     
!       do i=1,nop
!          dtbyarea=pt(i)%dt/pt(i)%cv
!          do j=1,nvar
!             pt(i)%qc(j)=pt(i)%qold(j)-alfa*dtbyarea*pt(i)%res(j)
!          enddo
!       enddo
      enddo
   elseif(timemode == 'lusgs')then
      call lusgs
   endif

   call residue
   write(99,'(i6,3e16.6)') iter, dlog10(fres),dlog10(fresi), cfl

   if(mod(iter,scrinterval) .eq. 0)then
      call tecplt 
      !call tecplt_del 
      call screen
   endif
enddo
close(99)

call tecplt 
call tecplt_del 


totaltime = etime(elapsed)
totaltime = totaltime/60.0
elapsed(1)= elapsed(1)/60.0
elapsed(2)= elapsed(1)/60.0
print *, 'Time: total=', totaltime, ' user=', elapsed(1), &
         ' system=', elapsed(2)

100 format(1x,i6,1x,4(f15.8,1x))



stop
end


subroutine screen
use param
use pri
use grid
implicit none
real(kind=dp) :: qc(nvar),con(nvar), cl, cd

integer(kind=i4)  :: is, i
real(kind=dp) :: q2, mach, ent 
 
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
do is=1,nop
    con(:)=pt(is)%qc(:)
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
enddo

100 format(1x,i6,1x,4(f15.6,1x))
write(*,9)('-',i=1,70)
9     format(70a)
write(*,10)m_inf,aoa_deg,Rey,CFL
10    format(' Mach =',f6.3,'    AOA =',f6.2, '  Rey = ',e15.4, &
       '   CFL = ',e15.4)
write(*,11)flux_type,ilimit,gridfile
11    format(' Flux = ',a7, '  Lim = ',i2,'     Grid= ',a30)
write(*,9)('-',i=1,70)
write(*,'(" Iterations        =",i12)')iter
write(*,'(" Global dt         =",e16.6)') dtglobal
write(*,'(" L2 residue        =",e16.6)')fres
write(*,'(" Linf residue      =",e16.6)')fresi
write(*,'(" Linf cell         =",i12)') iresi
write(*,*)
!write(*,'(" Cl, Cd            =",2f12.6)')cl,cd
write(*,*)
write(*,'(27x,"Min",8x,"Max")')          
write(*,'(" Density           =",2f12.6)')rmin, rmax
write(*,'(" Pressure          =",2f12.6)')pmin, pmax
write(*,'(" Mach number       =",2f12.6)')mmin, mmax
write(*,'(" x velocity        =",2f12.6)')umin, umax
write(*,'(" y velocity        =",2f12.6)')vmin, vmax
write(*,'(" Entropy           =",2f12.6)')emin, emax
write(*,*)
write(*,*)
write(*,*)
write(*,*)
call flush(6)

end


SUBROUTINE tecplt
use param
use grid
implicit none
integer(kind=i4) :: i,j,c, funit
real(kind=dp) :: con(nvar),qc(nvar)
real(kind=dp) :: q2, mach,entropy,mul,mutot,sutherland
real(kind=dp) :: ro , uo,vo,po,wt,x1,x2,y1,y2
character(len=15) :: ctype

if(maxval(cell(:)%nc2f) == 3) ctype="TRIANGLE"
if(maxval(cell(:)%nc2f) == 4) ctype="quadrilateral"

funit = 5
!call avg_c2v

OPEN(funit,file='tecplt.dat')
WRITE(funit,*) 'TITLE = "flo2d output" '

!WRITE(funit,*) 'VARIABLES="X","Y","rho","u","v","p" '
WRITE(funit,*) 'VARIABLES="X","Y","rho","u","v","p","M"'
!WRITE(funit,*) 'ZONE F=FEPOINT,ET=quadrilateral'
WRITE(funit,*) 'ZONE F=FEPOINT,ET=',trim(ctype)
WRITE(funit,*) 'N=',nop,',E=',noc
do i=1,nop

    ro=pt(i)%qp(1)
    uo=pt(i)%qp(2)
    vo=pt(i)%qp(3)
    po=pt(i)%qp(4)

    q2   = uo**2 + vo**2
    mach = dsqrt(q2*ro/(GAMMA*po))
    entropy=log10(po/ro**GAMMA/ent_inf)

    !WRITE(funit,*)pt(i)%x,pt(i)%y,r,u,v,p,mach
    WRITE(funit,100)pt(i)%x,pt(i)%y,ro,uo,vo,po,mach
      
END DO
DO i=1,noc
   !WRITE(funit,*)elem(1,i), elem(2,i), elem(3,i) 
   WRITE(funit,*)(cell(i)%c2v(j),j=1,cell(i)%nc2v) 
END DO

CLOSE(funit)

<<<<<<< HEAD
!DO i=1,noc
!   do j=1,cell(i)%nc2v
!      c=cell(i)%c2v(j)
!      x2 = pt(c)%x ; y2 = pt(c)%y
!      write(34,*)x2,y2
!   END DO
!      x1 = pt(cell(i)%c2v(1))%x ; y1 = pt(cell(i)%c2v(1))%y
!      write(34,*)x1,y1
!      write(34,*)
!END DO
=======
DO i=1,noc
   do j=1,cell(i)%nc2v
      c=cell(i)%c2v(j)
      x2 = pt(c)%x ; y2 = pt(c)%y
      write(34,*)x2,y2
   END DO
      x1 = pt(cell(i)%c2v(1))%x ; y1 = pt(cell(i)%c2v(1))%y
      write(34,*)x1,y1
      write(34,*)
END DO
>>>>>>> a688afb70276d1c80f99d1548c4045d94d229af4

100 format(1x,7(f15.6,2x))


END SUBROUTINE


SUBROUTINE tecplt_del
use param
use grid
implicit none
integer(kind=i4) :: i,j,c, funit
real(kind=dp) :: con(nvar),qc(nvar)
real(kind=dp) :: q2, mach,entropy,mul,mutot,sutherland
real(kind=dp) :: ro , uo,vo,po,wt,x1,x2,y1,y2

funit = 5
!call avg_Grad_c2v

OPEN(funit,file='tecplt_del.dat')
WRITE(funit,*) 'TITLE = "flo2d derivatives" '
WRITE(funit,*)'VARIABLES="X","Y","rhodx","rhody","udx","udy","vdx","vdy","pdx","pdy"'
WRITE(funit,*) 'ZONE F=FEPOINT,ET=quadrilateral'
WRITE(funit,*) 'N=',nop,',E=',noc
do i=1,nop
    WRITE(funit,100)pt(i)%x,pt(i)%y,pt(i)%grad(1,1),pt(i)%grad(2,1),pt(i)%grad(1,2),pt(i)%grad(2,2),pt(i)%grad(1,3),pt(i)%grad(2,3),pt(i)%grad(1,4),pt(i)%grad(2,4)
END DO
DO i=1,noc
   !WRITE(funit,*)elem(1,i), elem(2,i), elem(3,i) 
   WRITE(funit,*)(cell(i)%c2v(j),j=1,cell(i)%nc2v) 
END DO

CLOSE(funit)

DO i=1,noc
   do j=1,cell(i)%nc2v
      c=cell(i)%c2v(j)
      x2 = pt(c)%x ; y2 = pt(c)%y
      write(34,*)x2,y2
   END DO
      x1 = pt(cell(i)%c2v(1))%x ; y1 = pt(cell(i)%c2v(1))%y
      write(34,*)x1,y1
      write(34,*)
END DO

100 format(1x,7(f15.6,2x))


END SUBROUTINE
