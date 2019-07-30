subroutine cirs 
! Method taken from the J. Blazek  unstructured 2d solver code
use commons
use grid
use param
implicit none
integer(kind=i4) :: i, j, ie,itirs,nitirs,in,out 
!integer(kind=i4) :: ncontr(noc) 
!real(kind=dp),parameter :: epsirs=0.5_dp  
real(kind=dp)    :: epsirs
real(kind=dp)    :: OldRes(nvar,noc),ResItr(nvar,noc),den

OldRes=0.0_dp
!ncontr=0
nitirs=2

epsirs=dmax1(0.25_dp*( (cfl/cfl_max)**2-1.0_dp),0.0_dp)
!epsirs=0.5_dp

do i=1,noc
   do j=1,nvar
      OldRes(j,i)=cell(i)%res(j) 
   enddo
enddo

! Jacobi iterations
do itirs=1,nitirs

   ResItr=0.0_dp

   do ie=startFC,endFC
      in = fc(ie)%in
      out = fc(ie)%out
      ResItr(:,in)=ResItr(:,in)+cell(out)%res(:) 
      ResItr(:,out)=ResItr(:,out)+cell(in)%res(:) 
   enddo 

   do i=1,noc
      !den=1.0_dp/(1.0_dp+epsirs*dble(ncontr(i)))
      den=1.0_dp/(1.0_dp+epsirs*dble(cell(i)%nc2c))
      cell(i)%res(:)=(ResItr(:,i)*epsirs+OldRes(:,i))*den 
   enddo
   
enddo


end subroutine cirs 
