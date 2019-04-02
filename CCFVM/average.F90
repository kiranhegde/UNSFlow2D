
subroutine avg_c2v 
use grid 
use pri
use commons
implicit none
integer(kind=i4) :: i,j,c,k
real(kind=dp) :: con(nvar)
real(kind=dp) :: wt,cwt


do i=1,noc
      call con2prim(cell(i)%qc(1:nvar))
      do k=1,nvar
         cell(i)%qp(k)=prim(k)
      enddo
enddo

do i=1,nop

    wt=0.d0
    con(:)=0.d0
    do j=1,pt(i)%nv2c
       c=pt(i)%v2c(j)
       !cwt=1.d0/cell(c)%cv
       cwt=pt(i)%wt(j)
       wt=wt+cwt
       do k=1,nvar
          con(k)=con(k)+cell(c)%qp(k)*cwt 
       enddo
    enddo
    do k=1,nvar
       pt(i)%prim(k)=con(k)/wt    
    enddo
end do


end subroutine avg_c2v

! Gradients
subroutine avg_Grad_c2v 
use grid 
use pri
use commons
implicit none
integer(kind=i4) :: i,j,k,c
real(kind=dp) :: grad(ndim,nvar)
real(kind=dp) :: wt,cwt

grad=0.0_dp


do i=1,nop

    wt=0.d0
    do j=1,pt(i)%nv2c
       c=pt(i)%v2c(j)
       !cwt=1.d0/cell(c)%cv
       cwt=pt(i)%wt(j)
       wt=wt+cwt
       do k=1,nvar
          grad(1:ndim,k)=grad(1:ndim,k)+cell(c)%grad(1:ndim,k)*cwt 
       enddo
    enddo
    pt(i)%grad(1:ndim,1:nvar)=grad(1:ndim,1:nvar)/wt    

end do

end subroutine avg_Grad_c2v 
