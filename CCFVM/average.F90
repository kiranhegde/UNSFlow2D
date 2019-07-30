subroutine avg_c2v 
use grid 
use pri
use commons
implicit none
integer(kind=i4) :: i,j,c,k,p1,p2
real(kind=dp) :: con(npvar)
real(kind=dp) :: wt,cwt

do i=1,nop

    wt=0.d0
    con(:)=0.d0
    do j=1,pt(i)%nv2c
       c=pt(i)%v2c(j)
       !cwt=1.d0/cell(c)%cv
       cwt=pt(i)%wt(j)
       wt=wt+cwt
       do k=1,npvar
          con(k)=con(k)+cell(c)%qp(k)*cwt 
       enddo
    enddo
    do k=1,npvar
       pt(i)%prim(k)=con(k)/wt    
    enddo
end do

! face average
do i=1,nof
   p1=fc(i)%pt(1)   
   p2=fc(i)%pt(2)
   fc(i)%qp(:)=0.5d0*(pt(p1)%prim(:)+pt(p2)%prim(:))
end do


end subroutine avg_c2v

! Gradients
subroutine avg_Grad_c2v 
use grid 
use pri
use commons
implicit none
integer(kind=i4) :: i,j,k,c
real(kind=dp) :: grad(ndim,ngrad)
real(kind=dp) :: wt,cwt

grad=0.0_dp


do i=1,nop

    wt=0.d0
    do j=1,pt(i)%nv2c
       c=pt(i)%v2c(j)
       !cwt=1.d0/cell(c)%cv
       cwt=pt(i)%wt(j)
       wt=wt+cwt
       do k=1,ngrad
          grad(1:ndim,k)=grad(1:ndim,k)+cell(c)%grad(1:ndim,k)*cwt 
       enddo
    enddo
    pt(i)%grad(1:ndim,1:ngrad)=grad(1:ndim,1:ngrad)/wt    

end do

end subroutine avg_Grad_c2v 
