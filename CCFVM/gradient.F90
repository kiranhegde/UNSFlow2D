!======================================================================================
subroutine Gradient_LSQR
use grid 
use commons
use param,only:ilimit
implicit none
integer(kind=i4) :: i,j,k,c
real(kind=dp)    :: xc,yc,dx,dy,wt
real(kind=dp)    :: wx,wy     
real(kind=dp)    :: r11,r12,r22,alfa1,alfa2


do i=1,noc
   cell(i)%grad(:,:)=0.d0
enddo

do i=1,noc
   xc=cell(i)%cen(1)
   yc=cell(i)%cen(2)
   r11=cell(i)%r11
   r12=cell(i)%r12
   r22=cell(i)%r22

  do k=1,nvar
   do j=1,cell(i)%nc2c
      c=cell(i)%c2c(j)
      dx=cell(c)%cen(1)-xc
      dy=cell(c)%cen(2)-yc

      alfa1=dx/r11/r11
      alfa2=(dy-dx*r12/r11)/r22/r22
      wx=alfa1-alfa2*r12/r11
      wy=alfa2

      cell(i)%grad(1,k)=cell(i)%grad(1,k)+(cell(c)%qp(k)-cell(i)%qp(k))*wx
      cell(i)%grad(2,k)=cell(i)%grad(2,k)+(cell(c)%qp(k)-cell(i)%qp(k))*wy
    enddo
   enddo

enddo

if(ILIMIT==1) call limit

end subroutine Gradient_LSQR

!======================================================================================

!======================================================================================
subroutine Gradient_GG
use grid 
use commons
use param,only:ilimit

implicit none
integer(kind=i4) :: i,j,k,in,out,c,p1,p2,ps
real(kind=dp) :: var,dx,dy,xc,yc 
real(kind=dp) :: ql(nvar),qr(nvar)
real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,D_L0,check,dotprod

do i=1,noc
   cell(i)%grad(:,:)=0.d0
enddo

do k=1,noc

   do i=1,cell(k)%nc2v-1
      c=cell(k)%c2v(i)
      p1=cell(k)%c2v(i)
      p2=cell(k)%c2v(i+1)
      dx=pt(p2)%y-pt(p1)%y 
      dy=-(pt(p2)%x-pt(p1)%x)
      do j=1,nvar
         var=0.5d0*(pt(p1)%prim(j)+pt(p2)%prim(j))
         cell(k)%grad(1,j)=cell(k)%grad(1,j)+var*dx
         cell(k)%grad(2,j)=cell(k)%grad(2,j)+var*dy
      enddo
   enddo

enddo

do i=1,noc
   var=cell(i)%cv
   do j=1,nvar
   cell(i)%grad(1,j)=cell(i)%grad(1,j)/var
   cell(i)%grad(2,j)=cell(i)%grad(2,j)/var
   enddo
enddo

if(ILIMIT==1) call limit

end subroutine Gradient_GG

!======================================================================================

subroutine Gradient_GG_FC
use grid 
use commons
use param,only:ilimit
implicit none
integer(kind=i4) :: i,j,k,ie,in,out,p1,p2,c
real(kind=dp) :: q2, pv(nvar)
real(kind=dp)    :: x1,y1,x2,y2,dx,dy

real(kind=dp) :: ql(nvar),qr(nvar)
real(kind=dp) :: nr, dr, phi, q_min, q_max,kappa
real(kind=dp) :: TOL,alfa,D_L,D_L0


do i=1,noc
   cell(i)%grad(:,:)=0.d0
enddo

do i=1,nof
      in = fc(i)%in
      out = fc(i)%out
      p1=fc(i)%pt(1)       
      p2=fc(i)%pt(2)       

      do j=1,nvar    
          fc(i)%grad(1,j)=0.d0
          fc(i)%grad(2,j)=0.d0
          if(in/=0.and.out/=0) then
                  pv(j)=0.5d0*(pt(p1)%prim(j)+cell(out)%qp(j)) 
                  x1 = pt(p1)%x    ; y1 = pt(p1)%y
                  x2 = cell(out)%cen(1) ; y2 = cell(out)%cen(2)
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
         
                  pv(j)=0.5d0*(cell(out)%qp(j)+pt(p2)%prim(j)) 
                  x1 = cell(out)%cen(1) ; y1 = cell(out)%cen(2)
                  x2 = pt(p2)%x    ; y2 = pt(p2)%y
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
              
                  pv(j)=0.5d0*(pt(p2)%prim(j)+cell(in)%qp(j)) 
                  x1 = pt(p2)%x    ; y1 = pt(p2)%y
                  x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
              
                  pv(j)=0.5d0*(cell(in)%qp(j)+pt(p1)%prim(j)) 
                  x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
                  x2 = pt(p1)%x    ; y2 = pt(p1)%y
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
          endif
    
          if(in==0.and.out/=0) then
                  pv(j)=0.5d0*(pt(p1)%prim(j)+cell(out)%qp(j)) 
                  x1 = pt(p1)%x    ; y1 = pt(p1)%y
                  x2 = cell(out)%cen(1) ; y2 = cell(out)%cen(2)
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
         
                  pv(j)=0.5d0*(cell(out)%qp(j)+pt(p2)%prim(j)) 
                  x1 = cell(out)%cen(1) ; y1 = cell(out)%cen(2)
                  x2 = pt(p2)%x    ; y2 = pt(p2)%y
                  dx = y2-y1       ;  dy = -(x2-x1) 
                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
         
                  pv(j)=0.5d0*(pt(p2)%prim(j)+pt(p1)%prim(j)) 
                  x1 = pt(p2)%x ; y1 = pt(p2)%y
                  x2 = pt(p1)%x ; y2 = pt(p1)%y
                  dx = y2-y1    ;  dy = -(x2-x1) 
                  fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                  fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
          endif
    
          if(out==0.and.in/=0) then
                 pv(j)=0.5d0*(pt(p1)%prim(j)+pt(p2)%prim(j)) 
                 x1 = pt(p1)%x ; y1 = pt(p1)%y
                 x2 = pt(p2)%x ; y2 = pt(p2)%y
                 dx = y2-y1    ;  dy = -(x2-x1) 
                 fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                 fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
        
                 pv(j)=0.5d0*(pt(p2)%prim(j)+cell(in)%qp(j)) 
                 x1 = pt(p2)%x    ; y1 = pt(p2)%y
                 x2 = cell(in)%cen(1) ; y2 = cell(in)%cen(2)
                 dx = y2-y1       ;  dy = -(x2-x1) 
                 fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                 fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
             
                 pv(j)=0.5d0*(cell(in)%qp(j)+pt(p1)%prim(j)) 
                 x1 = cell(in)%cen(1) ; y1 = cell(in)%cen(2)
                 x2 = pt(p1)%x    ; y2 = pt(p1)%y
                 dx = y2-y1       ;  dy = -(x2-x1) 
                 fc(i)%grad(1,j)=fc(i)%grad(1,j)+pv(j)*dx
                 fc(i)%grad(2,j)=fc(i)%grad(2,j)+pv(j)*dy
          endif
    
          fc(i)%grad(1,j)=fc(i)%grad(1,j)/fc(i)%cov
          fc(i)%grad(2,j)=fc(i)%grad(2,j)/fc(i)%cov
    
          if(in/=0) then
          cell(in)%grad(1,j)=cell(in)%grad(1,j)+fc(i)%grad(1,j)*fc(i)%cov  
          cell(in)%grad(2,j)=cell(in)%grad(2,j)+fc(i)%grad(2,j)*fc(i)%cov  
          endif
    
          if(out/=0) then
          cell(out)%grad(1,j)=cell(out)%grad(1,j)+fc(i)%grad(1,j)*fc(i)%cov  
          cell(out)%grad(2,j)=cell(out)%grad(2,j)+fc(i)%grad(2,j)*fc(i)%cov  
          endif
      enddo 

enddo


do i=1,noc
   do j=1,nvar
   cell(i)%grad(1,j)= cell(i)%grad(1,j)/cell(i)%cov
   cell(i)%grad(2,j)= cell(i)%grad(2,j)/cell(i)%cov
   enddo
enddo

if(ILIMIT==1) call limit


!100 format(1x,4(f15.8,1x))
end subroutine Gradient_GG_FC

