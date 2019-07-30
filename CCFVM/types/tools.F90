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


subroutine pt2linemirror(p,a,px)
use data_type
implicit none
type (vector) :: a
type (coord)  :: p,px
real(kind=dp) :: dx,dy,r,x,y

    dx=a%e%x-a%s%x
    dy=a%e%y-a%s%y
     x=a%s%x-p%x
     y=a%s%y-p%y
     r=1.0_dp/(dx*dx+dy*dy)  
     
   px%x=p%x+2.0_dp*(x-x*dx*dx*r-y*dx*dy*r) 
   px%y=p%y+2.0_dp*(y-y*dy*dy*r-x*dx*dy*r) 

end subroutine pt2linemirror


end module tools

module strings
use data_type
implicit none

integer(kind=i4), parameter, private :: EOS = -1


contains

  subroutine int2str(str, iarr, n)
    !use bl_error_module
    use data_type
    character(len=*), intent(out) :: str
    integer(kind=i4), intent(in) :: n
    integer(kind=i4), intent(in) :: iarr(n)
    integer(kind=i4):: i

    if ( len(str) < n ) then
       print*,"int2str: str to large for iarr: size(iarr) = ", n
       !call bl_error("INT2STR: iarr to large for str: len = ", len(str))
    end if
    do i = 1, n
       if ( iarr(i) == EOS ) exit
       str(i:i) = char(iarr(i))
    end do

  end subroutine int2str

  !! Converts a Fortran string to an integer encoding.
  subroutine str2int(iarr, n, str)
    use data_type
    !use bl_error_module
    character(len=*), intent(in) :: str
    integer(kind=i4), intent(in) :: n
    integer(kind=i4) :: i, j
    integer(kind=i4), intent(out) :: iarr(n)

    if ( n <= len_trim(str) ) then
       print*,"str2int: str to large for iarr: size(iarr) = ", n
       !call bl_error("STR2INT: str to large for iarr: size(iarr) = ", n)
    end if

    iarr = 0
    j = 1
    do i = 1, len_trim(str)
       iarr(j) = ichar(str(i:i))
       j = j + 1
    end do

    iarr(j) = EOS

  end subroutine str2int

end  module strings
