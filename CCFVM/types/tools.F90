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
