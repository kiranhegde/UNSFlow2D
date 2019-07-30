module visc
use data_type
use commons
implicit none
! laminar paramters
! Rey = Reynolds number
! SCONST = Constant in Sutherland Law
character(len=24)  :: flow_type,grad_type
real(kind=dp) :: Rey, SCONST
real(kind=dp), parameter:: prandtl=0.72_dp, prandtl_turb=0.9_dp

! Parameters in Spallart-Allmaras model
real(kind=dp) :: Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1,  &
            Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2

end module visc

