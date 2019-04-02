module inf
use data_type
use commons
implicit none
! Freestream values
real(kind=dp) :: q_inf,m_inf, aoa, aoa_deg, u_inf, v_inf,  &
                 r_inf, p_inf, T_inf, T_infd, ent_inf,  &
                 conv_inf(nvar), qinf(nvar),a_inf, H_inf

! Vortex correction for farfield BC
integer(kind=i4)  :: vortex
real(kind=dp) :: xref, yref
end module inf
