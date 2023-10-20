
module structs
  use constants
  implicit none

  type control 
    integer :: herit_p
    character(len = 6 ) :: JJordTyp
    character(len = 6 ) :: which_corr_u_mat
    character(len = short_str ) :: path
    character(len = short_str ) :: dbgstr
    logical :: yes_inherit
    logical :: swap_1st_p
    logical :: yes_print
    logical :: yes_rad
    logical :: yes_skip
    real(long) :: skip_beg, skip_end
    real(long), dimension(ms) :: xplus
    integer, dimension(ms) :: jjord
    integer, dimension(me,ms) :: man_swp
    integer, dimension(me,ms) :: man_sgn
    logical :: yes_man_swp
    integer :: mxe, nchop ! number sample energy points, number open chan
  end type control

  type SDmat
    integer :: mxe
    integer :: nchop
    integer :: mlast, mfin
    real(long) :: Excit ! excitation energy w.r.t. IP1
    real(long), dimension(ms,ms) :: U
    real(long), dimension(ms) :: miu
    real(long), dimension(ms,max_in) :: DL
    real(long), dimension(ms,max_in) :: DV
    real(long), dimension(ms,max_in) :: Di
    real(long), dimension(max_in) :: Eistat ! w.r.t. IP1
    real(long), dimension(ms,6) :: jjls
    real(long), dimension(ms,3) :: LSJ
    real(long) :: JT
    integer, dimension(ms) :: chan_targcf
    character(len = short_str), dimension(ms) :: chan_ion_str, chan_eig_str
    real(long), dimension(ms,ms) :: Ug
  end type SDmat
  type fname
    integer :: unit_dbg = 14
    character(len = long_str) :: dbg = 'debug.log'
    integer :: unit_smooth_ctrl = 100
    character(len = long_str) :: smooth_ctrl = 'smooth.inp'
    integer :: unit_S_in = 123 
    character(len = long_str) :: S_in = 'ESUM-unf'
    integer :: unit_S_out_unf = 124
    character(len = long_str) :: S_out_unf = 'ESUM-unf-smooth'
    integer :: unit_S_out_for = 133
    character(len = long_str) :: S_out_for = 'miuang.out'
    integer :: unit_D_in = 125
    character(len = long_str) :: D_in = 'Dalpha-unf'
    integer :: unit_D_out_unf = 126
    character(len = long_str) :: D_out_unf = 'Dalpha-unf.smth'
    integer :: unit_D_out_for = 127
    character(len = long_str) :: D_out_for = 'Dalpha.txt'
    integer :: unit_inherit = 138
    character(len = long_str) :: inherit = 'inherit.txt'
    integer :: unit_herit = 139
    character(len = long_str) :: herit = 'herit.txt'
    integer :: unit_u_order = 132
    character(len = long_str) :: U_order = 'u-order.out'
    integer :: unit_u_sign = 134
    character(len = long_str) :: U_sign = 'sign.out'
  end type fname
end module structs
