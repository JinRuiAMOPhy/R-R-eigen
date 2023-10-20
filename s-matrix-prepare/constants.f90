module constants 
  implicit none
  integer, parameter :: long = kind(1.d0)
  integer, parameter :: me = 10000, ms = 500
  integer, parameter :: max_in = 20
  real(long), parameter :: pai = 3.141592653589793_long
  real(long), parameter :: ZERO  = 0.0_long
  real(long), parameter :: ONE   = 1.0_long
  real(long), parameter :: TWO   = 2.0_long
  real(long), parameter :: THREE = 3.0_long
  real(long), parameter :: FOUR  = 4.0_long
  real(long), parameter :: FIVE  = 5.0_long
  real(long), parameter :: SIX   = 6.0_long
  real(long), parameter :: SEVEN = 7.0_long
  real(long), parameter :: EIGHT = 8.0_long
  real(long), parameter :: NINE  = 9.0_long
  real(long), parameter :: TEN   = 10.0_long
  real(long), parameter :: THIRTY  = 30.0_long
  real(long), parameter :: HALF  = 0.5_long
  real(long), parameter :: TENTH  = 0.1_long
  real(long), parameter :: DPI  = 3.141592653589793_long
  real(long), parameter :: THIRD = ONE / THREE
  real(long), parameter :: TWOTHIRD = TWO / THREE
  real(long), parameter :: FOURTHIRD = FOUR / THREE
  integer, parameter :: short_str = 32
  integer, parameter :: max_col = 100
  integer, parameter :: longlong_str = 3000
  integer, parameter :: long_str = longlong_str / 2
  ! print control
  integer, parameter :: print_col_l_fmt = 12
  integer, parameter :: print_col_s_fmt = 5
  integer, parameter :: STDIN  = 5      ! unit number of standard input
  integer, parameter :: STDOUT = 6      ! unit number of standard output
  integer, parameter :: STDERR = 0      ! unit number of standard error
  integer, parameter :: funit_dbg = 888      ! unit number of debug file
  !
  ! physical constants from NIST
  !   ref) Mohr & Taylor, J. Phys. Chem. Ref. Data 28, 1713 (1999)
  !        Mohr & Taylor, Rev. Mod. Phys. 72, 351 (2000)
  !   E: electron, P: proton, and N: neutron
  !
  real(long), parameter :: RYDBERG=109737.3177_long    ! 1 Ryd.->1cm-1


end module constants 

