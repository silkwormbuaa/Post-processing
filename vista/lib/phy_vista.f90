! File    :   phy_vista.f90
! Time    :   2023/09/22 
! Author  :   Wencan WU 
! Version :   1.0
! Email   :   w.wu-3@tudelft.nl
! Desc    :   module of physics related functions

module phy

   implicit none

   contains

!-----------------------------------------------------------------------------
! >>> sutherland's law                                                  ( Nr. )
!-----------------------------------------------------------------------------
!
! Wencan Wu : w.wu-3@tudelft.nl
!
! History
!
! 2023/09/22  - created
!
! Desc
!
! Modules
!

   function sutherland( Ts ) result( mu )
!---- interface
      real(kind=8),   intent(in   ) :: Ts
      real(kind=8)                  :: mu

!---- private variables
      real(kind=8)                  :: T_ref, mu_ref, S_ref,                   &
                                       a1, b1

      T_ref  = 273.15
      mu_ref = 17.16e-6
      S_ref  = 110.4
      
      a1 = (T_ref + S_ref) / (Ts + S_ref)
      b1 = (Ts/T_ref) ** 1.5

      mu = mu_ref * a1 * b1
   
   end function sutherland




end module phy