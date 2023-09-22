! File    :   mth_vista.f90
! Time    :   2023/09/22 
! Author  :   Wencan WU 
! Version :   1.0
! Email   :   w.wu-3@tudelft.nl
! Desc    :   module of math related functions

module mth

   implicit none

   contains

!-----------------------------------------------------------------------------
! >>> bilinear interpolation                                        ( Nr. )
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
!     f(3)  # . . . . . . . . # f(4)
!    (x1,y2).                 . (x2,y2)
!           .                 .
!           .         f(x,y)  .
!           .            *    .
!           .                 .
!     f(1)  # . . . . . . . . # f(2)
!       (x1,y1)              (x2,y1)
!     
!

   function bilin_interp( x1, x2, y1, y2, f, x, y ) result( value )
!---- interface
      real(kind=8),   intent(in   ) :: x1,x2,y1,y2, f(4), x, y
      real(kind=8)                  :: value

!---- private variables
      real(kind=8)                  :: u, v   ! weight in x,y directions

!---- calculate weights

      u = (x - x1) / (x2 - x1)
      v = (y - y1) / (y2 - y1)

!---- calculate interpolated value

      value =  (1.0-u)*(1.0-v) * f(1)                                          &
             +      u *(1.0-v) * f(2)                                          &
             + (1.0-u)*     v  * f(3)                                          &
             +      u *     v  * f(4)

   end function bilin_interp


end module mth
