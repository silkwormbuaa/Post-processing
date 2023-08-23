! File    :   geo_vista.f90
! Time    :   2023/08/23 
! Author  :   Wencan WU 
! Version :   1.0
! Email   :   w.wu-3@tudelft.nl
! Desc    :   None

module geo

   implicit none

   contains

!-----------------------------------------------------------------------------
! >>> point_in_triangle                                                ( Nr. )
!-----------------------------------------------------------------------------
!
! Wencan Wu : w.wu-3@tudelft.nl
!
! History
!
! 2023/08/23  - created
!
! Desc
!   - use barycentric method to detect if a point is in a triangle
!   !! only applicable when the point and triangle are within the same plane 
!
! Modules
!

   function point_in_triangle(p,p1,p2,p3)  result(flag)
!---- interface
      real(kind=8)  , intent(in   ) :: p(3), p1(3), p2(3), p3(3)
      integer                       :: flag
   
!---- private variables
      real(kind=8)                  :: v0(3), v1(3), v2(3)
      real(kind=8)                  :: dot00, dot01, dot02, dot11, dot12
      real(kind=8)                  :: d, u, v

!---- check bounding box
!     [SH] must include tolerances for RAY_TRACER
!        , but this would make this check too expensive
!     if( p(1) < p1(1) .and. p(1) < p2(1) .and. p(1) < p3(1) )RETURN
!     if( p(1) > p1(1) .and. p(1) > p2(1) .and. p(1) > p3(1) )RETURN
!     if( p(2) < p1(2) .and. p(2) < p2(2) .and. p(2) < p3(2) )RETURN
!     if( p(2) > p1(2) .and. p(2) > p2(2) .and. p(2) > p3(2) )RETURN
!     if( p(3) < p1(3) .and. p(3) < p2(3) .and. p(3) < p3(3) )RETURN
!     if( p(3) > p1(3) .and. p(3) > p2(3) .and. p(3) > p3(3) )RETURN
   
!---- take p1 as origin
   
      v0 = p3 - p1
      v1 = p2 - p1
      v2 = p  - p1
   
!---- compute dot products
   
      dot00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3) ! v0 o v0
      dot01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3) ! v0 o v1
      dot02 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3) ! v0 o v2
      dot11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3) ! v1 o v1
      dot12 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) ! v1 o v2
   
!---- compute barycentric coordinates
   
      d = dot00 * dot11 - dot01 * dot01
      u = dot11 * dot02 - dot01 * dot12
      v = dot00 * dot12 - dot01 * dot02
   
!---- evaluate

      flag = 0

      if( d >= 0.0 )then
         if ( u >= 0.0 .and. v >= 0.0 .and. u + v <= d ) then 
            flag = 1
         end if

      elseif( d < 0.0 )then
         if ( u <= 0.0 .and. v <= 0.0 .and. u + v >= d ) then
            flag = 1
         end if
      end if

   end function point_in_triangle


!-----------------------------------------------------------------------------
! >>> ray_tracer_core                                                 ( Nr. )
!-----------------------------------------------------------------------------
!
! Wencan Wu : w.wu-3@tudelft.nl
!
! History
!
! 2023/08/23  - created
!
! Desc
!
! Modules
!

   subroutine ray_tracer_core     ( x,  x_ref, norm, v0, v1, v2,               &  
                                    signum, n_tri )

!---- interface
      real(kind=8)  , intent(in   ) :: x(3), x_ref(3)
      integer       , intent(in   ) :: n_tri
      real(kind=8)  , intent(in   ) :: norm(n_tri,3), v0(n_tri,3),             & 
                                       v1(n_tri,3),   v2(n_tri,3)
      real(kind=8)  , intent(inout) :: signum

!---- private variables
      integer                       :: i, n, cnt, flag
      real(kind=8)                  :: lam, v_n, ip(3), dx(3),                 &
                                       cover(3), tri_n(3),                     &
                                       tri_p1(3),tri_p2(3),  tri_p3(3),        &
                                       lam_list(1000), x_max(3), x_min(3)
      real*16                       :: tol, tri_cen(3)
      logical                       :: new
      real*16           , parameter :: C13 = real(1.0_16)/real(3.0_16)
      real*16           , parameter :: EPS = 2.0 * epsilon(1.0)

      dx = x - x_ref

      x_max = max( x_ref, x )
      x_min = min( x_ref, x )

      cnt = 0

      do i = 1, n_tri
      
         tri_n  = norm(i,:)
         tri_p1 = v0(i,:)
         tri_p2 = v1(i,:)
         tri_p3 = v2(i,:)

!   - compute overlap of triangle and doi

         cover = min( max( tri_p1, tri_p2, tri_p3 ), x_max )   &
               - max( min( tri_p1, tri_p2, tri_p3 ), x_min )

         if( all( cover >= 0.0 ) )then

            v_n = tri_n(1)*dx(1) + tri_n(2)*dx(2) + tri_n(3)*dx(3)

            if( v_n /= 0.0 )then ! not parallel to ray
      
               tri_cen = C13 * ( tri_p1 + tri_p2 + tri_p3 )
      
               lam = ( tri_n(1) * ( tri_cen(1) - x_ref(1) )      &
                     + tri_n(2) * ( tri_cen(2) - x_ref(2) )      &
                     + tri_n(3) * ( tri_cen(3) - x_ref(3) ) )
               lam = lam / v_n
      
               if( 0.0 <= lam .and. lam <= 1.0 )then
      
!   - check whether this intersection point was found before
      
                  tol = sum( max(abs(tri_cen), abs(x_ref)) ) * EPS
      
                  new = .true.
                  do n = 1, cnt
                     if( abs(lam-lam_list(n)) < tol )then
                        new = .false.
                        EXIT
                     end if
                  end do
      
!   - increase counter if triangle intersects with ray in new point
      
                  if( new )then
      
                     ip = x_ref + dx * lam

                     flag = point_in_triangle( ip, tri_p1, tri_p2, tri_p3 )

                     if( flag==1 )then
                        cnt = cnt + 1
                        if( cnt == 1001 )cnt = 999 ! buffer overflow
                        lam_list(cnt) = lam
                     end if ! found valid intersection
      
                  end if ! new
               end if ! intersection within ray range
            end if ! triangle not parallel to ray
         end if ! triangle within doi
      end do

      if( mod( cnt, 2 ) == 1 ) signum = -signum

      end subroutine ray_tracer_core


end module geo