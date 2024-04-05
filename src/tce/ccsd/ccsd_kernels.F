!
! Reference implementation written by Karol Kowalski
! Loops fully converted to DGEMM later.
!
      subroutine t2_p8(h1d,h2d,p3d,p4d,p5d,p6d,
     1                 t2sub,v2sub,r2sub,factor)
      IMPLICIT NONE
      integer, intent(in) :: h1d,h2d,p3d,p4d,p5d,p6d
      double precision, intent(in) :: factor
      double precision, intent(in) :: t2sub(h2d*h1d,p6d*p5d)
      double precision, intent(in) :: v2sub(p6d*p5d,p4d*p3d)
      double precision, intent(inout) :: r2sub(h2d*h1d,p4d*p3d)
#if 0
      r2sub = r2sub + factor * matmul(t2sub,v2sub)
#else
      call DGEMM('n','n',h2d*h1d,p4d*p3d,p6d*p5d,
     &           factor,t2sub,h2d*h1d,
     &                  v2sub,p6d*p5d,
     &           1.0d0, r2sub,h2d*h1d)
#endif
      end

      subroutine t2_p8_x(h1d,h2d,p3d,p4d,p5d,p6d,
     1                   t2sub,v2sub,r2sub,factor)
      IMPLICIT NONE
      integer, intent(in) :: h1d,h2d,p3d,p4d,p5d,p6d
      double precision, intent(in) :: factor
      double precision, intent(in) :: t2sub(h2d*h1d,p6d*p5d)
      double precision, intent(in) :: v2sub(p6d*p5d,p4d*p3d)
      double precision, intent(inout) :: r2sub(h2d*h1d,p4d*p3d)
#if 0
      r2sub = r2sub + factor * matmul(t2sub,v2sub)
#else
      call DGEMM('n','n',h2d*h1d,p4d*p3d,p6d*p5d,
     &           factor,t2sub,h2d*h1d,
     &                  v2sub,p6d*p5d,
     &           1.0d0, r2sub,h2d*h1d)
#endif
      end
