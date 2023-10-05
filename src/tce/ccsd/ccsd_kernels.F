!
! Reference implementation written by Karol Kowalski
!
      subroutine t2_p8(h1d,h2d,p3d,p4d,p5d,p6d,
     1                 t2sub,v2sub,r2sub,factor)
      IMPLICIT NONE
      integer h1d,h2d,p3d,p4d,p5d,p6d
      integer h1,h2,p3,p4,p5,p6
      double precision t2sub(h2d,h1d,p6d,p5d)
      double precision v2sub(p6d,p5d,p4d,p3d)
      double precision r2sub(h2d,h1d,p4d,p3d)
      double precision factor
      if ((p5d.lt.8).or.(p6d.lt.8)) then
!$omp  parallel do collapse(2)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,p3,p4,p5,p6)
        do p3=1,p3d
        do p4=1,p4d
        do h1=1,h1d
        do h2=1,h2d
        do p5=1,p5d
!$omp simd
        do p6=1,p6d
         r2sub(h2,h1,p4,p3)=r2sub(h2,h1,p4,p3)
     &   + factor*t2sub(h2,h1,p6,p5)*v2sub(p6,p5,p4,p3)
        enddo
!$omp end simd
        enddo
        enddo
        enddo
        enddo
        enddo
!$omp end parallel do
!
! All dims are at least 8, so more SIMD optimizations allowed.
!
      else ! inner loops at least 8
!$omp  parallel do collapse(2)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,p3,p4,p5,p6)
        do p3=1,p3d
        do p4=1,p4d
        do h1=1,h1d
        do h2=1,h2d
!dir$ loop count min(8), max(80), avg(32)
!dec$ unroll_and_jam = 8
        do p5=1,p5d
!dir$ loop count min(8), max(80), avg(32)
!dec$ unroll_and_jam = 8
!$omp simd
        do p6=1,p6d
         r2sub(h2,h1,p4,p3)=r2sub(h2,h1,p4,p3)
     &   + factor*t2sub(h2,h1,p6,p5)*v2sub(p6,p5,p4,p3)
      enddo
!$omp end simd
        enddo
        enddo
        enddo
        enddo
        enddo
!$omp end parallel do
      endif ! inner loops at least 8
      return
      end

      subroutine t2_p8_x(h1d,h2d,p3d,p4d,p5d,p6d,
     1                   t2sub,v2sub,r2sub,factor)
      IMPLICIT NONE
      integer, intent(in) :: h1d,h2d,p3d,p4d,p5d,p6d
      double precision, intent(in) :: factor
#if 0
      double precision, intent(in) :: t2sub(h2d,h1d,p6d,p5d)
      double precision, intent(in) :: v2sub(p6d,p5d,p4d,p3d)
      double precision, intent(inout) :: r2sub(h2d,h1d,p4d,p3d)
      integer :: h1,h2,p3,p4,p5,p6
      do p3=1,p3d
      do p4=1,p4d
      do h1=1,h1d
      do h2=1,h2d
      do p5=1,p5d
      do p6=1,p6d
       r2sub(h2,h1,p4,p3)=r2sub(h2,h1,p4,p3)
     & + factor*t2sub(h2,h1,p6,p5)*v2sub(p6,p5,p4,p3)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#elif 0
      double precision, intent(in) :: t2sub(h2d*h1d,p6d*p5d)
      double precision, intent(in) :: v2sub(p6d*p5d,p4d*p3d)
      double precision, intent(inout) :: r2sub(h2d*h1d,p4d*p3d)
      integer :: p34, h12, p56
      do p34=1,p3d*p4d
      do h12=1,h1d*h2d
      do p56=1,p5d*p6d
       r2sub(h12,p34) = r2sub(h12,p34)
     & + factor*t2sub(h12,p56)*v2sub(p56,p34)
      enddo
      enddo
      enddo
#else
      ! lda = h2d*h1d
      ! ldb = p6d*p5d
      ! ldc = h2d*h1d
      ! m = h2d*h1d
      ! n = p4d*p3d
      ! k = p6d*p5d
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
#endif
      return
      end
