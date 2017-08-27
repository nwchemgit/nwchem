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
!dir$ simd
        do p6=1,p6d
         r2sub(h2,h1,p4,p3)=r2sub(h2,h1,p4,p3)
     &   + factor*t2sub(h2,h1,p6,p5)*v2sub(p6,p5,p4,p3)
        enddo
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
!dir$ simd
        do p6=1,p6d
         r2sub(h2,h1,p4,p3)=r2sub(h2,h1,p4,p3)
     &   + factor*t2sub(h2,h1,p6,p5)*v2sub(p6,p5,p4,p3)
        enddo
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
      integer h1d,h2d,p3d,p4d,p5d,p6d
      integer h1,h2,p3,p4,p5,p6
      double precision t2sub(h2d,h1d,p6d,p5d)
      double precision v2sub(p6d,p5d,p4d,p3d)
      double precision r2sub(h2d,h1d,p4d,p3d)
      double precision factor
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
      return
      end
