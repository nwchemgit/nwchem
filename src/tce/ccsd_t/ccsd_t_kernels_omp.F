#define OMP_COLLAPSE_LEVEL 3
c
c These have been separated out from ccsd_t_singles_l.F and ccsd_t_doubles_l.F
c
      subroutine sd_t_s1_1(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h2d,h1d,p6d,p5d,p4d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p4=1,p4d
      do p5=1,p5d
      do p6=1,p6d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
#ifdef USE_OPENMP
!DIR$ IVDEP
#endif
      do h3=1,h3d
       triplesx(h3,h2,h1,p6,p5,p4)=triplesx(h3,h2,h1,p6,p5,p4)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_2(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p6d,p5d,p4d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p4=1,p4d
      do p5=1,p5d
      do p6=1,p6d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d ! interchanging h1 and h2 was a huge win with Cray+OpenMP
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
#ifdef USE_OPENMP
!DIR$ IVDEP
#endif
      do h3=1,h3d
       triplesx(h3,h1,h2,p6,p5,p4)=triplesx(h3,h1,h2,p6,p5,p4)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_3(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p6d,p5d,p4d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
!!!!dir$ prefetch t1sub
      do p4=1,p4d
      do p5=1,p5d
      do p6=1,p6d
#ifdef CRAYFTN
      do h1=1,h1d
#endif
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
#ifdef USE_OPENMP
!DIR$ IVDEP
#endif
#ifndef CRAYFTN
      do h1=1,h1d
#endif
       triplesx(h1,h3,h2,p6,p5,p4)=triplesx(h1,h3,h2,p6,p5,p4)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_4(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h2d,h1d,p6d,p4d,p5d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p5=1,p5d
!!!dir$ prefetch t1sub
      do p4=1,p4d
      do p6=1,p6d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
#ifdef USE_OPENMP
!DIR$ IVDEP
#endif
      do h3=1,h3d
       triplesx(h3,h2,h1,p6,p4,p5)=triplesx(h3,h2,h1,p6,p4,p5)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_5(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p6d,p4d,p5d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p5=1,p5d
!!!dir$ prefetch t1sub
      do p4=1,p4d
      do p6=1,p6d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
#ifdef USE_OPENMP
!DIR$ IVDEP
#endif
      do h3=1,h3d
       triplesx(h3,h1,h2,p6,p4,p5)=triplesx(h3,h1,h2,p6,p4,p5)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_6(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p6d,p4d,p5d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p5=1,p5d
!!!dir$ prefetch t1sub
      do p4=1,p4d
      do p6=1,p6d
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
#ifdef USE_OPENMP
!DIR$ IVDEP
#endif
      do h1=1,h1d
       triplesx(h1,h3,h2,p6,p4,p5)=triplesx(h1,h3,h2,p6,p4,p5)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_7(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h2d,h1d,p4d,p6d,p5d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p5=1,p5d
      do p6=1,p6d
!!!dir$ prefetch t1sub
      do p4=1,p4d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
#ifdef USE_OPENMP
!DIR$ IVDEP
#endif
      do h3=1,h3d
       triplesx(h3,h2,h1,p4,p6,p5)=triplesx(h3,h2,h1,p4,p6,p5)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_8(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p4d,p6d,p5d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p5=1,p5d
      do p6=1,p6d
!!!dir$ prefetch t1sub
      do p4=1,p4d
      do h2=1,h2d
!!!dir$ loop count min(8)
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
!DIR$ IVDEP
      do h3=1,h3d
       triplesx(h3,h1,h2,p4,p6,p5)=triplesx(h3,h1,h2,p4,p6,p5)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_s1_9(h3d,h2d,h1d,p6d,p5d,p4d,
     2                     triplesx,t1sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p4d,p6d,p5d)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p5=1,p5d
      do p6=1,p6d
!!!dir$ prefetch t1sub
      do p4=1,p4d
      do h2=1,h2d
!!!dir$ loop count min(8)
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dec$ vector always nontemporal
!dir$ simd
!DIR$ IVDEP
      do h1=1,h1d
       triplesx(h1,h3,h2,p4,p6,p5)=triplesx(h1,h3,h2,p4,p6,p5)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_1(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h2d,h1d,p6d,p5d,p4d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p5=1,p5d
      do p6=1,p6d
      do p4=1,p4d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h3,h2,h1,p6,p5,p4)=triplesx(h3,h2,h1,p6,p5,p4)
     1  -t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_2(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p6d,p5d,p4d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p4=1,p4d
      do p5=1,p5d
      do p6=1,p6d
      do h2=1,h2d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h3,h1,h2,p6,p5,p4)=triplesx(h3,h1,h2,p6,p5,p4)
     1   + t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_3(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d,h2d,p6d,p5d,p4d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p4=1,p4d
      do p5=1,p5d
      do p6=1,p6d
      do h2=1,h2d
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h1,h3,h2,p6,p5,p4)=triplesx(h1,h3,h2,p6,p5,p4)
     1  -t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_4(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h2d,h1d,p5d,p4d,p6d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
#ifdef __ppc64__
      do p4=1,p4d ! this loop used to be between p6 and p5
#endif
      do p5=1,p5d
      do h1=1,h1d
      do h2=1,h2d
#ifndef __ppc64__
      do p4=1,p4d ! this loop used to be between p6 and p5
#endif
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h3,h2,h1,p5,p4,p6)=triplesx(h3,h2,h1,p5,p4,p6)
     1  -t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_5(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p5d,p4d,p6d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6)
      do p6=1,p6d
      do p4=1,p4d
      do p5=1,p5d
      do h2=1,h2d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h3,h1,h2,p5,p4,p6)=triplesx(h3,h1,h2,p5,p4,p6)
     1   + t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_6(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d,h2d,p5d,p4d,p6d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do p4=1,p4d
      do p5=1,p5d
      do h2=1,h2d
      do h3=1,h3d
!dec$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h1=1,h1d
!dec$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h1,h3,h2,p5,p4,p6)=triplesx(h1,h3,h2,p5,p4,p6)
     1  -t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_7(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h2d,h1d,p5d,p6d,p4d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
#if !defined(CRAYFTN)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h2=1,h2d
      do h7=1,h7d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
#endif
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p4=1,p4d
      do p6=1,p6d
      do p5=1,p5d
#ifdef CRAYFTN
      do h2=1,h2d ! interchanging h1 and h2 was a huge win with Cray+OpenMP
      do h1=1,h1d
#else
      do h1=1,h1d
      do h2=1,h2d
#endif
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h3,h2,h1,p5,p6,p4)=triplesx(h3,h2,h1,p5,p6,p4)
#if !defined(CRAYFTN)
     1   + t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
#else
     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
#endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_8(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p5d,p6d,p4d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p4=1,p4d
      do p6=1,p6d
      do p5=1,p5d
      do h2=1,h2d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h3,h1,h2,p5,p6,p4)=triplesx(h3,h1,h2,p5,p6,p4)
     1  -t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d1_9(h3d,h2d,h1d,p6d,p5d,p4d,h7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d,h2d,p5d,p6d,p4d)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision v2tmp(h7d,h3d,h2d,p6d)
!$omp  parallel do collapse(3)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p6=1,p6d
      do h7=1,h7d
      do h2=1,h2d
!dec$ vector always nontemporal
!DIR$ IVDEP
      do h3=1,h3d
        v2tmp(h7,h3,h2,p6) = v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,h7)
      do p4=1,p4d
      do p6=1,p6d
      do p5=1,p5d
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h1=1,h1d
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do h7=1,h7d
       triplesx(h1,h3,h2,p5,p6,p4)=triplesx(h1,h3,h2,p5,p6,p4)
     1   + t2sub(h7,p4,p5,h1)*v2tmp(h7,h3,h2,p6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_1(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p6d,p5d,p4d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p4=1,p4d
      do p5=1,p5d
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h3,h2,h1,p6,p5,p4)=triplesx(h3,h2,h1,p6,p5,p4)
     1  -t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_2(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h1d,h3d,p6d,p5d,p4d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p6=1,p6d
      do h1=1,h1d
      do p4=1,p4d
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h2,h1,h3,p6,p5,p4)=triplesx(h2,h1,h3,p6,p5,p4)
     1  -t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_3(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p6d,p5d,p4d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p6=1,p6d
      do h1=1,h1d
      do p4=1,p4d
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h2,h3,h1,p6,p5,p4)=triplesx(h2,h3,h1,p6,p5,p4)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_4(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p6d,p4d,p5d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p4=1,p4d
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h3,h2,h1,p6,p4,p5)=triplesx(h3,h2,h1,p6,p4,p5)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_5(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h1d,h3d,p6d,p4d,p5d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p4=1,p4d
      do p6=1,p6d
      do h1=1,h1d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h2,h1,h3,p6,p4,p5)=triplesx(h2,h1,h3,p6,p4,p5)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_6(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p6d,p4d,p5d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p4=1,p4d
      do p6=1,p6d
      do h1=1,h1d
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h2,h3,h1,p6,p4,p5)=triplesx(h2,h3,h1,p6,p4,p5)
     1  -t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_7(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p4d,p6d,p5d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p6=1,p6d
      do p4=1,p4d
      do h1=1,h1d
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h3,h2,h1,p4,p6,p5)=triplesx(h3,h2,h1,p4,p6,p5)
     1  -t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_8(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h1d,h3d,p4d,p6d,p5d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p6=1,p6d
      do p4=1,p4d
!      do h3=1,h3d ! this loop order drops Edison from 80 to 50
      do h1=1,h1d
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h2=1,h2d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h2,h1,h3,p4,p6,p5)=triplesx(h2,h1,h3,p4,p6,p5)
     1  -t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c
      subroutine sd_t_d2_9(h3d,h2d,h1d,p6d,p5d,p4d,p7d,
     2                     triplesx,t2sub,v2sub)
      IMPLICIT NONE
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p4d,p6d,p5d)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
!$omp  parallel do collapse(OMP_COLLAPSE_LEVEL)
!$omp& default(shared) schedule(static)
!$omp& private(h1,h2,h3,p4,p5,p6,p7)
      do p5=1,p5d
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
      do p4=1,p4d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
      do h3=1,h3d
!!!dir$ loop count min(8)
!dec$ unroll_and_jam = 8
!dir$ simd
!DIR$ IVDEP
      do p7=1,p7d
       triplesx(h2,h3,h1,p4,p6,p5)=triplesx(h2,h3,h1,p4,p6,p5)
     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
