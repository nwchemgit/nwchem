#ifdef OFFLOAD_CODE
      subroutine offl_ccsd_t_transpt2(
#else
      subroutine ccsd_t_transpt2(
#endif
     C     t2sub,scratch,
     c     h1d,p5d,p4d,h7d)
      implicit none
      integer h1d,p5d,p4d,h7d
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision scratch(h1d,p5d,p4d,h7d)
c
      integer h1,p5,p4,h7
c
      integer pp5,hh1
!DIR$ ASSUME_ALIGNED t2sub: 64
!DIR$ ASSUME_ALIGNED scratch: 64
c
#define CHUNK 16
#ifdef OFFLOAD_CODE
!$omp target map(to:t2sub)
#endif
!$omp parallel do private(h7,p4,pp5,hh1,p5,h1) collapse(2)
      do h7=1,h7d
         do p4=1,p4d
            do pp5=1,p5d,CHUNK
               do hh1=1,h1d,CHUNK
                  do p5=pp5,min(pp5+CHUNK-1,p5d)
!DEC$ LOOP COUNT AVG=CHUNK
cc!deC$ SIMD
                     do h1=hh1,min(hh1+CHUNK-1,h1d)
                        scratch(h1,p5,p4,h7)=t2sub(h7,p4,p5,h1)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
#ifdef OFFLOAD_CODE
!$omp end target
#endif
      return
      end


#ifdef OFFLOAD_CODE
      subroutine offl_ccsd_t_transpt2_7124(
#else
      subroutine ccsd_t_transpt2_7124(
#endif
     t     t2sub,scratch,
     c    p4d,h1d,h2d,p7d)
      implicit none
      integer p4d,h1d,h2d,p7d
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision scratch(p7d,h2d,h1d,p4d)
c
      integer p4,h1,h2,p7
c
      integer hh2,pp7
!DIR$ ASSUME_ALIGNED t2sub: 64
!DIR$ ASSUME_ALIGNED scratch: 64
c
#define CHUNK 16
#ifdef OFFLOAD_CODE
!$omp target map(to:t2sub)
#endif
!$omp parallel do private(p4,h1,hh2,pp7,h2,p7) collapse(2)
      do p4=1,p4d
         do h1=1,h1d
            do hh2=1,h2d,CHUNK
               do pp7=1,p7d,CHUNK
                  do h2=hh2,min(hh2+CHUNK-1,h2d)
!DEC$ LOOP COUNT AVG=CHUNK
cc!deC$ SIMD
                     do p7=pp7,min(pp7+CHUNK-1,p7d)
                        scratch(p7,h2,h1,p4)=t2sub(p7,p4,h1,h2)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
#ifdef OFFLOAD_CODE
!$omp end target
#endif
      return
      end


#ifdef OFFLOAD_CODE
      subroutine offl_ccsd_t_transpt2_7326(
#else
      subroutine ccsd_t_transpt2_7326(
#endif
     t     v2sub,scratch,
     c    h3d,h2d,p6d,h7d)
      implicit none
      integer h3d,h2d,p6d,h7d
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision scratch(h7d,h3d,h2d,p6d)
c
      integer h3,h2,p6,h7
c
      integer hh3,hh7
!DIR$ ASSUME_ALIGNED v2sub: 64
!DIR$ ASSUME_ALIGNED scratch: 64
c
#ifdef OFFLOAD_CODE
!$omp target map(to:v2sub)
#endif
#define CHUNK 16
!$omp parallel do private(p6,h2,hh3,hh7,h3,h7) collapse(2)
      do p6=1,p6d
         do h2=1,h2d
            do hh3=1,h3d,CHUNK
               do hh7=1,h7d,CHUNK
                  do h3=hh3,min(hh3+CHUNK-1,h3d)
!DEC$ LOOP COUNT AVG=CHUNK
cc!deC$ SIMD
                     do h7=hh7,min(hh7+CHUNK-1,h7d)
                        scratch(h7,h3,h2,p6)=v2sub(h3,h2,p6,h7)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
#ifdef OFFLOAD_CODE
!$omp end target
#endif
      return
      end
#ifndef OFFLOAD_CODE
#ifdef USE_OFFLOAD
#define OFFLOAD_CODE
c
c     Compile source again for offload code
c
#include "ccsd_t_transpose.F"
#endif
#endif

