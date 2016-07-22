#define OMPCOLLAPSE 3
c#define INTEL_PRINT_TIMINGS 1

      subroutine sd_6ds_1(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_1'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3h2,h1,p6,p5,p4
      double precision triplesx(h3d*h2d,h1d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d*h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h3h2)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
!dec$ loop count max=40, min=8
      do h1=1,h1d
!dec$ loop count max=600, min=20
      do h3h2=1,h2d*h3d
       triplesx(h3h2,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3h2,h1,p6,p5-jstart+1,p4-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_2(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_2'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
!dec$ loop count max=40, min=8
      do h1=1,h1d
!dec$ loop count max=40, min=8
      do h3=1,h3d
       triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_3(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     2 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_3'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
!dec$ loop count max=40, min=8
      do h3=1,h3d
!dec$ loop count max=40, min=8
      do h1=1,h1d
       triplesx(h1,h3,h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h1,h3,h2,p6,p5-jstart+1,p4-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_4(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_4'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3h2,h1,p6,p5,p4
      double precision triplesx(h3d*h2d,h1d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d*h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h3h2)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
!dec$ loop count max=40, min=8
      do h1=1,h1d
!dec$ loop count max=600, min=20
      do h3h2=1,h3d*h2d
       triplesx(h3h2,h1,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h3h2,h1,p6,p4-jstart+1,p5-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_5(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_5'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
!dec$ loop count max=40, min=8
      do h1=1,h1d
!dec$ loop count max=40, min=8
      do h3=1,h3d
       triplesx(h3,h1,h2,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h3,h1,h2,p6,p4-jstart+1,p5-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_6(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_6'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p6d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
!dec$ loop count max=40, min=8
      do h3=1,h3d
!dec$ loop count max=40, min=8
      do h1=1,h1d
       triplesx(h1,h3,h2,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h1,h3,h2,p6,p4-jstart+1,p5-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_7(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_7'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3h2,h1,p6,p5,p4
      double precision triplesx(h3d*h2d,h1d,p4d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d*h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h3h2)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
!dec$ loop count max=40, min=20
      do h1=1,h1d
!dec$ loop count max=600, min=20
      do h3h2=1,h3d*h2d
       triplesx(h3h2,h1,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h3h2,h1,p4,p6-jstart+1,p5-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_8(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_8'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h3d,h1d,h2d,p4d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h2=1,h2d
!dec$ loop count max=40, min=20
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h3=1,h3d
       triplesx(h3,h1,h2,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h3,h1,h2,p4,p6-jstart+1,p5-istart+1)
     1   - t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_6ds_9(h3d,h2d,h1d,p6d,p5d,p4d,
     2               triplesx,t1sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_6ds_9'
      integer istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d
      integer h3,h2,h1,p6,p5,p4
      double precision triplesx(h1d,h3d,h2d,p4d,maxp5,maxp4)
      double precision t1sub(p4d,h1d)
      double precision v2sub(h3d,h2d,p6d,p5d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t1sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p6d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h2=1,h2d
!dec$ loop count max=40, min=20
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do h1=1,h1d
       triplesx(h1,h3,h2,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h1,h3,h2,p4,p6-jstart+1,p5-istart+1)
     1   + t1sub(p4,h1)*v2sub(h3,h2,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end



      subroutine sd_tx1_1(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_1'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d*h2d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d*h2d,p6d,h7d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h3h2,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h7=1,h7d
!dec$ loop count max=600, min=20
      do h3h2=1,h2d*h3d
       triplesx(h3h2,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3h2,h1,p6,p5-jstart+1,p4-istart+1)
     1   - t2sub(h7,p4,p5,h1)*v2sub(h3h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_2(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_2'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p6d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h2,h3,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h2=1,h2d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do h7=1,h7d
       triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h1,h2,p6,p5-jstart+1,p4-istart+1)
     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_3(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_3'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d*h2d,p6d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision scratch(h7d,h3d*h2d,p6d)
      external ccsd_t_transpt2_7326
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7326(v2sub,scratch,h3d,h2d,p6d,h7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,scratch)
!$omp&            private(p4,p5,p6,h1,h3h2,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
!dec$ loop count max=600, min=30
      do h3h2=1,h2d*h3d
!dec$ loop count max=40, min=8
      do h1=1,h1d
!dec$ loop count max=40, min=8
      do h7=1,h7d
       triplesx(h1,h3h2,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h1,h3h2,p6,p5-jstart+1,p4-istart+1)
     1    - t2sub(h7,p4,p5,h1)*scratch(h7,h3h2,p6)
c     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_4(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_4'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d*h2d,h1d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d*h2d,p6d,h7d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h3h2,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p6=istart,istop
      do p4=jstart,jstop
      do p5=1,p5d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h7=1,h7d
!dec$ loop count max=600, min=20
      do h3h2=1,h2d*h3d
       triplesx(h3h2,h1,p5,p4-jstart+1,p6-istart+1)=
     1     triplesx(h3h2,h1,p5,p4-jstart+1,p6-istart+1)
     1   - t2sub(h7,p4,p5,h1)*v2sub(h3h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_5(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_5'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision scratch(h7d,h3d,h2d,p6d)
      external ccsd_t_transpt2_7326
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7326(v2sub,scratch,h3d,h2d,p6d,h7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,scratch)
!$omp&            private(p4,p5,p6,h1,h2,h3,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p6=istart,istop
      do p4=jstart,jstop
      do p5=1,p5d
!dec$ loop count max=40, min=8
      do h2=1,h2d
!dec$ loop count max=40, min=8
      do h1=1,h1d
!dec$ loop count max=40, min=8
!dec$ unroll_and_jam (2)
      do h3=1,h3d
!dec$ loop count max=40, min=8
      do h7=1,h7d
       triplesx(h3,h1,h2,p5,p4-jstart+1,p6-istart+1)=
     1     triplesx(h3,h1,h2,p5,p4-jstart+1,p6-istart+1)
     1   + t2sub(h7,p4,p5,h1)*scratch(h7,h3,h2,p6)
c     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_6(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_6'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d*h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision scratch(h7d,h3d*h2d,p6d)
      external ccsd_t_transpt2_7326
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7326(v2sub,scratch,h3d,h2d,p6d,h7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,scratch)
!$omp&            private(p4,p5,p6,h1,h3h2,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p6=istart,istop
      do p4=jstart,jstop
      do p5=1,p5d
!dec$ loop count max=900, min=30
      do h3h2=1,h2d*h3d
!dec$ loop count max=40, min=8
      do h1=1,h1d
!dec$ loop count max=40, min=8
      do h7=1,h7d
       triplesx(h1,h3h2,p5,p4-jstart+1,p6-istart+1)=
     1     triplesx(h1,h3h2,p5,p4-jstart+1,p6-istart+1)
     1   - t2sub(h7,p4,p5,h1)*scratch(h7,h3h2,p6)
c     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_7(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_7'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d*h2d,h1d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d*h2d,p6d,h7d)
#include "ccsd_t_dts_kernels_omp_top.fh"
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,v2sub)
!$omp&            private(p4,p5,p6,h1,h3h2,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p6=jstart,jstop
      do p5=1,p5d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h7=1,h7d
!dec$ loop count max=600, min=20
      do h3h2=1,h2d*h3d
       triplesx(h3h2,h1,p5,p6-jstart+1,p4-istart+1)=
     1     triplesx(h3h2,h1,p5,p6-jstart+1,p4-istart+1)
     1   + t2sub(h7,p4,p5,h1)*v2sub(h3h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_8(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_8'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3,h2,h1,p6,p5,p4,h7
      double precision triplesx(h3d,h1d,h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision scratch(h7d,h3d,h2d,p6d)
      external ccsd_t_transpt2_7326
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7326(v2sub,scratch,h3d,h2d,p6d,h7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,scratch)
!$omp&            private(p4,p5,p6,h1,h2,h3,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p6=jstart,jstop
      do p5=1,p5d
      do h2=1,h2d
!dec$ loop count max=40, min=20
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do h7=1,h7d
       triplesx(h3,h1,h2,p5,p6-jstart+1,p4-istart+1)=
     1     triplesx(h3,h1,h2,p5,p6-jstart+1,p4-istart+1)
     1   - t2sub(h7,p4,p5,h1)*scratch(h7,h3,h2,p6)
c     1   - t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_tx1_9(h3d,h2d,h1d,p6d,p5d,p4d,
     1               h7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_tx1_9'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,h7d
      integer h3h2,h1,p6,p5,p4,h7
      double precision triplesx(h1d,h3d*h2d,p5d,maxp5,maxp4)
      double precision t2sub(h7d,p4d,p5d,h1d)
      double precision v2sub(h3d,h2d,p6d,h7d)
      double precision scratch(h7d,h3d*h2d,p6d)
      external ccsd_t_transpt2_7326
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7326(v2sub,scratch,h3d,h2d,p6d,h7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,t2sub,scratch)
!$omp&            private(p4,p5,p6,h1,h3h2,h7)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,h1d,h2d,h3d,h7d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p6=jstart,jstop
      do p5=1,p5d
!dec$ loop count max=900, min=30
      do h3h2=1,h2d*h3d
!dec$ loop count max=40, min=20
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h7=1,h7d
       triplesx(h1,h3h2,p5,p6-jstart+1,p4-istart+1)=
     1     triplesx(h1,h3h2,p5,p6-jstart+1,p4-istart+1)
     1   + t2sub(h7,p4,p5,h1)*scratch(h7,h3h2,p6)
c     1   + t2sub(h7,p4,p5,h1)*v2sub(h3,h2,p6,h7)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_1(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      character(len=*), parameter :: funcname = 'sd_ty1_1'
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(p7d,h2d,h1d,p4d)
      external ccsd_t_transpt2_7124
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7124(t2sub,scratch,p4d,h1d,h2d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h2=1,h2d
!dec$ loop count max=40, min=20
!dec$ unroll_and_jam(4)
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do p7=1,p7d
       triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h3,h2,h1,p6,p5-jstart+1,p4-istart+1)
     1   - scratch(p7,h2,h1,p4)*v2sub(p7,h3,p6,p5)
c     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_2(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_2'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2h1,p6,p5,p4,p7
      double precision triplesx(h2d*h1d,h3d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(h2d*h1d,p4d,p7d)
      external ccsd_t_transpt2
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2(t2sub,scratch,h2d,h1d,p4d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h2h1,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do p7=1,p7d
!dec$ loop count max=600, min=30
      do h2h1=1,h1d*h2d
       triplesx(h2h1,h3,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h2h1,h3,p6,p5-jstart+1,p4-istart+1)
     1   - scratch(h2h1,p4,p7)*v2sub(p7,h3,p6,p5)
c     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_3(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_3'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(p7d,h2d,h1d,p4d)
      external ccsd_t_transpt2_7124
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7124(t2sub,scratch,p4d,h1d,h2d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p4=istart,istop
      do p5=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
      do h3=1,h3d
!dec$ loop count max=40, min=20
!dec$ unroll_and_jam(4)
      do h2=1,h2d
!dec$ loop count max=40, min=20
      do p7=1,p7d
       triplesx(h2,h3,h1,p6,p5-jstart+1,p4-istart+1)=
     1     triplesx(h2,h3,h1,p6,p5-jstart+1,p4-istart+1)
     1   + scratch(p7,h2,h1,p4)*v2sub(p7,h3,p6,p5)
c     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_4(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_4'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(p7d,h2d,h1d,p4d)
      external ccsd_t_transpt2_7124
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7124(t2sub,scratch,p4d,h1d,h2d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h2=1,h2d
!dec$ loop count max=40, min=20
!dec$ unroll_and_jam(4)
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do p7=1,p7d
       triplesx(h3,h2,h1,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h3,h2,h1,p6,p4-jstart+1,p5-istart+1)
     1   + scratch(p7,h2,h1,p4)*v2sub(p7,h3,p6,p5)
c     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_5(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_5'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2h1,p6,p5,p4,p7
      double precision triplesx(h2d*h1d,h3d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(h2d*h1d,p4d,p7d)
      external ccsd_t_transpt2
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2(t2sub,scratch,h2d,h1d,p4d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h2h1,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h3=1,h3d
      do p7=1,p7d
      do h2h1=1,h2d*h1d
       triplesx(h2h1,h3,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h2h1,h3,p6,p4-jstart+1,p5-istart+1)
     1   + scratch(h2h1,p4,p7)*v2sub(p7,h3,p6,p5)
c     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_6(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_6'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p6d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(p7d,h2d,h1d,p4d)
      external ccsd_t_transpt2_7124
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7124(t2sub,scratch,p4d,h1d,h2d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p4=jstart,jstop
      do p6=1,p6d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h3=1,h3d
!dec$ loop count max=40, min=20
c!dec$ unroll_and_jam(4)
      do h2=1,h2d
!dec$ loop count max=40, min=20
      do p7=1,p7d
       triplesx(h2,h3,h1,p6,p4-jstart+1,p5-istart+1)=
     1     triplesx(h2,h3,h1,p6,p4-jstart+1,p5-istart+1)
     1   - scratch(p7,h2,h1,p4)*v2sub(p7,h3,p6,p5)
c     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_7(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_7'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h3d,h2d,h1d,p4d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(p7d,h2d,h1d,p4d)
      external ccsd_t_transpt2_7124
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7124(t2sub,scratch,p4d,h1d,h2d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h2=1,h2d
!dec$ loop count max=40, min=20
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do p7=1,p7d
       triplesx(h3,h2,h1,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h3,h2,h1,p4,p6-jstart+1,p5-istart+1)
     1   - scratch(p7,h2,h1,p4)*v2sub(p7,h3,p6,p5)
c     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_8(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_8'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2h1,p6,p5,p4,p7
      double precision triplesx(h2d*h1d,h3d,p4d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(h2d*h1d,p4d,p7d)
      external ccsd_t_transpt2
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2(t2sub,scratch,h2d,h1d,p4d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h2h1,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do p7=1,p7d
!dec$ loop count max=600, min=40
      do h2h1=1,h1d*h2d
       triplesx(h2h1,h3,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h2h1,h3,p4,p6-jstart+1,p5-istart+1)
     1   - scratch(h2h1,p4,p7)*v2sub(p7,h3,p6,p5)
c     1   - t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
c
      subroutine sd_ty1_9(h3d,h2d,h1d,p6d,p5d,p4d,
     1               p7d,
     2               triplesx,t2sub,v2sub,
     3 istart,istop,jstart,jstop,maxp4,maxp5)
      IMPLICIT NONE
      character(len=*), parameter :: funcname = 'sd_ty1_9'
      INTEGER istart,istop,jstart,jstop,maxp4,maxp5
      integer h3d,h2d,h1d,p6d,p5d,p4d,p7d
      integer h3,h2,h1,p6,p5,p4,p7
      double precision triplesx(h2d,h3d,h1d,p4d,maxp5,maxp4)
      double precision t2sub(p7d,p4d,h1d,h2d)
      double precision v2sub(p7d,h3d,p6d,p5d)
      double precision scratch(p7d,h2d,h1d,p4d)
      external ccsd_t_transpt2_7124
#include "ccsd_t_dts_kernels_omp_top.fh"
      call ccsd_t_transpt2_7124(t2sub,scratch,p4d,h1d,h2d,p7d)
#ifdef USE_OPENMP
!$omp parallel do default(none)
!$omp&            shared(triplesx,scratch,v2sub)
!$omp&            private(p4,p5,p6,p7,h1,h2,h3)
!$omp&            firstprivate(istart,istop,jstart,jstop)
!$omp&            firstprivate(p4d,p5d,p6d,p7d,h1d,h2d,h3d)
!$omp&            collapse(OMPCOLLAPSE)
!$omp&            schedule(static)
#endif
      do p5=istart,istop
      do p6=jstart,jstop
      do p4=1,p4d
      do h1=1,h1d
!dec$ loop count max=40, min=20
      do h3=1,h3d
!dec$ loop count max=40, min=20
      do h2=1,h2d
!dec$ loop count max=40, min=20
      do p7=1,p7d
       triplesx(h2,h3,h1,p4,p6-jstart+1,p5-istart+1)=
     1     triplesx(h2,h3,h1,p4,p6-jstart+1,p5-istart+1)
     1   + scratch(p7,h2,h1,p4)*v2sub(p7,h3,p6,p5)
c     1   + t2sub(p7,p4,h1,h2)*v2sub(p7,h3,p6,p5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#include "ccsd_t_dts_kernels_omp_bottom.fh"
      return
      end
