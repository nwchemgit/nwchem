C$Id$

c     Modifications (c) 2021 Advanced Micro Devices, Inc. All Rights Reserved.

#ifdef USE_OMP_TEAMS_DISTRIBUTE
#define TEAMS_DISTRIBUTE teams distribute
#else
#define TEAMS_DISTRIBUTE
#endif

#ifdef USE_OMP_SCHEDULE_STATIC_1
#define OMP_SCHEDULE schedule(static,1)
#else
#define OMP_SCHEDULE
#endif

#ifdef OFFLOAD_CODE
      subroutine offl_ccsd_t_dot(
#else
      subroutine ccsd_t_dot(
#endif
     &                      a_singles, a_doubles,
     &                      restricted,
     &                      h1b,h2b,h3b,p4b,p5b,p6b,
     &                      o_h1,o_h2,o_h3,
     &                      o_p4,o_p5,o_p6,
     &                      r_h1,r_h2,r_h3,
     &                      r_p4,r_p5,r_p6,
     &                      energy1,energy2)
C
C     $Id$
C
      implicit none
      integer h1b, h2b, h3b, p4b, p5b, p6b
      integer r_h1,r_h2,r_h3
      integer r_p4,r_p5,r_p6
      double precision o_h1(r_h1),o_h2(r_h2),o_h3(r_h3)
      double precision o_p4(r_p4),o_p5(r_p5),o_p6(r_p6)
      double precision a_singles(r_h3,r_h2,r_h1,r_p6,r_p5,r_p4)
      double precision a_doubles(r_h3,r_h2,r_h1,r_p6,r_p5,r_p4)
      logical restricted
      double precision energy1,energy2
c 
      integer h1, h2, h3, p4, p5, p6
      double precision factor,denom,sing,doub,denom_0
      double precision d_p4,d_p5,d_p6
      double precision d_h1,d_h2,d_h3
c
c factor = [ 1/36, 1/18, 1/12, 1/6, 1/4, 1/3, 1/2, 1, 2]
c

      if (restricted) then
        factor = 2.0d0
      else
        factor = 1.0d0
      endif
      if ((p4b.eq.p5b).and.(p5b.eq.p6b)) then
        factor = factor / 6.0d0
      else if ((p4b.eq.p5b).or.(p5b.eq.p6b)) then
        factor = factor / 2.0d0
      endif
      if ((h1b.eq.h2b).and.(h2b.eq.h3b)) then
        factor = factor / 6.0d0
      else if ((h1b.eq.h2b).or.(h2b.eq.h3b)) then
        factor = factor / 2.0d0
      endif

#ifdef TCE_OPENACC
!$acc  parallel loop collapse(3)
!$acc& private(p4,p5,p6,h1,h2,h3,sing,doub,denom)
!$acc& private(d_p4,d_p5,d_p6,d_h1,d_h2,d_h3)
!$acc& reduction(+:energy1) reduction(+:energy2)

#else
#ifdef OFFLOAD_CODE
#ifdef VERY_VERBOSE
      write (*,'(A)') "X"
#endif
!$omp target map(tofrom:o_p4) map(tofrom:o_p5)
!$omp& map(tofrom:o_p6) map(tofrom:o_h1)
!$omp& map(tofrom:o_h2) map(tofrom:o_h3)
!$omp& map(tofrom:energy1) map(tofrom:energy2)
c
!$omp TEAMS_DISTRIBUTE parallel do
!$omp& collapse(3) OMP_SCHEDULE
!$omp& default(shared)
!$omp& private(p4,p5,p6,h1,h2,h3,sing,doub,denom)
!$omp& private(d_p4,d_p5,d_p6,d_h1,d_h2,d_h3)
!$omp& reduction(+:energy1) reduction(+:energy2)
#else
#ifdef USE_OPENMP
!$omp  parallel do collapse(3)
!$omp& schedule(static)
!$omp& default(shared)
!$omp& private(p4,p5,p6,h1,h2,h3,sing,doub,denom)
!$omp& private(d_p4,d_p5,d_p6,d_h1,d_h2,d_h3)
!$omp& reduction(+:energy1) reduction(+:energy2)
#endif
#endif
#endif
      do p4 = 1, r_p4
       do p5 = 1, r_p5
        do p6 = 1, r_p6
         d_p4 = o_p4(p4)
         d_p5 = o_p5(p5)
         d_p6 = o_p6(p6)
         denom_0=-(d_p4+d_p5+d_p6)
         do h1 = 1, r_h1
          d_h1 = o_h1(h1)
          do h2 = 1, r_h2
           d_h2 = o_h2(h2)
           do h3 = 1, r_h3
            d_h3 = o_h3(h3)
            sing = a_singles(h3,h2,h1,p6,p5,p4)
            doub = a_doubles(h3,h2,h1,p6,p5,p4)
            denom = doub * factor / ( d_h1+d_h2+d_h3+denom_0)
            energy1 = energy1 + denom*doub
            energy2 = energy2 + denom*(doub+sing)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
#ifdef TCE_OPENACC
!$acc end parallel loop
#else
#ifdef OFFLOAD_CODE
!$omp end TEAMS_DISTRIBUTE parallel do
!$omp end target
#else
#ifdef USE_OPENMP
!$omp end parallel do
#endif
#endif
#endif
      return
      end
#ifndef OFFLOAD_CODE
#ifdef USE_OFFLOAD
#define OFFLOAD_CODE
c
c     Compile source again for offload code
c
#include "ccsd_t_dot.F"
#endif
#endif
