      subroutine hfmkr_1(R,PQ,IJK,NAB,Lr,Lr3)
c $Id: hfmkr_1.F,v 1.1 1994-08-12 00:43:34 d3e129 Exp $
      implicit none
c
      integer NAB, Lr, Lr3
      integer IJK(0:Lr,0:Lr,0:Lr)
      double precision PQ(NAB,3)
      double precision R(NAB,0:Lr,Lr3)
c
      integer m
      double precision x,y,z
#if defined(INT_TIMING)
#include "int_tim.fh"
#include "pstat.fh"      
#endif

#if defined(INT_TIMING)
      call pstat_on(hand_int(1))
#endif
c
      IJK(0,0,0)=1
      if(Lr.eq.0) goto 00001
      IJK(1,0,0)=2
      IJK(0,1,0)=3
      IJK(0,0,1)=4
c
00001 continue
      if(Lr.eq.0) goto 00002
      do 00100 m=1,NAB
        X=PQ(m,1)
        Y=PQ(m,2)
        Z=PQ(m,3)
        R(M,0,2)=X*R(m,1,1)
        R(M,0,3)=Y*R(m,1,1)
        R(M,0,4)=Z*R(m,1,1)
00100 continue
c
00002 continue
#if defined(INT_TIMING)
      call pstat_off(hand_int(1))
#endif
      end
