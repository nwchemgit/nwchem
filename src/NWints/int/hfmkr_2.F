      subroutine hfmkr_2(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id: hfmkr_2.F,v 1.1 1994-08-12 00:43:35 d3e129 Exp $
c
      integer NAB, Lr, Lr3
      integer IJK(0:Lr,0:Lr,0:Lr)
      double precision PQ(NAB,3)
      double precision R(NAB,0:Lr,Lr3)
c
      integer m
      double precision x,y,z
      double precision R21, R11
#if defined(INT_TIMING)
#include "int_tim.fh"
#include "pstat.fh"      
#endif

#if defined(INT_TIMING)
      call pstat_on(hand_int(2))
#endif
c
      IJK(0,0,0)=1
      IJK(1,0,0)=2
      IJK(0,1,0)=3
      IJK(0,0,1)=4
      IJK(2,0,0)=5
      IJK(1,1,0)=6
      IJK(1,0,1)=7
      IJK(0,2,0)=8
      IJK(0,1,1)=9
      IJK(0,0,2)=10
c
      do 00100 m=1,NAB
        X=PQ(m,1)
        Y=PQ(m,2)
        Z=PQ(m,3)
        R11=R(m,1,1)
        R21=R(m,2,1)
        R(M,0,2) =X*R11
        R(M,0,3) =Y*R11
        R(M,0,4) =Z*R11
        R(M,0,5) =X*X*R21+R11
        R(M,0,6) =X*Y*R21
        R(M,0,7) =X*Z*R21
        R(M,0,8) =Y*Y*R21+R11
        R(M,0,9) =Y*Z*R21
        R(M,0,10)=Z*Z*R21+R11
00100 continue
c
#if defined(INT_TIMING)
      call pstat_off(hand_int(2))
#endif
      end
