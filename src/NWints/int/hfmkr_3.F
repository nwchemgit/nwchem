      subroutine hfmkr_3(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id: hfmkr_3.F,v 1.2 1994-09-22 13:35:03 d3e129 Exp $
c
      integer NAB, Lr, Lr3
      integer IJK(0:Lr,0:Lr,0:Lr)
      double precision PQ(NAB,3)
      double precision R(NAB,0:Lr,Lr3)
c
      integer m
      double precision x,y,z
      double precision x_2, y_2, z_2
      double precision x_3, y_3, z_3
      double precision R11,R21,R31
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
      IJK(3,0,0)=11
      IJK(2,1,0)=12
      IJK(2,0,1)=13
      IJK(1,2,0)=14
      IJK(1,1,1)=15
      IJK(1,0,2)=16
      IJK(0,3,0)=17
      IJK(0,2,1)=18
      IJK(0,1,2)=19
      IJK(0,0,3)=20
c
      do 00100 m=1,NAB
        X=PQ(m,1)
        Y=PQ(m,2)
        Z=PQ(m,3)
        X_2=X*X
        Y_2=Y*Y
        Z_2=Z*Z
        X_3=X_2*X
        Y_3=Y_2*Y
        Z_3=Z_2*Z
        R11=R(m,1,1)
        R21=R(m,2,1)
        R31=R(m,3,1)
        R(M,0,2)=X*R11
        R(M,0,3)=Y*R11
        R(M,0,4)=Z*R11
        R(M,0,5)=X_2*R21+R11
        R(M,0,6)=X*Y*R21
        R(M,0,7)=X*Z*R21
        R(M,0,8)=Y_2*R21+R11
        R(M,0,9)=Y*Z*R21
        R(M,0,10)=Z_2*R21+R11
        R(M,0,11)=X_3*R31+3.0d00*X*R21
        R(M,0,12)=X_2*Y*R31+Y*R21
        R(M,0,13)=X_2*Z*R31+Z*R21
        R(M,0,14)=X*Y_2*R31+X*R21
        R(M,0,15)=X*Y*Z*R31
        R(M,0,16)=X*Z_2*R31+X*R21
        R(M,0,17)=Y_3*R31+3.0d00*Y*R21
        R(M,0,18)=Y_2*Z*R31+Z*R21
        R(M,0,19)=Y*Z_2*R31+Y*R21
        R(M,0,20)=Z_3*R31+3.0d00*Z*R21
00100 continue
c
      end
