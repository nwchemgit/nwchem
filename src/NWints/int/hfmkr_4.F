      subroutine hfmkr_4(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id: hfmkr_4.F,v 1.2 1994-09-22 13:35:04 d3e129 Exp $
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
      double precision x_4, y_4, z_4
      double precision R11,R21,R31,R41
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
      IJK(4,0,0)=21
      IJK(3,1,0)=22
      IJK(3,0,1)=23
      IJK(2,2,0)=24
      IJK(2,1,1)=25
      IJK(2,0,2)=26
      IJK(1,3,0)=27
      IJK(1,2,1)=28
      IJK(1,1,2)=29
      IJK(1,0,3)=30
      IJK(0,4,0)=31
      IJK(0,3,1)=32
      IJK(0,2,2)=33
      IJK(0,1,3)=34
      IJK(0,0,4)=35
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
        X_4=X_3*X
        Y_4=Y_3*Y
        Z_4=Z_3*Z
c
        R11=R(m,1,1)
        R21=R(m,2,1)
        R31=R(m,3,1)
        R41=R(m,4,1)
c
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
        R(M,0,21)=X_4*R41+6.0d00*X_2*R31+3.0d00*R21
        R(M,0,22)=X_3*Y*R41+3.0d00*X*Y*R31
        R(M,0,23)=X_3*Z*R41+3.0d00*X*Z*R31
        R(M,0,24)=X_2*Y_2*R41+(X_2+Y_2)*R31+R21
        R(M,0,25)=X_2*Y*Z*R41+Y*Z*R31
        R(M,0,26)=X_2*Z_2*R41+(X_2+Z_2)*R31+R21
        R(M,0,27)=X*Y_3*R41+3.0d00*X*Y*R31
        R(M,0,28)=X*Y_2*Z*R41+X*Z*R31
        R(M,0,29)=X*Y*Z_2*R41+X*Y*R31
        R(M,0,30)=X*Z_3*R41+3.0d00*X*Z*R31
        R(M,0,31)=Y_4*R41+6.0d00*Y_2*R31+3.0d00*R21
        R(M,0,32)=Y_3*Z*R41+3.0d00*Y*Z*R31
        R(M,0,33)=Y_2*Z_2*R41+(Y_2+Z_2)*R31+R21
        R(M,0,34)=Y*Z_3*R41+3.0d00*Y*Z*R31
        R(M,0,35)=Z_4*R41+6.0d00*Z_2*R31+3.0d00*R21
00100 continue
c
      end
