      Subroutine hfmkr(R,IJK,PQ,NAB,Lp,Lp3)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Dimension R(NAB,0:Lp,Lp3),IJK(0:Lp,0:Lp,0:Lp),PQ(NAB,3)

      Dimension Nxyz(3)
c
c Recursively define the auxiliary integrals needed for ERIs and NAIs.
c
c     Recursion Formula:
c
c        R        = PQx R            + (N-1) R
c         N,L,M,j        N-1,L,M,j+1          N-2,L,M,j+1
c
c        R        = PQy R            + (L-1) R
c         N,L,M,j        N,L-1,M,j+1          N,L-2,M,j+1
c
c        R        = PQz R            + (M-1) R
c         N,L,M,j        N,L,M-1,j+1          N,L,M-2,j+1
c
c******************************************************************************

#if defined(USE_SPECIAL_MAKER)
      if (Lp.le.1) then
        call hfmkr_1(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif(Lp.eq.2) then
        call hfmkr_2(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif(Lp.eq.3) then
        call hfmkr_3(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif(Lp.eq.4) then
        call hfmkr_4(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif(Lp.eq.5) then
        call hfmkr_5(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif(Lp.eq.6) then
        call hfmkr_6(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif(Lp.eq.7) then
        call hfmkr_7(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif(Lp.eq.8) then
        call hfmkr_8(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      elseif (Lp.eq.9) then
        call hfmkr_9(R,PQ,IJK,NAB,Lp,Lp3)
        goto 00099
      endif
#endif

#if defined(INT_DEBUG)
      write(6,*)'lp',lp,'      lp3',lp3,'     nab',nab
      write(6,*)' PQ coords start'
      call intintp(PQ,nab,3,'hfmkr pq')
      write(6,*)' PQ coords end'
      write(6,*)' r b4 recursion'
      do ii = 1,NAB
        do jj = 0,Lp
          do kk = 1,lp3
            if (abs(R(ii,jj,kk)).gt.1.0d-12) then
              write(6,*)'R (',ii,jj,kk,') =',R(ii,jj,kk)
            endif
          enddo
        enddo
      enddo
#endif
      
      n0 = 1

      IJK(0,0,0) = n0

      do 180 l = 1,Lp

       l2 = ((l+1)*(l+2))/2

       do 170 i = 1,l2

        call getNxyz(l,i,Nxyz)

        Ip = Nxyz(1)
        Jp = Nxyz(2)
        Kp = Nxyz(3)

        n0 = n0 + 1

        IJK(Ip,Jp,Kp) = n0

        if( Ip.gt.0 )then

         if( Ip.eq.1 )then

          n1 = IJK(Ip-1,Jp,Kp)

          do 115 j = 0,Lp-l
           do 110 m = 1,NAB
            R(m,j,n0) = PQ(m,1)*R(m,j+1,n1)
  110      continue
  115     continue

         else

          n1 = IJK(Ip-1,Jp,Kp)
          n2 = IJK(Ip-2,Jp,Kp)

          do 125 j = 0,Lp-l
           do 120 m = 1,NAB
            R(m,j,n0) = PQ(m,1)*R(m,j+1,n1) + (Ip-1)*R(m,j+1,n2)
  120      continue
  125     continue

         end if

        else if( Jp.gt.0 )then

         if( Jp.eq.1 )then

          n1 = IJK(Ip,Jp-1,Kp)

          do 135 j = 0,Lp-l
           do 130 m = 1,NAB
            R(m,j,n0) = PQ(m,2)*R(m,j+1,n1)
  130      continue
  135     continue

         else

          n1 = IJK(Ip,Jp-1,Kp)
          n2 = IJK(Ip,Jp-2,Kp)

          do 145 j = 0,Lp-l
           do 140 m = 1,NAB
            R(m,j,n0) = PQ(m,2)*R(m,j+1,n1) + (Jp-1)*R(m,j+1,n2)
  140      continue
  145     continue

         end if

        else if( Kp.gt.0 )then

         if( Kp.eq.1 )then

          n1 = IJK(Ip,Jp,Kp-1)

          do 155 j = 0,Lp-l
           do 150 m = 1,NAB
            R(m,j,n0) = PQ(m,3)*R(m,j+1,n1)
  150      continue
  155     continue

         else

          n1 = IJK(Ip,Jp,Kp-1)
          n2 = IJK(Ip,Jp,Kp-2)

          do 165 j = 0,Lp-l
           do 160 m = 1,NAB
            R(m,j,n0) = PQ(m,3)*R(m,j+1,n1) + (Kp-1)*R(m,j+1,n2)
  160      continue
  165     continue

         end if

        end if

  170  continue

  180 continue

00099 continue

#if defined(INT_DEBUG)
      write(6,*)' r after recursion'
      do ii = 1,NAB
        do jj = 0,Lp
          do kk = 1,lp3
            if (abs(R(ii,jj,kk)).gt.1.0d-12) then
              write(6,*)'R (',ii,jj,kk,') =',R(ii,jj,kk)
            endif
          enddo
        enddo
      enddo
#endif
      end


#if defined(USE_SPECIAL_MAKER)
      subroutine hfmkr_1(R,PQ,IJK,NAB,Lr,Lr3)
c $Id$
      implicit none
c
      integer NAB, Lr, Lr3
      integer IJK(0:Lr,0:Lr,0:Lr)
      double precision PQ(NAB,3)
      double precision R(NAB,0:Lr,Lr3)
c
      integer m
      double precision x,y,z
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
c
      end
      subroutine hfmkr_2(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
c
      integer NAB, Lr, Lr3
      integer IJK(0:Lr,0:Lr,0:Lr)
      double precision PQ(NAB,3)
      double precision R(NAB,0:Lr,Lr3)
c
      integer m
      double precision x,y,z
      double precision R21, R11
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
      end
      subroutine hfmkr_3(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
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
      subroutine hfmkr_4(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
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
      subroutine hfmkr_5(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
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
      double precision x_5, y_5, z_5
      double precision R11,R21,R31,R41,R51
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
      IJK(5,0,0)=36
      IJK(4,1,0)=37
      IJK(4,0,1)=38
      IJK(3,2,0)=39
      IJK(3,1,1)=40
      IJK(3,0,2)=41
      IJK(2,3,0)=42
      IJK(2,2,1)=43
      IJK(2,1,2)=44
      IJK(2,0,3)=45
      IJK(1,4,0)=46
      IJK(1,3,1)=47
      IJK(1,2,2)=48
      IJK(1,1,3)=49
      IJK(1,0,4)=50
      IJK(0,5,0)=51
      IJK(0,4,1)=52
      IJK(0,3,2)=53
      IJK(0,2,3)=54
      IJK(0,1,4)=55
      IJK(0,0,5)=56
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
        X_5=X_4*X
        Y_5=Y_4*Y
        Z_5=Z_4*Z
c
        R11=R(m,1,1)
        R21=R(m,2,1)
        R31=R(m,3,1)
        R41=R(m,4,1)
        R51=R(m,5,1)
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
        R(M,0,36)=X_5*R51+10.0d00* X_3*R41+15.0d00*X*R31
        R(M,0,37)=X_4*Y*R51+6.0d00*X_2*Y*R41+3.0d00*Y*R31
        R(M,0,38)=X_4*Z*R51+6.0d00*X_2*Z*R41+3.0d00*Z*R31
        R(M,0,39)=X_3*Y_2*R51+(X*(X_2+Y_2)+2.0d00*X*Y_2)*R41
     &      +3.0d00*X*R31
        R(M,0,40)=X_3*Y*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,41)=X_3*Z_2*R51+(X*(X_2+Z_2)+2.0d00*X*Z_2)*R41
     &      +3.0d00*X*R31
        R(M,0,42)=X_2*Y_3*R51+(3.0d00*X_2*Y+Y_3)*R41
     &      +3.0d00*Y*R31
        R(M,0,43)=X_2*Y_2*Z*R51+(X_2*Z+Y_2*Z)*R41
     &      +Z*R31
        R(M,0,44)=X_2*Y*Z_2*R51+(X_2*Y+Y*Z_2)*R41
     &      +Y*R31
        R(M,0,45)=X_2*Z_3*R51+(3.0d00*X_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,46)=X*Y_4*R51+6.0d00*X*Y_2*R41+3.0d00*X*R31
        R(M,0,47)=X*Y_3*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,48)=X*Y_2*Z_2*R51+X*(Y_2+Z_2)*R41+
     &      X*R31
        R(M,0,49)=X*Y*Z_3*R51+3.0d00*X*Y*Z*R41
        R(M,0,50)=X*Z_4*R51+6.0d00*X*Z_2*R41+3.0d00*X*R31
        R(M,0,51)=Y_5*R51+10.0d00* Y_3*R41+15.0d00*Y*R31
        R(M,0,52)=Y_4*Z*R51+6.0d00*Y_2*Z*R41+3.0d00*Z*R31
        R(M,0,53)=Y_3*Z_2*R51+(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R41
     &      +3.0d00*Y*R31
        R(M,0,54)=Y_2*Z_3*R51+(3.0d00*Y_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,55)=Y*Z_4*R51+6.0d00*Y*Z_2*R41+3.0d00*Y*R31
        R(M,0,56)=Z_5*R51+10.0d00* Z_3*R41+15.0d00*Z*R31
00100 continue
c
      end
      subroutine hfmkr_6(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
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
      double precision x_5, y_5, z_5
      double precision x_6, y_6, z_6
      double precision R11,R21,R31,R41,R51,R61
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
      IJK(5,0,0)=36
      IJK(4,1,0)=37
      IJK(4,0,1)=38
      IJK(3,2,0)=39
      IJK(3,1,1)=40
      IJK(3,0,2)=41
      IJK(2,3,0)=42
      IJK(2,2,1)=43
      IJK(2,1,2)=44
      IJK(2,0,3)=45
      IJK(1,4,0)=46
      IJK(1,3,1)=47
      IJK(1,2,2)=48
      IJK(1,1,3)=49
      IJK(1,0,4)=50
      IJK(0,5,0)=51
      IJK(0,4,1)=52
      IJK(0,3,2)=53
      IJK(0,2,3)=54
      IJK(0,1,4)=55
      IJK(0,0,5)=56
      IJK(6,0,0)=57
      IJK(5,1,0)=58
      IJK(5,0,1)=59
      IJK(4,2,0)=60
      IJK(4,1,1)=61
      IJK(4,0,2)=62
      IJK(3,3,0)=63
      IJK(3,2,1)=64
      IJK(3,1,2)=65
      IJK(3,0,3)=66
      IJK(2,4,0)=67
      IJK(2,3,1)=68
      IJK(2,2,2)=69
      IJK(2,1,3)=70
      IJK(2,0,4)=71
      IJK(1,5,0)=72
      IJK(1,4,1)=73
      IJK(1,3,2)=74
      IJK(1,2,3)=75
      IJK(1,1,4)=76
      IJK(1,0,5)=77
      IJK(0,6,0)=78
      IJK(0,5,1)=79
      IJK(0,4,2)=80
      IJK(0,3,3)=81
      IJK(0,2,4)=82
      IJK(0,1,5)=83
      IJK(0,0,6)=84
c
      do 00100 m=1,NAB
        X=PQ(m,1)
        Y=PQ(m,2)
        Z=PQ(m,3)
        R11=R(m,1,1)
        R(M,0,2)=X*R11
        R(M,0,3)=Y*R11
        R(M,0,4)=Z*R11
        X_2=X*X
        Y_2=Y*Y
        Z_2=Z*Z
        R21=R(m,2,1)
        R(M,0,5)=X_2*R21+R11
        R(M,0,6)=X*Y*R21
        R(M,0,7)=X*Z*R21
        R(M,0,8)=Y_2*R21+R11
        R(M,0,9)=Y*Z*R21
        R(M,0,10)=Z_2*R21+R11
        X_3=X_2*X
        Y_3=Y_2*Y
        Z_3=Z_2*Z
        R31=R(m,3,1)
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
        X_4=X_3*X
        Y_4=Y_3*Y
        Z_4=Z_3*Z
        R41=R(m,4,1)
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
        X_5=X_4*X
        Y_5=Y_4*Y
        Z_5=Z_4*Z
        R51=R(m,5,1)
        R(M,0,36)=X_5*R51+10.0d00* X_3*R41+15.0d00*X*R31
        R(M,0,37)=X_4*Y*R51+6.0d00*X_2*Y*R41+3.0d00*Y*R31
        R(M,0,38)=X_4*Z*R51+6.0d00*X_2*Z*R41+3.0d00*Z*R31
        R(M,0,39)=X_3*Y_2*R51+(X*(X_2+Y_2)+2.0d00*X*Y_2)*R41
     &      +3.0d00*X*R31
        R(M,0,40)=X_3*Y*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,41)=X_3*Z_2*R51+(X*(X_2+Z_2)+2.0d00*X*Z_2)*R41
     &      +3.0d00*X*R31
        R(M,0,42)=X_2*Y_3*R51+(3.0d00*X_2*Y+Y_3)*R41
     &      +3.0d00*Y*R31
        R(M,0,43)=X_2*Y_2*Z*R51+(X_2*Z+Y_2*Z)*R41
     &      +Z*R31
        R(M,0,44)=X_2*Y*Z_2*R51+(X_2*Y+Y*Z_2)*R41
     &      +Y*R31
        R(M,0,45)=X_2*Z_3*R51+(3.0d00*X_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,46)=X*Y_4*R51+6.0d00*X*Y_2*R41+3.0d00*X*R31
        R(M,0,47)=X*Y_3*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,48)=X*Y_2*Z_2*R51+X*(Y_2+Z_2)*R41+
     &      X*R31
        R(M,0,49)=X*Y*Z_3*R51+3.0d00*X*Y*Z*R41
        R(M,0,50)=X*Z_4*R51+6.0d00*X*Z_2*R41+3.0d00*X*R31
        R(M,0,51)=Y_5*R51+10.0d00* Y_3*R41+15.0d00*Y*R31
        R(M,0,52)=Y_4*Z*R51+6.0d00*Y_2*Z*R41+3.0d00*Z*R31
        R(M,0,53)=Y_3*Z_2*R51+(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R41
     &      +3.0d00*Y*R31
        R(M,0,54)=Y_2*Z_3*R51+(3.0d00*Y_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,55)=Y*Z_4*R51+6.0d00*Y*Z_2*R41+3.0d00*Y*R31
        R(M,0,56)=Z_5*R51+10.0d00* Z_3*R41+15.0d00*Z*R31
        X_6=X_5 * X
        Y_6=Y_5 * Y
        Z_6=Z_5 * Z
        R61=R(m,6,1)
        R(M,0,57)=X_6*R61+15.0d00*X_4*R51+45.0d00*X_2*R41
     &      +15.0d00*R31
        R(M,0,58)=X_5*Y*R61+10.0d00* X_3*Y*R51+15.0d00*X*Y*R41
        R(M,0,59)=X_5*Z*R61+10.0d00* X_3*Z*R51+15.0d00*X*Z*R41
        R(M,0,60)=X_4*Y_2*R61
     &      +(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+2.0d00*X*Y_2))*R51
     &      +(6.0d00*X_2+3.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,61)=X_4*Y*Z*R61+6.0d00*X_2*Y*Z*R51+3.0d00*Y*Z*R41
        R(M,0,62)=X_4*Z_2*R61
     &      +(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+2.0d00*X*Z_2))*R51
     &      +(6.0d00*X_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,63)=X_3*Y_3*R61
     &      +(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3)*R51+9.0d00*X*Y*R41
        R(M,0,64)=X_3*Y_2*Z*R61
     &      +(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z)*R51+3.0d00*X*Z*R41
        R(M,0,65)=X_3*Y*Z_2*R61
     &      +(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,66)=X_3*Z_3*R61
     &      +(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3)*R51+9.0d00*X*Z*R41
        R(M,0,67)=X_2*Y_4*R61+(Y_4+6.0d00*X_2*Y_2)*R51
     &      +(3.0d00*X_2+6.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,68)=X_2*Y_3*Z*R61+(3.0d00*X_2*Y*Z+Y_3*Z)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,69)=X_2*Y_2*Z_2*R61
     &      +(Y_2*Z_2+X_2*(Y_2+Z_2))*R51
     &      +(X_2+Y_2+Z_2)*R41+R31
        R(M,0,70)=X_2*Y*Z_3*R61+(3.0d00*X_2*Y*Z+Y*Z_3)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,71)=X_2*Z_4*R61+(Z_4+6.0d00*X_2*Z_2)*R51
     &      +(3.0d00*X_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,72)=X*Y_5*R61+10.0d00* X*Y_3*R51+15.0d00*X*Y*R41
        R(M,0,73)=X*Y_4*Z*R61+6.0d00*X*Y_2*Z*R51+3.0d00*X*Z*R41
        R(M,0,74)=X*Y_3*Z_2*R61
     &      +X*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,75)=X*Y_2*Z_3*R61+X*(3.0d00*Y_2*Z+Z_3)*R51
     &      +3.0d00*X*Z*R41
        R(M,0,76)=X*Y*Z_4*R61+6.0d00*X*Y*Z_2*R51+3.0d00*X*Y*R41
        R(M,0,77)=X*Z_5*R61+10.0d00* X*Z_3*R51+15.0d00*X*Z*R41
        R(M,0,78)=Y_6*R61+15.0d00*Y_4*R51+45.0d00*Y_2*R41
     &      +15.0d00*R31
        R(M,0,79)=Y_5*Z*R61+10.0d00* Y_3*Z*R51+15.0d00*Y*Z*R41
        R(M,0,80)=Y_4*Z_2*R61
     &      +(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))*R51
     &      +(6.0d00*Y_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,81)=Y_3*Z_3*R61
     &      +(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R51+9.0d00*Y*Z*R41
        R(M,0,82)=Y_2*Z_4*R61+(Z_4+6.0d00*Y_2*Z_2)*R51
     &      +(3.0d00*Y_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,83)=Y*Z_5*R61+10.0d00* Y*Z_3*R51+15.0d00*Y*Z*R41
        R(M,0,84)=Z_6*R61+15.0d00*Z_4*R51+45.0d00*Z_2*R41
     &      +15.0d00*R31
00100 continue
c
      end
      subroutine hfmkr_7(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
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
      double precision x_5, y_5, z_5
      double precision x_6, y_6, z_6
      double precision x_7, y_7, z_7
      double precision R11,R21,R31,R41,R51,R61,R71
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
      IJK(5,0,0)=36
      IJK(4,1,0)=37
      IJK(4,0,1)=38
      IJK(3,2,0)=39
      IJK(3,1,1)=40
      IJK(3,0,2)=41
      IJK(2,3,0)=42
      IJK(2,2,1)=43
      IJK(2,1,2)=44
      IJK(2,0,3)=45
      IJK(1,4,0)=46
      IJK(1,3,1)=47
      IJK(1,2,2)=48
      IJK(1,1,3)=49
      IJK(1,0,4)=50
      IJK(0,5,0)=51
      IJK(0,4,1)=52
      IJK(0,3,2)=53
      IJK(0,2,3)=54
      IJK(0,1,4)=55
      IJK(0,0,5)=56
      IJK(6,0,0)=57
      IJK(5,1,0)=58
      IJK(5,0,1)=59
      IJK(4,2,0)=60
      IJK(4,1,1)=61
      IJK(4,0,2)=62
      IJK(3,3,0)=63
      IJK(3,2,1)=64
      IJK(3,1,2)=65
      IJK(3,0,3)=66
      IJK(2,4,0)=67
      IJK(2,3,1)=68
      IJK(2,2,2)=69
      IJK(2,1,3)=70
      IJK(2,0,4)=71
      IJK(1,5,0)=72
      IJK(1,4,1)=73
      IJK(1,3,2)=74
      IJK(1,2,3)=75
      IJK(1,1,4)=76
      IJK(1,0,5)=77
      IJK(0,6,0)=78
      IJK(0,5,1)=79
      IJK(0,4,2)=80
      IJK(0,3,3)=81
      IJK(0,2,4)=82
      IJK(0,1,5)=83
      IJK(0,0,6)=84
      IJK(7,0,0)=85
      IJK(6,1,0)=86
      IJK(6,0,1)=87
      IJK(5,2,0)=88
      IJK(5,1,1)=89
      IJK(5,0,2)=90
      IJK(4,3,0)=91
      IJK(4,2,1)=92
      IJK(4,1,2)=93
      IJK(4,0,3)=94
      IJK(3,4,0)=95
      IJK(3,3,1)=96
      IJK(3,2,2)=97
      IJK(3,1,3)=98
      IJK(3,0,4)=99
      IJK(2,5,0)=100
      IJK(2,4,1)=101
      IJK(2,3,2)=102
      IJK(2,2,3)=103
      IJK(2,1,4)=104
      IJK(2,0,5)=105
      IJK(1,6,0)=106
      IJK(1,5,1)=107
      IJK(1,4,2)=108
      IJK(1,3,3)=109
      IJK(1,2,4)=110
      IJK(1,1,5)=111
      IJK(1,0,6)=112
      IJK(0,7,0)=113
      IJK(0,6,1)=114
      IJK(0,5,2)=115
      IJK(0,4,3)=116
      IJK(0,3,4)=117
      IJK(0,2,5)=118
      IJK(0,1,6)=119
      IJK(0,0,7)=120
c
      do 00100 m=1,NAB
        X=PQ(m,1)
        Y=PQ(m,2)
        Z=PQ(m,3)
        R11=R(m,1,1)
        R(M,0,2)=X*R11
        R(M,0,3)=Y*R11
        R(M,0,4)=Z*R11
        X_2=X*X
        Y_2=Y*Y
        Z_2=Z*Z
        R21=R(m,2,1)
        R(M,0,5)=X_2*R21+R11
        R(M,0,6)=X*Y*R21
        R(M,0,7)=X*Z*R21
        R(M,0,8)=Y_2*R21+R11
        R(M,0,9)=Y*Z*R21
        R(M,0,10)=Z_2*R21+R11
        X_3=X_2*X
        Y_3=Y_2*Y
        Z_3=Z_2*Z
        R31=R(m,3,1)
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
        X_4=X_3*X
        Y_4=Y_3*Y
        Z_4=Z_3*Z
        R41=R(m,4,1)
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
        X_5=X_4*X
        Y_5=Y_4*Y
        Z_5=Z_4*Z
        R51=R(m,5,1)
        R(M,0,36)=X_5*R51+10.0d00* X_3*R41+15.0d00*X*R31
        R(M,0,37)=X_4*Y*R51+6.0d00*X_2*Y*R41+3.0d00*Y*R31
        R(M,0,38)=X_4*Z*R51+6.0d00*X_2*Z*R41+3.0d00*Z*R31
        R(M,0,39)=X_3*Y_2*R51+(X*(X_2+Y_2)+2.0d00*X*Y_2)*R41
     &      +3.0d00*X*R31
        R(M,0,40)=X_3*Y*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,41)=X_3*Z_2*R51+(X*(X_2+Z_2)+2.0d00*X*Z_2)*R41
     &      +3.0d00*X*R31
        R(M,0,42)=X_2*Y_3*R51+(3.0d00*X_2*Y+Y_3)*R41
     &      +3.0d00*Y*R31
        R(M,0,43)=X_2*Y_2*Z*R51+(X_2*Z+Y_2*Z)*R41
     &      +Z*R31
        R(M,0,44)=X_2*Y*Z_2*R51+(X_2*Y+Y*Z_2)*R41
     &      +Y*R31
        R(M,0,45)=X_2*Z_3*R51+(3.0d00*X_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,46)=X*Y_4*R51+6.0d00*X*Y_2*R41+3.0d00*X*R31
        R(M,0,47)=X*Y_3*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,48)=X*Y_2*Z_2*R51+X*(Y_2+Z_2)*R41+
     &      X*R31
        R(M,0,49)=X*Y*Z_3*R51+3.0d00*X*Y*Z*R41
        R(M,0,50)=X*Z_4*R51+6.0d00*X*Z_2*R41+3.0d00*X*R31
        R(M,0,51)=Y_5*R51+10.0d00* Y_3*R41+15.0d00*Y*R31
        R(M,0,52)=Y_4*Z*R51+6.0d00*Y_2*Z*R41+3.0d00*Z*R31
        R(M,0,53)=Y_3*Z_2*R51+(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R41
     &      +3.0d00*Y*R31
        R(M,0,54)=Y_2*Z_3*R51+(3.0d00*Y_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,55)=Y*Z_4*R51+6.0d00*Y*Z_2*R41+3.0d00*Y*R31
        R(M,0,56)=Z_5*R51+10.0d00* Z_3*R41+15.0d00*Z*R31
        X_6=X_5 * X
        Y_6=Y_5 * Y
        Z_6=Z_5 * Z
        R61=R(m,6,1)
        R(M,0,57)=X_6*R61+15.0d00*X_4*R51+45.0d00*X_2*R41
     &      +15.0d00*R31
        R(M,0,58)=X_5*Y*R61+10.0d00* X_3*Y*R51+15.0d00*X*Y*R41
        R(M,0,59)=X_5*Z*R61+10.0d00* X_3*Z*R51+15.0d00*X*Z*R41
        R(M,0,60)=X_4*Y_2*R61
     &      +(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+2.0d00*X*Y_2))*R51
     &      +(6.0d00*X_2+3.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,61)=X_4*Y*Z*R61+6.0d00*X_2*Y*Z*R51+3.0d00*Y*Z*R41
        R(M,0,62)=X_4*Z_2*R61
     &      +(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+2.0d00*X*Z_2))*R51
     &      +(6.0d00*X_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,63)=X_3*Y_3*R61
     &      +(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3)*R51+9.0d00*X*Y*R41
        R(M,0,64)=X_3*Y_2*Z*R61
     &      +(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z)*R51+3.0d00*X*Z*R41
        R(M,0,65)=X_3*Y*Z_2*R61
     &      +(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,66)=X_3*Z_3*R61
     &      +(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3)*R51+9.0d00*X*Z*R41
        R(M,0,67)=X_2*Y_4*R61+(Y_4+6.0d00*X_2*Y_2)*R51
     &      +(3.0d00*X_2+6.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,68)=X_2*Y_3*Z*R61+(3.0d00*X_2*Y*Z+Y_3*Z)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,69)=X_2*Y_2*Z_2*R61
     &      +(Y_2*Z_2+X_2*(Y_2+Z_2))*R51
     &      +(X_2+Y_2+Z_2)*R41+R31
        R(M,0,70)=X_2*Y*Z_3*R61+(3.0d00*X_2*Y*Z+Y*Z_3)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,71)=X_2*Z_4*R61+(Z_4+6.0d00*X_2*Z_2)*R51
     &      +(3.0d00*X_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,72)=X*Y_5*R61+10.0d00* X*Y_3*R51+15.0d00*X*Y*R41
        R(M,0,73)=X*Y_4*Z*R61+6.0d00*X*Y_2*Z*R51+3.0d00*X*Z*R41
        R(M,0,74)=X*Y_3*Z_2*R61
     &      +X*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,75)=X*Y_2*Z_3*R61+X*(3.0d00*Y_2*Z+Z_3)*R51
     &      +3.0d00*X*Z*R41
        R(M,0,76)=X*Y*Z_4*R61+6.0d00*X*Y*Z_2*R51+3.0d00*X*Y*R41
        R(M,0,77)=X*Z_5*R61+10.0d00* X*Z_3*R51+15.0d00*X*Z*R41
        R(M,0,78)=Y_6*R61+15.0d00*Y_4*R51+45.0d00*Y_2*R41
     &      +15.0d00*R31
        R(M,0,79)=Y_5*Z*R61+10.0d00* Y_3*Z*R51+15.0d00*Y*Z*R41
        R(M,0,80)=Y_4*Z_2*R61
     &      +(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))*R51
     &      +(6.0d00*Y_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,81)=Y_3*Z_3*R61
     &      +(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R51+9.0d00*Y*Z*R41
        R(M,0,82)=Y_2*Z_4*R61+(Z_4+6.0d00*Y_2*Z_2)*R51
     &      +(3.0d00*Y_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,83)=Y*Z_5*R61+10.0d00* Y*Z_3*R51+15.0d00*Y*Z*R41
        R(M,0,84)=Z_6*R61+15.0d00*Z_4*R51+45.0d00*Z_2*R41
     &      +15.0d00*R31
        X_7=X_6 * X
        Y_7=Y_6 * Y
        Z_7=Z_6 * Z
        R71=R(m,7,1)
        R(M,0,85)=X_7*R71+21.0d00*X_5*R61+105.0d00*X_3*R51
     &      +105.0d00*X*R41
        R(M,0,86)=X_6*Y*R71+15.0d00*X_4*Y*R61+45.0d00*X_2*Y*R51
     &      +15.0d00*Y*R41
        R(M,0,87)=X_6*Z*R71+15.0d00*X_4*Z*R61+45.0d00*X_2*Z*R51
     &      +15.0d00*Z*R41
        R(M,0,88)=X_5*Y_2*R71
     &      +(X*(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+
     &      2.0d00*X*Y_2))+4.0d00*X_3*Y_2)*R61
     &      +(X*(6.0d00*X_2+3.0d00*Y_2)+
     &      4.0d00*X*(X_2+Y_2)+8.0d00*X*Y_2)*R51
     &      +15.0d00*X*R41
        R(M,0,89)=X_5*Y*Z*R71+10.0d00* X_3*Y*Z*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,90)=X_5*Z_2*R71
     &      +(X*(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+
     &      2.0d00*X*Z_2))+4.0d00*X_3*Z_2)*R61
     &      +(X*(6.0d00*X_2+3.0d00*Z_2)+
     &      4.0d00*X*(X_2+Z_2)+8.0d00*X*Z_2)*R51
     &      +15.0d00*X*R41
        R(M,0,91)=X_4*Y_3*R71
     &      +(3.0d00*X_2*Y_3+X*(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3))*R61
     &      +(18.0d00*X_2*Y+3.0d00*Y_3)*R51+9.0d00*Y*R41
        R(M,0,92)=X_4*Y_2*Z*R71
     &      +(3.0d00*X_2*Y_2*Z+X*(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z))*R61
     &      +(6.0d00*X_2*Z+3.0d00*Y_2*Z)*R51+3.0d00*Z*R41
        R(M,0,93)=X_4*Y*Z_2*R71
     &      +(X*(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)+3.0d00*X_2*Y*Z_2)*R61
     &      +(6.0d00*X_2*Y+3.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,94)=X_4*Z_3*R71
     &      +(3.0d00*X_2*Z_3+X*(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3))*R61
     &      +(18.0d00*X_2*Z+3.0d00*Z_3)*R51+9.0d00*Z*R41
        R(M,0,95)=X_3*Y_4*R71
     &      +(X*(Y_4+6.0d00*X_2*Y_2)+2.0d00*X*Y_4)*R61
     &      +(X*(3.0d00*X_2+6.0d00*Y_2)+12.0d00*X*Y_2)*R51+9.0d00*X*R41
        R(M,0,96)=X_3*Y_3*Z*R71
     &      +(X*(3.0d00*X_2*Y*Z+Y_3*Z)+2.0d00*X*Y_3*Z)*R61
     &      +9.0d00*X*Y*Z*R51
        R(M,0,97)=X_3*Y_2*Z_2*R71
     &      +(X*(Y_2*Z_2+X_2*(Y_2+Z_2))+2.0d00*X*Y_2*Z_2)*R61
     &      +(X*(X_2+Y_2+Z_2)+2.0d00*X*(Y_2+Z_2))*R51+
     &      3.0d00*X*R41
        R(M,0,98)=X_3*Y*Z_3*R71
     &      +(X*(3.0d00*X_2*Y*Z+Y*Z_3)+2.0d00*X*Y*Z_3)*R61
     &      +9.0d00*X*Y*Z*R51
        R(M,0,99)=X_3*Z_4*R71
     &      +(X*(Z_4+6.0d00*X_2*Z_2)+2.0d00*X*Z_4)*R61
     &      +(X*(3.0d00*X_2+6.0d00*Z_2)+12.0d00*X*Z_2)*R51+9.0d00*X*R41
        R(M,0,100)=X_2*Y_5*R71+(Y_5+10.0d00* X_2*Y_3)*R61
     &      +(10.0d00*Y_3+15.0d00*X_2*Y)*R51+15.0d00*Y*R41
        R(M,0,101)=X_2*Y_4*Z*R71+
     &      (Y_4*Z+6.0d00*X_2*Y_2*Z)*R61
     &      +(3.0d00*X_2*Z+6.0d00*Y_2*Z)*R51+3.0d00*Z*R41
        R(M,0,102)=X_2*Y_3*Z_2*R71
     &      +(X_2*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+Y_3*Z_2)*R61
     &      +(3.0d00*X_2*Y+Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,103)=X_2*Y_2*Z_3*R71
     &      +(Y_2*Z_3+X_2*(3.0d00*Y_2*Z+Z_3))*R61
     &      +(3.0d00*X_2*Z+3.0d00*Y_2*Z+Z_3)*R51+3.0d00*Z*R41
        R(M,0,104)=X_2*Y*Z_4*R71+
     &      (6.0d00*X_2*Y*Z_2+Y*Z_4)*R61
     &      +(3.0d00*X_2*Y+6.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,105)=X_2*Z_5*R71+(Z_5+10.0d00* X_2*Z_3)*R61
     &      +(10.0d00*Z_3+15.0d00*X_2*Z)*R51+15.0d00*Z*R41
        R(M,0,106)=X*Y_6*R71+15.0d00*X*Y_4*R61+
     &      45.0d00*X*Y_2*R51+15.0d00*X*R41
        R(M,0,107)=X*Y_5*Z*R71+10.0d00* X*Y_3*Z*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,108)=X*Y_4*Z_2*R71
     &      +X*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))*R61
     &      +X*(6.0d00*Y_2+3.0d00*Z_2)*R51+3.0d00*X*R41
        R(M,0,109)=X*Y_3*Z_3*R71
     &      +X*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R61+9.0d00*X*Y*Z*R51
        R(M,0,110)=X*Y_2*Z_4*R71+X*(Z_4+6.0d00*Y_2*Z_2)*R61
     &      +X*(3.0d00*Y_2+6.0d00*Z_2)*R51+3.0d00*X*R41
        R(M,0,111)=X*Y*Z_5*R71+10.0d00* X*Y*Z_3*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,112)=X*Z_6*R71+15.0d00*X*Z_4*R61+
     &      45.0d00*X*Z_2*R51
     &      +15.0d00*X*R41
        R(M,0,113)=Y_7*R71+21.0d00*Y_5*R61+105.0d00*Y_3*R51
     &      +105.0d00*Y*R41
        R(M,0,114)=Y_6*Z*R71+15.0d00*Y_4*Z*R61+
     &      45.0d00*Y_2*Z*R51
     &      +15.0d00*Z*R41
        R(M,0,115)=Y_5*Z_2*R71
     &      +(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2)*R61
     &      +(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)*R51
     &      +15.0d00*Y*R41
        R(M,0,116)=Y_4*Z_3*R71
     &      +(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3))*R61
     &      +(18.0d00*Y_2*Z+3.0d00*Z_3)*R51+9.0d00*Z*R41
        R(M,0,117)=Y_3*Z_4*R71
     &      +(Y*(Z_4+6.0d00*Y_2*Z_2)+2.0d00*Y*Z_4)*R61
     &      +(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)*R51+9.0d00*Y*R41
        R(M,0,118)=Y_2*Z_5*R71+(Z_5+10.0d00* Y_2*Z_3)*R61
     &      +(10.0d00*Z_3+15.0d00*Y_2*Z)*R51+15.0d00*Z*R41
        R(M,0,119)=Y*Z_6*R71+15.0d00*Y*Z_4*R61+
     &      45.0d00*Y*Z_2*R51
     &      +15.0d00*Y*R41
        R(M,0,120)=Z_7*R71+21.0d00*Z_5*R61+105.0d00*Z_3*R51
     &      +105.0d00*Z*R41
00100 continue
c
      end
      subroutine hfmkr_8(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
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
      double precision x_5, y_5, z_5
      double precision x_6, y_6, z_6
      double precision x_7, y_7, z_7
      double precision x_8, y_8, z_8
      double precision R11,R21,R31,R41,R51,R61,R71,R81
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
      IJK(5,0,0)=36
      IJK(4,1,0)=37
      IJK(4,0,1)=38
      IJK(3,2,0)=39
      IJK(3,1,1)=40
      IJK(3,0,2)=41
      IJK(2,3,0)=42
      IJK(2,2,1)=43
      IJK(2,1,2)=44
      IJK(2,0,3)=45
      IJK(1,4,0)=46
      IJK(1,3,1)=47
      IJK(1,2,2)=48
      IJK(1,1,3)=49
      IJK(1,0,4)=50
      IJK(0,5,0)=51
      IJK(0,4,1)=52
      IJK(0,3,2)=53
      IJK(0,2,3)=54
      IJK(0,1,4)=55
      IJK(0,0,5)=56
      IJK(6,0,0)=57
      IJK(5,1,0)=58
      IJK(5,0,1)=59
      IJK(4,2,0)=60
      IJK(4,1,1)=61
      IJK(4,0,2)=62
      IJK(3,3,0)=63
      IJK(3,2,1)=64
      IJK(3,1,2)=65
      IJK(3,0,3)=66
      IJK(2,4,0)=67
      IJK(2,3,1)=68
      IJK(2,2,2)=69
      IJK(2,1,3)=70
      IJK(2,0,4)=71
      IJK(1,5,0)=72
      IJK(1,4,1)=73
      IJK(1,3,2)=74
      IJK(1,2,3)=75
      IJK(1,1,4)=76
      IJK(1,0,5)=77
      IJK(0,6,0)=78
      IJK(0,5,1)=79
      IJK(0,4,2)=80
      IJK(0,3,3)=81
      IJK(0,2,4)=82
      IJK(0,1,5)=83
      IJK(0,0,6)=84
      IJK(7,0,0)=85
      IJK(6,1,0)=86
      IJK(6,0,1)=87
      IJK(5,2,0)=88
      IJK(5,1,1)=89
      IJK(5,0,2)=90
      IJK(4,3,0)=91
      IJK(4,2,1)=92
      IJK(4,1,2)=93
      IJK(4,0,3)=94
      IJK(3,4,0)=95
      IJK(3,3,1)=96
      IJK(3,2,2)=97
      IJK(3,1,3)=98
      IJK(3,0,4)=99
      IJK(2,5,0)=100
      IJK(2,4,1)=101
      IJK(2,3,2)=102
      IJK(2,2,3)=103
      IJK(2,1,4)=104
      IJK(2,0,5)=105
      IJK(1,6,0)=106
      IJK(1,5,1)=107
      IJK(1,4,2)=108
      IJK(1,3,3)=109
      IJK(1,2,4)=110
      IJK(1,1,5)=111
      IJK(1,0,6)=112
      IJK(0,7,0)=113
      IJK(0,6,1)=114
      IJK(0,5,2)=115
      IJK(0,4,3)=116
      IJK(0,3,4)=117
      IJK(0,2,5)=118
      IJK(0,1,6)=119
      IJK(0,0,7)=120
      IJK(8,0,0)=121
      IJK(7,1,0)=122
      IJK(7,0,1)=123
      IJK(6,2,0)=124
      IJK(6,1,1)=125
      IJK(6,0,2)=126
      IJK(5,3,0)=127
      IJK(5,2,1)=128
      IJK(5,1,2)=129
      IJK(5,0,3)=130
      IJK(4,4,0)=131
      IJK(4,3,1)=132
      IJK(4,2,2)=133
      IJK(4,1,3)=134
      IJK(4,0,4)=135
      IJK(3,5,0)=136
      IJK(3,4,1)=137
      IJK(3,3,2)=138
      IJK(3,2,3)=139
      IJK(3,1,4)=140
      IJK(3,0,5)=141
      IJK(2,6,0)=142
      IJK(2,5,1)=143
      IJK(2,4,2)=144
      IJK(2,3,3)=145
      IJK(2,2,4)=146
      IJK(2,1,5)=147
      IJK(2,0,6)=148
      IJK(1,7,0)=149
      IJK(1,6,1)=150
      IJK(1,5,2)=151
      IJK(1,4,3)=152
      IJK(1,3,4)=153
      IJK(1,2,5)=154
      IJK(1,1,6)=155
      IJK(1,0,7)=156
      IJK(0,8,0)=157
      IJK(0,7,1)=158
      IJK(0,6,2)=159
      IJK(0,5,3)=160
      IJK(0,4,4)=161
      IJK(0,3,5)=162
      IJK(0,2,6)=163
      IJK(0,1,7)=164
      IJK(0,0,8)=165
c
      do 00100 m=1,NAB
        X=PQ(m,1)
        Y=PQ(m,2)
        Z=PQ(m,3)
        R11=R(m,1,1)
        R(M,0,2)=X*R11
        R(M,0,3)=Y*R11
        R(M,0,4)=Z*R11
        X_2=X*X
        Y_2=Y*Y
        Z_2=Z*Z
        R21=R(m,2,1)
        R(M,0,5)=X_2*R21+R11
        R(M,0,6)=X*Y*R21
        R(M,0,7)=X*Z*R21
        R(M,0,8)=Y_2*R21+R11
        R(M,0,9)=Y*Z*R21
        R(M,0,10)=Z_2*R21+R11
        X_3=X_2*X
        Y_3=Y_2*Y
        Z_3=Z_2*Z
        R31=R(m,3,1)
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
        X_4=X_3*X
        Y_4=Y_3*Y
        Z_4=Z_3*Z
        R41=R(m,4,1)
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
        X_5=X_4*X
        Y_5=Y_4*Y
        Z_5=Z_4*Z
        R51=R(m,5,1)
        R(M,0,36)=X_5*R51+10.0d00* X_3*R41+15.0d00*X*R31
        R(M,0,37)=X_4*Y*R51+6.0d00*X_2*Y*R41+3.0d00*Y*R31
        R(M,0,38)=X_4*Z*R51+6.0d00*X_2*Z*R41+3.0d00*Z*R31
        R(M,0,39)=X_3*Y_2*R51+(X*(X_2+Y_2)+2.0d00*X*Y_2)*R41
     &      +3.0d00*X*R31
        R(M,0,40)=X_3*Y*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,41)=X_3*Z_2*R51+(X*(X_2+Z_2)+2.0d00*X*Z_2)*R41
     &      +3.0d00*X*R31
        R(M,0,42)=X_2*Y_3*R51+(3.0d00*X_2*Y+Y_3)*R41
     &      +3.0d00*Y*R31
        R(M,0,43)=X_2*Y_2*Z*R51+(X_2*Z+Y_2*Z)*R41
     &      +Z*R31
        R(M,0,44)=X_2*Y*Z_2*R51+(X_2*Y+Y*Z_2)*R41
     &      +Y*R31
        R(M,0,45)=X_2*Z_3*R51+(3.0d00*X_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,46)=X*Y_4*R51+6.0d00*X*Y_2*R41+3.0d00*X*R31
        R(M,0,47)=X*Y_3*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,48)=X*Y_2*Z_2*R51+X*(Y_2+Z_2)*R41+
     &      X*R31
        R(M,0,49)=X*Y*Z_3*R51+3.0d00*X*Y*Z*R41
        R(M,0,50)=X*Z_4*R51+6.0d00*X*Z_2*R41+3.0d00*X*R31
        R(M,0,51)=Y_5*R51+10.0d00* Y_3*R41+15.0d00*Y*R31
        R(M,0,52)=Y_4*Z*R51+6.0d00*Y_2*Z*R41+3.0d00*Z*R31
        R(M,0,53)=Y_3*Z_2*R51+(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R41
     &      +3.0d00*Y*R31
        R(M,0,54)=Y_2*Z_3*R51+(3.0d00*Y_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,55)=Y*Z_4*R51+6.0d00*Y*Z_2*R41+3.0d00*Y*R31
        R(M,0,56)=Z_5*R51+10.0d00* Z_3*R41+15.0d00*Z*R31
        X_6=X_5 * X
        Y_6=Y_5 * Y
        Z_6=Z_5 * Z
        R61=R(m,6,1)
        R(M,0,57)=X_6*R61+15.0d00*X_4*R51+45.0d00*X_2*R41
     &      +15.0d00*R31
        R(M,0,58)=X_5*Y*R61+10.0d00* X_3*Y*R51+15.0d00*X*Y*R41
        R(M,0,59)=X_5*Z*R61+10.0d00* X_3*Z*R51+15.0d00*X*Z*R41
        R(M,0,60)=X_4*Y_2*R61
     &      +(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+2.0d00*X*Y_2))*R51
     &      +(6.0d00*X_2+3.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,61)=X_4*Y*Z*R61+6.0d00*X_2*Y*Z*R51+3.0d00*Y*Z*R41
        R(M,0,62)=X_4*Z_2*R61
     &      +(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+2.0d00*X*Z_2))*R51
     &      +(6.0d00*X_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,63)=X_3*Y_3*R61
     &      +(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3)*R51+9.0d00*X*Y*R41
        R(M,0,64)=X_3*Y_2*Z*R61
     &      +(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z)*R51+3.0d00*X*Z*R41
        R(M,0,65)=X_3*Y*Z_2*R61
     &      +(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,66)=X_3*Z_3*R61
     &      +(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3)*R51+9.0d00*X*Z*R41
        R(M,0,67)=X_2*Y_4*R61+(Y_4+6.0d00*X_2*Y_2)*R51
     &      +(3.0d00*X_2+6.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,68)=X_2*Y_3*Z*R61+(3.0d00*X_2*Y*Z+Y_3*Z)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,69)=X_2*Y_2*Z_2*R61
     &      +(Y_2*Z_2+X_2*(Y_2+Z_2))*R51
     &      +(X_2+Y_2+Z_2)*R41+R31
        R(M,0,70)=X_2*Y*Z_3*R61+(3.0d00*X_2*Y*Z+Y*Z_3)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,71)=X_2*Z_4*R61+(Z_4+6.0d00*X_2*Z_2)*R51
     &      +(3.0d00*X_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,72)=X*Y_5*R61+10.0d00* X*Y_3*R51+15.0d00*X*Y*R41
        R(M,0,73)=X*Y_4*Z*R61+6.0d00*X*Y_2*Z*R51+3.0d00*X*Z*R41
        R(M,0,74)=X*Y_3*Z_2*R61
     &      +X*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,75)=X*Y_2*Z_3*R61+X*(3.0d00*Y_2*Z+Z_3)*R51
     &      +3.0d00*X*Z*R41
        R(M,0,76)=X*Y*Z_4*R61+6.0d00*X*Y*Z_2*R51+3.0d00*X*Y*R41
        R(M,0,77)=X*Z_5*R61+10.0d00* X*Z_3*R51+15.0d00*X*Z*R41
        R(M,0,78)=Y_6*R61+15.0d00*Y_4*R51+45.0d00*Y_2*R41
     &      +15.0d00*R31
        R(M,0,79)=Y_5*Z*R61+10.0d00* Y_3*Z*R51+15.0d00*Y*Z*R41
        R(M,0,80)=Y_4*Z_2*R61
     &      +(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))*R51
     &      +(6.0d00*Y_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,81)=Y_3*Z_3*R61
     &      +(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R51+9.0d00*Y*Z*R41
        R(M,0,82)=Y_2*Z_4*R61+(Z_4+6.0d00*Y_2*Z_2)*R51
     &      +(3.0d00*Y_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,83)=Y*Z_5*R61+10.0d00* Y*Z_3*R51+15.0d00*Y*Z*R41
        R(M,0,84)=Z_6*R61+15.0d00*Z_4*R51+45.0d00*Z_2*R41
     &      +15.0d00*R31
        X_7=X_6 * X
        Y_7=Y_6 * Y
        Z_7=Z_6 * Z
        R71=R(m,7,1)
        R(M,0,85)=X_7*R71+21.0d00*X_5*R61+105.0d00*X_3*R51
     &      +105.0d00*X*R41
        R(M,0,86)=X_6*Y*R71+15.0d00*X_4*Y*R61+45.0d00*X_2*Y*R51
     &      +15.0d00*Y*R41
        R(M,0,87)=X_6*Z*R71+15.0d00*X_4*Z*R61+45.0d00*X_2*Z*R51
     &      +15.0d00*Z*R41
        R(M,0,88)=X_5*Y_2*R71
     &      +(X*(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+
     &      2.0d00*X*Y_2))+4.0d00*X_3*Y_2)*R61
     &      +(X*(6.0d00*X_2+3.0d00*Y_2)+
     &      4.0d00*X*(X_2+Y_2)+8.0d00*X*Y_2)*R51
     &      +15.0d00*X*R41
        R(M,0,89)=X_5*Y*Z*R71+10.0d00* X_3*Y*Z*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,90)=X_5*Z_2*R71
     &      +(X*(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+
     &      2.0d00*X*Z_2))+4.0d00*X_3*Z_2)*R61
     &      +(X*(6.0d00*X_2+3.0d00*Z_2)+
     &      4.0d00*X*(X_2+Z_2)+8.0d00*X*Z_2)*R51
     &      +15.0d00*X*R41
        R(M,0,91)=X_4*Y_3*R71
     &      +(3.0d00*X_2*Y_3+X*(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3))*R61
     &      +(18.0d00*X_2*Y+3.0d00*Y_3)*R51+9.0d00*Y*R41
        R(M,0,92)=X_4*Y_2*Z*R71
     &      +(3.0d00*X_2*Y_2*Z+X*(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z))*R61
     &      +(6.0d00*X_2*Z+3.0d00*Y_2*Z)*R51+3.0d00*Z*R41
        R(M,0,93)=X_4*Y*Z_2*R71
     &      +(X*(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)+3.0d00*X_2*Y*Z_2)*R61
     &      +(6.0d00*X_2*Y+3.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,94)=X_4*Z_3*R71
     &      +(3.0d00*X_2*Z_3+X*(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3))*R61
     &      +(18.0d00*X_2*Z+3.0d00*Z_3)*R51+9.0d00*Z*R41
        R(M,0,95)=X_3*Y_4*R71
     &      +(X*(Y_4+6.0d00*X_2*Y_2)+2.0d00*X*Y_4)*R61
     &      +(X*(3.0d00*X_2+6.0d00*Y_2)+12.0d00*X*Y_2)*R51+9.0d00*X*R41
        R(M,0,96)=X_3*Y_3*Z*R71
     &      +(X*(3.0d00*X_2*Y*Z+Y_3*Z)+2.0d00*X*Y_3*Z)*R61
     &      +9.0d00*X*Y*Z*R51
        R(M,0,97)=X_3*Y_2*Z_2*R71
     &      +(X*(Y_2*Z_2+X_2*(Y_2+Z_2))+2.0d00*X*Y_2*Z_2)*R61
     &      +(X*(X_2+Y_2+Z_2)+2.0d00*X*(Y_2+Z_2))*R51+
     &      3.0d00*X*R41
        R(M,0,98)=X_3*Y*Z_3*R71
     &      +(X*(3.0d00*X_2*Y*Z+Y*Z_3)+2.0d00*X*Y*Z_3)*R61
     &      +9.0d00*X*Y*Z*R51
        R(M,0,99)=X_3*Z_4*R71
     &      +(X*(Z_4+6.0d00*X_2*Z_2)+2.0d00*X*Z_4)*R61
     &      +(X*(3.0d00*X_2+6.0d00*Z_2)+12.0d00*X*Z_2)*R51+9.0d00*X*R41
        R(M,0,100)=X_2*Y_5*R71+(Y_5+10.0d00* X_2*Y_3)*R61
     &      +(10.0d00*Y_3+15.0d00*X_2*Y)*R51+15.0d00*Y*R41
        R(M,0,101)=X_2*Y_4*Z*R71+
     &      (Y_4*Z+6.0d00*X_2*Y_2*Z)*R61
     &      +(3.0d00*X_2*Z+6.0d00*Y_2*Z)*R51+3.0d00*Z*R41
        R(M,0,102)=X_2*Y_3*Z_2*R71
     &      +(X_2*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+Y_3*Z_2)*R61
     &      +(3.0d00*X_2*Y+Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,103)=X_2*Y_2*Z_3*R71
     &      +(Y_2*Z_3+X_2*(3.0d00*Y_2*Z+Z_3))*R61
     &      +(3.0d00*X_2*Z+3.0d00*Y_2*Z+Z_3)*R51+3.0d00*Z*R41
        R(M,0,104)=X_2*Y*Z_4*R71+
     &      (6.0d00*X_2*Y*Z_2+Y*Z_4)*R61
     &      +(3.0d00*X_2*Y+6.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,105)=X_2*Z_5*R71+(Z_5+10.0d00* X_2*Z_3)*R61
     &      +(10.0d00*Z_3+15.0d00*X_2*Z)*R51+15.0d00*Z*R41
        R(M,0,106)=X*Y_6*R71+15.0d00*X*Y_4*R61+
     &      45.0d00*X*Y_2*R51+15.0d00*X*R41
        R(M,0,107)=X*Y_5*Z*R71+10.0d00* X*Y_3*Z*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,108)=X*Y_4*Z_2*R71
     &      +X*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))*R61
     &      +X*(6.0d00*Y_2+3.0d00*Z_2)*R51+3.0d00*X*R41
        R(M,0,109)=X*Y_3*Z_3*R71
     &      +X*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R61+9.0d00*X*Y*Z*R51
        R(M,0,110)=X*Y_2*Z_4*R71+X*(Z_4+6.0d00*Y_2*Z_2)*R61
     &      +X*(3.0d00*Y_2+6.0d00*Z_2)*R51+3.0d00*X*R41
        R(M,0,111)=X*Y*Z_5*R71+10.0d00* X*Y*Z_3*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,112)=X*Z_6*R71+15.0d00*X*Z_4*R61+
     &      45.0d00*X*Z_2*R51
     &      +15.0d00*X*R41
        R(M,0,113)=Y_7*R71+21.0d00*Y_5*R61+105.0d00*Y_3*R51
     &      +105.0d00*Y*R41
        R(M,0,114)=Y_6*Z*R71+15.0d00*Y_4*Z*R61+
     &      45.0d00*Y_2*Z*R51
     &      +15.0d00*Z*R41
        R(M,0,115)=Y_5*Z_2*R71
     &      +(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2)*R61
     &      +(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)*R51
     &      +15.0d00*Y*R41
        R(M,0,116)=Y_4*Z_3*R71
     &      +(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3))*R61
     &      +(18.0d00*Y_2*Z+3.0d00*Z_3)*R51+9.0d00*Z*R41
        R(M,0,117)=Y_3*Z_4*R71
     &      +(Y*(Z_4+6.0d00*Y_2*Z_2)+2.0d00*Y*Z_4)*R61
     &      +(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)*R51+9.0d00*Y*R41
        R(M,0,118)=Y_2*Z_5*R71+(Z_5+10.0d00* Y_2*Z_3)*R61
     &      +(10.0d00*Z_3+15.0d00*Y_2*Z)*R51+15.0d00*Z*R41
        R(M,0,119)=Y*Z_6*R71+15.0d00*Y*Z_4*R61+
     &      45.0d00*Y*Z_2*R51
     &      +15.0d00*Y*R41
        R(M,0,120)=Z_7*R71+21.0d00*Z_5*R61+105.0d00*Z_3*R51
     &      +105.0d00*Z*R41
        X_8=X_7 * X
        Y_8=Y_7 * Y
        Z_8=Z_7 * Z
        R81=R(m,8,1)
        R(M,0,121)=X_8*R81+28.0d00*X_6*R71+210.0d00*X_4*R61
     &      +420.0d00*X_2*R51+105.0d00*R41
        R(M,0,122)=X_7*Y*R81+21.0d00*X_5*Y*R71
     &      +105.0d00*X_3*Y*R61+105.0d00*X*Y*R51
        R(M,0,123)=X_7*Z*R81+21.0d00*X_5*Z*R71
     &      +105.0d00*X_3*Z*R61+105.0d00*X*Z*R51
        R(M,0,124)=X_6*Y_2*R81
     &      +(5.0d00*X_4*Y_2+X*(X*(3.0d00*X_2*Y_2+
     &      X*(X*(X_2+Y_2)+2.0d00*X*Y_2))+4.0d00*X_3*Y_2))*R71
     &      +(X*(X*(6.0d00*X_2+3.0d00*Y_2)+4.0d00*X*(X_2+Y_2)+
     &      8.0d00*X*Y_2)+
     &      5.0d00*X*(X*(X_2+Y_2)+2.0d00*X*Y_2)+15.0d00*X_2*Y_2)*R61
     &      +(45.0d00*X_2+15.0d00*Y_2)*R51+15.0d00*R41
        R(M,0,125)=X_6*Y*Z*R81+15.0d00*X_4*Y*Z*R71
     &      +45.0d00*X_2*Y*Z*R61+15.0d00*Y*Z*R51
        R(M,0,126)=X_6*Z_2*R81
     &      +(5.0d00*X_4*Z_2+X*(X*(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+
     &      2.0d00*X*Z_2))+4.0d00*X_3*Z_2))*R71
     &      +(X*(X*(6.0d00*X_2+3.0d00*Z_2)+4.0d00*X*(X_2+Z_2)+
     &      8.0d00*X*Z_2)+
     &      5.0d00*X*(X*(X_2+Z_2)+2.0d00*X*Z_2)+15.0d00*X_2*Z_2)*R61
     &      +(45.0d00*X_2+15.0d00*Z_2)*R51+15.0d00*R41
        R(M,0,127)=X_5*Y_3*R81
     &      +(X*(3.0d00*X_2*Y_3+X*(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3))+
     &      4.0d00*X_3*Y_3)*R71
     &      +(X*(18.0d00*X_2*Y+3.0d00*Y_3)+4.0d00*X*(3.0d00*X_2*Y+Y_3)+
     &      8.0d00*X*Y_3)*R61+45.0d00*X*Y*R51
        R(M,0,128)=X_5*Y_2*Z*R81
     &      +(X*(3.0d00*X_2*Y_2*Z+X*(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z))+
     &      4.0d00*X_3*Y_2*Z)*R71
     &      +(X*(6.0d00*X_2*Z+3.0d00*Y_2*Z)+4.0d00*X*(X_2*Z+Y_2*Z)+
     &      8.0d00*X*Y_2*Z)*R61+15.0d00*X*Z*R51
        R(M,0,129)=X_5*Y*Z_2*R81
     &      +(X*(X*(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)+3.0d00*X_2*Y*Z_2)+
     &      4.0d00*X_3*Y*Z_2)*R71
     &      +(X*(6.0d00*X_2*Y+3.0d00*Y*Z_2)+4.0d00*X*(X_2*Y+Y*Z_2)+
     &      8.0d00*X*Y*Z_2)*R61+15.0d00*X*Y*R51
        R(M,0,130)=X_5*Z_3*R81
     &      +(X*(3.0d00*X_2*Z_3+X*(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3))+
     &      4.0d00*X_3*Z_3)*R71
     &      +(X*(18.0d00*X_2*Z+3.0d00*Z_3)+4.0d00*X*(3.0d00*X_2*Z+Z_3)+
     &      8.0d00*X*Z_3)*R61+45.0d00*X*Z*R51
        R(M,0,131)=X_4*Y_4*R81
     &      +(3.0d00*X_2*Y_4+X*(X*(Y_4+6.0d00*X_2*Y_2)+
     &      2.0d00*X*Y_4))*R71
     &      +(X*(X*(3.0d00*X_2+6.0d00*Y_2)+12.0d00*X*Y_2)+
     &      18.0d00*X_2*Y_2+3.0d00*Y_4)*
     &      R61+(18.0d00*X_2+18.0d00*Y_2)*R51+9.0d00*R41
        R(M,0,132)=X_4*Y_3*Z*R81
     &      +(3.0d00*X_2*Y_3*Z+X*(X*(3.0d00*X_2*Y*Z+Y_3*Z)+
     &      2.0d00*X*Y_3*Z))*
     &      R71+(18.0d00*X_2*Y*Z+3.0d00*Y_3*Z)*R61+9.0d00*Y*Z*R51
        R(M,0,133)=X_4*Y_2*Z_2*R81
     &      +(3.0d00*X_2*Y_2*Z_2+X*(X*(Y_2*Z_2+X_2*(Y_2+Z_2))+
     &      2.0d00*X*Y_2*Z_2))*R71
     &      +(X*(X*(X_2+Y_2+Z_2)+2.0d00*X*(Y_2+Z_2))+
     &      3.0d00*X_2*(Y_2+Z_2)+3.0d00*Y_2*Z_2)*R61
     &      +(6.0d00*X_2+3.0d00*Y_2+3.0d00*Z_2)*R51+3.0d00*R41
        R(M,0,134)=X_4*Y*Z_3*R81
     &      +(X*(X*(3.0d00*X_2*Y*Z+Y*Z_3)+2.0d00*X*Y*Z_3)+
     &      3.0d00*X_2*Y*Z_3)*
     &      R71+(18.0d00*X_2*Y*Z+3.0d00*Y*Z_3)*R61+9.0d00*Y*Z*R51
        R(M,0,135)=X_4*Z_4*R81
     &      +(3.0d00*X_2*Z_4+X*(X*(Z_4+6.0d00*X_2*Z_2)+
     &      2.0d00*X*Z_4))*R71
     &      +(X*(X*(3.0d00*X_2+6.0d00*Z_2)+12.0d00*X*Z_2)+
     &      18.0d00*X_2*Z_2+3.0d00*Z_4)*
     &      R61+(18.0d00*X_2+18.0d00*Z_2)*R51+9.0d00*R41
        R(M,0,136)=X_3*Y_5*R81
     &      +(X*(Y_5+10.0d00* X_2*Y_3)+2.0d00*X*Y_5)*R71
     &      +(X*(10.0d00*Y_3+15.0d00*X_2*Y)+20.0d00*X*Y_3)*R61+
     &      45.0d00*X*Y*R51
        R(M,0,137)=X_3*Y_4*Z*R81
     &      +(X*(Y_4*Z+6.0d00*X_2*Y_2*Z)+2.0d00*X*Y_4*Z)*R71
     &      +(X*(3.0d00*X_2*Z+6.0d00*Y_2*Z)+12.0d00*X*Y_2*Z)*R61+
     &      9.0d00*X*Z*R51
        R(M,0,138)=X_3*Y_3*Z_2*R81
     &      +(X*(X_2*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+Y_3*Z_2)+
     &      2.0d00*X*Y_3*Z_2)*R71
     &      +(X*(3.0d00*X_2*Y+Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+
     &      2.0d00*X*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))*R61+9.0d00*X*Y*R51
        R(M,0,139)=X_3*Y_2*Z_3*R81
     &      +(X*(Y_2*Z_3+X_2*(3.0d00*Y_2*Z+Z_3))+2.0d00*X*Y_2*Z_3)*
     &      R71+(X*(3.0d00*X_2*Z+3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*X*(3.0d00*Y_2*Z+Z_3))*
     &      R61+9.0d00*X*Z*R51
        R(M,0,140)=X_3*Y*Z_4*R81
     &      +(X*(6.0d00*X_2*Y*Z_2+Y*Z_4)+2.0d00*X*Y*Z_4)*R71
     &      +(X*(3.0d00*X_2*Y+6.0d00*Y*Z_2)+12.0d00*X*Y*Z_2)*R61+
     &      9.0d00*X*Y*R51
        R(M,0,141)=X_3*Z_5*R81
     &      +(X*(Z_5+10.0d00* X_2*Z_3)+2.0d00*X*Z_5)*R71
     &      +(X*(10.0d00*Z_3+15.0d00*X_2*Z)+20.0d00*X*Z_3)*R61+
     &      45.0d00*X*Z*R51
        R(M,0,142)=X_2*Y_6*R81+(Y_6+15.0d00*X_2*Y_4)*R71
     &      +(45.0d00*X_2*Y_2+15.0d00*Y_4)*R61+
     &      (15.0d00*X_2+45.0d00*Y_2)*R51+15.0d00*R41
        R(M,0,143)=X_2*Y_5*Z*R81+(Y_5*Z+10.0d00* X_2*Y_3*Z)*
     &      R71+(10.0d00*Y_3*Z+15.0d00*X_2*Y*Z)*R61+15.0d00*Y*Z*R51
        R(M,0,144)=X_2*Y_4*Z_2*R81
     &      +(Y_4*Z_2+X_2*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2)))*R71+
     &      (X_2*(6.0d00*Y_2+3.0d00*Z_2)+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+
     &      3.0d00*Y_2*Z_2)*R61+(3.0d00*X_2+6.0d00*Y_2+3.0d00*Z_2)*R51
     &      +3.0d00*R41
        R(M,0,145)=X_2*Y_3*Z_3*R81
     &      +(X_2*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)+Y_3*Z_3)*R71
     &      +(9.0d00*X_2*Y*Z+Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R61
     &      +9.0d00*Y*Z*R51
        R(M,0,146)=X_2*Y_2*Z_4*R81
     &      +(Y_2*Z_4+X_2*(Z_4+6.0d00*Y_2*Z_2))*R71
     &      +(X_2*(3.0d00*Y_2+6.0d00*Z_2)+6.0d00*Y_2*Z_2+Z_4)*R61
     &      +(3.0d00*X_2+3.0d00*Y_2+6.0d00*Z_2)*R51+3.0d00*R41
        R(M,0,147)=X_2*Y*Z_5*R81+(10.0d00*X_2*Y*Z_3+Y*Z_5)*
     &      R71+(15.0d00*X_2*Y*Z+10.0d00* Y*Z_3)*R61+15.0d00*Y*Z*R51
        R(M,0,148)=X_2*Z_6*R81+(Z_6+15.0d00*X_2*Z_4)*R71
     &      +(45.0d00*X_2*Z_2+15.0d00*Z_4)*R61+(15.0d00*X_2+
     &      45.0d00*Z_2)*
     &      R51+15.0d00*R41
        R(M,0,149)=X*Y_7*R81+21.0d00*X*Y_5*R71
     &      +105.0d00*X*Y_3*R61+105.0d00*X*Y*R51
        R(M,0,150)=X*Y_6*Z*R81+15.0d00*X*Y_4*Z*R71
     &      +45.0d00*X*Y_2*Z*R61+15.0d00*X*Z*R51
        R(M,0,151)=X*Y_5*Z_2*R81
     &      +X*(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))+
     &      4.0d00*Y_3*Z_2)*R71
     &      +X*(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)*R61
     &      +15.0d00*X*Y*R51
        R(M,0,152)=X*Y_4*Z_3*R81
     &      +X*(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*Y*Z_3))*R71
     &      +X*(18.0d00*Y_2*Z+3.0d00*Z_3)*R61+9.0d00*X*Z*R51
        R(M,0,153)=X*Y_3*Z_4*R81
     &      +X*(Y*(Z_4+6.0d00*Y_2*Z_2)+2.0d00*Y*Z_4)*R71
     &      +X*(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)*R61+
     &      9.0d00*X*Y*R51
        R(M,0,154)=X*Y_2*Z_5*R81+X*(Z_5+10.0d00* Y_2*Z_3)*R71
     &      +X*(10.0d00*Z_3+15.0d00*Y_2*Z)*R61+15.0d00*X*Z*R51
        R(M,0,155)=X*Y*Z_6*R81+15.0d00*X*Y*Z_4*R71
     &      +45.0d00*X*Y*Z_2*R61+15.0d00*X*Y*R51
        R(M,0,156)=X*Z_7*R81+21.0d00*X*Z_5*R71
     &      +105.0d00*X*Z_3*R61+105.0d00*X*Z*R51
        R(M,0,157)=Y_8*R81+28.0d00*Y_6*R71+210.0d00*Y_4*R61
     &      +420.0d00*Y_2*R51+105.0d00*R41
        R(M,0,158)=Y_7*Z*R81+21.0d00*Y_5*Z*R71
     &      +105.0d00*Y_3*Z*R61+105.0d00*Y*Z*R51
        R(M,0,159)=Y_6*Z_2*R81
     &      +(5.0d00*Y_4*Z_2+Y*(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2))*R71
     &      +(Y*(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)+
     &      5.0d00*Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+15.0d00*Y_2*Z_2)*R61
     &      +(45.0d00*Y_2+15.0d00*Z_2)*R51+15.0d00*R41
        R(M,0,160)=Y_5*Z_3*R81
     &      +(Y*(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3))+
     &      4.0d00*Y_3*Z_3)*R71
     &      +(Y*(18.0d00*Y_2*Z+3.0d00*Z_3)+4.0d00*Y*(3.0d00*Y_2*Z+Z_3)+
     &      8.0d00*Y*Z_3)*
     &      R61+45.0d00*Y*Z*R51
        R(M,0,161)=Y_4*Z_4*R81
     &      +(3.0d00*Y_2*Z_4+Y*(Y*(Z_4+6.0d00*Y_2*Z_2)+
     &      2.0d00*Y*Z_4))*R71
     &      +(Y*(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)+
     &      18.0d00*Y_2*Z_2+3.0d00*Z_4)*
     &      R61+(18.0d00*Y_2+18.0d00*Z_2)*R51+9.0d00*R41
        R(M,0,162)=Y_3*Z_5*R81
     &      +(Y*(Z_5+10.0d00* Y_2*Z_3)+2.0d00*Y*Z_5)*R71
     &      +(Y*(10.0d00*Z_3+15.0d00*Y_2*Z)+20.0d00*Y*Z_3)*R61+
     &      45.0d00*Y*Z*R51
        R(M,0,163)=Y_2*Z_6*R81+(Z_6+15.0d00*Y_2*Z_4)*R71
     &      +(45.0d00*Y_2*Z_2+15.0d00*Z_4)*R61+
     &      (15.0d00*Y_2+45.0d00*Z_2)*R51
     &      +15.0d00*R41
        R(M,0,164)=Y*Z_7*R81+21.0d00*Y*Z_5*R71
     &      +105.0d00*Y*Z_3*R61+105.0d00*Y*Z*R51
        R(M,0,165)=Z_8*R81+28.0d00*Z_6*R71+210.0d00*Z_4*R61
     &      +420.0d00*Z_2*R51+105.0d00*R41
00100 continue
c
      end
      subroutine hfmkr_9(R,PQ,IJK,NAB,Lr,Lr3)
      implicit none
c $Id$
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
      double precision x_5, y_5, z_5
      double precision x_6, y_6, z_6
      double precision x_7, y_7, z_7
      double precision x_8, y_8, z_8
      double precision R11,R21,R31,R41,R51,R61,R71,R81,R91
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
      IJK(5,0,0)=36
      IJK(4,1,0)=37
      IJK(4,0,1)=38
      IJK(3,2,0)=39
      IJK(3,1,1)=40
      IJK(3,0,2)=41
      IJK(2,3,0)=42
      IJK(2,2,1)=43
      IJK(2,1,2)=44
      IJK(2,0,3)=45
      IJK(1,4,0)=46
      IJK(1,3,1)=47
      IJK(1,2,2)=48
      IJK(1,1,3)=49
      IJK(1,0,4)=50
      IJK(0,5,0)=51
      IJK(0,4,1)=52
      IJK(0,3,2)=53
      IJK(0,2,3)=54
      IJK(0,1,4)=55
      IJK(0,0,5)=56
      IJK(6,0,0)=57
      IJK(5,1,0)=58
      IJK(5,0,1)=59
      IJK(4,2,0)=60
      IJK(4,1,1)=61
      IJK(4,0,2)=62
      IJK(3,3,0)=63
      IJK(3,2,1)=64
      IJK(3,1,2)=65
      IJK(3,0,3)=66
      IJK(2,4,0)=67
      IJK(2,3,1)=68
      IJK(2,2,2)=69
      IJK(2,1,3)=70
      IJK(2,0,4)=71
      IJK(1,5,0)=72
      IJK(1,4,1)=73
      IJK(1,3,2)=74
      IJK(1,2,3)=75
      IJK(1,1,4)=76
      IJK(1,0,5)=77
      IJK(0,6,0)=78
      IJK(0,5,1)=79
      IJK(0,4,2)=80
      IJK(0,3,3)=81
      IJK(0,2,4)=82
      IJK(0,1,5)=83
      IJK(0,0,6)=84
      IJK(7,0,0)=85
      IJK(6,1,0)=86
      IJK(6,0,1)=87
      IJK(5,2,0)=88
      IJK(5,1,1)=89
      IJK(5,0,2)=90
      IJK(4,3,0)=91
      IJK(4,2,1)=92
      IJK(4,1,2)=93
      IJK(4,0,3)=94
      IJK(3,4,0)=95
      IJK(3,3,1)=96
      IJK(3,2,2)=97
      IJK(3,1,3)=98
      IJK(3,0,4)=99
      IJK(2,5,0)=100
      IJK(2,4,1)=101
      IJK(2,3,2)=102
      IJK(2,2,3)=103
      IJK(2,1,4)=104
      IJK(2,0,5)=105
      IJK(1,6,0)=106
      IJK(1,5,1)=107
      IJK(1,4,2)=108
      IJK(1,3,3)=109
      IJK(1,2,4)=110
      IJK(1,1,5)=111
      IJK(1,0,6)=112
      IJK(0,7,0)=113
      IJK(0,6,1)=114
      IJK(0,5,2)=115
      IJK(0,4,3)=116
      IJK(0,3,4)=117
      IJK(0,2,5)=118
      IJK(0,1,6)=119
      IJK(0,0,7)=120
      IJK(8,0,0)=121
      IJK(7,1,0)=122
      IJK(7,0,1)=123
      IJK(6,2,0)=124
      IJK(6,1,1)=125
      IJK(6,0,2)=126
      IJK(5,3,0)=127
      IJK(5,2,1)=128
      IJK(5,1,2)=129
      IJK(5,0,3)=130
      IJK(4,4,0)=131
      IJK(4,3,1)=132
      IJK(4,2,2)=133
      IJK(4,1,3)=134
      IJK(4,0,4)=135
      IJK(3,5,0)=136
      IJK(3,4,1)=137
      IJK(3,3,2)=138
      IJK(3,2,3)=139
      IJK(3,1,4)=140
      IJK(3,0,5)=141
      IJK(2,6,0)=142
      IJK(2,5,1)=143
      IJK(2,4,2)=144
      IJK(2,3,3)=145
      IJK(2,2,4)=146
      IJK(2,1,5)=147
      IJK(2,0,6)=148
      IJK(1,7,0)=149
      IJK(1,6,1)=150
      IJK(1,5,2)=151
      IJK(1,4,3)=152
      IJK(1,3,4)=153
      IJK(1,2,5)=154
      IJK(1,1,6)=155
      IJK(1,0,7)=156
      IJK(0,8,0)=157
      IJK(0,7,1)=158
      IJK(0,6,2)=159
      IJK(0,5,3)=160
      IJK(0,4,4)=161
      IJK(0,3,5)=162
      IJK(0,2,6)=163
      IJK(0,1,7)=164
      IJK(0,0,8)=165
      IJK(9,0,0)=166
      IJK(8,1,0)=167
      IJK(8,0,1)=168
      IJK(7,2,0)=169
      IJK(7,1,1)=170
      IJK(7,0,2)=171
      IJK(6,3,0)=172
      IJK(6,2,1)=173
      IJK(6,1,2)=174
      IJK(6,0,3)=175
      IJK(5,4,0)=176
      IJK(5,3,1)=177
      IJK(5,2,2)=178
      IJK(5,1,3)=179
      IJK(5,0,4)=180
      IJK(4,5,0)=181
      IJK(4,4,1)=182
      IJK(4,3,2)=183
      IJK(4,2,3)=184
      IJK(4,1,4)=185
      IJK(4,0,5)=186
      IJK(3,6,0)=187
      IJK(3,5,1)=188
      IJK(3,4,2)=189
      IJK(3,3,3)=190
      IJK(3,2,4)=191
      IJK(3,1,5)=192
      IJK(3,0,6)=193
      IJK(2,7,0)=194
      IJK(2,6,1)=195
      IJK(2,5,2)=196
      IJK(2,4,3)=197
      IJK(2,3,4)=198
      IJK(2,2,5)=199
      IJK(2,1,6)=200
      IJK(2,0,7)=201
      IJK(1,8,0)=202
      IJK(1,7,1)=203
      IJK(1,6,2)=204
      IJK(1,5,3)=205
      IJK(1,4,4)=206
      IJK(1,3,5)=207
      IJK(1,2,6)=208
      IJK(1,1,7)=209
      IJK(1,0,8)=210
      IJK(0,9,0)=211
      IJK(0,8,1)=212
      IJK(0,7,2)=213
      IJK(0,6,3)=214
      IJK(0,5,4)=215
      IJK(0,4,5)=216
      IJK(0,3,6)=217
      IJK(0,2,7)=218
      IJK(0,1,8)=219
      IJK(0,0,9)=220
c
      do 00100 m=1,NAB
        X=PQ(m,1)
        Y=PQ(m,2)
        Z=PQ(m,3)
        R11=R(m,1,1)
        R(M,0,2)=X*R11
        R(M,0,3)=Y*R11
        R(M,0,4)=Z*R11
        X_2=X*X
        Y_2=Y*Y
        Z_2=Z*Z
        R21=R(m,2,1)
        R(M,0,5)=X_2*R21+R11
        R(M,0,6)=X*Y*R21
        R(M,0,7)=X*Z*R21
        R(M,0,8)=Y_2*R21+R11
        R(M,0,9)=Y*Z*R21
        R(M,0,10)=Z_2*R21+R11
        X_3=X_2*X
        Y_3=Y_2*Y
        Z_3=Z_2*Z
        R31=R(m,3,1)
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
        X_4=X_3*X
        Y_4=Y_3*Y
        Z_4=Z_3*Z
        R41=R(m,4,1)
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
        X_5=X_4*X
        Y_5=Y_4*Y
        Z_5=Z_4*Z
        R51=R(m,5,1)
        R(M,0,36)=X_5*R51+10.0d00* X_3*R41+15.0d00*X*R31
        R(M,0,37)=X_4*Y*R51+6.0d00*X_2*Y*R41+3.0d00*Y*R31
        R(M,0,38)=X_4*Z*R51+6.0d00*X_2*Z*R41+3.0d00*Z*R31
        R(M,0,39)=X_3*Y_2*R51+(X*(X_2+Y_2)+2.0d00*X*Y_2)*R41
     &      +3.0d00*X*R31
        R(M,0,40)=X_3*Y*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,41)=X_3*Z_2*R51+(X*(X_2+Z_2)+2.0d00*X*Z_2)*R41
     &      +3.0d00*X*R31
        R(M,0,42)=X_2*Y_3*R51+(3.0d00*X_2*Y+Y_3)*R41
     &      +3.0d00*Y*R31
        R(M,0,43)=X_2*Y_2*Z*R51+(X_2*Z+Y_2*Z)*R41
     &      +Z*R31
        R(M,0,44)=X_2*Y*Z_2*R51+(X_2*Y+Y*Z_2)*R41
     &      +Y*R31
        R(M,0,45)=X_2*Z_3*R51+(3.0d00*X_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,46)=X*Y_4*R51+6.0d00*X*Y_2*R41+3.0d00*X*R31
        R(M,0,47)=X*Y_3*Z*R51+3.0d00*X*Y*Z*R41
        R(M,0,48)=X*Y_2*Z_2*R51+X*(Y_2+Z_2)*R41+
     &      X*R31
        R(M,0,49)=X*Y*Z_3*R51+3.0d00*X*Y*Z*R41
        R(M,0,50)=X*Z_4*R51+6.0d00*X*Z_2*R41+3.0d00*X*R31
        R(M,0,51)=Y_5*R51+10.0d00* Y_3*R41+15.0d00*Y*R31
        R(M,0,52)=Y_4*Z*R51+6.0d00*Y_2*Z*R41+3.0d00*Z*R31
        R(M,0,53)=Y_3*Z_2*R51+(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R41
     &      +3.0d00*Y*R31
        R(M,0,54)=Y_2*Z_3*R51+(3.0d00*Y_2*Z+Z_3)*R41
     &      +3.0d00*Z*R31
        R(M,0,55)=Y*Z_4*R51+6.0d00*Y*Z_2*R41+3.0d00*Y*R31
        R(M,0,56)=Z_5*R51+10.0d00* Z_3*R41+15.0d00*Z*R31
        X_6=X_5 * X
        Y_6=Y_5 * Y
        Z_6=Z_5 * Z
        R61=R(m,6,1)
        R(M,0,57)=X_6*R61+15.0d00*X_4*R51+45.0d00*X_2*R41
     &      +15.0d00*R31
        R(M,0,58)=X_5*Y*R61+10.0d00* X_3*Y*R51+15.0d00*X*Y*R41
        R(M,0,59)=X_5*Z*R61+10.0d00* X_3*Z*R51+15.0d00*X*Z*R41
        R(M,0,60)=X_4*Y_2*R61
     &      +(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+2.0d00*X*Y_2))*R51
     &      +(6.0d00*X_2+3.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,61)=X_4*Y*Z*R61+6.0d00*X_2*Y*Z*R51+3.0d00*Y*Z*R41
        R(M,0,62)=X_4*Z_2*R61
     &      +(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+2.0d00*X*Z_2))*R51
     &      +(6.0d00*X_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,63)=X_3*Y_3*R61
     &      +(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3)*R51+9.0d00*X*Y*R41
        R(M,0,64)=X_3*Y_2*Z*R61
     &      +(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z)*R51+3.0d00*X*Z*R41
        R(M,0,65)=X_3*Y*Z_2*R61
     &      +(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,66)=X_3*Z_3*R61
     &      +(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3)*R51+9.0d00*X*Z*R41
        R(M,0,67)=X_2*Y_4*R61+(Y_4+6.0d00*X_2*Y_2)*R51
     &      +(3.0d00*X_2+6.0d00*Y_2)*R41+3.0d00*R31
        R(M,0,68)=X_2*Y_3*Z*R61+(3.0d00*X_2*Y*Z+Y_3*Z)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,69)=X_2*Y_2*Z_2*R61
     &      +(Y_2*Z_2+X_2*(Y_2+Z_2))*R51
     &      +(X_2+Y_2+Z_2)*R41+R31
        R(M,0,70)=X_2*Y*Z_3*R61+(3.0d00*X_2*Y*Z+Y*Z_3)*R51
     &      +3.0d00*Y*Z*R41
        R(M,0,71)=X_2*Z_4*R61+(Z_4+6.0d00*X_2*Z_2)*R51
     &      +(3.0d00*X_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,72)=X*Y_5*R61+10.0d00* X*Y_3*R51+15.0d00*X*Y*R41
        R(M,0,73)=X*Y_4*Z*R61+6.0d00*X*Y_2*Z*R51+3.0d00*X*Z*R41
        R(M,0,74)=X*Y_3*Z_2*R61
     &      +X*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R51+3.0d00*X*Y*R41
        R(M,0,75)=X*Y_2*Z_3*R61+X*(3.0d00*Y_2*Z+Z_3)*R51
     &      +3.0d00*X*Z*R41
        R(M,0,76)=X*Y*Z_4*R61+6.0d00*X*Y*Z_2*R51+3.0d00*X*Y*R41
        R(M,0,77)=X*Z_5*R61+10.0d00* X*Z_3*R51+15.0d00*X*Z*R41
        R(M,0,78)=Y_6*R61+15.0d00*Y_4*R51+45.0d00*Y_2*R41
     &      +15.0d00*R31
        R(M,0,79)=Y_5*Z*R61+10.0d00* Y_3*Z*R51+15.0d00*Y*Z*R41
        R(M,0,80)=Y_4*Z_2*R61
     &      +(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))*R51
     &      +(6.0d00*Y_2+3.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,81)=Y_3*Z_3*R61
     &      +(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R51+9.0d00*Y*Z*R41
        R(M,0,82)=Y_2*Z_4*R61+(Z_4+6.0d00*Y_2*Z_2)*R51
     &      +(3.0d00*Y_2+6.0d00*Z_2)*R41+3.0d00*R31
        R(M,0,83)=Y*Z_5*R61+10.0d00* Y*Z_3*R51+15.0d00*Y*Z*R41
        R(M,0,84)=Z_6*R61+15.0d00*Z_4*R51+45.0d00*Z_2*R41
     &      +15.0d00*R31
        X_7=X_6 * X
        Y_7=Y_6 * Y
        Z_7=Z_6 * Z
        R71=R(m,7,1)
        R(M,0,85)=X_7*R71+21.0d00*X_5*R61+105.0d00*X_3*R51
     &      +105.0d00*X*R41
        R(M,0,86)=X_6*Y*R71+15.0d00*X_4*Y*R61+45.0d00*X_2*Y*R51
     &      +15.0d00*Y*R41
        R(M,0,87)=X_6*Z*R71+15.0d00*X_4*Z*R61+45.0d00*X_2*Z*R51
     &      +15.0d00*Z*R41
        R(M,0,88)=X_5*Y_2*R71
     &      +(X*(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+
     &      2.0d00*X*Y_2))+4.0d00*X_3*Y_2)*R61
     &      +(X*(6.0d00*X_2+3.0d00*Y_2)+
     &      4.0d00*X*(X_2+Y_2)+8.0d00*X*Y_2)*R51
     &      +15.0d00*X*R41
        R(M,0,89)=X_5*Y*Z*R71+10.0d00* X_3*Y*Z*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,90)=X_5*Z_2*R71
     &      +(X*(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+
     &      2.0d00*X*Z_2))+4.0d00*X_3*Z_2)*R61
     &      +(X*(6.0d00*X_2+3.0d00*Z_2)+
     &      4.0d00*X*(X_2+Z_2)+8.0d00*X*Z_2)*R51
     &      +15.0d00*X*R41
        R(M,0,91)=X_4*Y_3*R71
     &      +(3.0d00*X_2*Y_3+X*(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3))*R61
     &      +(18.0d00*X_2*Y+3.0d00*Y_3)*R51+9.0d00*Y*R41
        R(M,0,92)=X_4*Y_2*Z*R71
     &      +(3.0d00*X_2*Y_2*Z+X*(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z))*R61
     &      +(6.0d00*X_2*Z+3.0d00*Y_2*Z)*R51+3.0d00*Z*R41
        R(M,0,93)=X_4*Y*Z_2*R71
     &      +(X*(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)+3.0d00*X_2*Y*Z_2)*R61
     &      +(6.0d00*X_2*Y+3.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,94)=X_4*Z_3*R71
     &      +(3.0d00*X_2*Z_3+X*(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3))*R61
     &      +(18.0d00*X_2*Z+3.0d00*Z_3)*R51+9.0d00*Z*R41
        R(M,0,95)=X_3*Y_4*R71
     &      +(X*(Y_4+6.0d00*X_2*Y_2)+2.0d00*X*Y_4)*R61
     &      +(X*(3.0d00*X_2+6.0d00*Y_2)+12.0d00*X*Y_2)*R51+9.0d00*X*R41
        R(M,0,96)=X_3*Y_3*Z*R71
     &      +(X*(3.0d00*X_2*Y*Z+Y_3*Z)+2.0d00*X*Y_3*Z)*R61
     &      +9.0d00*X*Y*Z*R51
        R(M,0,97)=X_3*Y_2*Z_2*R71
     &      +(X*(Y_2*Z_2+X_2*(Y_2+Z_2))+2.0d00*X*Y_2*Z_2)*R61
     &      +(X*(X_2+Y_2+Z_2)+2.0d00*X*(Y_2+Z_2))*R51+
     &      3.0d00*X*R41
        R(M,0,98)=X_3*Y*Z_3*R71
     &      +(X*(3.0d00*X_2*Y*Z+Y*Z_3)+2.0d00*X*Y*Z_3)*R61
     &      +9.0d00*X*Y*Z*R51
        R(M,0,99)=X_3*Z_4*R71
     &      +(X*(Z_4+6.0d00*X_2*Z_2)+2.0d00*X*Z_4)*R61
     &      +(X*(3.0d00*X_2+6.0d00*Z_2)+12.0d00*X*Z_2)*R51+9.0d00*X*R41
        R(M,0,100)=X_2*Y_5*R71+(Y_5+10.0d00* X_2*Y_3)*R61
     &      +(10.0d00*Y_3+15.0d00*X_2*Y)*R51+15.0d00*Y*R41
        R(M,0,101)=X_2*Y_4*Z*R71+
     &      (Y_4*Z+6.0d00*X_2*Y_2*Z)*R61
     &      +(3.0d00*X_2*Z+6.0d00*Y_2*Z)*R51+3.0d00*Z*R41
        R(M,0,102)=X_2*Y_3*Z_2*R71
     &      +(X_2*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+Y_3*Z_2)*R61
     &      +(3.0d00*X_2*Y+Y*(Y_2+Z_2)+2.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,103)=X_2*Y_2*Z_3*R71
     &      +(Y_2*Z_3+X_2*(3.0d00*Y_2*Z+Z_3))*R61
     &      +(3.0d00*X_2*Z+3.0d00*Y_2*Z+Z_3)*R51+3.0d00*Z*R41
        R(M,0,104)=X_2*Y*Z_4*R71+
     &      (6.0d00*X_2*Y*Z_2+Y*Z_4)*R61
     &      +(3.0d00*X_2*Y+6.0d00*Y*Z_2)*R51+3.0d00*Y*R41
        R(M,0,105)=X_2*Z_5*R71+(Z_5+10.0d00* X_2*Z_3)*R61
     &      +(10.0d00*Z_3+15.0d00*X_2*Z)*R51+15.0d00*Z*R41
        R(M,0,106)=X*Y_6*R71+15.0d00*X*Y_4*R61+
     &      45.0d00*X*Y_2*R51+15.0d00*X*R41
        R(M,0,107)=X*Y_5*Z*R71+10.0d00* X*Y_3*Z*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,108)=X*Y_4*Z_2*R71
     &      +X*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))*R61
     &      +X*(6.0d00*Y_2+3.0d00*Z_2)*R51+3.0d00*X*R41
        R(M,0,109)=X*Y_3*Z_3*R71
     &      +X*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R61+9.0d00*X*Y*Z*R51
        R(M,0,110)=X*Y_2*Z_4*R71+X*(Z_4+6.0d00*Y_2*Z_2)*R61
     &      +X*(3.0d00*Y_2+6.0d00*Z_2)*R51+3.0d00*X*R41
        R(M,0,111)=X*Y*Z_5*R71+10.0d00* X*Y*Z_3*R61
     &      +15.0d00*X*Y*Z*R51
        R(M,0,112)=X*Z_6*R71+15.0d00*X*Z_4*R61+
     &      45.0d00*X*Z_2*R51
     &      +15.0d00*X*R41
        R(M,0,113)=Y_7*R71+21.0d00*Y_5*R61+105.0d00*Y_3*R51
     &      +105.0d00*Y*R41
        R(M,0,114)=Y_6*Z*R71+15.0d00*Y_4*Z*R61+
     &      45.0d00*Y_2*Z*R51
     &      +15.0d00*Z*R41
        R(M,0,115)=Y_5*Z_2*R71
     &      +(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2)*R61
     &      +(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)*R51
     &      +15.0d00*Y*R41
        R(M,0,116)=Y_4*Z_3*R71
     &      +(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3))*R61
     &      +(18.0d00*Y_2*Z+3.0d00*Z_3)*R51+9.0d00*Z*R41
        R(M,0,117)=Y_3*Z_4*R71
     &      +(Y*(Z_4+6.0d00*Y_2*Z_2)+2.0d00*Y*Z_4)*R61
     &      +(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)*R51+9.0d00*Y*R41
        R(M,0,118)=Y_2*Z_5*R71+(Z_5+10.0d00* Y_2*Z_3)*R61
     &      +(10.0d00*Z_3+15.0d00*Y_2*Z)*R51+15.0d00*Z*R41
        R(M,0,119)=Y*Z_6*R71+15.0d00*Y*Z_4*R61+
     &      45.0d00*Y*Z_2*R51
     &      +15.0d00*Y*R41
        R(M,0,120)=Z_7*R71+21.0d00*Z_5*R61+105.0d00*Z_3*R51
     &      +105.0d00*Z*R41
        X_8=X_7 * X
        Y_8=Y_7 * Y
        Z_8=Z_7 * Z
        R81=R(m,8,1)
        R(M,0,121)=X_8*R81+28.0d00*X_6*R71+210.0d00*X_4*R61
     &      +420.0d00*X_2*R51+105.0d00*R41
        R(M,0,122)=X_7*Y*R81+21.0d00*X_5*Y*R71
     &      +105.0d00*X_3*Y*R61+105.0d00*X*Y*R51
        R(M,0,123)=X_7*Z*R81+21.0d00*X_5*Z*R71
     &      +105.0d00*X_3*Z*R61+105.0d00*X*Z*R51
        R(M,0,124)=X_6*Y_2*R81
     &      +(5.0d00*X_4*Y_2+X*(X*(3.0d00*X_2*Y_2+
     &      X*(X*(X_2+Y_2)+2.0d00*X*Y_2))+4.0d00*X_3*Y_2))*R71
     &      +(X*(X*(6.0d00*X_2+3.0d00*Y_2)+4.0d00*X*(X_2+Y_2)+
     &      8.0d00*X*Y_2)+
     &      5.0d00*X*(X*(X_2+Y_2)+2.0d00*X*Y_2)+15.0d00*X_2*Y_2)*R61
     &      +(45.0d00*X_2+15.0d00*Y_2)*R51+15.0d00*R41
        R(M,0,125)=X_6*Y*Z*R81+15.0d00*X_4*Y*Z*R71
     &      +45.0d00*X_2*Y*Z*R61+15.0d00*Y*Z*R51
        R(M,0,126)=X_6*Z_2*R81
     &      +(5.0d00*X_4*Z_2+X*(X*(3.0d00*X_2*Z_2+X*(X*(X_2+Z_2)+
     &      2.0d00*X*Z_2))+4.0d00*X_3*Z_2))*R71
     &      +(X*(X*(6.0d00*X_2+3.0d00*Z_2)+4.0d00*X*(X_2+Z_2)+
     &      8.0d00*X*Z_2)+
     &      5.0d00*X*(X*(X_2+Z_2)+2.0d00*X*Z_2)+15.0d00*X_2*Z_2)*R61
     &      +(45.0d00*X_2+15.0d00*Z_2)*R51+15.0d00*R41
        R(M,0,127)=X_5*Y_3*R81
     &      +(X*(3.0d00*X_2*Y_3+X*(X*(3.0d00*X_2*Y+Y_3)+2.0d00*X*Y_3))+
     &      4.0d00*X_3*Y_3)*R71
     &      +(X*(18.0d00*X_2*Y+3.0d00*Y_3)+4.0d00*X*(3.0d00*X_2*Y+Y_3)+
     &      8.0d00*X*Y_3)*R61+45.0d00*X*Y*R51
        R(M,0,128)=X_5*Y_2*Z*R81
     &      +(X*(3.0d00*X_2*Y_2*Z+X*(X*(X_2*Z+Y_2*Z)+2.0d00*X*Y_2*Z))+
     &      4.0d00*X_3*Y_2*Z)*R71
     &      +(X*(6.0d00*X_2*Z+3.0d00*Y_2*Z)+4.0d00*X*(X_2*Z+Y_2*Z)+
     &      8.0d00*X*Y_2*Z)*R61+15.0d00*X*Z*R51
        R(M,0,129)=X_5*Y*Z_2*R81
     &      +(X*(X*(X*(X_2*Y+Y*Z_2)+2.0d00*X*Y*Z_2)+3.0d00*X_2*Y*Z_2)+
     &      4.0d00*X_3*Y*Z_2)*R71
     &      +(X*(6.0d00*X_2*Y+3.0d00*Y*Z_2)+4.0d00*X*(X_2*Y+Y*Z_2)+
     &      8.0d00*X*Y*Z_2)*R61+15.0d00*X*Y*R51
        R(M,0,130)=X_5*Z_3*R81
     &      +(X*(3.0d00*X_2*Z_3+X*(X*(3.0d00*X_2*Z+Z_3)+2.0d00*X*Z_3))+
     &      4.0d00*X_3*Z_3)*R71
     &      +(X*(18.0d00*X_2*Z+3.0d00*Z_3)+4.0d00*X*(3.0d00*X_2*Z+Z_3)+
     &      8.0d00*X*Z_3)*R61+45.0d00*X*Z*R51
        R(M,0,131)=X_4*Y_4*R81
     &      +(3.0d00*X_2*Y_4+X*(X*(Y_4+6.0d00*X_2*Y_2)+
     &      2.0d00*X*Y_4))*R71
     &      +(X*(X*(3.0d00*X_2+6.0d00*Y_2)+12.0d00*X*Y_2)+
     &      18.0d00*X_2*Y_2+3.0d00*Y_4)*
     &      R61+(18.0d00*X_2+18.0d00*Y_2)*R51+9.0d00*R41
        R(M,0,132)=X_4*Y_3*Z*R81
     &      +(3.0d00*X_2*Y_3*Z+X*(X*(3.0d00*X_2*Y*Z+Y_3*Z)+
     &      2.0d00*X*Y_3*Z))*
     &      R71+(18.0d00*X_2*Y*Z+3.0d00*Y_3*Z)*R61+9.0d00*Y*Z*R51
        R(M,0,133)=X_4*Y_2*Z_2*R81
     &      +(3.0d00*X_2*Y_2*Z_2+X*(X*(Y_2*Z_2+X_2*(Y_2+Z_2))+
     &      2.0d00*X*Y_2*Z_2))*R71
     &      +(X*(X*(X_2+Y_2+Z_2)+2.0d00*X*(Y_2+Z_2))+
     &      3.0d00*X_2*(Y_2+Z_2)+3.0d00*Y_2*Z_2)*R61
     &      +(6.0d00*X_2+3.0d00*Y_2+3.0d00*Z_2)*R51+3.0d00*R41
        R(M,0,134)=X_4*Y*Z_3*R81
     &      +(X*(X*(3.0d00*X_2*Y*Z+Y*Z_3)+2.0d00*X*Y*Z_3)+
     &      3.0d00*X_2*Y*Z_3)*
     &      R71+(18.0d00*X_2*Y*Z+3.0d00*Y*Z_3)*R61+9.0d00*Y*Z*R51
        R(M,0,135)=X_4*Z_4*R81
     &      +(3.0d00*X_2*Z_4+X*(X*(Z_4+6.0d00*X_2*Z_2)+
     &      2.0d00*X*Z_4))*R71
     &      +(X*(X*(3.0d00*X_2+6.0d00*Z_2)+12.0d00*X*Z_2)+
     &      18.0d00*X_2*Z_2+3.0d00*Z_4)*
     &      R61+(18.0d00*X_2+18.0d00*Z_2)*R51+9.0d00*R41
        R(M,0,136)=X_3*Y_5*R81
     &      +(X*(Y_5+10.0d00* X_2*Y_3)+2.0d00*X*Y_5)*R71
     &      +(X*(10.0d00*Y_3+15.0d00*X_2*Y)+20.0d00*X*Y_3)*R61+
     &      45.0d00*X*Y*R51
        R(M,0,137)=X_3*Y_4*Z*R81
     &      +(X*(Y_4*Z+6.0d00*X_2*Y_2*Z)+2.0d00*X*Y_4*Z)*R71
     &      +(X*(3.0d00*X_2*Z+6.0d00*Y_2*Z)+12.0d00*X*Y_2*Z)*R61+
     &      9.0d00*X*Z*R51
        R(M,0,138)=X_3*Y_3*Z_2*R81
     &      +(X*(X_2*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+Y_3*Z_2)+
     &      2.0d00*X*Y_3*Z_2)*R71
     &      +(X*(3.0d00*X_2*Y+Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+
     &      2.0d00*X*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))*R61+9.0d00*X*Y*R51
        R(M,0,139)=X_3*Y_2*Z_3*R81
     &      +(X*(Y_2*Z_3+X_2*(3.0d00*Y_2*Z+Z_3))+2.0d00*X*Y_2*Z_3)*
     &      R71+(X*(3.0d00*X_2*Z+3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*X*(3.0d00*Y_2*Z+Z_3))*
     &      R61+9.0d00*X*Z*R51
        R(M,0,140)=X_3*Y*Z_4*R81
     &      +(X*(6.0d00*X_2*Y*Z_2+Y*Z_4)+2.0d00*X*Y*Z_4)*R71
     &      +(X*(3.0d00*X_2*Y+6.0d00*Y*Z_2)+12.0d00*X*Y*Z_2)*R61+
     &      9.0d00*X*Y*R51
        R(M,0,141)=X_3*Z_5*R81
     &      +(X*(Z_5+10.0d00* X_2*Z_3)+2.0d00*X*Z_5)*R71
     &      +(X*(10.0d00*Z_3+15.0d00*X_2*Z)+20.0d00*X*Z_3)*R61+
     &      45.0d00*X*Z*R51
        R(M,0,142)=X_2*Y_6*R81+(Y_6+15.0d00*X_2*Y_4)*R71
     &      +(45.0d00*X_2*Y_2+15.0d00*Y_4)*R61+
     &      (15.0d00*X_2+45.0d00*Y_2)*R51+15.0d00*R41
        R(M,0,143)=X_2*Y_5*Z*R81+(Y_5*Z+10.0d00* X_2*Y_3*Z)*
     &      R71+(10.0d00*Y_3*Z+15.0d00*X_2*Y*Z)*R61+15.0d00*Y*Z*R51
        R(M,0,144)=X_2*Y_4*Z_2*R81
     &      +(Y_4*Z_2+X_2*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2)))*R71+
     &      (X_2*(6.0d00*Y_2+3.0d00*Z_2)+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+
     &      3.0d00*Y_2*Z_2)*R61+(3.0d00*X_2+6.0d00*Y_2+3.0d00*Z_2)*R51
     &      +3.0d00*R41
        R(M,0,145)=X_2*Y_3*Z_3*R81
     &      +(X_2*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)+Y_3*Z_3)*R71
     &      +(9.0d00*X_2*Y*Z+Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)*R61
     &      +9.0d00*Y*Z*R51
        R(M,0,146)=X_2*Y_2*Z_4*R81
     &      +(Y_2*Z_4+X_2*(Z_4+6.0d00*Y_2*Z_2))*R71
     &      +(X_2*(3.0d00*Y_2+6.0d00*Z_2)+6.0d00*Y_2*Z_2+Z_4)*R61
     &      +(3.0d00*X_2+3.0d00*Y_2+6.0d00*Z_2)*R51+3.0d00*R41
        R(M,0,147)=X_2*Y*Z_5*R81+(10.0d00*X_2*Y*Z_3+Y*Z_5)*
     &      R71+(15.0d00*X_2*Y*Z+10.0d00* Y*Z_3)*R61+15.0d00*Y*Z*R51
        R(M,0,148)=X_2*Z_6*R81+(Z_6+15.0d00*X_2*Z_4)*R71
     &      +(45.0d00*X_2*Z_2+15.0d00*Z_4)*R61+(15.0d00*X_2+
     &      45.0d00*Z_2)*
     &      R51+15.0d00*R41
        R(M,0,149)=X*Y_7*R81+21.0d00*X*Y_5*R71
     &      +105.0d00*X*Y_3*R61+105.0d00*X*Y*R51
        R(M,0,150)=X*Y_6*Z*R81+15.0d00*X*Y_4*Z*R71
     &      +45.0d00*X*Y_2*Z*R61+15.0d00*X*Z*R51
        R(M,0,151)=X*Y_5*Z_2*R81
     &      +X*(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))+
     &      4.0d00*Y_3*Z_2)*R71
     &      +X*(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)*R61
     &      +15.0d00*X*Y*R51
        R(M,0,152)=X*Y_4*Z_3*R81
     &      +X*(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*Y*Z_3))*R71
     &      +X*(18.0d00*Y_2*Z+3.0d00*Z_3)*R61+9.0d00*X*Z*R51
        R(M,0,153)=X*Y_3*Z_4*R81
     &      +X*(Y*(Z_4+6.0d00*Y_2*Z_2)+2.0d00*Y*Z_4)*R71
     &      +X*(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)*R61+
     &      9.0d00*X*Y*R51
        R(M,0,154)=X*Y_2*Z_5*R81+X*(Z_5+10.0d00* Y_2*Z_3)*R71
     &      +X*(10.0d00*Z_3+15.0d00*Y_2*Z)*R61+15.0d00*X*Z*R51
        R(M,0,155)=X*Y*Z_6*R81+15.0d00*X*Y*Z_4*R71
     &      +45.0d00*X*Y*Z_2*R61+15.0d00*X*Y*R51
        R(M,0,156)=X*Z_7*R81+21.0d00*X*Z_5*R71
     &      +105.0d00*X*Z_3*R61+105.0d00*X*Z*R51
        R(M,0,157)=Y_8*R81+28.0d00*Y_6*R71+210.0d00*Y_4*R61
     &      +420.0d00*Y_2*R51+105.0d00*R41
        R(M,0,158)=Y_7*Z*R81+21.0d00*Y_5*Z*R71
     &      +105.0d00*Y_3*Z*R61+105.0d00*Y*Z*R51
        R(M,0,159)=Y_6*Z_2*R81
     &      +(5.0d00*Y_4*Z_2+Y*(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2))*R71
     &      +(Y*(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)+
     &      5.0d00*Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+15.0d00*Y_2*Z_2)*R61
     &      +(45.0d00*Y_2+15.0d00*Z_2)*R51+15.0d00*R41
        R(M,0,160)=Y_5*Z_3*R81
     &      +(Y*(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3))+
     &      4.0d00*Y_3*Z_3)*R71
     &      +(Y*(18.0d00*Y_2*Z+3.0d00*Z_3)+4.0d00*Y*(3.0d00*Y_2*Z+Z_3)+
     &      8.0d00*Y*Z_3)*
     &      R61+45.0d00*Y*Z*R51
        R(M,0,161)=Y_4*Z_4*R81
     &      +(3.0d00*Y_2*Z_4+Y*(Y*(Z_4+6.0d00*Y_2*Z_2)+
     &      2.0d00*Y*Z_4))*R71
     &      +(Y*(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)+
     &      18.0d00*Y_2*Z_2+3.0d00*Z_4)*
     &      R61+(18.0d00*Y_2+18.0d00*Z_2)*R51+9.0d00*R41
        R(M,0,162)=Y_3*Z_5*R81
     &      +(Y*(Z_5+10.0d00* Y_2*Z_3)+2.0d00*Y*Z_5)*R71
     &      +(Y*(10.0d00*Z_3+15.0d00*Y_2*Z)+20.0d00*Y*Z_3)*R61+
     &      45.0d00*Y*Z*R51
        R(M,0,163)=Y_2*Z_6*R81+(Z_6+15.0d00*Y_2*Z_4)*R71
     &      +(45.0d00*Y_2*Z_2+15.0d00*Z_4)*R61+
     &      (15.0d00*Y_2+45.0d00*Z_2)*R51
     &      +15.0d00*R41
        R(M,0,164)=Y*Z_7*R81+21.0d00*Y*Z_5*R71
     &      +105.0d00*Y*Z_3*R61+105.0d00*Y*Z*R51
        R(M,0,165)=Z_8*R81+28.0d00*Z_6*R71+210.0d00*Z_4*R61
     &      +420.0d00*Z_2*R51+105.0d00*R41
        R91=R(m,9,1)
        R(M,0,166)=X_5*X_4*R91+36.0d00*X_7*R81+378.0d00*X_5*R71
     &      +1260.0d00*X_3*R61+945.0d00*X*R51
        R(M,0,167)=X_8*Y*R91+28.0d00*X_6*Y*R81
     &      +210.0d00*X_4*Y*R71+420.0d00*X_2*Y*R61+105.0d00*Y*R51
        R(M,0,168)=X_8*Z*R91+28.0d00*X_6*Z*R81
     &      +210.0d00*X_4*Z*R71+420.0d00*X_2*Z*R61+105.0d00*Z*R51
        R(M,0,169)=X_7*Y_2*R91
     &      +(6.0d00*X_5*Y_2+X*(5.0d00*X_4*Y_2+X*(X*(3.0d00*
     &      X_2*Y_2+X*(X*(X_2+
     &      Y_2)+2.0d00*X*Y_2))+4.0d00*X_3*Y_2)))*R81
     &      +(6.0d00*X*(3.0d00*X_2*Y_2+X*(X*(X_2+Y_2)+2.0d00*X*Y_2))+
     &      X*(X*(X*(6.0d00*X_2+3.0d00*Y_2)+4.0d00*X*(X_2+Y_2)+
     &      8.0d00*X*Y_2)+
     &      5.0d00*X*(X*(X_2+Y_2)+2.0d00*X*Y_2)+15.0d00*X_2*Y_2)+
     &      24.0d00*X_3*Y_2)*
     &      R71+(X*(45.0d00*X_2+15.0d00*Y_2)+24.0d00*X*(X_2+Y_2)+
     &      48.0d00*X*Y_2+
     &      6.0d00*X*(6.0d00*X_2+3.0d00*Y_2))*R61+105.0d00*X*R51
        R(M,0,170)=X_7*Y*Z*R91+21.0d00*X_5*Y*Z*R81
     &      +105.0d00*X_3*Y*Z*R71+105.0d00*X*Y*Z*R61
        R(M,0,171)=X_7*Z_2*R91
     &      +(6.0d00*X_5*Z_2+X*(5.0d00*X_4*Z_2+X*(X*(3.0d00*X_2*Z_2+
     &      X*(X*(X_2+
     &      Z_2)+2.0d00*X*Z_2))+4.0d00*X_3*Z_2)))*R81+(6.0d00*X*
     &      (3.0d00*X_2*Z_2
     &      +X*(X*(X_2+Z_2)+2.0d00*X*Z_2))+X*(X*(X*(6.0d00*X_2+
     &      3.0d00*Z_2)+4.0d00*X*
     &      (X_2+Z_2)+8.0d00*X*Z_2)+5.0d00*X*(X*(X_2+Z_2)+
     &      2.0d00*X*Z_2)+
     &      15.0d00*X_2*Z_2)+24.0d00*X_3*Z_2)*R71+(X*(45.0d00*X_2+
     &      15.0d00*Z_2)+
     &      24.0d00*X*(X_2+Z_2)+48.0d00*X*Z_2+6.0d00*X*(6.0d00*X_2+
     &      3.0d00*Z_2))*R61
     &      +105.0d00*X*R51
        R(M,0,172)=X_6*Y_3*R91
     &      +(5.0d00*X_4*Y_3+X*(X*(3.0d00*X_2*Y_3+
     &      X*(X*(3.0d00*X_2*Y+Y_3)+
     &      2.0d00*X*Y_3))+4.0d00*X_3*Y_3))*R81+(X*(X*(18.0d00*X_2*Y+
     &      3.0d00*Y_3)+4.0d00*X*(3.0d00*X_2*Y+Y_3)+8.0d00*X*Y_3)+
     &      5.0d00*X*(X*(3.0d00*X_2*Y+
     &      Y_3)+2.0d00*X*Y_3)+15.0d00*X_2*Y_3)*R71+(135.0d00*X_2*Y+
     &      15.0d00*Y_3)*R61+45.0d00*Y*R51
        R(M,0,173)=X_6*Y_2*Z*R91
     &      +(5.0d00*X_4*Y_2*Z+X*(X*(3.0d00*X_2*Y_2*Z+
     &      X*(X*(X_2*Z+Y_2*Z)+
     &      2.0d00*X*Y_2*Z))+4.0d00*X_3*Y_2*Z))*R81+(X*(X*
     &      (6.0d00*X_2*Z+
     &      3.0d00*Y_2*Z)+4.0d00*X*(X_2*Z+Y_2*Z)+
     &      8.0d00*X*Y_2*Z)+5.0d00*X*(X*(X_2*Z+
     &      Y_2*Z)+2.0d00*X*Y_2*Z)+15.0d00*X_2*Y_2*Z)*R71
     &      +(45.0d00*X_2*Z+15.0d00*Y_2*Z)*R61+15.0d00*Z*R51
        R(M,0,174)=X_6*Y*Z_2*R91
     &      +(5.0d00*X_4*Y*Z_2+X*(X*(X*(X*(X_2*Y+Y*Z_2)+
     &      2.0d00*X*Y*Z_2)+
     &      3.0d00*X_2*Y*Z_2)+4.0d00*X_3*Y*Z_2))*R81+
     &      (5.0d00*X*(X*(X_2*Y+
     &      Y*Z_2)+2.0d00*X*Y*Z_2)+X*(X*(6.0d00*X_2*Y+
     &      3.0d00*Y*Z_2)+4.0d00*X*(X_2*Y+
     &      Y*Z_2)+8.0d00*X*Y*Z_2)+15.0d00*X_2*Y*Z_2)*R71+
     &      (45.0d00*X_2*Y+
     &      15.0d00*Y*Z_2)*R61+15.0d00*Y*R51
        R(M,0,175)=X_6*Z_3*R91
     &      +(5.0d00*X_4*Z_3+X*(X*(3.0d00*X_2*Z_3+X*(X*
     &      (3.0d00*X_2*Z+Z_3)+
     &      2.0d00*X*Z_3))+4.0d00*X_3*Z_3))*R81+(X*(X*(18.0d00*X_2*Z+
     &      3.0d00*Z_3)+4.0d00*X*(3.0d00*X_2*Z+Z_3)+8.0d00*X*Z_3)+
     &      5.0d00*X*(X*(3.0d00*X_2*Z+
     &      Z_3)+2.0d00*X*Z_3)+15.0d00*X_2*Z_3)*R71+(135.0d00*X_2*Z+
     &      15.0d00*Z_3)*R61+45.0d00*Z*R51
        R(M,0,176)=X_5*Y_4*R91
     &      +(X*(3.0d00*X_2*Y_4+X*(X*(Y_4+6.0d00*X_2*Y_2)+
     &      2.0d00*X*Y_4))+
     &      4.0d00*X_3*Y_4)*R81+(X*(X*(X*(3.0d00*X_2+
     &      6.0d00*Y_2)+12.0d00*X*Y_2)+
     &      18.0d00*X_2*Y_2+3.0d00*Y_4)+4.0d00*X*(Y_4+
     &      6.0d00*X_2*Y_2)+8.0d00*X*Y_4)*
     &      R71+(X*(18.0d00*X_2+18.0d00*Y_2)+
     &      4.0d00*X*(3.0d00*X_2+6.0d00*Y_2)+
     &      48.0d00*X*Y_2)*R61+45.0d00*X*R51
        R(M,0,177)=X_5*Y_3*Z*R91
     &      +(X*(3.0d00*X_2*Y_3*Z+X*(X*(3.0d00*X_2*Y*Z+Y_3*Z)+
     &      2.0d00*X*Y_3*Z))+
     &      4.0d00*X_3*Y_3*Z)*R81+(X*(18.0d00*X_2*Y*Z+3.0d00*Y_3*Z)+
     &      4.0d00*X*(3.0d00*X_2*Y*Z+Y_3*Z)+8.0d00*X*Y_3*Z)*R71+
     &      45.0d00*X*Y*Z*R61
        R(M,0,178)=X_5*Y_2*Z_2*R91
     &      +(X*(3.0d00*X_2*Y_2*Z_2+X*(X*(Y_2*Z_2+X_2*(Y_2+Z_2))+
     &      2.0d00*X*Y_2*Z_2))+4.0d00*X_3*Y_2*Z_2)*R81+(X*(X*(X*(X_2+
     &      Y_2+Z_2)+2.0d00*X*(Y_2+Z_2))+3.0d00*X_2*(Y_2+Z_2)+
     &      3.0d00*Y_2*
     &      Z_2)+4.0d00*X*(Y_2*Z_2+X_2*(Y_2+Z_2))+8.0d00*X*Y_2*Z_2)*
     &      R71+(X*(6.0d00*X_2+3.0d00*Y_2+3.0d00*Z_2)+
     &      4.0d00*X*(X_2+Y_2+Z_2)+
     &      8.0d00*X*(Y_2+Z_2))*R61+15.0d00*X*R51
        R(M,0,179)=X_5*Y*Z_3*R91
     &      +(X*(X*(X*(3.0d00*X_2*Y*Z+Y*Z_3)+2.0d00*X*Y*Z_3)+
     &      3.0d00*X_2*Y*Z_3)+
     &      4.0d00*X_3*Y*Z_3)*R81+(X*(18.0d00*X_2*Y*Z+
     &      3.0d00*Y*Z_3)+4.0d00*X*
     &      (3.0d00*X_2*Y*Z+Y*Z_3)+8.0d00*X*Y*Z_3)*R71+
     &      45.0d00*X*Y*Z*R61
        R(M,0,180)=X_5*Z_4*R91
     &      +(X*(3.0d00*X_2*Z_4+X*(X*(Z_4+6.0d00*X_2*Z_2)+
     &      2.0d00*X*Z_4))+
     &      4.0d00*X_3*Z_4)*R81+(X*(X*(X*(3.0d00*X_2+6.0d00*Z_2)+
     &      12.0d00*X*Z_2)+
     &      18.0d00*X_2*Z_2+3.0d00*Z_4)+4.0d00*X*(Z_4+6.0d00*X_2*Z_2)+
     &      8.0d00*X*Z_4)*
     &      R71+(X*(18.0d00*X_2+18.0d00*Z_2)+4.0d00*X*(3.0d00*X_2+
     &      6.0d00*Z_2)+48.0d00*X*Z_2)*R61+45.0d00*X*R51
        R(M,0,181)=X_4*Y_5*R91
     &      +(3.0d00*X_2*Y_5+X*(X*(Y_5+10.0d00* X_2*Y_3)+
     &      2.0d00*X*Y_5))*R81
     &      +(X*(X*(10.0d00*Y_3+15.0d00*X_2*Y)+20.0d00*X*Y_3)+
     &      30.0d00*X_2*Y_3+
     &      3.0d00*Y_5)*R71+(90.0d00*X_2*Y+30.0d00*Y_3)*R61+
     &      45.0d00*Y*R51
        R(M,0,182)=X_4*Y_4*Z*R91+(3.0d00*X_2*Y_4*Z+X*(X*(Y_4*Z+
     &      6.0d00*X_2*Y_2*Z)+2.0d00*X*Y_4*Z))*R81+(X*(X*(
     &      3.0d00*X_2*Z+
     &      6.0d00*Y_2*Z)+12.0d00*X*Y_2*Z)+18.0d00*X_2*Y_2*Z+
     &      3.0d00*Y_4*Z)*R71
     &      +(18.0d00*X_2*Z+18.0d00*Y_2*Z)*R61+9.0d00*Z*R51
        R(M,0,183)=X_4*Y_3*Z_2*R91
     &      +(3.0d00*X_2*Y_3*Z_2+X*(X*(X_2*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+
     &      Y_3*Z_2)+2.0d00*X*Y_3*Z_2))*R81+(3.0d00*X_2*(Y*(Y_2+
     &      Z_2)+2.0d00*Y*Z_2)+X*(X*(3.0d00*X_2*Y+Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2)+
     &      2.0d00*X*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))+3.0d00*Y_3*Z_2)*R71+
     &      (18.0d00*X_2*Y+3.0d00*Y*(Y_2+Z_2)+6.0d00*Y*Z_2)*R61+
     &      9.0d00*Y*R51
        R(M,0,184)=X_4*Y_2*Z_3*R91
     &      +(3.0d00*X_2*Y_2*Z_3+X*(X*(Y_2*Z_3+X_2*(3.0d00*Y_2*Z+Z_3))+
     &      2.0d00*X*Y_2*Z_3))*R81+(X*(X*(3.0d00*X_2*Z+
     &      3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*X*(3.0d00*Y_2*Z+Z_3))+
     &      3.0d00*X_2*(3.0d00*Y_2*Z+Z_3)+3.0d00*Y_2*Z_3)*
     &      R71+(9.0d00*Y_2*Z+3.0d00*Z_3+18.0d00*X_2*Z)*R61+
     &      9.0d00*Z*R51
        R(M,0,185)=X_4*Y*Z_4*R91
     &      +(X*(X*(6.0d00*X_2*Y*Z_2+Y*Z_4)+2.0d00*X*Y*Z_4)+
     &      3.0d00*X_2*Y*Z_4)*
     &      R81+(X*(X*(3.0d00*X_2*Y+6.0d00*Y*Z_2)+12.0d00*X*Y*Z_2)+
     &      18.0d00*X_2*Y*Z_2+3.0d00*Y*Z_4)*R71
     &      +(18.0d00*X_2*Y+18.0d00*Y*Z_2)*R61+9.0d00*Y*R51
        R(M,0,186)=X_4*Z_5*R91
     &      +(3.0d00*X_2*Z_5+X*(X*(Z_5+10.0d00* X_2*Z_3)+
     &      2.0d00*X*Z_5))*R81
     &      +(X*(X*(10.0d00*Z_3+15.0d00*X_2*Z)+20.0d00*X*Z_3)+
     &      30.0d00*X_2*Z_3+
     &      3.0d00*Z_5)*R71+(90.0d00*X_2*Z+30.0d00*Z_3)*R61+
     &      45.0d00*Z*R51
        R(M,0,187)=X_3*Y_6*R91
     &      +(X*(Y_6+15.0d00*X_2*Y_4)+2.0d00*X*Y_6)*R81
     &      +(X*(45.0d00*X_2*Y_2+15.0d00*Y_4)+30.0d00*X*Y_4)*R71
     &      +(X*(15.0d00*X_2+45.0d00*Y_2)+90.0d00*X*Y_2)*R61+
     &      45.0d00*X*R51
        R(M,0,188)=X_3*Y_5*Z*R91
     &      +(X*(Y_5*Z+10.0d00* X_2*Y_3*Z)+2.0d00*X*Y_5*Z)*R81
     &      +(X*(10.0d00*Y_3*Z+15.0d00*X_2*Y*Z)+20.0d00*X*Y_3*Z)*R71
     &      +45.0d00*X*Y*Z*R61
        R(M,0,189)=X_3*Y_4*Z_2*R91
     &      +(X*(Y_4*Z_2+X_2*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2)))+2.0d00*X*Y_4*Z_2)*R81
     &      +(X*(X_2*(6.0d00*Y_2+
     &      3.0d00*Z_2)+Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+
     &      3.0d00*Y_2*Z_2)+2.0d00*X*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2)))*R71+(X*(3.0d00*X_2+6.0d00*Y_2+3.0d00*Z_2)+
     &      2.0d00*X*(6.0d00*Y_2+3.0d00*Z_2))*R61+9.0d00*X*R51
        R(M,0,190)=X_3*Y_3*Z_3*R91
     &      +(X*(X_2*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)+Y_3*Z_3)+
     &      2.0d00*X*Y_3*Z_3)*R81
     &      +(X*(9.0d00*X_2*Y*Z+Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)+
     &      2.0d00*X*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3))*R71+
     &      27.0d00*X*Y*Z*R61
        R(M,0,191)=X_3*Y_2*Z_4*R91
     &      +(X*(Y_2*Z_4+X_2*(Z_4+6.0d00*Y_2*Z_2))+2.0d00*X*Y_2*Z_4)*
     &      R81+(X*(X_2*(3.0d00*Y_2+6.0d00*Z_2)+6.0d00*Y_2*Z_2+Z_4)+
     &      2.0d00*X*(Z_4+6.0d00*Y_2*Z_2))*R71
     &      +(X*(3.0d00*X_2+3.0d00*Y_2+6.0d00*Z_2)+
     &      2.0d00*X*(3.0d00*Y_2+6.0d00*Z_2))*R61
     &      +9.0d00*X*R51
        R(M,0,192)=X_3*Y*Z_5*R91
     &      +(X*(10.0d00*X_2*Y*Z_3+Y*Z_5)+2.0d00*X*Y*Z_5)*R81
     &      +(X*(15.0d00*X_2*Y*Z+10.0d00* Y*Z_3)+20.0d00*X*Y*Z_3)*R71
     &      +45.0d00*X*Y*Z*R61
        R(M,0,193)=X_3*Z_6*R91
     &      +(X*(Z_6+15.0d00*X_2*Z_4)+2.0d00*X*Z_6)*R81
     &      +(X*(45.0d00*X_2*Z_2+15.0d00*Z_4)+30.0d00*X*Z_4)*R71
     &      +(X*(15.0d00*X_2+45.0d00*Z_2)+90.0d00*X*Z_2)*R61+
     &      45.0d00*X*R51
        R(M,0,194)=X_2*Y_7*R91+(Y_7+21.0d00*X_2*Y_5)*R81
     &      +(105.0d00*X_2*Y_3+21.0d00*Y_5)*R71
     &      +(105.0d00*X_2*Y+105.0d00*Y_3)*R61+105.0d00*Y*R51
        R(M,0,195)=X_2*Y_6*Z*R91+(Y_6*Z+15.0d00*X_2*Y_4*Z)*
     &      R81+(45.0d00*X_2*Y_2*Z+15.0d00*Y_4*Z)*R71
     &      +(15.0d00*X_2*Z+45.0d00*Y_2*Z)*R61+15.0d00*Z*R51
        R(M,0,196)=X_2*Y_5*Z_2*R91
     &      +(Y_5*Z_2+X_2*(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2))*R81+(Y*(3.0d00*Y_2*Z_2+
     &      Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))+X_2*(Y*(6.0d00*Y_2+
     &      3.0d00*Z_2)+
     &      4.0d00*Y*(Y_2+Z_2)+8.0d00*Y*Z_2)+4.0d00*Y_3*Z_2)*R71
     &      +(15.0d00*X_2*Y+4.0d00*Y*(Y_2+Z_2)+8.0d00*Y*Z_2+
     &      Y*(6.0d00*Y_2+3.0d00*Z_2))*
     &      R61+15.0d00*Y*R51
        R(M,0,197)=X_2*Y_4*Z_3*R91
     &      +(Y_4*Z_3+X_2*(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*Y*Z_3)))*R81+(X_2*(18.0d00*Y_2*Z+3.0d00*Z_3)+
     &      Y*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)+3.0d00*Y_2*Z_3)*R71
     &      +(18.0d00*Y_2*Z+3.0d00*Z_3+9.0d00*X_2*Z)*R61+9.0d00*Z*R51
        R(M,0,198)=X_2*Y_3*Z_4*R91
     &      +(X_2*(Y*(Z_4+6.0d00*Y_2*Z_2)+2.0d00*Y*Z_4)+Y_3*Z_4)*R81
     &      +(X_2*(Y*(3.0d00*Y_2+6.0d00*Z_2)+
     &      12.0d00*Y*Z_2)+Y*(Z_4+6.0d00*Y_2*Z_2)+
     &      2.0d00*Y*Z_4)*R71
     &      +(9.0d00*X_2*Y+Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)*R61
     &      +9.0d00*Y*R51
        R(M,0,199)=X_2*Y_2*Z_5*R91
     &      +(Y_2*Z_5+X_2*(Z_5+10.0d00* Y_2*Z_3))*R81
     &      +(X_2*(10.0d00*Z_3+15.0d00*Y_2*Z)+10.0d00* Y_2*Z_3+Z_5)*R71
     &      +(15.0d00*Y_2*Z+10.0d00* Z_3+15.0d00*X_2*Z)*R61+
     &      15.0d00*Z*R51
        R(M,0,200)=X_2*Y*Z_6*R91+
     &      (15.0d00*X_2*Y*Z_4+Y*Z_6)*R81
     &      +(45.0d00*X_2*Y*Z_2+15.0d00*Y*Z_4)*R71
     &      +(15.0d00*X_2*Y+45.0d00*Y*Z_2)*R61+15.0d00*Y*R51
        R(M,0,201)=X_2*Z_7*R91+(Z_7+21.0d00*X_2*Z_5)*R81
     &      +(105.0d00*X_2*Z_3+21.0d00*Z_5)*R71
     &      +(105.0d00*X_2*Z+105.0d00*Z_3)*R61+105.0d00*Z*R51
        R(M,0,202)=X*Y_8*R91+28.0d00*X*Y_6*R81
     &      +210.0d00*X*Y_4*R71+420.0d00*X*Y_2*R61+105.0d00*X*R51
        R(M,0,203)=X*Y_7*Z*R91+21.0d00*X*Y_5*Z*R81
     &      +105.0d00*X*Y_3*Z*R71+105.0d00*X*Y*Z*R61
        R(M,0,204)=X*Y_6*Z_2*R91
     &      +X*(5.0d00*Y_4*Z_2+Y*(Y*(3.0d00*Y_2*Z_2+Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2))*R81
     &      +X*(Y*(Y*(6.0d00*Y_2+3.0d00*Z_2)+4.0d00*Y*(Y_2+Z_2)+
     &      8.0d00*Y*Z_2)+
     &      5.0d00*Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2)+15.0d00*Y_2*Z_2)*R71
     &      +X*(45.0d00*Y_2+15.0d00*Z_2)*R61+15.0d00*X*R51
        R(M,0,205)=X*Y_5*Z_3*R91
     &      +X*(Y*(3.0d00*Y_2*Z_3+Y*(Y*(3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*Y*Z_3))+4.0d00*Y_3*Z_3)*R81
     &      +X*(Y*(18.0d00*Y_2*Z+3.0d00*Z_3)+
     &      4.0d00*Y*(3.0d00*Y_2*Z+Z_3)+8.0d00*Y*Z_3)*R71+
     &      45.0d00*X*Y*Z*R61
        R(M,0,206)=X*Y_4*Z_4*R91
     &      +X*(3.0d00*Y_2*Z_4+Y*(Y*(Z_4+6.0d00*Y_2*Z_2)+
     &      2.0d00*Y*Z_4))*
     &      R81+X*(Y*(Y*(3.0d00*Y_2+6.0d00*Z_2)+12.0d00*Y*Z_2)+
     &      18.0d00*Y_2*Z_2+
     &      3.0d00*Z_4)*R71+X*(18.0d00*Y_2+18.0d00*Z_2)*R61+
     &      9.0d00*X*R51
        R(M,0,207)=X*Y_3*Z_5*R91
     &      +X*(Y*(Z_5+10.0d00* Y_2*Z_3)+2.0d00*Y*Z_5)*R81
     &      +X*(Y*(10.0d00*Z_3+15.0d00*Y_2*Z)+20.0d00*Y*Z_3)*R71
     &      +45.0d00*X*Y*Z*R61
        R(M,0,208)=X*Y_2*Z_6*R91+X*(Z_6+15.0d00*Y_2*Z_4)*R81
     &      +X*(45.0d00*Y_2*Z_2+15.0d00*Z_4)*R71
     &      +X*(15.0d00*Y_2+45.0d00*Z_2)*R61+15.0d00*X*R51
        R(M,0,209)=X*Y*Z_7*R91+21.0d00*X*Y*Z_5*R81
     &      +105.0d00*X*Y*Z_3*R71+105.0d00*X*Y*Z*R61
        R(M,0,210)=X*Z_8*R91+28.0d00*X*Z_6*R81
     &      +210.0d00*X*Z_4*R71+420.0d00*X*Z_2*R61+105.0d00*X*R51
        R(M,0,211)=Y_5*Y_4*R91+36.0d00*Y_7*R81+378.0d00*Y_5*R71
     &      +1260.0d00*Y_3*R61+945.0d00*Y*R51
        R(M,0,212)=Y_8*Z*R91+28.0d00*Y_6*Z*R81
     &      +210.0d00*Y_4*Z*R71+420.0d00*Y_2*Z*R61+105.0d00*Z*R51
        R(M,0,213)=Y_7*Z_2*R91
     &      +(6.0d00*Y_5*Z_2+Y*(5.0d00*Y_4*Z_2+Y*(Y*(3.0d00*Y_2*Z_2+
     &      Y*(Y*(Y_2+
     &      Z_2)+2.0d00*Y*Z_2))+4.0d00*Y_3*Z_2)))*R81+(6.0d00*Y*
     &      (3.0d00*Y_2*Z_2
     &      +Y*(Y*(Y_2+Z_2)+2.0d00*Y*Z_2))+Y*(Y*(Y*(6.0d00*Y_2+
     &      3.0d00*Z_2)+4.0d00*Y*
     &      (Y_2+Z_2)+8.0d00*Y*Z_2)+5.0d00*Y*(Y*(Y_2+Z_2)+
     &      2.0d00*Y*Z_2)+15.0d00*
     &      Y_2*Z_2)+24.0d00*Y_3*Z_2)*R71
     &      +(Y*(45.0d00*Y_2+15.0d00*Z_2)+24.0d00*Y*(Y_2+Z_2)+
     &      48.0d00*Y*Z_2+
     &      6.0d00*Y*(6.0d00*Y_2+3.0d00*Z_2))*R61+105.0d00*Y*R51
        R(M,0,214)=Y_6*Z_3*R91
     &      +(5.0d00*Y_4*Z_3+Y*(Y*(3.0d00*Y_2*Z_3+Y*(Y*
     &      (3.0d00*Y_2*Z+Z_3)+
     &      2.0d00*Y*Z_3))+4.0d00*Y_3*Z_3))*R81
     &      +(Y*(Y*(18.0d00*Y_2*Z+3.0d00*Z_3)+4.0d00*Y*(3.0d00*Y_2*Z+
     &      Z_3)+8.0d00*Y*Z_3)+
     &      5.0d00*Y*(Y*(3.0d00*Y_2*Z+Z_3)+2.0d00*Y*Z_3)+
     &      15.0d00*Y_2*Z_3)*R71
     &      +(135.0d00*Y_2*Z+15.0d00*Z_3)*R61+45.0d00*Z*R51
        R(M,0,215)=Y_5*Z_4*R91
     &      +(Y*(3.0d00*Y_2*Z_4+Y*(Y*(Z_4+6.0d00*Y_2*Z_2)+
     &      2.0d00*Y*Z_4))+
     &      4.0d00*Y_3*Z_4)*R81+(Y*(Y*(Y*(3.0d00*Y_2+6.0d00*Z_2)+
     &      12.0d00*Y*Z_2)+18.0d00*Y_2*Z_2+3.0d00*Z_4)+4.0d00*Y*(Z_4+
     &      6.0d00*Y_2*Z_2)+
     &      8.0d00*Y*Z_4)*R71
     &      +(Y*(18.0d00*Y_2+18.0d00*Z_2)+4.0d00*Y*(3.0d00*Y_2+
     &      6.0d00*Z_2)+48.0d00*Y*Z_2)*
     &      R61+45.0d00*Y*R51
        R(M,0,216)=Y_4*Z_5*R91
     &      +(3.0d00*Y_2*Z_5+Y*(Y*(Z_5+10.0d00* Y_2*Z_3)+
     &      2.0d00*Y*Z_5))*R81
     &      +(Y*(Y*(10.0d00*Z_3+15.0d00*Y_2*Z)+20.0d00*Y*Z_3)+
     &      30.0d00*Y_2*Z_3+3.0d00*Z_5)
     &      *R71+(90.0d00*Y_2*Z+30.0d00*Z_3)*R61+45.0d00*Z*R51
        R(M,0,217)=Y_3*Z_6*R91
     &      +(Y*(Z_6+15.0d00*Y_2*Z_4)+2.0d00*Y*Z_6)*R81
     &      +(Y*(45.0d00*Y_2*Z_2+15.0d00*Z_4)+30.0d00*Y*Z_4)*R71
     &      +(Y*(15.0d00*Y_2+45.0d00*Z_2)+90.0d00*Y*Z_2)*R61+
     &      45.0d00*Y*R51
        R(M,0,218)=Y_2*Z_7*R91+(Z_7+21.0d00*Y_2*Z_5)*R81
     &      +(105.0d00*Y_2*Z_3+21.0d00*Z_5)*R71
     &      +(105.0d00*Y_2*Z+105.0d00*Z_3)*R61+105.0d00*Z*R51
        R(M,0,219)=Y*Z_8*R91+28.0d00*Y*Z_6*R81
     &      +210.0d00*Y*Z_4*R71+420.0d00*Y*Z_2*R61+105.0d00*Y*R51
        R(M,0,220)=Z_5*Z_4*R91+36.0d00*Z_7*R81+378.0d00*Z_5*R71
     &      +1260.0d00*Z_3*R61+945.0d00*Z*R51
00100 continue
c
      end
#endif
