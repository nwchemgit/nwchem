      Subroutine hfmkr(R,IJK,PQ,NAB,Lp,Lp3)
c $Id: hfmkr.F,v 1.6 1996-10-11 10:13:01 d3e129 Exp $

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
