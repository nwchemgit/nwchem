      Subroutine hf2oi(E,Sab,Nint,NPP,La,Lb,Li,canAB)
c $Id: hf2oi.f,v 1.5 1998-05-05 21:49:04 d3e129 Exp $

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> 2-Center Overlap Integrals

      Dimension Sab(Nint)

      Dimension Nxyz(3)
c
c Compute the 2-ctr overlap integrals.
c
c     formula:
c
c             Ia,Ib   Ja,Jb   Ka,Kb
c     Sab = Ex      Ey      Ez
c             0       0       0
c
c******************************************************************************

c Initialize the block of NAIs.

      do 10 nn = 1,Nint
       Sab(nn) = 0.D0
   10 continue

c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2

c Loop over shell components.

      nn = 0

      do 50 ma = 1,La2

c Define the angular momentum indices for shell "A".

       call getNxyz(La,ma,Nxyz)

       Ia = Nxyz(1)
       Ja = Nxyz(2)
       Ka = Nxyz(3)

       if( canAB )then
        mb_limit = ma
       else
        mb_limit = Lb2
       end if

       do 40 mb = 1,mb_limit

c Define the angular momentum indices for shell "B".

        call getNxyz(Lb,mb,Nxyz)

        Ib = Nxyz(1)
        Jb = Nxyz(2)
        Kb = Nxyz(3)

        nn = nn + 1

        do 20 mp = 1,NPP
*         e_tmp = E(1,mp,0,Ia,Ib)*
*     &           E(2,mp,0,Ja,Jb)*
*     &           E(3,mp,0,Ka,Kb)
*         Sab(nn) = Sab(nn) + e_tmp
*         write(6,10000)
*     &         ' int=',nn,' mp =',mp,
*     &         ' val = ',e_tmp,' integral=',Sab(nn)

         Sab(nn) = Sab(nn) + E(1,mp,0,Ia,Ib)*
     &                       E(2,mp,0,Ja,Jb)*
     &                       E(3,mp,0,Ka,Kb)
   20   continue

   40  continue

   50 continue
*10000 format(a,i4,a,i4,a,1pd20.10,a,1pd20.10)
      end
      Subroutine hf2oi_gc(E,Sab,SabP,SabH,Nint,
     &    NCA,NCB,NPP,
     &    La,Lb,Li,gct_a,gct_b,canAB)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB

c--> Hermite Linear Expansion Coefficients

      Dimension E(3,NPP,0:((La+Li)+(Lb+Li)),0:(La+Li),0:(Lb+Li))

c--> 2-Center Overlap Integrals

      Dimension Sab(Nint*nca*ncb)
      double precision SabP(NPP,Nint)
      double precision SabH(NPP,Nint,NCA)
c--> general contraction matrices
      double precision gct_a(NPP,NCA) ! [output] general contraction coefs for A multiply
      double precision gct_b(NCB,NPP) ! [output] general contraction coefs for B multiply

      Dimension Nxyz(3)
*      double precision SabC
*      double precision sscp(50000)
c
c Compute the 2-ctr overlap integrals.
c
c     formula:
c
c             Ia,Ib   Ja,Jb   Ka,Kb
c     Sab = Ex      Ey      Ez
c             0       0       0
c
c******************************************************************************

c Initialize the block of NAIs.

      call dfill(Nint*nca*ncb,0.0d00,Sab,1)
      call dfill(NPP*Nint,0.0d00,SabP,1)
      call dfill(NPP*Nint*NCA,0.0d00,SabH,1)
*      do 10 nn = 1,(Nint*nca*ncb)
*       Sab(nn) = 0.D0
*   10 continue

c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2

      if (nint.ne.(la2*lb2)) call errquit
     &    ('hf2oi_gc: nint/la2/lb2 error',(nint-(la2*lb2)))
c Loop over shell components.

      nn = 0

      do 50 ma = 1,La2

c Define the angular momentum indices for shell "A".

       call getNxyz(La,ma,Nxyz)

       Ia = Nxyz(1)
       Ja = Nxyz(2)
       Ka = Nxyz(3)

       if( canAB )then
        mb_limit = ma
       else
        mb_limit = Lb2
       end if

       do 40 mb = 1,mb_limit

c Define the angular momentum indices for shell "B".

        call getNxyz(Lb,mb,Nxyz)

        Ib = Nxyz(1)
        Jb = Nxyz(2)
        Kb = Nxyz(3)

        nn = nn + 1

        do 20 mp = 1,NPP
         SabP(mp,nn) = SabP(mp,nn) + E(1,mp,0,Ia,Ib)*
     &                               E(2,mp,0,Ja,Jb)*
     &                               E(3,mp,0,Ka,Kb)
   20   continue

   40  continue

   50 continue

*      icount = 0
*      do ica = 1,nca
*        do icb = 1,ncb
*          write(6,*)' simple integrals for ica = ',ica,'and icb = ',icb
*          do nn = 1,Nint
*            SabC = 0.0d00
*            do mp = 1,NPP
*              e_tmp = SabP(mp,nn)*gct_a(mp,ica)*gct_b(icb,mp)
*              SabC = SabC + e_tmp
*              write(6,10000)
*     &              ' int=',nn,' mp =',mp,
*     &              ' val = ',e_tmp,' integral=',SabC
*            enddo
*            icount = icount + 1
*            write(79,*)' simple ',ica,icb,nn,SabC,icount
*          enddo
*        enddo
*      enddo
            
c take primitives and half transformed multiplied by A general contraction matrix

      do 00100 ica = 1,NCA
        do 00200 nn = 1,Nint
          do 00300 mp = 1,NPP
            SabH(mp,nn,ica) = SabP(mp,nn)*gct_a(mp,ica)
00300     continue
00200   continue
00100 continue

c* broke ica/iii/icb norm=16 sp
c* broke ica/icb/iii norm=16 sp
c* broke iii/ica/icb norm=12 sp
c* broke iii/icb/ica norm=13.4349
c* broke icb/iii/ica norm=17.5349
c* try icb/ica/iii
      nn = 0
*      iwiw = 1
      do 00400 ica = 1,NCA
        do 00500 icb = 1,NCB
          do 00600 iii = 1,Nint
*            if (iwiw.eq.0)
*     &            write(6,*)'complex integrals for ica = ',
*     &            ica,'and icb = ',icb
*            iwiw = iwiw + 1
            nn = nn + 1
            Sab(nn) = 0.0d00
            do 00700 mp = 1,NPP
              e_tmp = SabH(mp,iii,ica)*gct_b(icb,mp)
              Sab(nn) = Sab(nn) + e_tmp
*              write(6,10000)
*     &              ' int=',nn,' mp =',mp,
*     &              ' val = ',e_tmp,' integral=',Sab(nn)
00700       continue
*            write(80,*)' complex ',ica,icb,iii,Sab(nn),nn
00600     continue
00500   continue
00400 continue
10000 format(a,i4,a,i4,a,1pd20.10,a,1pd20.10)
c
c      copy integrals
c
      call dcopy((nint*nca*ncb),Sab,1,SabH,1)
      call hf1_tran_shift(Sab,SabH,(nca*ncb),la,lb,nca,ncb)
      end
