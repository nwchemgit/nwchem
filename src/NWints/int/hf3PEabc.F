      subroutine hf3PEabc(Pabc, E, R0, IJK,
     &    Nabc, La, Lb, Lc, Lg, Lg3, nint, lim_nint)
c $Id$
      implicit none
#include "errquit.fh"
c
      integer La, Lb, Lc, Lg, Lg3, Nabc, nint, lim_nint
      integer IJK(0:Lg,0:Lg,0:Lg)
      double precision E(Nabc,3,0:Lg,0:La,0:Lb,0:Lc)
      double precision R0(Nabc,Lg3)
      double precision Pabc(*)
      double precision pabc_int
c
      integer Nxyz(3)
      integer La2, Lb2, Lc2, ica, icb, icc
      integer Ia, Ja, Ka
      integer Ib, Jb, Kb
      integer Ic, Jc, Kc
      integer Ip, Jp, Kp
      integer np, mp
c
c Define the number of shell components on each center.

      La2 = ((La+1)*(La+2))/2
      Lb2 = ((Lb+1)*(Lb+2))/2
      Lc2 = ((Lc+1)*(Lc+2))/2
      nint = la2*lb2*lc2
      if (nint.gt.lim_nint) then
        write(6,*)' la = ',la
        write(6,*)' lb = ',lb
        write(6,*)' lc = ',lc
        write(6,*)' nint     = ',nint
        write(6,*)' lim_nint = ',lim_nint
        call errquit('hf3PEabc: something hosed',911, INT_ERR)
      endif

c Loop over shell components.
c
      nint = 0   
c
      do 00100 ica = 1,La2

        call getNxyz(La,ica,Nxyz)
        Ia = Nxyz(1)
        Ja = Nxyz(2)
        Ka = Nxyz(3)
        
        do 00200 icb = 1,Lb2
          
          call getNxyz(Lb,icb,Nxyz)
          Ib = Nxyz(1)
          Jb = Nxyz(2)
          Kb = Nxyz(3)
          
          do 00300 icc = 1, Lc2

            call getNxyz(Lc,icc,Nxyz)
            Ic = Nxyz(1)
            Jc = Nxyz(2)
            Kc = Nxyz(3)

            if (nint.eq.lim_nint) then
              write(6,*)' nint =   ',nint
              write(6,*)' lim_nint ',lim_nint
              call errquit
     &            ('hf3PEabc: buffer too small ',911, MEM_ERR)
            endif

            nint = nint + 1


            pabc_int = 0.0d00
            do 00400 Ip = 0, (Ia + Ib + Ic)
              do 00500 Jp = 0, (Ja + Jb +Jc)
                do 00600 Kp = 0, (Ka + Kb +Kc)

                  np = IJK(Ip, Jp, Kp)
                  
                  do 00700 mp = 1, Nabc
                    pabc_int = pabc_int + 
     &                  E(mp,1,Ip,Ia,Ib,Ic) *
     &                  E(mp,2,Jp,Ja,Jb,Jc) *
     &                  E(mp,3,Kp,Ka,Kb,Kc) * R0(mp,np)
                      
00700             continue
00600           continue
00500         continue
00400       continue
            Pabc(nint) = pabc_int
00300     continue
00200   continue
00100 continue
c
      end
