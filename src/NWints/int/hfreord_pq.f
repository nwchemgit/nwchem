      subroutine hfreord_pq(eri,scr,szp,szq)
c $Id: hfreord_pq.f,v 1.1 1995-03-28 20:12:21 d3e129 Exp $ 
c
c routine to reorder <Q,P> integrals to <P,Q> 
c
      implicit none
c

      integer szp  ! [input] size of P integral block (la2*lb2)
      integer szq  ! [input] size of Q integral block (lc2*ld2)
      double precision scr(szp*szq)  ! [input]  <Q|P> integrals
      double precision eri(szp*szq)  ! [output] <P|Q> integrals
c
      integer ito, ifrom, i, j
c
c  reorder p and q integrals 
c
      ito = 0
      do 00100 i = 1,szp
        ifrom = i
        do 00200 j = 1,szq
          ito = ito + 1
          eri(ito) = scr(ifrom)
          ifrom = ifrom + szp
00200   continue
00100 continue
      end
