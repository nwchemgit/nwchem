       integer function moints_numgr( basis, blocklen )
C$Id: moints_gsh1.F,v 1.2 1995-02-02 23:18:47 d3g681 Exp $
       implicit none
#include "bas.fh"
       integer basis, blocklen
       integer nsh, grlen, ngr, ishlen, ibfhi, ibflo, ish

       if (.not.bas_numcont(basis,nsh))
     $      call errquit('moints_numshgr: cannot get info',0)

       ngr = 1
       grlen = 0
       do ish=1,nsh
         if (.not. bas_cn2bfr(basis,ish,ibflo,ibfhi))
     $        call errquit('moints_numshgr: bas_cn2bfr',ish)
         ishlen = ibfhi - ibflo + 1
         if ((grlen+ishlen).gt.blocklen) then
           grlen = ishlen
           ngr = ngr + 1
         else
           grlen = grlen + ishlen
         endif
       enddo
       moints_numgr = ngr
       return
       end

           






       subroutine moints_grmap( basis, blocklen, ngr, grmap )
       implicit none
#include "bas.fh"
       integer basis, blocklen, ngr, grmap(4,ngr)
       integer nsh, grlen, grct, ishlen, ibfhi, ibflo, ish

       if (.not.bas_numcont(basis,nsh))
     $      call errquit('moints_numshgr: cannot get info',0)

       grlen = 0
       grct = 0
       do ish=1,nsh
         if (.not. bas_cn2bfr(basis,ish,ibflo,ibfhi))
     $        call errquit('moints_numshgr: bas_cn2bfr',ish)
         ishlen = ibfhi - ibflo + 1
         if ((grlen.eq.0).or.(grlen+ishlen).gt.blocklen) then
           grct = grct + 1
           grlen = ishlen
           grmap(1,grct) = ish
           grmap(2,grct) = ish
           grmap(3,grct) = ibflo
           grmap(4,grct) = ibfhi
         else
           grlen = grlen + ishlen
           grmap(2,grct) = ish
           grmap(4,grct) = ibfhi
         endif
       enddo

#ifdef DEBUG
       do ish=1,ngr
         write(6,912) ish,grmap(1,ish),grmap(2,ish),
     $                    grmap(3,ish),grmap(4,ish)
 912     format(i4,'  shlo:shhi ',2i4,'  bflo:bfhi ',2i4)
       enddo
#endif
       return
       end






      integer function gr_len( ngr, grmap, ish, ilo, ihi )
      implicit none
      integer ngr
      integer grmap(4,ngr)
      integer ish, ilo, ihi

      gr_len = grmap(4,ish) - grmap(3,ish) + 1
      ilo = grmap(3,ish)
      ihi = grmap(4,ish)
      return
      end




      subroutine gr_sh_idx( ngr, grmap, ish, ilo, ihi )
      implicit none
      integer ngr
      integer grmap(4,ngr)
      integer ish, ilo, ihi

      ilo = grmap(1,ish)
      ihi = grmap(2,ish)
      return
      end








      subroutine moints_shorder( basis, nsh, maxshlen, map )
      implicit none
#include "global.fh"
#include "bas.fh"
      integer basis
      integer nsh, maxshlen
      integer map(nsh)
      integer shlen, ii, ish, ibflo, ibfhi
      logical status

      ii = 0
      do shlen=maxshlen,1,-1
        do ish=1,nsh
          status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
          if ((ibfhi - ibflo + 1).eq.shlen) then
            ii = ii + 1
            map(ii) = ish
          endif
        enddo
      enddo
#ifdef DEBUG
      if (ga_nodeid().eq.0) then
        write(6,900)
 900    format(/,'Shell map:')
        write(6,901) (map(ii),ii=1,nsh)
 901    format(16i5)
        write(6,*)
      endif
#endif
      return
      end

        
