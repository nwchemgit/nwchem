       integer function moints_numgr( basis, blocklen )
C$Id: moints_gsh1.F,v 1.4 1995-10-19 00:24:01 d3g681 Exp $
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
*#ifdef DEBUG
      if (ga_nodeid().eq.0) then
        write(6,900)
 900    format(/,'Shell map:')
        write(6,901) (map(ii),ii=1,nsh)
 901    format(16i5)
        write(6,*)
      endif
*#endif
      return
      end

        




c
c  Memory and buffering arithmetic for 
c  semi-direct algorithm
c
c    Do a dummy run over outer shell loops
c    and calculate memory mapping of shell pairs
c    into global buffer.
c    Keep track of filling and flushing the buffer
c    to IO
c
      subroutine moints_locmap( basis, nsh, nbf, tol2e, shmap,
     $                          g_buf, gloc, nnbf, rloc,
     $                          iopass, rlen )
      implicit none
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, nsh, nbf                    ! [input]  Usual meaning
      double precision tol2e                     ! [input]  Tolerance
      integer shmap(nsh)                         ! [input]  Map loop index -> shell index
      integer g_buf                              ! [input]  GA buffer
      integer gloc(nbf,nbf)                      ! [output] Map basis -> GA memory
      integer nnbf                               ! [output] Screened basf pairs <= (nbf*(nbf+1))/2
      integer rloc(*)                            ! [output] Scatter array (IO record -> upper triangle)
      integer iopass(*)                          ! [output] IO pass for each task
      integer rlen(*)                            ! [output] Record length for IO write
c
      integer ish0, jsh0, ish, jsh
      integer ibflo, ibfhi, jbflo, jbfhi
      integer i, j, ijlen, incr
      integer ntask, ipass
      integer gtype, dim1, dim2
      logical status
      
      call ga_inquire(g_buf,gtype,dim1,dim2)
      nnbf = 0
      incr =  0
      ntask = 0
      ipass = 1
      do ish0=1,nsh
        do jsh0=1,ish0
          ish = max(shmap(ish0),shmap(jsh0))
          jsh = min(shmap(ish0),shmap(jsh0))
          if ((schwarz_shell(ish,jsh)*schwarz_max()).ge.tol2e) then
            ntask = ntask + 1
            status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
            status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
            ijlen = (ibfhi - ibflo + 1)*(jbfhi - jbflo + 1)
            if ((incr+ijlen).gt.dim1) then
              rlen(ipass) = incr
              ipass = ipass + 1
              incr = 0
            endif
            iopass(ntask) = ipass
            do i=ibflo,ibfhi
              do j=jbflo,jbfhi
                if ((ish.ne.jsh).or.((ish.eq.jsh).and.(i.ge.j))) then
                  incr = incr + 1
                  gloc(i,j) = incr
                  gloc(j,i) = incr
                  nnbf = nnbf + 1
                  rloc(nnbf) = (max(i,j)-1)*nbf+min(i,j)
                endif
              enddo
            enddo
          endif
        enddo
      enddo
      rlen(ipass) = incr
      do i=1,ga_nnodes()
        iopass(ntask+i) = ipass + 1
      enddo		

      IF (GA_NODEID().EQ.0) THEN
c$$$      WRITE(6,901)
c$$$ 901  FORMAT(/,' LOCATION MAP: ',/)
c$$$      DO I=1,NBF
c$$$        WRITE(6,900) (GLOC(I,J),J=1,NBF)
c$$$      ENDDO
c$$$      WRITE(6,902)
c$$$ 902  FORMAT(/,'REVERSE LOCATION MAP:',/)
c$$$      WRITE(6,900) (RLOC(I),I=1,NNBF)
c$$$      WRITE(6,*)
c$$$        WRITE(6,903)
c$$$ 903    FORMAT(/,'IO PASS:',/)
c$$$        WRITE(6,900) (IOPASS(I),I=1,NTASK+GA_NNODES())
c$$$        WRITE(6,*)
c$$$      WRITE(6,904)
c$$$ 904  FORMAT(/,'RECORD LENGTHS:',/)
c$$$      WRITE(6,900) (RLEN(I),I=1,IPASS)
      ENDIF

 900  FORMAT(16I4)
      return
      end


          
          
          
          
