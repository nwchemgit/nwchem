       integer function moints_numgr( basis, blocklen )
C$Id$
       implicit none
#include "errquit.fh"
#include "bas.fh"
       integer basis, blocklen
       integer nsh, grlen, ngr, ishlen, ibfhi, ibflo, ish

       if (.not.bas_numcont(basis,nsh))
     $      call errquit('moints_numshgr: cannot get info',0, BASIS_ERR)

       ngr = 1
       grlen = 0
       do ish=1,nsh
         if (.not. bas_cn2bfr(basis,ish,ibflo,ibfhi))
     $        call errquit('moints_numshgr: bas_cn2bfr',ish, BASIS_ERR)
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
#include "errquit.fh"
#include "bas.fh"
       integer basis, blocklen, ngr, grmap(4,ngr)
       integer nsh, grlen, grct, ishlen, ibfhi, ibflo, ish

       if (.not.bas_numcont(basis,nsh))
     $      call errquit('moints_numshgr: cannot get info',0, BASIS_ERR)

       grlen = 0
       grct = 0
       do ish=1,nsh
         if (.not. bas_cn2bfr(basis,ish,ibflo,ibfhi))
     $        call errquit('moints_numshgr: bas_cn2bfr',ish, BASIS_ERR)
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







c
c  Reorder and group shells in optimal
c  order for AO integral blocking and
c  dynamic load balancing
c
c
      subroutine moints_shorder( basis, nsh, nbf, blen,
     $                           ngr, glo, ghi, shmap,
     $                           bfmap, rbfmap )
     $                  
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "bas.fh"
      integer basis
      integer nsh, nbf
      integer blen                         ! [input]  Blocking length (upper bound on length of groups)
      integer ngr                          ! [output] Number of groups formed
      integer glo(*)                       ! [output] Low indices of groups
      integer ghi(*)                       ! [output] High indices of groups
      integer shmap(nsh)                   ! [output] Shell new  --> orig
      integer bfmap(nbf)                   ! [output] Basis new  --> orig
      integer rbfmap(nbf)                  ! [output] Basis orig --> new
c      
      integer ice, itype, nprm, ncnt, spch, wght, maxwght
      integer ii, j, k, i, ish, jsh, ibflo, ibfhi
      integer tlen, ilen, glen, nsht
      integer l_t, k_t
      logical status, more, oreorder
      integer lwght(-2:6)
      data lwght/5,3,1,2,4,6,7,8,9/
      data oreorder/.true./                 ! toggle shell reordering

c
c  Ordering by shell type (the default)
c
      if (oreorder) then
        call ifill(nsh,0,bfmap,1)
        do jsh=1,nsh
          maxwght = 0
          do ish=1,nsh
            status = bas_cn2uce( basis, ish, ice )
            status = bas_continfo(basis, ish, itype, nprm, ncnt, spch )
            wght = ice + 50*(nprm + 50*(ncnt + 50*lwght(itype)))
            if ((bfmap(ish).eq.0).and.(wght.gt.maxwght)) then
              maxwght = wght
              ii = ish
            endif
          enddo
          bfmap(ii) = maxwght
          shmap(jsh)  = ii
        enddo
c
c  Types of shells
c  Large groups of similar shells may be split
c  to fit into a block.
c
        if (.not.ma_push_get(mt_int,(4*nsh),'gsh tmp',l_t,k_t))
     $       call errquit('moints_shorder: cannot allocate tmp',(4*nsh),
     &       MA_ERR)
        wght = 0
        nsht = 0
        tlen = 0
        do ii=1,nsh
          status = bas_cn2bfr(basis,shmap(ii),ibflo,ibfhi)
          ilen = ibfhi - ibflo + 1
          if ((bfmap(shmap(ii)).ne.wght).or.((tlen+ilen).gt.blen)) then
            nsht = nsht + 1
            wght = bfmap(shmap(ii))
            int_mb(k_t+4*(nsht-1)) = ii
            int_mb(k_t+4*(nsht-1)+2) = 0
            tlen = 0
          endif
          tlen = tlen + ilen
          int_mb(k_t+4*(nsht-1)+1) = ii
          int_mb(k_t+4*(nsht-1)+3) = tlen
        enddo
c     
c  Group together types subject to block 
c  length constraint. Obtain final shell 
c  ordering based on contigious groups of types.
c
        call ifill(nsh,0,bfmap,1)
        j = 0
        ngr = 0
 101    continue
        glen = 0
        ngr = ngr + 1
        glo(ngr) = 0
 102    continue
        more = .false.
        do ii=0,nsht-1
          if (int_mb(k_t+4*ii+2).eq.0) then
            ilen = int_mb(k_t+4*ii+3)
            if ((glen+ilen).le.blen) then
              int_mb(k_t+4*ii+2) = 1
              glen = glen + ilen
              more = .true.
              if (glo(ngr).eq.0) glo(ngr) = j + 1
              do k=int_mb(k_t+4*ii),int_mb(k_t+4*ii+1)
                j = j + 1
                bfmap(j) = shmap(k)
              enddo
              ghi(ngr) = j
            endif
          endif
        enddo
        if (more) goto 102
        do ii=0,nsht-1
          more = more.or.(int_mb(k_t+4*ii+2).eq.0)
        enddo
        if (more) goto 101
        call icopy(nsh,bfmap,1,shmap,1)
c     
c  Debug print
c     
        if (util_print('shell order',print_debug).and.
     $       (ga_nodeid().eq.0)) then
          write(6,841) nsht
 841      format(/,'Number of shell types:',i5)
          do ii=0,nsht-1
            write(6,842) int_mb(k_t+4*ii),int_mb(k_t+4*ii+1),
     $           int_mb(k_t+4*ii+3)
 842        format(5x,3i5) 
          enddo
        endif
        if (.not.ma_pop_stack(l_t))
     $       call errquit('moints_shorder: cannot pop stack',0,
     &       MA_ERR)
      else
c
c  No reordering of shells -- still need grouping
c
        do ish=1,nsh
          shmap(ish) = ish
        enddo
        ngr = 1
        glo(ngr) = 0
        glen = 0
        do ish=1,nsh
          if (glo(ngr).eq.0) glo(ngr) = ish
          status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
          ilen = ibfhi - ibflo + 1
          if ((glen+ilen).le.blen) then
            glen = glen + ilen
          else
            glen = ilen
            ngr = ngr + 1
            glo(ngr) = ish
          endif
          ghi(ngr) = ish
        enddo
      endif
c
c  Basis functions maps
c
      ii = 0
      do ish=1,nsh
        status = bas_cn2bfr(basis,shmap(ish),ibflo,ibfhi)
        do i=ibflo,ibfhi
          ii = ii + 1
          bfmap(ii) = i         ! new idx   -->  orig idx
          rbfmap(i) = ii        ! orig idx  -->  new idx
        enddo
      enddo
c
c  More debug print
c
      if (util_print('shell order',print_debug).and.
     $   (ga_nodeid().eq.0)) then
        write(6,550) ngr
 550    format(/,'Number of groups:',i5,/)
        do k=1,ngr              ! Example of how to use group info
          write(6,552) k
 552      format(/,2x,'Group:',i2,/,
     $             2x,'idx',3x,'sh',3x,'l',4x,'ct',3x,'pr',3x,'ce',/,
     $             2x,29('-'))
          do ii=glo(k),ghi(k)
            ish = shmap(ii)
            status = bas_continfo(basis, ish, itype, nprm, ncnt, spch )
            status = bas_cn2uce( basis, ish, ice )
            write(6,553) ii, ish, itype, ncnt, nprm, ice
 553        format(6i5)
          enddo
        enddo
        write(6,*)
      endif
      if (util_print('shell map',print_debug).and.
     $   (ga_nodeid().eq.0)) then
        write(6,900)
 900    format(/,'Shell map:')
        write(6,901) (shmap(ii),ii=1,nsh)
 901    format(16i5)
        write(6,902)
 902    format(/,'Basis function maps:')
        write(6,901) (bfmap(i),i=1,nbf)
        write(6,903)
 903    format(/,'Reverse Basis function maps:')
        write(6,901) (rbfmap(i),i=1,nbf)
        write(6,*)
      endif

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
     $                          ioreclen, gloc, nnbf, rloc,
     $                          nioflush, ionext, rlen )
      implicit none
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, nsh, nbf                    ! [input]  Usual meaning
      double precision tol2e                     ! [input]  Tolerance
      integer shmap(nsh)                         ! [input]  Map loop index -> shell index
      integer ioreclen                              ! [input]  I/O reclen (double words)
      integer gloc(nbf,nbf)                      ! [output] Map basis -> GA memory
      integer nnbf                               ! [output] Screened basf pairs <= (nbf*(nbf+1))/2
      integer rloc(*)                            ! [output] Scatter array (IO record -> upper triangle)
      integer nioflush                           ! [output] Number of IO flushes
      integer ionext(*)                          ! [output] Next values for IO flush
      integer rlen(*)                            ! [output] Record length for IO flush
      integer ish0, jsh0, ish, jsh
      integer ibflo, ibfhi, jbflo, jbfhi
      integer i, j, ijlen, incr
      integer ntask
      logical status
      
      nnbf = 0
      incr =  0
      ntask = 0
      nioflush = 0
      do ish0=1,nsh
        do jsh0=1,ish0
          ish = max(shmap(ish0),shmap(jsh0))
          jsh = min(shmap(ish0),shmap(jsh0))
          if ((schwarz_shell(ish,jsh)*schwarz_max()).ge.tol2e) then
            ntask = ntask + 1
            status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
            status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
            ijlen = (ibfhi - ibflo + 1)*(jbfhi - jbflo + 1)
            if ((incr+ijlen).gt.ioreclen) then
              nioflush = nioflush + 1
              rlen(nioflush) = incr
              ionext(nioflush) = ntask - 1
              incr = 0
            endif
C            iopass(ntask) = ipass
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
      if (incr.gt.0) then
        nioflush = nioflush + 1
        rlen(nioflush) = incr
        ionext(nioflush) = ntask
      endif
C      do i=1,ga_nnodes()
C        iopass(ntask+i) = ipass + 1
C      enddo		

      IF (GA_NODEID().EQ.0) THEN
         write(6,*) ' vrange nnbf ',nnbf
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
c$$$        WRITE(6,900) (IONEXT(I),I=1,NIOFLUSH)
c$$$        WRITE(6,*)
CCc$$$      WRITE(6,904)
c$$$ 904  FORMAT(/,'RECORD LENGTHS:',/)
c$$$      WRITE(6,900) (RLEN(I),I=1,IPASS)
      ENDIF

 900  FORMAT(16I4)
      return
      end


          
          
          
          
