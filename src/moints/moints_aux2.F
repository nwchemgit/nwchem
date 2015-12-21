#ifdef OLD_AOBLOCK
c
c  General interface to block AO integral routine
c  Branch on whether to use 
c          - labels/shell blocking
c  or
c          - quartet at a time
c  
c
      subroutine moints_gblk_old( basis, ish, jsh, kshlo, kshhi,
     $                            lshlo, lshhi, schw_ij, tol2e,
     $                            erilen, eri, scrlen, iscr,
     $                            ibflo, ibfhi, jbflo, jbfhi,
     $                            kblo, kbhi, lblo, lbhi,
     $                            ssbb, osym, oblk )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      logical osym
      logical oblk
c
      integer mxshq, imem
      integer k_tmp, l_tmp
      integer k_iv, k_jv, k_kv, k_lv
      integer k_il, k_jl, k_kl, k_ll
      integer k_qv, l_qv
c
c
c
      if (.not.(oblk)) then
        call moints_gblk_sq_old( basis, ish, jsh, kshlo, kshhi,
     $                       lshlo, lshhi, schw_ij, tol2e,
     $                       erilen, eri, scrlen, iscr,
     $                       ibflo, ibfhi, jbflo, jbfhi,
     $                       kblo, kbhi, lblo, lbhi,
     $                       ssbb, osym )
      else
**         call errquit('BLOCKING INTERFACE', 0)
        mxshq = (kshhi - kshlo + 1)*(lshhi - lshlo + 1)
        imem = 4*mxshq + 4*erilen
        if (.not. ma_push_get(MT_INT, imem, 'gblk tmp',l_tmp, k_tmp))
     $     call errquit(
     $       'moints_gblk: no memory (stack) for labels',imem, MA_ERR)
        k_iv = k_tmp
        k_jv = k_tmp + mxshq
        k_kv = k_tmp + 2*mxshq
        k_lv = k_tmp + 3*mxshq
        k_il = k_tmp + 4*mxshq
        k_jl = k_tmp + 4*mxshq + erilen
        k_kl = k_tmp + 4*mxshq + 2*erilen
        k_ll = k_tmp + 4*mxshq + 3*erilen
        if (.not.ma_push_get(MT_DBL, mxshq, 'gblk q',l_qv, k_qv))
     $     call errquit('moints_gblk: no memory for sym factors',imem,
     &       MA_ERR)
        call moints_gblk_mq_old( basis, ish, jsh, kshlo, kshhi,
     $                       lshlo, lshhi, schw_ij, tol2e,
     $                       erilen, eri, scrlen, iscr,
     $                       ibflo, ibfhi, jbflo, jbfhi,
     $                       kblo, kbhi, lblo, lbhi,
     $                       ssbb, osym, mxshq, dbl_mb(k_qv),
     $                       int_mb(k_iv), int_mb(k_jv),
     $                       int_mb(k_kv), int_mb(k_lv),
     $                       int_mb(k_il), int_mb(k_jl),
     $                       int_mb(k_kl), int_mb(k_ll) )
        if (.not.ma_pop_stack(l_qv))
     $    call errquit('moints_gblk: failed to pop',0, MA_ERR)
        if (.not.ma_pop_stack(l_tmp))
     $    call errquit('moints_gblk: failed to pop',0, MA_ERR)
      endif
c
c
c
      return
      end

        



      





c
c  This uses block NWints interface
c  with multiple shell quartets and labels
c
c
      subroutine moints_gblk_mq_old( basis, ish, jsh, kshlo, kshhi,
     $                           lshlo, lshhi, schw_ij, tol2e,
     $                           erilen, eri, scrlen, iscr,
     $                           ibflo, ibfhi, jbflo, jbfhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb, osym, mxshq, qv4, 
     $                           ishv, jshv, kshv, lshv,
     $                           ilab, jlab, klab, llab )
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      logical osym
      integer mxshq
      double precision qv4(mxshq)
      integer ishv(mxshq), jshv(mxshq), lshv(mxshq), kshv(mxshq)
      integer ilab(erilen), jlab(erilen), klab(erilen), llab(erilen)
c
c
c
      double precision beffect, q4
      integer nshq, neri, neritotal, ieri, nblksiz
      integer ii, jj, kk, ll, ksh, lsh
      integer ilen, jlen, klen
      logical more
      logical odoit
      logical sym_shell_quartet
      external sym_shell_quartet
c
      logical intb_2e4c, intb_init4c
      external intb_2e4c, intb_init4c
c
      neritotal = 0
      nblksiz = (lbhi-lblo+1)*(kbhi-kblo+1)*
     $          (jbfhi-jbflo+1)*(ibfhi-ibflo+1)
      call dfill(nblksiz, 0.d0, ssbb, 1 )
c
c  Fill out shell arrays with interacting quartet labels
c        
      q4 = 1.d0
      nshq = 0
      do ksh=kshlo,kshhi
        do lsh=lshlo,lshhi
          if (schwarz_shell(ksh,lsh)*schw_ij.ge.tol2e) then
            odoit = .true.
            if (osym) odoit = sym_shell_quartet(basis, ish, jsh,
     $                                          ksh, lsh, q4)
            if (odoit) then
              nshq = nshq + 1
              kshv(nshq) = ksh
              lshv(nshq) = lsh
              qv4(nshq)  = q4
            endif
          endif
        enddo
      enddo
      call ifill( nshq, ish, ishv, 1 )
      call ifill( nshq, jsh, jshv, 1 )
c
c  Prepare texas for shell block
c
      if (.not. intb_init4c( basis, ishv, jshv, basis, kshv, lshv,
     $     nshq, qv4, osym, 
     $     scrlen, iscr, erilen, beffect ))
     $   call errquit('moints_gblk: intb_init4c failed',0,
     &       INT_ERR)
c
c  Get some batch of integrals
c
 100  more = intb_2e4c( basis, ishv, jshv, basis, kshv, lshv,
     $                  nshq, qv4, osym, tol2e, .false.,
     $                  ilab, jlab, klab, llab,
     $                  eri, erilen, neri, scrlen, iscr )
c
c  Unpack labels and assign to integrals
c
      if (neri.gt.0) then
        neritotal = neritotal + neri
        do ieri=1,neri
          ii = ilab(ieri)
          jj = jlab(ieri)
          kk = klab(ieri)
          ll = llab(ieri)
          ssbb(ll,kk,jj,ii) = eri(ieri)
        enddo
      endif
      if (more) goto 100
c
c  Reflect diagonal shell block
c
      klen = kbhi - kblo + 1
      ilen = ibfhi - ibflo + 1
      jlen = jbfhi - jbflo + 1
      if ((kshlo.eq.lshlo).and.(osym))
     $   call eriblktr(ilen, jlen, klen, ssbb )

c$$$      LL = (LBHI - LBLO + 1)
c$$$      KK = (KBHI - KBLO + 1)
c$$$      WRITE(6,881) IBFLO, JBFLO
c$$$ 881  FORMAT(/,'+++',2I5)
c$$$      CALL MOINTS_MATPRINT(LL,KK,SSBB(LBLO,KBLO,JBFLO,IBFLO))
      return
      end







c
c  This is using older NWints interface
c  with one shell quartet at a time
c
c
      subroutine moints_gblk_sq_old( basis, ish, jsh, kshlo, kshhi,
     $                           lshlo, lshhi, schw_ij, tol2e,
     $                           erilen, eri, scrlen, iscr, 
     $                           ibflo, ibfhi, jbflo, jbfhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb, osym )
C     $Id$
      implicit none
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      double precision q4
      logical sym_shell_quartet, osym, odoit
      external sym_shell_quartet
c     
      integer ilen, jlen
      integer ksh, lsh, kbflo, kbfhi, lbflo, lbfhi, ltop
      integer klen, llen, kblen, lblen, bsize
      logical status


      q4 = 1.0d0                ! If not using symmetry
      ilen = ibfhi - ibflo + 1
      jlen = jbfhi - jbflo + 1
      kblen = kbhi - kblo + 1
      lblen = lbhi - lblo + 1
      bsize = kblen*lblen*ilen*jlen
      call dfill(bsize,0.d0,ssbb,1)
      do ksh=kshlo,kshhi
         status = bas_cn2bfr(basis,ksh,kbflo,kbfhi)
         klen = kbfhi - kbflo + 1
         ltop = lshhi
         if (kshlo.eq.lshlo) ltop = ksh
         do lsh=lshlo,ltop
            odoit = (schwarz_shell(ksh,lsh)*schw_ij) .ge. tol2e
            if (odoit .and. osym) then
               odoit = odoit.and.
     $         sym_shell_quartet(basis, ish, jsh, ksh, lsh, q4)
            endif
            if (odoit) then
               status = bas_cn2bfr(basis,lsh,lbflo,lbfhi)
               llen = lbfhi - lbflo + 1
               call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                       scrlen, iscr, erilen, eri )
               call eri2blk( ilen, jlen, klen, llen, eri,
     $              ssbb(lbflo,kbflo,jbflo,ibflo), lblen, kblen, q4 )
            endif
         enddo
      enddo
      if (kshlo.eq.lshlo) call eriblktr(ilen, jlen, kblen, ssbb )

c$$$      WRITE(6,881) IBFLO, JBFLO
c$$$ 881  FORMAT(/,'+++',2I5)
c$$$      CALL MOINTS_MATPRINT(LBLEN,KBLEN,SSBB(LBLO,KBLO,JBFLO,IBFLO))
      return
      end
c
c
c
c
c
      subroutine eri2blk( ilen, jlen, klen, llen, eri, blk,
     $                    lblen, kblen, q4)
      implicit none
      integer ilen,jlen,klen,llen,lblen,kblen
      double precision blk(lblen,kblen,jlen,ilen)
      double precision eri(llen,klen,jlen,ilen)
      double precision q4
      integer k,l,i,j

      do i=1,ilen
        do j=1,jlen
          do k=1,klen
            do l=1,llen
              blk(l,k,j,i) = eri(l,k,j,i)*q4
            enddo
          enddo
        enddo
      enddo

      return
      end
#endif









c
c   ======================================================================
c
c                     NEW REORDERED GROUPED SHELL VERSIONS
c
c   ======================================================================
c
c
c
c  General interface to block AO integral routine
c  Branch on whether to use 
c          - labels/shell blocking
c  or
c          - quartet at a time
c  
c
      subroutine moints_gblk( basis, ish, jsh, kshlo, kshhi,
     $                        lshlo, lshhi, shmap, rbfmap,
     $                        schw_ij, tol2e, osym, oblk,
     $                        erilen, eri, scrlen, iscr,
     $                        ibflo, ibfhi, jbflo, jbfhi,
     $                        kblo, kbhi, lblo, lbhi,
     $                        ssbb )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer shmap(*), rbfmap(*)
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      logical osym, oblk
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
c
      integer mxshq, imem
      integer k_tmp, l_tmp
      integer k_iv, k_jv, k_kv, k_lv
      integer k_il, k_jl, k_kl, k_ll
      integer k_qv, l_qv
      logical status
      integer blkid
c
      logical moints_gblk_fromdisk, moints_gblk_todisk
      external moints_gblk_fromdisk, moints_gblk_todisk
c
c  Check to see if we can retrieve from disk
c
      blkid = kshlo*1000 + lshlo
      if (moints_gblk_fromdisk( blkid, ish, jsh, kshlo, lshlo, 
     $                          ibflo, ibfhi, jbflo, jbfhi, 
     $                          kblo, kbhi, lblo, lbhi,
     $                          ssbb )) then
        return
      endif

c
c  Otherwise compute directly
c
      if (oblk) then
        mxshq = (kshhi - kshlo + 1)*(lshhi - lshlo + 1)
        imem = 4*mxshq + 4*erilen
        if (.not. ma_push_get(MT_INT, imem, 'gblk tmp',l_tmp, k_tmp))
     $     call errquit('moints_gblk: no memory for labels',imem,
     &       MA_ERR)
        k_iv = k_tmp
        k_jv = k_iv + mxshq
        k_kv = k_iv + 2*mxshq
        k_lv = k_iv + 3*mxshq
        k_il = k_iv + 4*mxshq
        k_jl = k_il + erilen
        k_kl = k_il + 2*erilen
        k_ll = k_il + 3*erilen
        if (.not.ma_push_get(MT_DBL, mxshq, 'gblk q',l_qv, k_qv))
     $     call errquit('moints_gblk: no memory for sym factors',imem,
     &       MA_ERR)
        call moints_gblk_mq( basis, ish, jsh, kshlo, kshhi,
     $                       lshlo, lshhi, shmap, rbfmap,
     $                       schw_ij, tol2e,
     $                       erilen, eri, scrlen, iscr,
     $                       ibflo, ibfhi, jbflo, jbfhi,
     $                       kblo, kbhi, lblo, lbhi,
     $                       ssbb, osym, mxshq, dbl_mb(k_qv),
     $                       int_mb(k_iv), int_mb(k_jv),
     $                       int_mb(k_kv), int_mb(k_lv),
     $                       int_mb(k_il), int_mb(k_jl),
     $                       int_mb(k_kl), int_mb(k_ll) )
        if (.not.ma_pop_stack(l_qv))
     $    call errquit('moints_gblk: failed to pop',0, MA_ERR)
        if (.not.ma_pop_stack(l_tmp))
     $    call errquit('moints_gblk: failed to pop',0, MA_ERR)
      else
        call moints_gblk_sq( basis, ish, jsh, kshlo, kshhi,
     $                       lshlo, lshhi, shmap, rbfmap,
     $                       schw_ij, tol2e, osym,
     $                       erilen, eri, scrlen, iscr,
     $                       ibflo, ibfhi, jbflo, jbfhi,
     $                       kblo, kbhi, lblo, lbhi,
     $                       ssbb )
      endif
c
c  Cache this block to disk if required
c
      status = moints_gblk_todisk( blkid, ish, jsh, 
     $                             kshlo, lshlo,
     $                             ibflo, ibfhi, jbflo, jbfhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             ssbb )
c
c
c
      return
      end

        


      subroutine moints_gblk_mq( basis, ish, jsh, kshlo, kshhi,
     $                           lshlo, lshhi, shmap, rbfmap,
     $                           schw_ij, tol2e,
     $                           erilen, eri, scrlen, iscr,
     $                           ibflo, ibfhi, jbflo, jbfhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb, osym, mxshq, qv4, 
     $                           ishv, jshv, kshv, lshv,
     $                           ilab, jlab, klab, llab )
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer shmap(*), rbfmap(*)
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      logical osym
      integer mxshq
      double precision qv4(mxshq)
      integer ishv(mxshq), jshv(mxshq), lshv(mxshq), kshv(mxshq)
      integer ilab(erilen), jlab(erilen), klab(erilen), llab(erilen)
c
c
c
      double precision beffect, q4
      integer nshq,  neri, neritotal, ieri, nblksiz
      integer ii, jj, kk, ll, ksh, lsh
      integer ilen, jlen, klen, ltop
      logical more, odoit
      logical sym_shell_quartet
      external sym_shell_quartet
c
      logical intb_2e4c, intb_init4c
      external intb_2e4c, intb_init4c
c
      neritotal = 0
      nblksiz = (lbhi-lblo+1)*(kbhi-kblo+1)*
     $          (jbfhi-jbflo+1)*(ibfhi-ibflo+1)
      call dfill(nblksiz, 0.d0, ssbb, 1 )
c
c  Fill out index arrays with interacting quartet labels
c        
      q4 = 1.d0
      nshq = 0
      do kk=kshlo,kshhi
        ksh = shmap(kk)
        ltop = lshhi
        if (kshlo.eq.lshlo) ltop = kk
        do ll=lshlo,ltop
          lsh = shmap(ll)
          if (schwarz_shell(ksh,lsh)*schw_ij.ge.tol2e) then
            odoit = .true.
            if (osym) odoit = sym_shell_quartet(basis, ish, jsh,
     $                                          max(ksh,lsh), 
     $                                          min(ksh,lsh), q4 )
            if (odoit) then
              nshq = nshq + 1
              kshv(nshq) = ksh
              lshv(nshq) = lsh
              qv4(nshq)  = q4
            endif
          endif
        enddo
      enddo
      call ifill( nshq, ish, ishv, 1 )
      call ifill( nshq, jsh, jshv, 1 )
c
c  Prepare texas for shell block
c
      if (.not. intb_init4c( basis, ishv, jshv, basis, kshv, lshv,
     $                       nshq, qv4, osym, scrlen, iscr, 
     $                       erilen, beffect ))
     $   call errquit('moints_gblk: intb_init4c failed',0, INT_ERR)
c
c  Loop through batches of integrals until exhausted
c
 100  more = intb_2e4c( basis, ishv, jshv, basis, kshv, lshv,
     $                  nshq, qv4, osym, TOL2E, .false.,
     $                  ilab, jlab, klab, llab,
     $                  eri, erilen, neri, scrlen, iscr )
c
c  Unpack labels and assign to integrals
c
      if (neri.gt.0) then
        neritotal = neritotal + neri
        do ieri=1,neri
          ii = ilab(ieri)
          jj = jlab(ieri)
          kk = rbfmap(klab(ieri))
          ll = rbfmap(llab(ieri))
          ssbb(ll,kk,jj,ii) = eri(ieri)
*          write(6,7) ilab(ieri),jlab(ieri),klab(ieri),llab(ieri),
*     $         eri(ieri)
* 7        format(4i5,f12.6)
        enddo
      endif
      if (more) goto 100
c
c  Reflect diagonal shell block
c
      klen = kbhi  - kblo + 1
      ilen = ibfhi - ibflo + 1
      jlen = jbfhi - jbflo + 1
      if (kshlo.eq.lshlo) call eriblktr(ilen, jlen, klen, ssbb )
c
c
c
c$$$      LL = (LBHI - LBLO + 1)
c$$$      KK = (KBHI - KBLO + 1)
c$$$      WRITE(6,881) IBFLO, JBFLO
c$$$ 881  FORMAT(/,'===',2I5)
c$$$      CALL MOINTS_MATPRINT(LL,KK,SSBB(LBLO,KBLO,JBFLO,IBFLO))
      return
      end







c
c  This computes one shell quartet at a time
c
      subroutine moints_gblk_sq( basis, ish, jsh, kshlo, kshhi,
     $                           lshlo, lshhi, shmap, rbfmap,
     $                           schw_ij, tol2e, osym,
     $                           erilen, eri, scrlen, iscr, 
     $                           ibflo, ibfhi, jbflo, jbfhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           ssbb )
C     $Id$
      implicit none
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer shmap(*), rbfmap(*)
      logical osym
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
c     
      integer ilen, jlen, kblen, lblen
      integer ksh, lsh, kbflo, kbfhi, lbflo, lbfhi, ltop
      integer kk, ll, bsize
      logical status, odoit
      double precision q4
      logical sym_shell_quartet
      external sym_shell_quartet
c
      q4 = 1.0d0
      ilen = ibfhi - ibflo + 1
      jlen = jbfhi - jbflo + 1
      kblen = kbhi - kblo + 1
      lblen = lbhi - lblo + 1
      bsize = kblen*lblen*ilen*jlen
      call dfill(bsize,0.d0,ssbb,1)
      do kk=kshlo,kshhi
        ksh = shmap(kk)
        status = bas_cn2bfr(basis, ksh, kbflo, kbfhi)
        ltop = lshhi
        if (kshlo.eq.lshlo) ltop = kk
        do ll=lshlo,ltop
          lsh = shmap(ll)
          if ((schwarz_shell(ksh,lsh)*schw_ij).ge.tol2e) then
            odoit = .true.
            if (osym) odoit = sym_shell_quartet(basis, ish, jsh,
     $                                         max(ksh,lsh),
     $                                         min(ksh,lsh), q4 )
            if (odoit) then
              status = bas_cn2bfr(basis, lsh, lbflo, lbfhi)
              call int_2e4c( basis, ish, jsh, basis, ksh, lsh,
     $                       scrlen, iscr, erilen, eri )            
              call eri2blkx(rbfmap, ilen, jlen,
     $                      kblo,  kbhi,  lblo,  lbhi,
     $                      kbflo, kbfhi, lbflo, lbfhi,
     $                      eri, q4, ssbb )
            endif
          endif
        enddo
      enddo
      if (kshlo.eq.lshlo) call eriblktr(ilen, jlen, kblen, ssbb )

c$$$      WRITE(6,881) ISH, JSH
c$$$ 881  FORMAT(/,'+++',2I5)
c$$$      CALL MOINTS_MATPRINT(LBLEN,KBLEN,SSBB(LBLO,KBLO,JBFLO,IBFLO))
      

      return
      end

      subroutine moints_trf1_new( nbf, qlo, qhi, mo1_lo, mo1_hi,
     $     ilen, jlen, klo, khi, llo, lhi,
     $     ssbb, ct, ssni, hlp )
      implicit none
      integer nbf, qlo, qhi, mo1_lo, mo1_hi
      integer ilen, jlen, klo, khi, llo, lhi
      double precision ssbb(llo:lhi,klo:khi,jlen,ilen)
      double precision ct(qlo:qhi,nbf)
      double precision ssni(nbf,jlen,ilen,mo1_lo:mo1_hi)
      double precision hlp(mo1_lo:mo1_hi,*)
c     
      double precision s
      double precision tol
      data tol/1.d-12/
c     
      integer i, j, k, l, o
      integer kk,ll,mm
      double precision dabsmax
      external dabsmax
c     
      kk=khi-klo+1
      ll=lhi-llo+1
      mm=mo1_hi-mo1_lo+1
      do i = 1, ilen
         do j = 1, jlen

            if (klo.ne.llo) then
               call dgemm('n','n',mm,kk,ll,
     A              1d0,ct(mo1_lo,llo),qhi-qlo+1,
     B              ssbb(llo,klo,j,i),ll,
     C              0d0,hlp(mo1_lo,1),mm)
cc     C              0d0,hlp(mo1_lo,klo),mm)
               do o = mo1_lo, mo1_hi
                  do k = klo, khi
                     ssni(k,j,i,o) = ssni(k,j,i,o) + hlp(o,k-klo+1)
c                     ssni(k,j,i,o) = ssni(k,j,i,o) + hlp(o,k)
                  enddo
               enddo
            endif
c     
            call dgemm('n','t',mm,ll,kk,
     A           1d0,ct(mo1_lo,klo),qhi-qlo+1,
     B           ssbb(llo,klo,j,i),ll,
     C           0d0,hlp(mo1_lo,1),mm)
c     C           0d0,hlp(mo1_lo,llo),mm)
            do o = mo1_lo, mo1_hi
               do l = llo, lhi
                  ssni(l,j,i,o) = ssni(l,j,i,o) + hlp(o,l-llo+1)
c                  ssni(l,j,i,o) = ssni(l,j,i,o) + hlp(o,l)
               enddo
            enddo
         enddo
      enddo
c     
      
c$$$  nmo1 = mo1_hi - mo1_lo + 1
c$$$  llen = lhi - llo + 1
c$$$  klen = khi - klo + 1
c$$$  kjilen = klen*jlen*ilen
c$$$  ljilen = llen*jlen*ilen
c$$$  
c$$$  call dgemm( 't', 'n', kjilen, nmo1, llen, scale, ssbb, llen,
c$$$  $            c(llo,mo1_lo), nbf, 0.d0, hlp, kjilen )
c$$$  call moints_1idxpush( nbf, ilen, jlen, klo, khi, nmo1,
c$$$  $                      hlp, ssni )
c$$$  call dgemm( 't', 'n', ljilen, nmo1, klen, scale, ssbbt, klen,
c$$$  $            c(klo,mo1_lo), nbf, 0.d0, hlp, ljilen )
c$$$  call moints_1idxpush( nbf, ilen, jlen, llo, lhi, nmo1,
c$$$  $                      hlp, ssni )

      return
      end












      
      subroutine moints_trf1( nbf, qlo, qhi, mo1_lo, mo1_hi,
     $     ilen, jlen, klo, khi, llo, lhi,
     $     scale, ssbb, ssbbt, c, ssni, hlp )
      implicit none
      integer nbf, qlo, qhi, mo1_lo, mo1_hi
      integer ilen, jlen, klo, khi, llo, lhi
      double precision scale
      double precision ssbb(llo:lhi,klo:khi,jlen,ilen)
      double precision ssbbt(klo:khi,llo:lhi,jlen,ilen)
      double precision c(nbf,qlo:qhi)
      double precision ssni(nbf,jlen,ilen,mo1_lo:mo1_hi)
      double precision hlp(*)
c     
      integer nmo1, llen, klen, kjilen, ljilen
c     
      nmo1 = mo1_hi - mo1_lo + 1
      llen = lhi - llo + 1
      klen = khi - klo + 1
      kjilen = klen*jlen*ilen
      ljilen = llen*jlen*ilen
c     
      call dgemm( 't', 'n', kjilen, nmo1, llen, scale, ssbb, llen,
     $            c(llo,mo1_lo), nbf, 0.d0, hlp, kjilen )
      call moints_1idxpush( nbf, ilen, jlen, klo, khi, nmo1,
     $                      hlp, ssni )
      call dgemm( 't', 'n', ljilen, nmo1, klen, scale, ssbbt, klen,
     $            c(klo,mo1_lo), nbf, 0.d0, hlp, ljilen )
      call moints_1idxpush( nbf, ilen, jlen, llo, lhi, nmo1,
     $                      hlp, ssni )
c   
      return
      end











c
c  Pack shell quartet ERI into bigger blocks
c  with remapped K,L indices of a contigious range
c
      subroutine eri2blkx( rbfmap, ilen, jlen,
     $                     kblo, kbhi, lblo, lbhi,
     $                     kbflo, kbfhi, lbflo, lbfhi, 
     $                     eri, q4, blk )
      implicit none
      integer rbfmap(*)                                          ! [input] basis map orig -> new
      integer ilen, jlen                                         ! [input] ij shell lengths
      integer kblo, kbhi                                         ! [input] k range in new order
      integer lblo, lbhi                                         ! [input] l range in new order
      integer kbflo, kbfhi                                       ! [input] original basis indices
      integer lbflo, lbfhi                                       ! [input] original basis indices
      double precision q4                                        ! [input] symmetry scale
      double precision eri(lbflo:lbfhi,kbflo:kbfhi,jlen,ilen)    ! [input] eri shell quartet
      double precision blk(lblo:lbhi,kblo:kbhi,jlen,ilen)        ! [output] eri block
      integer k,l,i,j,kk,ll


      do i=1,ilen
        do j=1,jlen
          do k=kbflo,kbfhi
            kk = rbfmap(k)
            do l=lbflo,lbfhi
              ll = rbfmap(l)
              blk(ll,kk,j,i) = eri(l,k,j,i)*q4
            enddo
          enddo
        enddo
      enddo
      
      return
      end










      subroutine eriblktr(ilen, jlen, kblen, ssbb )
      implicit none
      integer ilen,jlen,kblen
      double precision ssbb(kblen,kblen,jlen,ilen)
      integer i,j,k,l

      do i=1,ilen
        do j=1,jlen
          do k=1,kblen
            do l=1,k-1
              ssbb(k,l,j,i) = ssbb(l,k,j,i)
            enddo
          enddo
        enddo
      enddo
      return
      end










      subroutine moints_1idxpush( nbf, ilen, jlen, klo, khi, nmo1,
     $                            x, ssni )
      implicit none
      integer nbf, ilen, jlen, klo, khi, nmo1
      double precision x(klo:khi,jlen,ilen,nmo1)
      double precision ssni(nbf,jlen,ilen,nmo1)
      integer a,i,j,k
      
      do a=1,nmo1
        do i=1,ilen
          do j=1,jlen
            do k=klo,khi
              ssni(k,j,i,a) = ssni(k,j,i,a) + x(k,j,i,a)
            enddo
          enddo
        enddo
      enddo
      return
      end













      subroutine moints_blktr( ilen, jlen, klen, llen,
     $                         ssbb, ssbbt )
      implicit none
      integer ilen, jlen, klen, llen
      double precision ssbb(llen,klen,jlen,ilen)
      double precision ssbbt(klen,llen,jlen,ilen)
      integer i,j,k,l

      do i=1,ilen
        do j=1,jlen
          do l=1,llen
            do k=1,klen
              ssbbt(k,l,j,i) = ssbb(l,k,j,i)
            enddo
          enddo
        enddo
      enddo
      return
      end








c
c Coulomb 2nd index transform
c
      subroutine moints_trf2J( nbf, qlo, qhi, ostart, olo, ohi, 
     $                         ilo, ihi, jlo, jhi, ssni, h1, h2, 
     $                         c, g_coul )
      implicit none
#include "global.fh"
      integer nbf, qlo, qhi
      integer ostart, olo, ohi, ilo, ihi, jlo, jhi
      double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
      double precision h1(nbf,ilo:ihi)
      double precision h2(jlo:jhi,ilo:ihi)
      double precision c(nbf,qlo:qhi)
      integer g_coul
c
      integer nni, ijlo, ijhi, ilen, jlen, ijlen
      integer ab, aa, bb, a, b, i, aoff, ofroz
      integer olb,ain
c
      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      nni = ilen*nbf
      ijlen = ilen*jlen
#ifdef BLOCK_TRANSF
      ijlo = (ilo-1)*nbf + 1
      ijhi = ihi*nbf
#endif
cedo       do a=olo,ohi
cstaggering to avoid congestion in ga_acc
      olb=ohi-olo+1
      do ain=ga_nodeid()+1,ga_nodeid()+olb
        a=mod(ain,olb)+olo
        do b=ostart,a
          call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,a),
     $               nbf,c(1,b),nbf,0.d0,h2,ijlen)
#ifndef NOCOMMS
          aa = a - ofroz
          bb = b - ofroz
          ab = (aa*(aa-1))/2 + bb - aoff
#ifdef BLOCK_TRANSF
          call dfill(nni,0.d0,h1,1)
          do i=ilo,ihi
            jtop = jhi
            if (jhi.eq.ihi) jtop = i
            do j=jlo,jtop
              h1(j,i) = h2(j,i)
            enddo
          enddo
          call ga_acc(g_coul,ijlo,ijhi,ab,ab,h1,nni,1.d0)
#else
          do i=ilo,ihi
            ijlo = (i-1)*nbf + jlo
            ijhi = (i-1)*nbf + jhi
            if (jhi.eq.ihi) ijhi = (i-1)*nbf + i
            call ga_acc(g_coul,ijlo,ijhi,ab,ab,h2(jlo,i),1,1.d0)
          enddo
#endif
#endif
        enddo
      enddo
      return
      end







      subroutine moints_sym_zero(nbf, vlo, vhi, syms, i, j, x)
      implicit none
      integer nbf
      integer vlo, vhi
      integer syms(nbf)
      integer i, j
      double precision x(vlo:vhi, vlo:vhi)
c
      integer symij
      integer symijp
      integer p, q
#include "bitops.fh"
c
C$$$      WRITE(6,*) (SYMS(P),P=1,NBF)
C$$$      WRITE(6,*) ' MOINTS_SYM_ZERO ', I, J
c
      symij = ieor(syms(i),syms(j))
      do p = vlo, vhi
         symijp = ieor(symij,syms(p))
         do q = vlo, vhi
            if (ieor(symijp,syms(q)) .ne. 0) then
               x(q,p) = 0.0d0
               x(p,q) = 0.0d0
            endif
         enddo
      enddo
c
      end








      subroutine moints_Ktrf34( nbf, qlo, qhi, ostart, olo, ohi, 
     $                          vlo, vhi, otrp, c, tmp, osym, syms, 
     $                          g_exch )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nbf, qlo, qhi
      integer ostart, olo, ohi, vlo, vhi
      logical otrp, osym
      double precision c(nbf,qlo:qhi)
      double precision tmp(nbf,qlo:qhi)
      integer syms(nbf)
      integer g_exch
c
      integer nvir, clo, chi, rlo, rhi, my_id
      integer ofroz, a1, a2, aa1, aa2, aa, aoff, k_local, ld
#ifdef BAD_GACCESS
      integer l_local
#endif

      ofroz = ostart - 1
      nvir = vhi - vlo + 1
      my_id = ga_nodeid()
      call ga_distribution(g_exch, my_id, rlo, rhi, clo, chi )
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      do a1=olo,ohi
        aa1 = a1 - ofroz
        do a2=ostart,a1
          aa2 = a2 - ofroz
          aa = (aa1*(aa1-1))/2 + aa2 - aoff
          if ((aa.ge.clo).and.(aa.le.chi)) then
#ifdef BAD_GACCESS
      ld=rhi-rlo+1
      if(.not.ma_push_get(MT_DBL,ld,
     $  'scratch buff2x', l_local, k_local)) call
     $  errquit('moints2x: pushget failed',0,0)
            call ga_get(g_exch,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
#else
            call ga_access(g_exch,rlo,rhi,aa,aa,k_local,ld)
#endif
            call moints_Ktrf34a( otrp, nbf, qlo, qhi, vlo, vhi, c,
     $                           dbl_mb(k_local), tmp )
            if (osym) call moints_sym_zero(nbf, vlo, vhi, syms,
     $                                     a1, a2, dbl_mb(k_local))
#ifdef BAD_GACCESS
        call ga_put(g_exch,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('moints2x: popstack failed',0,0)
#else
            call ga_release(g_exch,rlo,rhi,aa,aa)
#endif
          endif
        enddo
      enddo
      return
      end



c
c Note: May require transpose the half-transformed matrix
c       to conform to convention for storage
c
      subroutine moints_Ktrf34a(otrp, nbf, qlo, qhi, 
     $                          vlo, vhi, c, x, tmp)
      implicit none
      logical otrp
      integer nbf, qlo, qhi, vlo, vhi
      double precision c(nbf,qlo:qhi)
      double precision x(vlo:vhi,vlo:vhi)
      double precision tmp(vlo:vhi,nbf)
      integer nv
      character*1 ttt
      
      ttt = 'n'
      if (otrp) ttt = 't'
      nv = vhi - vlo + 1
CRK temporary bug fix
CRK   call dgemm('t',ttt,nv,nbf,nbf,1.d0,c(1,vlo),nbf,x,nbf,
CRK  $           0.d0,tmp,nv)
CRK   call dgemm('n','n',nv,nv,nbf,1.d0,tmp,nv,c(1,vlo),nbf,
CRK  $           0.d0,x,nv)
      call dgemm('t',ttt,nv,nbf,nbf,1.d0,c(1,vlo),nbf,x,nbf,
     $           0.d0,tmp,nv)
      call dgemm('n','n',nv,nv,nbf,1.d0,tmp,nv,c(1,vlo),nbf,
     $           0.d0,x,nv)

      return
      end









c
c
c
      subroutine moints_Jtrf34( nbf, qlo, qhi, ostart, olo, ohi, 
     $                          vlo, vhi, c, tmp, osym, syms, g_coul )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nbf, qlo, qhi
      integer ostart, olo, ohi, vlo, vhi
      logical osym
      integer syms(nbf)
      double precision c(nbf,qlo:qhi)
      double precision tmp(nbf,qlo:qhi)
      integer g_coul
c
      integer clo, chi, rlo, rhi, my_id
      integer ofroz, aoff, a1, a2, aa1, aa2, aa, k_local, ld
#ifdef BAD_GACCESS
      integer l_local
#endif
c
      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      my_id = ga_nodeid()
      call ga_distribution(g_coul, my_id, rlo, rhi, clo, chi )
      do a1=olo,ohi
        aa1 = a1 - ofroz
        do a2=ostart,a1
          aa2 = a2 - ofroz
          aa = (aa1*(aa1-1))/2 + aa2 - aoff
          if ((aa.ge.clo).and.(aa.le.chi)) then
#ifdef BAD_GACCESS
      ld=rhi-rlo+1
      if(.not.ma_push_get(MT_DBL,ld,
     $  'scratch buff2x', l_local, k_local)) call
     $  errquit('moints2x: pushget failed',0,0)
            call ga_get(g_coul,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
#else
            call ga_access(g_coul,rlo,rhi,aa,aa,k_local,ld)
#endif
            call moints_Jtrf34a( nbf, qlo, qhi, vlo, vhi, c,
     $                           dbl_mb(k_local), tmp )
            if (osym) call moints_sym_zero(nbf, vlo, vhi, syms,
     $                                     a1, a2, dbl_mb(k_local))
#ifdef BAD_GACCESS
        call ga_put(g_coul,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('moints2x: popstack failed',0,0)
#else
            call ga_release(g_coul,rlo,rhi,aa,aa)
#endif
          endif
        enddo
      enddo
      return
      end





      subroutine moints_Jtrf34a( nbf, qlo, qhi, vlo, vhi, c, x, tmp )
      implicit none
      integer nbf, qlo, qhi, vlo, vhi
      double precision c(nbf,qlo:qhi)
      double precision x(nbf,nbf)
      double precision tmp(nbf,nbf)
      integer nv

      nv = vhi - vlo + 1
      call upper2square( nbf, x, tmp )
      call dgemm('t', 'n', nv, nbf, nbf, 1.d0, c(1,vlo), nbf, tmp, 
     $           nbf, 0.d0, x, nv)
      call dgemm('n', 'n', nv, nv, nbf, 1.d0, x, nv, c(1,vlo), nbf,
     $           0.d0, tmp, nv)
      call dcopy((nv*nv),tmp,1,x,1)

      return
      end







      subroutine print_shell_block(basis, eri, ish, jsh, ksh, lsh)
      implicit none
#include "bas.fh"
      integer basis
      double precision eri(*)
      integer ish, jsh, ksh, lsh
      logical status
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
c
      status = bas_cn2bfr(basis,ish,ilo,ihi)
      status = bas_cn2bfr(basis,jsh,jlo,jhi)
      status = bas_cn2bfr(basis,ksh,klo,khi)
      status = bas_cn2bfr(basis,lsh,llo,lhi)
c
      call psb2(eri,ilo,ihi,jlo,jhi,klo,khi,llo,lhi)
c
      end
      subroutine psb2(eri,ilo,ihi,jlo,jhi,klo,khi,llo,lhi)
      implicit none
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
      integer i,j,k,l
c
      do i = ilo,ihi
         do j= jlo,jhi
            do k=klo,khi
               do l=llo,lhi
                  write(6,1) i,j,k,l,eri(l,k,j,i)
 1                format(1x,4i5,1p,d20.12)
               enddo
            enddo
         enddo
      enddo
c
      end
      





      subroutine upper2square( n, a, b )
      integer n
      double precision a(n,n)
      double precision b(n,n)
      integer i,j

      do i=1,n
        do j=1,i
          b(j,i) = a(j,i)
          b(i,j) = a(j,i)
        enddo
      enddo
      return
      end









       logical function ga_check_JKblocked(g_a,n1,n2,jlo,jhi)
C$Id$
       implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
       integer g_a, n1, n2, jlo, jhi
       integer n1t, n2n2
       integer gtype, dim1, dim2
       integer ilo, ihi, my_id, jstat
       logical status

       n1t = (n1*(n1+1))/2
       n2n2 = n2*n2
       my_id = ga_nodeid()
       call ga_inquire(g_a,gtype,dim1,dim2)
       status = (dim1.eq.n2n2).and.(dim2.eq.n1t)
       call ga_distribution(g_a,my_id,ilo,ihi,jlo,jhi)
       if ((jlo.ne.0).and.(jhi.ne.-1)) then
         status = status.and.((ilo.eq.1).and.(ihi.eq.n2n2).and.
     $                        (jlo.ge.1).and.(jhi.le.n1t))
       endif

       call ga_sync()
       jstat = 0
       if (status) jstat = 1
! igop('*') detects the presence of one or more stat=0 inputs
#if NWCHEM_USE_IGOP_PROD
       call ga_igop(Msg_ChkJK,jstat,1,'*')
#else
       call ga_igop(Msg_ChkJK,jstat,1,'min')
#endif
#ifdef DEBUG       
       write(6,901) my_id,jstat,dim1,dim2,n2n2,n1t,
     $      ilo,ihi,n2n2,jlo,jhi,n1t,status
 901   format('ga_check_moblocked: ID:',i3,5x,i3,'==',i3,'?',5x,
     $        'dim:  [',i3,',',i3,']',2x,'(',i3,',',i3,')',5x,
     $        'col_block:',i3,'-',i3,1x,'(',i3,')',5x,
     $        'row_block:',i3,'-',i3,1x,'(',i3,')',5x,
     $        'status:',l1)
#endif

       ga_check_JKblocked = (jstat.ne.0)
       return
       end




c$$$       integer function ga_create_Jblocked(noper,nbf,label)
c$$$C$Id$
c$$$       implicit none
c$$$#include "errquit.fh"
c$$$#include "global.fh"
c$$$#include "mafdecls.fh"
c$$$       integer nbf,noper
c$$$       character*(*) label
c$$$       integer g_a
c$$$       integer row_dist, col_dist
c$$$
c$$$       col_dist = 1
c$$$       row_dist = nbf*nbf
c$$$       if (.not.ga_create(MT_DBL,(nbf*nbf),noper,label,
c$$$     $                    row_dist,col_dist,g_a))
c$$$     $      call errquit('ga_create_Jblocked: cannot allocate',0)
c$$$
c$$$       ga_create_Jblocked = g_a
c$$$
c$$$       return
c$$$       end
c$$$





c$$$       integer function ga_create_Kblocked(noper,nbf,nvir,label)
c$$$       implicit none
c$$$#include "errquit.fh"
c$$$#include "global.fh"
c$$$#include "mafdecls.fh"
c$$$       integer nbf,noper,nvir
c$$$       character*(*) label
c$$$       integer g_a
c$$$       integer row_dist, col_dist
c$$$
c$$$       col_dist = 1
c$$$       row_dist = nbf*nvir
c$$$       if (.not.ga_create(MT_DBL,(nbf*nvir),noper,label,
c$$$     $                    row_dist,col_dist,g_a))
c$$$     $      call errquit('ga_create_Kblocked: cannot allocate',0)
c$$$
c$$$       ga_create_Kblocked = g_a
c$$$
c$$$       return
c$$$       end






      integer function ga_create_JKblocked(noper,n1,n2,label)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer noper,n1,n2
      character*(*) label
      integer g_a
      integer row_dist, col_dist

      col_dist = 1
      row_dist = n1*n2
      if (.not.ga_create(MT_DBL,(n1*n2),noper,label,
     $                    row_dist,col_dist,g_a)) then
         write(6,*) ' ga_create_JKblocked: ', label
         call util_flush(6)
         call errquit('ga_create_JKblocked: cannot allocate',0,
     &       GA_ERR)
      endif
      ga_create_JKblocked = g_a
      
      return
      end









      subroutine row_exch( n, m, map, x, y )
      integer n
      integer m
      integer map(n)
      double precision x(n,m)
      double precision y(n,m)

      do i=1,m
        do j=1,n
          y(map(j),i) = x(j,i)
        enddo
      enddo
      return
      end






c
c  Required by moints_trp ... will remove later
c
      subroutine moints2xv_Ktrf4a( vlo, vhi, nbf, c, x, tmp )
      implicit none
      integer vlo, vhi, nbf
      double precision c(nbf,nbf)
      double precision x(vlo:vhi,nbf)
      double precision tmp(vlo:vhi,vlo:vhi)
      integer nvir
      
      nvir = vhi - vlo + 1
      call dgemm('t', 't', nvir, nvir, nbf, 1.d0, c(1,vlo), nbf,
     $            x, nvir, 0.d0, tmp, nvir)
      return
      end
 
 






      subroutine summaxmin( arg, ysum, ymax, ymin )
      implicit none
      double precision arg, ysum, ymax, ymin
      ysum = ysum + arg
      ymax = max(ymax,arg)
      ymin = min(ymin,arg)
      end

