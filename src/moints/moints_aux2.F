      subroutine moints_gblk( basis, ish, jsh, kshlo, kshhi,
     $     lshlo, lshhi, schw_ij, tol2e,
     $     max2e, eri, mem2, iscr, ilen,
     $     jlen, kblo, kbhi, lblo, lbhi, ssbb, osym )
C     $Id: moints_aux2.F,v 1.7 1995-10-26 06:26:22 rg240 Exp $
      implicit none
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer max2e, mem2
      double precision schw_ij, tol2e, eri(*), iscr(*)
      integer ilen, jlen, kblo, kbhi, lblo, lbhi
      double precision ssbb(lblo:lbhi,kblo:kbhi,jlen,ilen)
      double precision q4
      logical sym_shell_quartet, osym, odoit
      external sym_shell_quartet
c     
      integer ksh, lsh, kbflo, kbfhi, lbflo, lbfhi, ltop
      integer klen, llen, kblen, lblen, bsize
      logical status

      q4 = 1.0d0                ! If not using symmetry
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
*               write(6,1) ish, jsh, ksh, lsh, odoit, q4
* 1             format(1x,4i5,2x,l,2x,f6.2)
            endif
            if (odoit) then
               status = bas_cn2bfr(basis,lsh,lbflo,lbfhi)
               llen = lbfhi - lbflo + 1
               call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                       mem2, iscr, max2e, eri )
c     
*               call print_shell_block(basis,eri, ish, jsh, ksh, lsh)
c     
               call eri2blk( ilen, jlen, klen, llen, eri,
     $              ssbb(lbflo,kbflo,1,1), lblen, kblen, q4 )
            endif
         enddo
      enddo
      if (kshlo.eq.lshlo) call eriblktr(ilen, jlen, kblen, ssbb )
      return
      end





      
      subroutine moints_trf1( nbf, molo, mohi, ilen, jlen,
     $                          klo, khi, llo, lhi, scale,
     $                          ssbb, ssbbt, c, ssni, hlp )
      implicit none
      integer nbf, molo, mohi
      integer ilen, jlen, klo, khi, llo, lhi
      double precision scale
      double precision ssbb(llo:lhi,klo:khi,jlen,ilen)
      double precision ssbbt(klo:khi,llo:lhi,jlen,ilen)
      double precision c(nbf,nbf)
      double precision ssni(nbf,jlen,ilen,molo:mohi)
      double precision hlp(*)
c
      integer nmo, llen, klen, kjilen, ljilen
c
c
c      
      nmo = mohi - molo + 1
      llen = lhi - llo + 1
      klen = khi - klo + 1
      kjilen = klen*jlen*ilen
      ljilen = llen*jlen*ilen
      call dgemm( 't', 'n', kjilen, nmo, llen, scale, ssbb, llen,
     $            c(llo,molo), nbf, 0.d0, hlp, kjilen )
      call moints_1idxpush( nbf, ilen, jlen, klo, khi, nmo,
     $                      hlp, ssni )
      call dgemm( 't', 'n', ljilen, nmo, klen, scale, ssbbt, klen,
     $            c(klo,molo), nbf, 0.d0, hlp, ljilen )
      call moints_1idxpush( nbf, ilen, jlen, llo, lhi, nmo,
     $                      hlp, ssni )

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






      subroutine moints_1idxpush( nbf, ilen, jlen, klo, khi, nmo,
     $                            x, ssni )
      implicit none
      integer nbf, ilen, jlen, klo, khi, nmo
      double precision x(klo:khi,jlen,ilen,nmo)
      double precision ssni(nbf,jlen,ilen,nmo)
      integer a,i,j,k
      
      do a=1,nmo
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
      subroutine moints_trf2J( nbf, ostart, olo, ohi, ilo, ihi,
     $                         jlo, jhi, ssni, h1, h2, c, g_coul )
      implicit none
      integer nbf, ostart, olo, ohi, ilo, ihi, jlo, jhi
      double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
      double precision h1(nbf,ilo:ihi)
      double precision h2(jlo:jhi,ilo:ihi)
      double precision c(nbf,nbf)
      integer g_coul
c
      integer nni, ijlo, ijhi, ilen, jlen, ijlen
      integer ab, aa, bb, a, b, i, aoff, ofroz
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
      do a=olo,ohi
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








      subroutine moints_Ktrf34( g_exch, otrp, ostart, olo, ohi,
     $                          vlo, vhi, nbf, c, tmp, osym, syms )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_exch
      logical otrp
      integer nbf, ostart, olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)

      integer nvir, clo, chi, rlo, rhi, my_id
      integer ofroz, a1, a2, aa1, aa2, aa, aoff, k_local, ld

      logical osym
      integer syms(nbf)

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
            call ga_access(g_exch,rlo,rhi,aa,aa,k_local,ld)
            call moints_Ktrf34a( otrp, nbf, c, dbl_mb(k_local), tmp )
            if (osym) call moints_sym_zero(nbf, 1, nbf, syms, a1, a2,
     $                                     dbl_mb(k_local))
            call ga_release(g_exch,rlo,rhi,aa,aa)
          endif
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
      integer ieor
      intrinsic ieor
c
      write(6,*) ' moints_sym_zero ', i, j
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




c
c Note: May require transpose the half-transformed matrix
c       to conform to convention for storage
c
      subroutine moints_Ktrf34a( otrp, nbf, c, x, tmp )
      implicit none
      logical otrp
      integer nbf
      double precision c(nbf,nbf)
      double precision x(nbf,nbf)
      double precision tmp(nbf,nbf)
      character*1 ttt
      
      ttt = 'n'
      if (otrp) ttt = 't'
      call dgemm('t',ttt,nbf,nbf,nbf,1.d0,c,nbf,x,nbf,
     $           0.d0,tmp,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.d0,tmp,nbf,c,nbf,
     $           0.d0,x,nbf)

      return
      end









c
c
c
      subroutine moints_Jtrf34( g_coul, ostart, olo, ohi, vlo, vhi,
     $                          nbf, c, tmp, osym, syms )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_coul
      integer ostart, olo, ohi, vlo, vhi, nbf
      logical osym
      integer syms(nbf)
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)

      integer clo, chi, rlo, rhi, my_id
      integer ofroz, aoff, a1, a2, aa1, aa2, aa, k_local, ld

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
            call ga_access(g_coul,rlo,rhi,aa,aa,k_local,ld)
            call moints_Jtrf34a( nbf, c, dbl_mb(k_local), tmp )
            if (osym) call moints_sym_zero(nbf, 1, nbf, syms, a1, a2,
     $                                     dbl_mb(k_local))
            call ga_release(g_coul,rlo,rhi,aa,aa)
          endif
        enddo
      enddo
      return
      end





      subroutine moints_Jtrf34a( nbf, c, x, tmp )
      implicit none
      integer nbf
      double precision c(nbf,nbf)
      double precision x(nbf,nbf)
      double precision tmp(nbf,nbf)

      integer i, j

      do i=1,nbf
        do j=1,i
          tmp(j,i) = x(j,i)
          tmp(i,j) = x(j,i)
        enddo
      enddo

      call dgemm('t','n',nbf,nbf,nbf,1.d0,c,nbf,tmp,nbf,
     $           0.d0,x,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.d0,x,nbf,c,nbf,
     $           0.d0,tmp,nbf)

      do i=1,nbf
        do j=1,nbf
          x(j,i) = tmp(j,i)
        enddo
      enddo

      return
      end







      subroutine print_shell_block(basis, eri, ish, jsh, ksh, lsh)
      implicit none
#include "bas.fh"
      integer basis
      double precision eri
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
        do j=1,i-1
          b(j,i) = a(j,i)
          b(i,j) = a(j,i)
        enddo
      enddo
      return
      end

