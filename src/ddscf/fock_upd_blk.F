      subroutine fock_upd_blk( nfock, vfock, ilo, ihi,
     $     jlo, jhi, fac, buf, tmp )
c     
c     $Id$
c     
      implicit none
#include "global.fh"
c     
c     Accumulate Fock contribution of atom block
c     into array of handles
c     
      integer nfock, vfock(nfock)
      double precision fac(nfock)
      double precision buf(nfock,*)
      double precision tmp(*)
c     
      integer ilo, ihi, jlo, jhi
      integer ilen, jlen
      integer ii, i, j, ji
c     
      jlen = jhi - jlo + 1
      ilen = ihi - ilo + 1
c     
      if ((ilen.gt.0).and.(jlen.gt.0)) then
         if (nfock.gt.1) then
            do ii=1,nfock
               do j=1,jlen
                  ji = (j-1)*ilen
                  do i=1,ilen
                     tmp(ji+i) = buf(ii,ji+i)
                  enddo
               enddo
               call ga_acc( vfock(ii), ilo, ihi, jlo, jhi, tmp,
     $              ilen, fac(ii) )
            enddo
         else
            call ga_acc( vfock(1), ilo, ihi, jlo, jhi, buf,
     $           ilen, fac(1) )
         endif
      endif
c     
      end
      subroutine upd_atom_block(g_array, basis, iat, jat, buf)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
c
      integer g_array, basis, iat, jat
      double precision buf(*)
c
      integer ilo, ihi, jlo, jhi, idim, jdim
c
c     add atom block buf info of the matrix g_array (over basis functions)
c
      if (.not. bas_ce2bfr(basis, iat, ilo, ihi))
     $     call errquit('upd_atom_block: ce2bfr failed', 0, BASIS_ERR)
      if (.not. bas_ce2bfr(basis, jat, jlo, jhi))
     $     call errquit('upd_atom_block: ce2bfr failed', 0, BASIS_ERR)
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
c     
      if (idim.gt.0 .and. jdim.gt.0) then
         call ga_acc(g_array, ilo, ihi, jlo, jhi, buf, idim,
     $        1.0d0)
      endif
*      write(6,*) ' uab: ', iat, jat
*      call output(buf, 1, idim, 1, jdim, idim, jdim, 1)
c
      end
