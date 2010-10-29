      subroutine xc_atom_block(array, nbf,
     .     basis, iat, jat, buf)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
c
      integer basis, iat, jat
      integer nbf
      double precision array(nbf,*)
      double precision buf(*)
c
      integer ilo, ihi, jlo, jhi, idim, jdim
      integer i,j,icount
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
         do j=jlo,jhi
            icount=1+(j-jlo)*idim
            call daxpy(idim,1d0,buf(icount),1,array(ilo,j),1)
         enddo
      endif
c
      end
      subroutine xc_atom_blockd(array, nbf,
     .     basis, iat, jat, buf)
c     fills only lower triangle
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
c
      integer basis, iat, jat
      integer nbf
      double precision array(*)
      double precision buf(*)
c
      integer ilo, ihi, jlo, jhi, idim, jdim
      integer i,j,icount
      integer ij,nnn
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
         if(ilo.ne.jlo) then
            do j=jlo,jhi
               icount=1+(j-jlo)*idim
               ij=((j-1)*(2*(nbf+1)-j)+1)/2+ilo-j+1
               call daxpy(idim,1d0,buf(icount),1,array(ij),1)
            enddo
         else
c do only lower tr
            nnn=idim
            do j=jlo,jhi
               icount=1+(j-jlo)*idim+j-jlo
               ij=((j-1)*(2*(nbf+1)-j)+1)/2+
     +              1
               call daxpy(nnn,1d0,buf(icount),1,array(ij),1)
               nnn=nnn-1
            enddo
         endif
      endif
c
      end
