      subroutine fock_get_blk( narray, varray, ilo, ihi,
     $     jlo, jhi, buf, tmp)
c     
c     $Id$
c     
      implicit none
#include "global.fh"
c     
c     get the block of the matrix g_array (over basis functions)
c     result into buf(1:idim, 1:jdim) as a dense array. 
c     
      integer narray, varray(narray)
      double precision buf(narray,*)
      double precision tmp(*)
      integer ilo, ihi, jlo, jhi
c
      integer iarray, jarray, i, j, ji
      integer ilen, jlen
c     
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
c     
      if ((ilen.gt.0).and.(jlen.gt.0)) then
         if (narray.gt.1) then
            do iarray=1,narray
c     
c     determine if this array is the same as one already fetched
c     ... if so just copy the data
c     
               do jarray = 1, iarray-1
                  if (varray(jarray) .eq. varray(iarray)) then
                     do j=1,jlen
                        ji = (j-1)*ilen
                        do i=1,ilen
                           buf(iarray,ji+i) = buf(jarray,ji+i)
                        enddo
                     enddo
                     goto 10
                  endif
               enddo
c
c     Array not a duplicate ... got to get it
c
               call ga_get(varray(iarray), ilo, ihi, jlo, jhi,
     $              tmp, ilen)
               do j=1,jlen
                  ji = (j-1)*ilen
                  do i=1,ilen
                     buf(iarray,ji+i) = tmp(ji+i)
                  enddo
               enddo
c
 10            continue
            enddo
         else
            call ga_get(varray(1), ilo, ihi, jlo, jhi, buf, ilen )
         endif
      endif
c     
      end
      subroutine get_atom_block(g_array, basis, iat, jat, buf,
     $     idim, jdim)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
c
      integer g_array, basis, iat, jat, idim, jdim
      double precision buf(*)
c
      integer ilo, ihi, jlo, jhi
c
c     get the block of the matrix g_array (over basis functions)
c     corresponding to the atoms iat and jat, putting the
c     result into buf(1:idim, 1:jdim) as a dense array. 
c     idim, jdim are returned
c
      if (.not. bas_ce2bfr(basis, iat, ilo, ihi))
     $     call errquit('get_atom_block: ce2bfr failed', 0, BASIS_ERR)
      if (.not. bas_ce2bfr(basis, jat, jlo, jhi))
     $     call errquit('get_atom_block: ce2bfr failed', 0, BASIS_ERR)
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
c
      if (idim.gt.0 .and. jdim.gt.0) then
         call ga_get(g_array, ilo, ihi, jlo, jhi, buf, idim)
      endif
*      write(6,*) ' gab: ', iat, jat
*      call output(buf, 1, idim, 1, jdim, idim, jdim, 1)
c
      end
