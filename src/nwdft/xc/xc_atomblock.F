      subroutine xc_atom_block(array, nbf,
     .     basis, iat, jat, buf)
C$Id: xc_atomblock.F,v 1.2 2003-10-17 22:55:53 carlfahl Exp $
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
!      write(6,*) ' uabentry: ', ilo, ihi,jlo,jhi
!      call output(array, 1, nbf, 1, nbf, nbf, nbf, 1)
c     
      if (idim.gt.0 .and. jdim.gt.0) then
!         call ga_acc(g_array, ilo, ihi, jlo, jhi, buf, idim,
!     $        1.0d0)
         do j=jlo,jhi
            icount=1+(j-jlo)*idim
            call daxpy(idim,1d0,buf(icount),1,array(ilo,j),1)
!            do i=ilo,ihi
!               array(i,j)=array(i,j)+buf(icount+i-ilo)
!!               write(0,*) 'ij',i,j,array(i,j),buf(icount+i-ilo)
!            enddo
         enddo
      endif
!      write(6,*) ' uab: ', ilo, ihi,jlo,jhi
!      call output(buf, 1, idim, 1, jdim, idim, jdim, 1)
!      call output(array, 1, nbf, 1, nbf, nbf, nbf, 1)
c
      end
