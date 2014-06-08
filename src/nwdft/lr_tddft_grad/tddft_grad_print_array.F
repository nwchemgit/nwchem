c
c     Utility routine to print out a global array in a sensible way.
c     HvD 9/2006
c
      subroutine tddft_grad_print_array(ipol,nroot,g_array,scale)
c
c     It is assumed that the global array is stored in a 3-dimensional
c     array. In my code the roots are the first dimension and I want
c     the matrices printed a root at a time. The other dimensions 
c     are enquired from the GA it-self.
c
c     This is a collective operation to ensure the integrity of the
c     printed data.
c
      implicit none

#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
c     Input:
c
      integer ipol          ! =1 (restricted), =2 (unrestricted)
      integer nroot         ! the number of roots to be treated
      integer g_array(ipol) ! the global array we want to print
c
      double precision scale ! scale factor useful for comparing 
                             ! restricted and unrestricted results
c
c     Local:
c
      integer ilo(3),ihi(3,2),type,ndim,ld(3)
      integer nsize,lbuf,ibuf
      integer ip,ir,ielm
      logical oroot
c
      call ga_sync
c
      oroot = 0.eq.ga_nodeid()
c
      ilo(1) = 1
      ilo(2) = 1
      ilo(3) = 1
      nsize  = 0
      ndim   = ga_ndim(g_array(1))
      do ip = 1, ipol
        call nga_inquire(g_array(ip),type,ndim,ihi(1,ip))
        if (ndim.eq.2) then
          nsize = max(nsize,ihi(1,ip)*ihi(2,ip))
        else if (ndim.eq.3) then
          nsize = max(nsize,ihi(2,ip)*ihi(3,ip))
        else
          call errquit("tddft_grad_print_array: funny dimensions", 0, 
     +         UERR)
        endif
      enddo
      if (.not.MA_push_get(MT_DBL,nsize,"print buffer",lbuf,ibuf)) call
     +  errquit('tddft_grad_print_array: failed to allocate buffer',
     +  0,UERR)
      if (ndim.eq.2) then
        do ip = 1, ipol
          ld(1)     = ihi(1,ip)
          ld(2)     = ihi(2,ip)
          if (oroot) then
            call nga_get(g_array(ip),ilo,ihi(1,ip),dbl_mb(ibuf),ld)
            do ielm = 0, nsize-1
              dbl_mb(ibuf+ielm) = dbl_mb(ibuf+ielm)*scale
            enddo
            call print_matrix(dbl_mb(ibuf),ihi(1,ip),ihi(2,ip),ld(1))
          endif
        enddo
      else if (ndim.eq.3) then
        if (oroot.and.nroot.gt.1) then
          write(LuOut,"('- Printing matrices for ',i2,' roots.')")nroot
        endif
        do ir = 1, nroot
          if (oroot.and.nroot.gt.1) then
            write(LuOut,"('- Printing matrix for root ',i2)")ir
          endif
          do ip = 1, ipol
            ilo(1)    = ir
            ihi(1,ip) = ir
            ld(1)     = 1
            ld(2)     = ihi(2,ip)
            ld(3)     = ihi(3,ip)
            if (oroot) then
              call nga_get(g_array(ip),ilo,ihi(1,ip),dbl_mb(ibuf),ld)
              do ielm = 0, nsize-1
                dbl_mb(ibuf+ielm) = dbl_mb(ibuf+ielm)*scale
              enddo
              call print_matrix(dbl_mb(ibuf),ihi(2,ip),ihi(3,ip),ld(2))
            endif
          enddo
        enddo
      endif
      if (.not.MA_pop_stack(lbuf)) call errquit(
     +  'tddft_grad_print_array: failed to deallocate buffer',
     +  0,UERR)
      call ga_sync
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine print_matrix(a,nrow,ncol,ndim)
      implicit none
c
c     $Id$
c
c     Print the matrix A(1:nrow,1:ncol).
c
c     Written by Huub van Dam, September 2006
c
#include "stdio.fh"
c
c     Input:
c
      integer nrow
      integer ncol
      integer ndim
      double precision a(ndim,ncol)
c
c     Local:
c
      integer i,j,jlo,jhi
      integer mxcol
      parameter (mxcol=10)
c
      if (nrow.gt.0) then
        jlo = 1
        jhi = 0
        do while (jlo.le.ncol) 
          jhi = jhi+mxcol
          jhi = min(jhi,ncol)
          write(luout,100)(j,j=jlo,jhi)
          do i = 1, nrow
            write(luout,110)i,(a(i,j),j=jlo,jhi)
          enddo
          jlo = jhi + 1
        enddo
      endif
 100  format(6x,10(5x,i3,4x))
 110  format(i5,1x,10f12.6)
      end
c
c-----------------------------------------------------------------------
c
      subroutine print_matrix_patch(a,irow,jrow,icol,jcol,ndim)
      implicit none
c
c     $Id$
c
c     Print the matrix A(irow:jrow,icol:jcol).
c
c     Written by Huub van Dam, September 2006
c
#include "stdio.fh"
#include "errquit.fh"
c
c     Input:
c
      integer irow,jrow
      integer icol,jcol
      integer ndim
      double precision a(ndim,jcol)
c
c     Local:
c
      integer i,j,jlo,jhi
      integer mxcol
      parameter (mxcol=10)
c
      if (jrow.gt.ndim) then
        call errquit("print_matrix_patch: jrow out of range",0, UERR)
      endif
      if (jrow-irow+1.gt.0) then
        jlo = icol
        jhi = icol-1
        do while (jlo.le.jcol) 
          jhi = jhi+mxcol
          jhi = min(jhi,jcol)
          write(luout,100)(j,j=jlo,jhi)
          do i = irow, jrow
            write(luout,110)i,(a(i,j),j=jlo,jhi)
          enddo
          jlo = jhi + 1
        enddo
      endif
 100  format(6x,10(5x,i3,4x))
 110  format(i5,1x,10f12.6)
      end
