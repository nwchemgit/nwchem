      subroutine ga_scale_lh(g_a,vector)
c $Id$
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
      double precision vector(*)
c
c     Do full matrix times diagonal matrix with the diagonal
c     matrix on the left hand side stored as a vector
c     This boils down to a row wise scaling of the g_a
c    
c     do i=1,nbf
c        do j=1,nbf
c           g_a(i,j)=g_a(i,j)*vector(i)
c        enddo
c     enddo
c
      integer ma_type, dim1, dim2, n
      integer i, j, ilo, ihi, jlo, jhi
      double precision number
      double complex znumber
c
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.eq.mt_dbl) then
      else if (ma_type.eq.mt_dcpl) then
      else
        write(0,*) ' ma_type ',ma_type,' mt_dbl  ',mt_dbl
        write(0,*) ' ma_type ',ma_type,' mt_dcpl ',mt_dcpl
        call errquit
     $     ('ga_trace_diagonal: array is not real nor complex ', g_a,0)
      endif
      n = dim1
c
c     Extract and sum the diags local to each process
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         if (ma_type.eq.MT_DBL) then
            do i = ilo,ihi
               do j = jlo,jhi
                  call ga_get(g_a, i, i, j, j, number, 1)
                  number = number*vector(i)
                  call ga_put(g_a,i, i, j, j, number, 1)
               enddo
            enddo
         else if (ma_type.eq.MT_DCPL) then
            do i = ilo,ihi
               do j = jlo,jhi
                  call ga_get(g_a, i, i, j, j, znumber, 1)
                  znumber = znumber*vector(i)
                  call ga_put(g_a,i, i, j, j, znumber, 1)
               enddo
            enddo
         endif
      endif
c
      call ga_sync()
c
      end
c
      subroutine ga_scale_rh(g_a,vector)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
      double precision vector(*)
c
c     Do diagonal matrix times full matrix with the diagonal
c     matrix on the right hand side stored as a vector
c     This boils down to a column wise scaling of the g_a
c    
c     do i=1,nbf
c        do j=1,nbf
c           g_a(i,j)=g_a(i,j)*vector(j)
c        enddo
c     enddo
c
      integer ma_type, dim1, dim2, n
      integer i, j, ilo, ihi, jlo, jhi
      double precision number
c
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dbl .or. dim1.ne.dim2) call errquit
     $     ('ga_trace_diagonal: array is not square/real', g_a,0)
      n = dim1
c
c
c     Extract and sum the diags local to each process
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         do i = ilo,ihi
            do j = jlo,jhi
               call ga_get(g_a, i, i, j, j, number, 1)
               number = number*vector(j)
               call ga_put(g_a,i, i, j, j, number, 1)
            enddo
         enddo
      endif
c
      call ga_sync()
c
      end
c
      subroutine ga_scale_rlh(g_a,vector1,vector2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
      double precision vector1(*),vector2(*)
c
c     Do diagonal matrix times full matrix times diagonal matrix 
c     with the diagonal matrices on left and right hand side stored
c     as a vector. This boils down to a column plus row wise 
c     scaling of the g_a, effectively element wise
c    
c     do i=1,nbf
c        do j=1,nbf
c           g_a(i,j)=g_a(i,j)*vector1(i)*vector2(j)
c        enddo
c     enddo
c
      integer ma_type, dim1, dim2, n
      integer i,j, ilo, ihi, jlo, jhi
      double precision number
c
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dbl .or. dim1.ne.dim2) call errquit
     $     ('ga_trace_diagonal: array is not square/real', g_a,0)
      n = dim1
c
c
c     Extract and sum the diags local to each process
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         do i = ilo,ihi
            do j = jlo,jhi
               call ga_get(g_a, i, i, j, j, number, 1)
               number = number*vector1(i)*vector2(j)
               call ga_put(g_a,i, i, j, j, number, 1)
            enddo
         enddo
      endif
c
      call ga_sync()
c
      end
c
      subroutine ga_scale_divEpp(g_a,vector)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
      double precision vector(*)
c
c     Do diagonal matrix times full matrix times diagonal matrix 
c     with the diagonal matrices on left and right hand side stored
c     as a vector. This boils down to a column plus row wise 
c     scaling of the g_a, effectively element wise
c    
c     do i=1,nbf
c        do j=1,nbf
c           g_a(i,j)=g_a(i,j)/(vector(i)+vector(j))
c        enddo
c     enddo
c
      integer ma_type, dim1, dim2, n
      integer i,j, ilo, ihi, jlo, jhi
      double precision number
c
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dbl .or. dim1.ne.dim2) call errquit
     $     ('ga_trace_diagonal: array is not square/real', g_a,0)
      n = dim1
c
c
c     Extract and sum the diags local to each process
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         do i = ilo,ihi
            do j = jlo,jhi
               call ga_get(g_a, i, i, j, j, number, 1)
               number = number/(vector(i)+vector(j))
               call ga_put(g_a,i, i, j, j, number, 1)
            enddo
         enddo
      endif
c
      call ga_sync()
c
      end
c
      subroutine ga_copy_todiag(g_a,vector)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
      double precision vector(*)
c
c     Put the elements of the vector on the diagonal of the matrix in g_a
c    
c     do i=1,nbf
c        do i=1,nbf
c           g_a(i,i)=vector(i)
c        enddo
c     enddo
c
      integer ma_type, dim1, dim2, n
      integer i, ilo, ihi, jlo, jhi
c
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dbl .or. dim1.ne.dim2) call errquit
     $     ('ga_trace_diagonal: array is not square/real', g_a,0)
      n = dim1
c
c
c     Put elements on diagonal of g_a
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         ilo = max(ilo,jlo)
         ihi = min(ihi,jhi)
         do i = ilo,ihi
            call ga_put(g_a,i, i, i, i, vector(i), 1)
         enddo
      endif
c
      call ga_sync()
c
      end
c
      subroutine ga_add_todiag(g_a,value)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
      double precision value
c
c     Add constant to the diagonal of the matrix in g_a
c    
c     do i=1,nbf
c        do i=1,nbf
c           g_a(i,i)=g_a(i,i)+value
c        enddo
c     enddo
c
      integer ma_type, dim1, dim2, n
      integer i, ilo, ihi, jlo, jhi
      double precision number
c
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dbl .or. dim1.ne.dim2) call errquit
     $     ('ga_trace_diagonal: array is not square/real', g_a,0)
      n = dim1
c
c
c     Put elements on diagonal of g_a
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         ilo = max(ilo,jlo)
         ihi = min(ihi,jhi)
         do i = ilo,ihi
               call ga_get(g_a, i, i, i, i, number, 1)
               number = number+value
            call ga_put(g_a,i, i, i, i, number, 1)
         enddo
      endif
c
      call ga_sync()
c
      end
