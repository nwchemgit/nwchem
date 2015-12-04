      subroutine ga_reorder(g_a, orow, rmap, ocol, cmap)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_a               ! [input] Handle to array to be reordered
      logical orow              ! [input] If true reorder the rows
      integer rmap(*)         ! [input] if (orow) a(i,j) -> a(rmap(i),j)
      logical ocol              ! [input] If true reorder the columns
      integer cmap(*)         ! [input] if (ocol) a(i,j) -> a(i,cmap(j))
c
      integer g_d, i, j, l_v, k_v, l_vv, k_vv, ma_type, dim1, dim2, jj
c
      call ga_inquire(g_a, ma_type, dim1, dim2)
c     
      if (.not. ma_alloc_get(mt_dbl, dim1, 'gareo', l_v, k_v))
     $     call errquit('ga_reorder: could not allocate column', dim1,
     &       MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, dim1, 'gareo', l_vv, k_vv))
     $     call errquit('ga_reorder: could not allocate column2', dim1,
     &       MA_ERR)
c
*      write(6,*) ' BEFORE '
*      call util_flush(6)
*      call ga_print(g_a)
c
      call ga_sync
      if (.not. ga_duplicate(g_a, g_d, 'ga_reorder')) call errquit
     $     ('ga_reorder: duplicate failed', 0, GA_ERR)
      call ga_copy(g_a, g_d)
c
      do j = 1+ga_nodeid(), dim2, ga_nnodes()
         if (orow) then
            call ga_get(g_d, 1, dim1, j, j, dbl_mb(k_v), dim1)
            do i = 1, dim1
               dbl_mb(k_vv+rmap(i)-1) = dbl_mb(k_v+i-1)
            enddo
         else
            call ga_get(g_d, 1, dim1, j, j, dbl_mb(k_vv), dim1)
         endif
         jj = j
         if (ocol) jj = cmap(j)
         call ga_put(g_a, 1, dim1, jj, jj, dbl_mb(k_vv), dim1)
      enddo
      if (.not. ma_free_heap(l_vv)) call errquit('ga_reo: ma?',0,
     &       MA_ERR)
      if (.not. ma_free_heap(l_v))  call errquit('ga_reo: ma2?',0,
     &       MA_ERR)
      call ga_sync
      if (.not. ga_destroy(g_d)) call errquit('ga_reo: ga_destroy?',0,
     &       GA_ERR)
c
*      write(6,*) ' AFTER '
*      call util_flush(6)
*      call ga_print(g_a)
c
      end
      subroutine nga_reorder(g_a, orow, rmap, ocol, cmap)
      implicit none
#include "errquit.fh"
c
c This is basically just an extension of ga_reorder and is not very
c generic at this point.  As a matter of fact, it assumes (and tests)
c that the dimension is 3 and that you only want to reorder the last
c two indices.  This can be made more general after I test this version.
c Also, I am wasting a lot of memory by duplicating the whole ga.  This
c will need to be optimized in the future.
c
#include "global.fh"
#include "mafdecls.fh"
      integer g_a               ! [input] Handle to array to be reordered
      logical orow              ! [input] If true reorder the rows
      integer rmap(*)         ! [input] if (orow) a(i,j) -> a(rmap(i),j)
      logical ocol              ! [input] If true reorder the columns
      integer cmap(*)         ! [input] if (ocol) a(i,j) -> a(i,cmap(j))
c
      integer g_d, i, j, l_v, k_v, l_vv, k_vv, ma_type 
      integer dim0, dim1, dim2, jj
      integer ndim, dims(3), lo(3), hi(3), ld(2)
c
      ndim = ga_ndim(g_a)
      if (ndim.ne.3)
     $     call errquit('nga_reorder: must have 3 dimensions',ndim,
     &       GA_ERR)
      call nga_inquire(g_a, ma_type, ndim, dims)
      dim1 = dims(2)
      dim2 = dims(3)
c     
      if (.not. ma_alloc_get(mt_dbl, dim1, 'gareo', l_v, k_v))
     $     call errquit('ga_reorder: could not allocate column', dim1,
     &       GA_ERR)
      if (.not. ma_alloc_get(mt_dbl, dim1, 'gareo', l_vv, k_vv))
     $     call errquit('ga_reorder: could not allocate column2', dim1,
     &       GA_ERR)
c
*      write(6,*) ' BEFORE '
*      call util_flush(6)
*      call ga_print(g_a)
c
      call ga_sync
      if (.not. ga_duplicate(g_a, g_d, 'ga_reorder')) call errquit
     $     ('ga_reorder: duplicate failed', 0, GA_ERR)
      call ga_copy(g_a, g_d)
c
      ld(1) = 1
      lo(2) = 1
      hi(2) = dim1
      ld(2) = dim1
      do dim0 = 1, dims(1)
      lo(1) = dim0
      hi(1) = dim0
      do j = 1+ga_nodeid(), dim2, ga_nnodes()
         lo(3) = j
         hi(3) = j
         if (orow) then
            call nga_get(g_d, lo, hi, dbl_mb(k_v), ld)
            do i = 1, dim1
               dbl_mb(k_vv+rmap(i)-1) = dbl_mb(k_v+i-1)
            enddo
         else
            call nga_get(g_d, lo, hi, dbl_mb(k_vv), ld)
         endif
         jj = j
         if (ocol) jj = cmap(j)
         lo(3) = jj
         hi(3) = jj
         call nga_put(g_a, lo, hi, dbl_mb(k_vv), ld)
      enddo
      enddo
      if (.not. ma_free_heap(l_vv)) call errquit('ga_reo: ma?',0,
     &       MA_ERR)
      if (.not. ma_free_heap(l_v))  call errquit('ga_reo: ma2?',0,
     &       MA_ERR)
      call ga_sync
      if (.not. ga_destroy(g_d)) call errquit('ga_reo: ga_destroy?',0,
     &       GA_ERR)
c
*      write(6,*) ' AFTER '
*      call util_flush(6)
*      call ga_print(g_a)
c
      end
      subroutine matrix_reorder(dim1, dim2, a, orow, rmap, ocol, cmap)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer dim1, dim2
      double precision a(dim1,dim2) ! [input/output] Array to be reordered
      logical orow              ! [input] If true reorder the rows
      integer rmap(*)           ! [input] if (orow) a(i,j) -> a(rmap(i),j)
      logical ocol              ! [input] If true reorder the columns
      integer cmap(*)           ! [input] if (ocol) a(i,j) -> a(i,cmap(j))
c
      integer i, j, l_v, k_v, l_vv, k_vv,jj
      integer l_d, k_d
c
      if (.not. ma_alloc_get(mt_dbl, dim1*dim2, 'mareo', l_d, k_d))
     $     call errquit('ga_reorder: could not allocate dup', dim1*dim2,
     &       MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, dim1, 'mareo', l_v, k_v))
     $     call errquit('ga_reorder: could not allocate column', dim1,
     &       MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, dim1, 'mareo', l_vv, k_vv))
     $     call errquit('ga_reorder: could not allocate column2', dim1,
     &       MA_ERR)
c
      call dcopy(dim1*dim2, a, 1, dbl_mb(k_d), 1)
c
*      write(6,*) ' BEFORE '
*      call output(a, 1, dim1, 1, dim2, dim1, dim2, 1)
c
      do j = 1, dim2
         if (orow) then
            call dcopy(dim1, dbl_mb(k_d+(j-1)*dim1), 1, dbl_mb(k_v), 1)
            do i = 1, dim1
               dbl_mb(k_vv+rmap(i)-1) = dbl_mb(k_v+i-1)
            enddo
         else
            call dcopy(dim1, dbl_mb(k_d+(j-1)*dim1), 1, dbl_mb(k_vv), 1)
         endif
         jj = j
         if (ocol) jj = cmap(j)
         call dcopy(dim1, dbl_mb(k_vv), 1, a(1,jj), 1)
      enddo
c
*      write(6,*) ' AFTER '
*      call output(a, 1, dim1, 1, dim2, dim1, dim2, 1)
c
      if (.not. ma_free_heap(l_vv)) call errquit('ma_reo: ma?',0,
     &       MA_ERR)
      if (.not. ma_free_heap(l_v))  call errquit('ma_reo: ma2?',0,
     &       MA_ERR)
      if (.not. ma_free_heap(l_d)) call errquit('ma_reo: ma?',0,
     &       MA_ERR)
c
      end

