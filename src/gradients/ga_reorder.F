      subroutine ga_reorder(g_a, orow, rmap, ocol, cmap)
      implicit none
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
     $     call errquit('ga_reorder: could not allocate column', dim1)
      if (.not. ma_alloc_get(mt_dbl, dim1, 'gareo', l_vv, k_vv))
     $     call errquit('ga_reorder: could not allocate column2', dim1)
c
*      write(6,*) ' BEFORE '
*      call util_flush(6)
*      call ga_print(g_a)
c
      call ga_sync
      if (.not. ga_duplicate(g_a, g_d, 'ga_reorder')) call errquit
     $     ('ga_reorder: duplicate failed', 0)
      call ga_copy(g_a, g_d)
c
      do j = 1+ga_nodeid(), dim2, ga_nnodes()
         if (orow) then
            call ga_get(g_d, 1, dim1, j, j, dbl_mb(k_v), 1)
            do i = 1, dim1
               dbl_mb(k_vv+rmap(i)-1) = dbl_mb(k_v+i-1)
            enddo
         else
            call ga_get(g_d, 1, dim1, j, j, dbl_mb(k_vv), 1)
         endif
         jj = j
         if (ocol) jj = cmap(j)
         call ga_put(g_a, 1, dim1, jj, jj, dbl_mb(k_vv), 1)
      enddo
      if (.not. ma_free_heap(l_vv)) call errquit('ga_reo: ma?',0)
      if (.not. ma_free_heap(l_v))  call errquit('ga_reo: ma2?',0)
      call ga_sync
      if (.not. ga_destroy(g_d)) call errquit('ga_reo: ga_destroy?',0)
c
*      write(6,*) ' AFTER '
*      call util_flush(6)
*      call ga_print(g_a)
c
      end
      subroutine matrix_reorder(dim1, dim2, a, orow, rmap, ocol, cmap)
      implicit none
#include "mafdecls.fh"
      integer dim1, dim2
      double precision a(dim1,dim2) ! [input/output] Array to be reordered
      logical orow              ! [input] If true reorder the rows
      integer rmap(*)           ! [input] if (orow) a(i,j) -> a(rmap(i),j)
      logical ocol              ! [input] If true reorder the columns
      integer cmap(*)           ! [input] if (ocol) a(i,j) -> a(i,cmap(j))
c
      integer i, j, l_v, k_v, l_vv, k_vv, dim1, dim2, jj
      integer l_d, k_d
c
      if (.not. ma_alloc_get(mt_dbl, dim1*dim2, 'mareo', l_d, k_d))
     $     call errquit('ga_reorder: could not allocate dup', dim1*dim2)
      if (.not. ma_alloc_get(mt_dbl, dim1, 'mareo', l_v, k_v))
     $     call errquit('ga_reorder: could not allocate column', dim1)
      if (.not. ma_alloc_get(mt_dbl, dim1, 'mareo', l_vv, k_vv))
     $     call errquit('ga_reorder: could not allocate column2', dim1)
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
      if (.not. ma_free_heap(l_vv)) call errquit('ma_reo: ma?',0)
      if (.not. ma_free_heap(l_v))  call errquit('ma_reo: ma2?',0)
      if (.not. ma_free_heap(l_d)) call errquit('ma_reo: ma?',0)
c
      end

