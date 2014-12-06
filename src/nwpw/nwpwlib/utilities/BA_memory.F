
*     ******************************************
*     *                                        *
*     *             BA_alloc_get               *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save.

      logical function BA_alloc_get(datatype,nelem,name,memhandle,indx)
      implicit none
#include "mafdecls.fh"
      integer datatype,nelem
      character*(*) name
      MAPOINTER memhandle,indx

      MAPOINTER omp_memhandle,omp_indx
      logical   omp_ok
      common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok

!$OMP SINGLE
c      omp_ok = MA_set_error_print(.true.)
c      omp_ok = MA_set_auto_verify(.true.)
c      omp_ok = MA_set_hard_fail(.true.)
      omp_ok = MA_alloc_get(datatype,nelem,name,omp_memhandle,omp_indx)
c      omp_ok = MA_verify_allocator_stuff()
!$OMP END SINGLE
      memhandle    = omp_memhandle
      indx         = omp_indx
      BA_alloc_get = omp_ok
      return
      end

*     ******************************************
*     *                                        *
*     *             BA_free_heap               *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_free_heap to make it thread save.

      logical function BA_free_heap(memhandle)
      implicit none
#include "mafdecls.fh"
      MAPOINTER memhandle

      MAPOINTER omp_memhandle,omp_indx
      logical   omp_ok
      common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok

!$OMP SINGLE
c      omp_ok = MA_set_auto_verify(.true.)
c      omp_ok = MA_set_hard_fail(.true.)
      omp_ok = MA_free_heap(memhandle)
c      omp_ok = MA_verify_allocator_stuff()
!$OMP END SINGLE
      BA_free_heap = omp_ok
      return
      end

*     ******************************************
*     *                                        *
*     *             BA_push_get                *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_push_get to make it thread save.

      logical function BA_push_get(datatype,nelem,namee,memhandle,indx)
      implicit none
#include "mafdecls.fh"
      integer datatype,nelem
      character*(*) namee
      MAPOINTER memhandle,indx

      MAPOINTER omp_memhandle,omp_indx
      logical   omp_ok
      common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
!$OMP BARRIER
!$OMP SINGLE
c      omp_ok = MA_set_auto_verify(.true.)
c      omp_ok = MA_set_hard_fail(.true.)
      omp_ok = MA_push_get(datatype,nelem,namee,omp_memhandle,omp_indx)
c      omp_ok = MA_verify_allocator_stuff()
!$OMP END SINGLE
      memhandle   = omp_memhandle
      indx        = omp_indx
      BA_push_get = omp_ok
      return
      end

*     ******************************************
*     *                                        *
*     *             BA_pop_stack               *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_pop_stack to make it thread save.

      logical function BA_pop_stack(memhandle)
      implicit none
#include "mafdecls.fh"
      MAPOINTER memhandle

      MAPOINTER omp_memhandle,omp_indx
      logical   omp_ok
      common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
!$OMP SINGLE
c      omp_ok = MA_set_auto_verify(.true.)
c      omp_ok = MA_set_hard_fail(.true.)
      omp_ok = MA_pop_stack(memhandle)
c      omp_ok = MA_verify_allocator_stuff()
!$OMP END SINGLE
      BA_pop_stack = omp_ok
      return
      end

*     ******************************************
*     *                                        *
*     *             BA_get_index               *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_get_index to make it thread save.

      logical function BA_get_index(memhandle,indx)
      implicit none
#include "mafdecls.fh"
      MAPOINTER memhandle,indx

      MAPOINTER omp_memhandle,omp_indx
      logical   omp_ok
      common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok

!$OMP SINGLE
c      omp_ok = MA_set_auto_verify(.true.)
      omp_ok = MA_get_index(memhandle,omp_indx)
c      omp_ok = MA_verify_allocator_stuff()
!$OMP END SINGLE
      indx         = omp_indx
      BA_get_index = omp_ok
      return
      end


*     ******************************************
*     *                                        *
*     *             BA_set_auto_verify         *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_set_auto_verify to make it thread save.

      logical function BA_set_auto_verify(ok)
      implicit none
#include "mafdecls.fh"
      logical ok

      MAPOINTER omp_memhandle,omp_indx
      logical   omp_ok
      common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
!$OMP SINGLE
      omp_ok = MA_set_auto_verify(ok)
!$OMP END SINGLE
      BA_set_auto_verify = omp_ok
      return
      end

