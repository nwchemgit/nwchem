
*     ******************************************
*     *                                        *
*     *             BA_alloc_get               *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save.

      logical function BA_alloc_get(datatype,nelem,namee,memhandle,indx)
      implicit none
#include "mafdecls.fh"
      integer datatype,nelem
      character*(*) namee
      integer   memhandle,indx

      !MAPOINTER omp_memhandle,omp_indx
      !logical   omp_ok
      !common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
      integer tid,Parallel_threadid
      external Parallel_threadid
      logical ok
      tid = Parallel_threadid()

!$OMP SINGLE 
      ok = MA_alloc_get(datatype,nelem,namee,memhandle,indx)
      !write(*,*) tid," BA_alloc_get memhandle=",memhandle
!$OMP END SINGLE copyprivate(memhandle,indx,ok)

      BA_alloc_get = ok
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
      integer   memhandle

      !MAPOINTER omp_memhandle,omp_indx
      !logical   omp_ok
      !common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
      logical ok

!$OMP SINGLE
      ok = MA_free_heap(memhandle)
!$OMP END SINGLE copyprivate(ok)

      BA_free_heap = ok
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
      integer memhandle,indx

      !MAPOINTER omp_memhandle,omp_indx
      !logical   omp_ok
      !common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
      logical ok

!$OMP SINGLE
      ok = MA_push_get(datatype,nelem,namee,memhandle,indx)
!$OMP END SINGLE copyprivate(memhandle,indx,ok)
      BA_push_get = ok
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
      integer   memhandle

      !MAPOINTER omp_memhandle,omp_indx
      !logical   omp_ok
      !common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok

      logical ok
!$OMP SINGLE
      ok = MA_pop_stack(memhandle)
!$OMP END SINGLE copyprivate(ok)

      BA_pop_stack = ok
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
      integer   memhandle,indx

      !MAPOINTER omp_memhandle,omp_indx
      !logical   omp_ok
      !common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
      logical ok

!$OMP SINGLE
      ok = MA_get_index(memhandle,indx)
!$OMP END SINGLE copyprivate(ok,indx)

      BA_get_index = ok
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

      !MAPOINTER omp_memhandle,omp_indx
      !logical   omp_ok
      !common / omp_memhandle_cb / omp_memhandle,omp_indx,omp_ok
      logical ok2

!$OMP SINGLE
      ok2 = MA_set_auto_verify(ok)
!$OMP END SINGLE copyprivate(ok2)

      BA_set_auto_verify = ok2
      return
      end

