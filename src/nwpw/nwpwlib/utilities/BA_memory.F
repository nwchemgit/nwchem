
*     ******************************************
*     *                                        *
*     *             BA_alloc_get               *
*     *                                        *
*     ******************************************

*   This is a wrapper around BA_alloc_get to make it thread save.

      logical function BA_alloc_get(datatype,nelem,name,memhandle,index)
      implicit none
#include "mafdecls.fh"
      integer datatype,nelem
      character*(*) name
      MAPOINTER memhandle,index

!$OMP SINGLE
      BA_alloc_get = MA_alloc_get(datatype,nelem,name,memhandle,index)
!$OMP END SINGLE
      return
      end

*     ******************************************
*     *                                        *
*     *             BA_free_heap               *
*     *                                        *
*     ******************************************

*   This is a wrapper around BA_free_heap to make it thread save.

      logical function BA_free_heap(memhandle)
      implicit none
#include "mafdecls.fh"
      MAPOINTER memhandle

!$OMP SINGLE
      BA_free_heap = MA_free_heap(memhandle)
!$OMP END SINGLE
      return
      end

*     ******************************************
*     *                                        *
*     *             BA_push_get                *
*     *                                        *
*     ******************************************

*   This is a wrapper around BA_push_get to make it thread save.

      logical function BA_push_get(datatype,nelem,name,memhandle,index)
      implicit none
#include "mafdecls.fh"
      integer datatype,nelem
      character*(*) name
      MAPOINTER memhandle,index

!$OMP SINGLE
      BA_push_get = MA_push_get(datatype,nelem,name,memhandle,index)
!$OMP END SINGLE
      return
      end

*     ******************************************
*     *                                        *
*     *             BA_pop_stack               *
*     *                                        *
*     ******************************************

*   This is a wrapper around BA_pop_stack to make it thread save.

      logical function BA_pop_stack(memhandle)
      implicit none
#include "mafdecls.fh"
      integer*8 memhandle

!$OMP SINGLE
      BA_pop_stack = MA_pop_stack(memhandle)
!$OMP END SINGLE
      return
      end
