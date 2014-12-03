
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

!$OMP SINGLE
      BA_alloc_get = MA_alloc_get(datatype,nelem,name,memhandle,indx)
!$OMP END SINGLE
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

*   This is a wrapper around MA_push_get to make it thread save.

      logical function BA_push_get(datatype,nelem,name,memhandle,indx)
      implicit none
#include "mafdecls.fh"
      integer datatype,nelem
      character*(*) name
      MAPOINTER memhandle,indx

!$OMP SINGLE
      BA_push_get = MA_push_get(datatype,nelem,name,memhandle,indx)
!$OMP END SINGLE
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

!$OMP SINGLE
      BA_pop_stack = MA_pop_stack(memhandle)
!$OMP END SINGLE
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

!$OMP SINGLE
      BA_get_index = MA_get_index(memhandle,indx)
!$OMP END SINGLE
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

!$OMP SINGLE
      BA_set_auto_verify = MA_set_auto_verify(ok)
!$OMP END SINGLE
      return
      end

