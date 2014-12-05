
*     ******************************************
*     *                                        *
*     *             btdb_get                   *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save.

      logical function btdb_get(rtdb,namee,datatype,nelem,array)
      implicit none
#include "rtdb.fh"
      integer rtdb,datatype,nelem
      character*(*) namee
      integer array(*)

      logical ok

!$OMP CRITICAL
      ok = rtdb_get(rtdb,namee,datatype,nelem,array)
!$OMP END CRITICAL
      btdb_get = ok
      return
      end


*     ******************************************
*     *                                        *
*     *             btdb_put                   *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save.

      logical function btdb_put(rtdb,namee,datatype,nelem,array)
      implicit none
#include "rtdb.fh"
      integer rtdb,datatype,nelem
      character*(*) namee
      integer array(*)
      logical ok

!$OMP CRITICAL
      ok = rtdb_put(rtdb,namee,datatype,nelem,array)
!$OMP END CRITICAL
      btdb_put = ok
      return
      end


*     ******************************************
*     *                                        *
*     *             btdb_cget                  *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save.

      logical function btdb_cget(rtdb,namee,nelem,array)
      implicit none
#include "rtdb.fh"
      integer rtdb,datatype,nelem
      character*(*) namee
      character*(*) array
      logical ok

!$OMP CRITICAL
      ok = rtdb_cget(rtdb,namee,nelem,array)
!$OMP END CRITICAL
      btdb_cget = ok
      return
      end


*     ******************************************
*     *                                        *
*     *             btdb_cput                  *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save.

      logical function btdb_cput(rtdb,namee,nelem,array)
      implicit none
#include "rtdb.fh"
      integer rtdb,datatype,nelem
      character*(*) namee
      character*(*) array
      logical ok

!$OMP CRITICAL
      ok = rtdb_cput(rtdb,namee,nelem,array)
!$OMP END CRITICAL
      btdb_cput = ok
      return
      end


