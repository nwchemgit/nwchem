
*     ******************************************
*     *                                        *
*     *             btdb_get                   *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save. 
* Right now every thread calls rtdb_get so that array is set for each thread.

      logical function btdb_get(rtdb,namee,datatype,nelem,array)
      implicit none
      integer rtdb,datatype,nelem
      character*(*) namee
      integer array(*)

#include "rtdb.fh"

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
* Only one thread writes to the rtdb.

      logical function btdb_put(rtdb,namee,datatype,nelem,array)
      implicit none
      integer rtdb,datatype,nelem
      character*(*) namee
      integer array(*)

#include "rtdb.fh"

      logical omp_ok
      common / omp_memhandle_cb1 / omp_ok

!$OMP SINGLE
      omp_ok = rtdb_put(rtdb,namee,datatype,nelem,array)
!$OMP END SINGLE
      btdb_put = omp_ok
      return
      end


*     ******************************************
*     *                                        *
*     *             btdb_cget                  *
*     *                                        *
*     ******************************************

*   This is a wrapper around MA_alloc_get to make it thread save.
* Right now every thread calls rtdb_get so that array is set for each thread.

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
* Only one thread writes to the rtdb.

      logical function btdb_cput(rtdb,namee,nelem,array)
      implicit none
      integer rtdb,datatype,nelem
      character*(*) namee
      character*(*) array

#include "rtdb.fh"

      logical omp_ok
      common / omp_memhandle_cb1 / omp_ok

!$OMP SINGLE
      omp_ok = rtdb_cput(rtdb,namee,nelem,array)
!$OMP END SINGLE
      btdb_cput = omp_ok
      return
      end
