
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

!$OMP SINGLE
      btdb_get = rtdb_get(rtdb,namee,datatype,nelem,array)
!$OMP END SINGLE
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

!$OMP SINGLE
      btdb_put = rtdb_put(rtdb,namee,datatype,nelem,array)
!$OMP END SINGLE
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

!$OMP SINGLE
      btdb_cget = rtdb_cget(rtdb,namee,nelem,array)
!$OMP END SINGLE
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

!$OMP SINGLE
      btdb_cput = rtdb_cput(rtdb,namee,nelem,array)
!$OMP END SINGLE
      return
      end


