
*     ******************************************
*     *                                        *
*     *             beom_create               *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_create to make it thread safe.

      logical function beom_create(geom,namee)
      implicit none
      integer geom
      character*(*) namee

#include "geom.fh"

      integer omp_geom
      logical omp_ok
      common / omp_memhandle_geom / omp_geom,omp_ok

!$OMP SINGLE
      omp_ok = geom_create(omp_geom,namee)
!$OMP END SINGLE
      geom        = omp_geom
      beom_create = omp_ok

      return
      end


*     ******************************************
*     *                                        *
*     *             beom_rtdb_load             *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_rtdb_load to make it thread safe.

      logical function beom_rtdb_load(rtdb,geom,namee)
      implicit none
      integer rtdb,geom
      character*(*) namee

#include "geom.fh"

      integer omp_geom
      logical omp_ok
      common / omp_memhandle_geom / omp_geom,omp_ok

!$OMP SINGLE
      omp_ok = geom_rtdb_load(rtdb,geom,namee)
!$OMP END SINGLE
      beom_rtdb_load = omp_ok
      return
      end


*     ******************************************
*     *                                        *
*     *             beom_destroy               *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_destroy to make it thread safe.
      
      logical function beom_destroy(geom)
      implicit none
      integer geom

#include "geom.fh"

      integer omp_geom
      logical omp_ok
      common / omp_memhandle_geom / omp_geom,omp_ok

!$OMP SINGLE
      omp_ok = geom_destroy(geom)
!$OMP END SINGLE
      beom_destroy = omp_ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_cent_set              *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_set to make it thread safe.
      
      logical function beom_cent_set(geom, icent, t, c, q)
      implicit none
      integer geom
      integer icent
      character*16 t
      double precision c(3)
      double precision q

#include "geom.fh"

      integer omp_geom
      logical omp_ok
      common / omp_memhandle_geom / omp_geom,omp_ok

!$OMP SINGLE
      omp_ok = geom_cent_set(geom, icent, t, c, q)
!$OMP END SINGLE
      beom_cent_set = omp_ok
      return
      end




