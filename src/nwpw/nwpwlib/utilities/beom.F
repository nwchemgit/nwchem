
*     ******************************************
*     *                                        *
*     *             beom_create                *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_create to make it thread safe.

      logical function beom_create(geom,namee)
      implicit none
      integer geom
      character*(*) namee

#include "geom.fh"

      !integer omp_geom
      !logical omp_ok
      !common / omp_memhandle_geom / omp_geom,omp_ok
      logical ok

!$OMP SINGLE
      ok = geom_create(geom,namee)
!$OMP END SINGLE copyprivate(geom,ok)
      beom_create = ok

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

      !integer omp_geom
      !logical omp_ok
      !common / omp_memhandle_geom / omp_geom,omp_ok
      logical ok

!$OMP SINGLE
      ok = geom_rtdb_load(rtdb,geom,namee)
!$OMP END SINGLE copyprivate(ok)
      beom_rtdb_load = ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_rtdb_store            *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_rtdb_store to make it thread safe.

      logical function beom_rtdb_store(rtdb,geom,namee)
      implicit none
      integer rtdb,geom
      character*(*) namee

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_rtdb_store(rtdb,geom,namee)
!$OMP END SINGLE copyprivate(ok)
      beom_rtdb_store = ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_rtdb_delete           *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_rtdb_store to make it thread safe.

      logical function beom_rtdb_delete(rtdb,namee)
      implicit none
      integer rtdb
      character*(*) namee

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_rtdb_delete(rtdb,namee)
!$OMP END SINGLE copyprivate(ok)

      beom_rtdb_delete = ok
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

      !integer omp_geom
      !logical omp_ok
      !common / omp_memhandle_geom / omp_geom,omp_ok
      logical ok

!$OMP SINGLE
      ok = geom_destroy(geom)
!$OMP END SINGLE copyprivate(ok)
      beom_destroy = ok
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

      !integer omp_geom
      !logical omp_ok
      !common / omp_memhandle_geom / omp_geom,omp_ok
      logical ok

!$OMP SINGLE
      ok = geom_cent_set(geom, icent, t, c, q)
!$OMP END SINGLE copyprivate(ok)
      beom_cent_set = ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_ncent_set             *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_ncent_set to make it thread safe.
      
      logical function beom_ncent_set(geom, ncent)
      implicit none
      integer geom
      integer ncent


#include "geom.fh"

      !integer omp_geom
      !logical omp_ok
      !common / omp_memhandle_geom / omp_geom,omp_ok
      logical ok

!$OMP SINGLE
      ok = geom_ncent_set(geom, ncent)
!$OMP END SINGLE copyprivate(ok)
      beom_ncent_set = ok
      return
      end




*     ******************************************
*     *                                        *
*     *             beom_cent_get              *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_cent_get(geom,icent,t,c,q)
      implicit none
      integer geom
      integer icent
      character*16 t
      double precision c(3)
      double precision q

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_cent_get(geom, icent, t, c, q)
!$OMP END SINGLE copyprivate(t,c,q,ok)

      beom_cent_get = ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_cent_tag              *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_cent_tag(geom,icent,t)
      implicit none
      integer geom
      integer icent
      character*16 t

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_cent_tag(geom,icent,t)
!$OMP END SINGLE copyprivate(t,ok)

      beom_cent_tag = ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_tag_to_element        *
*     *                                        *
*     ******************************************
      
*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_tag_to_element(t,symbol,element,atn,
     S     symbol_size,element_size)
      implicit none
      character*16 t
      integer symbol_size,element_size
      integer symbol_mxsize,element_mxsize
      parameter(symbol_mxsize=2,element_mxsize=16)
      character*(symbol_mxsize) symbol
      character*(element_mxsize) element
      integer atn

#include "geom.fh"

      logical ok
      if(symbol_mxsize.ne.symbol_size.and.
     E     element_mxsize.ne.element_mxsize) then
         write(6,*) ' beom_tag_to_element: wrong args'
         beom_tag_to_element = .false.
         return
      endif
!$OMP SINGLE
      ok = geom_tag_to_element(t,symbol,element,atn)
!$OMP END SINGLE copyprivate(symbol,element,atn,ok)

      beom_tag_to_element = ok
      return
      end



*     ******************************************
*     *                                        *
*     *             beom_cart_get              *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_cart_get(geom,ncent,t,c,q)
      implicit none
      integer geom
      integer ncent
      character*16 t(ncent)
      double precision c(3,ncent)
      double precision q(ncent)

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_cart_get(geom,ncent,t,c,q)
!$OMP END SINGLE copyprivate(ncent,t,c,q,ok)

      beom_cart_get = ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_cart_set              *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_cart_set(geom,ncent,t,c,q)
      implicit none
      integer geom
      integer ncent
      character*16 t(ncent)
      double precision c(3,ncent)
      double precision q(ncent)

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_cart_set(geom,ncent,t,c,q)
!$OMP END SINGLE copyprivate(ok)

      beom_cart_set = ok
      return
      end


*     ******************************************
*     *                                        *
*     *             beom_cart_coords_get       *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_cart_coords_get(geom,c)
      implicit none
      integer geom
      double precision c(3,*)

#include "geom.fh"

      beom_cart_coords_get = geom_cart_coords_get(geom,c)
      return
      end





*     ******************************************
*     *                                        *
*     *             beom_masses_get            *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_masses_get(geom,ncent,masses)
      implicit none
      integer geom
      integer ncent
      double precision masses(ncent)

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_masses_get(geom,ncent,masses)
!$OMP END SINGLE copyprivate(ncent,masses,ok)

      beom_masses_get = ok
      return
      end

*     ******************************************
*     *                                        *
*     *             beom_mass_get              *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_mass_get(geom,icent,mass)
      implicit none
      integer geom
      integer icent
      double precision mass

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_mass_get(geom,icent,mass)
!$OMP END SINGLE copyprivate(mass,ok)

      beom_mass_get = ok
      return
      end





*     ******************************************
*     *                                        *
*     *             beom_masses_set            *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_masses_set to make it thread safe.

      logical function beom_masses_set(geom,ncent,masses)
      implicit none
      integer geom
      integer ncent
      double precision masses(ncent)

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_masses_set(geom,ncent,masses)
!$OMP END SINGLE copyprivate(ok)

      beom_masses_set = ok
      return
      end


*     ******************************************
*     *                                        *
*     *             beom_vel_get               *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_cent_get to make it thread safe.

      logical function beom_vel_get(geom,vel)
      implicit none
      integer geom
      double precision vel(3,*)

#include "geom.fh"

      logical ok

!$OMP CRITICAL
      ok = geom_vel_get(geom,vel)
!$OMP END CRITICAL

      beom_vel_get = ok
      return
      end


*     ******************************************
*     *                                        *
*     *             beom_vel_set               *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_vel_set to make it thread safe.

      logical function beom_vel_set(geom,vel)
      implicit none
      integer geom
      double precision vel(3,*)

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_vel_set(geom,vel)
!$OMP END SINGLE copyprivate(ok)

      beom_vel_set = ok
      return
      end



*     ******************************************
*     *                                        *
*     *             beom_ncent                 *
*     *                                        *
*     ******************************************

*   This is a wrapper around geom_ncent to make it thread safe.

      logical function beom_ncent(geom,ncent)
      implicit none
      integer geom
      integer ncent

#include "geom.fh"

      logical ok

!$OMP SINGLE
      ok = geom_ncent(geom,ncent)
!$OMP END SINGLE copyprivate(ncent,ok)

      beom_ncent = ok
      return
      end




