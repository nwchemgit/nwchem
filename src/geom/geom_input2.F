      logical function geom_zmt_get_varinfo(geom,zvarname,zvarsign)
*
* $Id$
*
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
      character*(*) zvarname(*)
      double precision zvarsign(*)
      logical geom_check_handle
      external geom_check_handle
c
      integer i
      geom_zmt_get_varinfo = geom_check_handle(geom, 'geom_zm_getvinfo')
      if (.not. geom_zmt_get_varinfo) return
c
      do i = 1, zmt_nzvar(geom)
         zvarname(i) = zmt_varname(i,geom)
         zvarsign(i) = zmt_varsign(i,geom)
      enddo
c
      end
      logical function geom_zmt_get_source(geom,source)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
      character*(*) source
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_get_source = geom_check_handle(geom, 'geom_zmt_source')
      if (.not. geom_zmt_get_source) return
c
      source = zmt_source(geom)
c
      end
      logical function geom_zmt_put_data(geom, nizmat, izmat,
     $     nzfrz, izfrz, izfrz_val, nzvar, zvarname, zvarsign, source)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom, nizmat, izmat(*), nzfrz, izfrz(*), nzvar
      character*(*) zvarname(*), source
      double precision zvarsign(*), izfrz_val(*)
c
      integer i
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_put_data = geom_check_handle(geom, 'geom_zm_put_data')
      if (.not. geom_zmt_put_data) return
c
      zmt_nizmat(geom) = nizmat
      do i = 1, nizmat
         zmt_izmat(i,geom) = izmat(i)
      enddo
      zmt_nzfrz(geom)  = nzfrz
      do i = 1, nzfrz
         zmt_izfrz(i,geom) = izfrz(i)
         zmt_izfrz_val(i,geom) = izfrz_val(i)
      enddo
      zmt_nzvar(geom) = nzvar
      do i = 1, nzvar
         zmt_varname(i,geom) = zvarname(i)
         zmt_varsign(i,geom) = zvarsign(i)
      enddo
      zmt_source(geom) = source
c
      end
      logical function geom_zmt_get_nzvar(geom,nzvar)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom
      integer nzvar            
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_get_nzvar = geom_check_handle(geom, 'geom_zm_get_nzvar')
      if (.not. geom_zmt_get_nzvar) return
c
      nzvar = zmt_nzvar(geom)
c
      end
      logical function geom_zmt_get_nizmat(geom,nizmat)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom
      integer nizmat            
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_get_nizmat = geom_check_handle(geom, 'geom_zm_g_ni')
      if (.not. geom_zmt_get_nizmat) return
c
      nizmat = zmt_nizmat(geom)
c
      end
      logical function geom_zmt_get_izmat(geom,
     1                                    izmat,nizmat)         
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom
      integer nizmat             
      integer izmat             
      integer i
      dimension izmat(*)
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_get_izmat = geom_check_handle(geom, 'geom_zm_g_i')
      if (.not. geom_zmt_get_izmat) return
c
      do i = 1, nizmat
         izmat(i) = zmt_izmat(i,geom)
      end do
c
      end
      logical function geom_zmt_get_izfrz(geom,izfrz)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom
      integer izfrz
      integer i
      dimension izfrz(*)
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_get_izfrz = geom_check_handle(geom, 'geom_us_g_iz')
      if (.not. geom_zmt_get_izfrz) return
c
      do i = 1, zmt_nzfrz(geom)
         izfrz(i) = zmt_izfrz(i,geom)
      enddo
c
      end
      logical function geom_zmt_get_izfrz_val(geom,izfrz_val)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom
      double precision izfrz_val(*)
      integer i
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_get_izfrz_val = geom_check_handle(geom, 'geom_us_g_iz')
      if (.not. geom_zmt_get_izfrz_val) return
c
      do i = 1, zmt_nzfrz(geom)
         izfrz_val(i) = zmt_izfrz_val(i,geom)
      enddo
c
      end
      logical function geom_zmt_get_nzfrz(geom,nzfrz)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom
      integer nzfrz 
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_get_nzfrz = geom_check_handle(geom, 'geom_us_g_nz')
      if (.not. geom_zmt_get_nzfrz) return
c
      nzfrz = zmt_nzfrz(geom)
c
      end
