      logical function geom_zmt_put_nizmat(geom,nizmat)
      implicit none
#include "geomP.fh"
c
      integer geom
      integer nizmat
      logical geom_check_handle
      external geom_check_handle
c
      geom_zmt_put_nizmat = geom_check_handle(geom, 'geom_zm_p_ni')
      if (.not. geom_zmt_put_nizmat) return
c
      zmt_nizmat(geom) = nizmat
c
      geom_zmt_put_nizmat=.true.
      return
      end
      logical function geom_zmt_get_nizmat(geom,nizmat)
      implicit none
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
c
      nizmat = zmt_nizmat(geom)
c
      geom_zmt_get_nizmat=.true.
      return
      end
      logical function geom_zmt_put_izmat(geom,
     1                                    izmat,nizmat)          
      implicit none
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
      geom_zmt_put_izmat = geom_check_handle(geom, 'geom_zm_p_i')
      if (.not. geom_zmt_put_izmat) return
c
c
      do i = 1, nizmat
         zmt_izmat(i,geom) = izmat(i)
      end do
c
      geom_zmt_put_izmat=.true.
      return
      end
      logical function geom_zmt_get_izmat(geom,
     1                                    izmat,nizmat)         
      implicit none
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
      geom_zmt_get_izmat=.true.
      return
      end
      logical function geom_usr_put_nizmat(geom,nizmat)
      implicit none
#include "geomP.fh"
c
      integer geom
      integer nizmat
      logical geom_check_handle
      external geom_check_handle
c
      geom_usr_put_nizmat = geom_check_handle(geom, 'geom_us_p_ni')
      if (.not. geom_usr_put_nizmat) return
c
c
      usr_nizmat(geom) = nizmat
c
      geom_usr_put_nizmat=.true.
      return
      end
      logical function geom_usr_get_nizmat(geom,nizmat)
      implicit none
#include "geomP.fh"
c
      integer geom
      integer nizmat            
      logical geom_check_handle
      external geom_check_handle
c
      geom_usr_get_nizmat = geom_check_handle(geom, 'geom_us_g_ni')
      if (.not. geom_usr_get_nizmat) return
c
c
      nizmat = usr_nizmat(geom)
c
      geom_usr_get_nizmat=.true.
      return
      end
      logical function geom_usr_put_izmat(geom,
     1                                    izmat,nizmat)
      implicit none
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
      geom_usr_put_izmat = geom_check_handle(geom, 'geom_us_p_i')
      if (.not. geom_usr_put_izmat) return
c
c
      do i = 1, nizmat
         usr_izmat(i,geom) = izmat(i)
      end do
c
      geom_usr_put_izmat=.true.
      return
      end
      logical function geom_usr_get_izmat(geom,
     1                                    izmat,nizmat)
      implicit none
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
      geom_usr_get_izmat = geom_check_handle(geom, 'geom_us_g_i')
      if (.not. geom_usr_get_izmat) return
c
c
      do i = 1, nizmat
         izmat(i) = usr_izmat(i,geom)
      end do
c
      geom_usr_get_izmat=.true.
      return
      end
      logical function geom_usr_put_izfrz(geom,
     1                                    izfrz,nzfrz)
      implicit none
#include "geomP.fh"
c
      integer geom
      integer nzfrz 
      integer izfrz
      integer i
      dimension izfrz(*)
      logical geom_check_handle
      external geom_check_handle
c
      geom_usr_put_izfrz = geom_check_handle(geom, 'geom_us_p_iz')
      if (.not. geom_usr_put_izfrz) return
c
      do i = 1, max_coor
         usr_izfrz(i,geom) = 0 
      enddo
c
      if(nzfrz.gt.0) then
         do i = 1, nzfrz
            usr_izfrz(i,geom) = izfrz(i)
         enddo
      endif
c
      geom_usr_put_izfrz=.true.
      return
      end
      logical function geom_usr_get_izfrz(geom,
     1                                    izfrz,nzfrz)
      implicit none
#include "geomP.fh"
c
      integer geom
      integer nzfrz
      integer izfrz
      integer i
      dimension izfrz(*)
      logical geom_check_handle
      external geom_check_handle
c
      geom_usr_get_izfrz = geom_check_handle(geom, 'geom_us_g_iz')
      if (.not. geom_usr_get_izfrz) return
c
      do i = 1, max_coor
         izfrz(i) = 0 
      enddo
c
      if(nzfrz.gt.0) then
         do i = 1, nzfrz
            izfrz(i) = usr_izfrz(i,geom)
         enddo
      endif
c
      geom_usr_get_izfrz=.true.
      return
      end
      logical function geom_usr_put_nzfrz(geom,nzfrz)
      implicit none
#include "geomP.fh"
c
      integer geom
      integer nzfrz
      logical geom_check_handle
      external geom_check_handle
c
      geom_usr_put_nzfrz = geom_check_handle(geom, 'geom_us_p_nz')
      if (.not. geom_usr_put_nzfrz) return
c
c
      usr_nzfrz(geom) = nzfrz
c
      geom_usr_put_nzfrz=.true.
      return
      end
      logical function geom_usr_get_nzfrz(geom,nzfrz)
      implicit none
#include "geomP.fh"
c
      integer geom
      integer nzfrz 
      logical geom_check_handle
      external geom_check_handle
c
      geom_usr_get_nzfrz = geom_check_handle(geom, 'geom_us_g_nz')
      if (.not. geom_usr_get_nzfrz) return
c
c
      nzfrz = usr_nzfrz(geom)
c
      geom_usr_get_nzfrz=.true.
      return
      end
      logical function geom_zmt_put_zmat(geom,zmat,nzmat)          
      implicit none
c     
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer geom
      integer nzmat
      integer i
      double precision  zmt_zmat
      common/hnd_zmtsav/zmt_zmat(mxcoor)
      double precision      zmat(*)
c
      do i = 1, nzmat
         zmt_zmat(i) = zmat(i)
      end do
c
      geom_zmt_put_zmat=.true.
      return
      end
      logical function geom_zmt_get_zmat(geom,zmat,nzmat)         
      implicit none
c     
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer geom
      integer nzmat
      integer i
      double precision  zmt_zmat
      common/hnd_zmtsav/zmt_zmat(mxcoor)
      double precision      zmat(*)
c
      do i = 1, nzmat
         zmat(i) = zmt_zmat(i)
      end do
c
      geom_zmt_get_zmat=.true.
      return
      end
