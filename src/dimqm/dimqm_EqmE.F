      subroutine dimqm_EqmE(rtdb, g_dens, geom, basis,
     $                      fld, xyz)
c  
c      Author: Justin Moore
c 
c      Called from: dimqm_main.F dimqm_addop.F dimqm_indDipoles.F
c                   dimqm_Dens2Dipoles.F
c
c      Calls: hnd_elfcon.F
c
c      Subroutine to calculate the electric field due to the QM
c      electrons
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
#include "geom.fh"
c    Input Variables
      integer rtdb
      integer g_dens
      integer geom
      integer basis
      double precision fld(3, nDIM)
      double precision xyz(3, nDIM)
      double precision scrn(nDIM)
      integer i, j, n
      logical stat
c
      integer nNuc
      integer l_qmxyz, k_qmxyz
c      double precision dimfactors(nDIM)
      double precision temp
      double precision pref
      double precision ox, oy, oz
      double precision ac, aq, aa

c      ac = 1.0d0
c      aq = 1.0d0
c      aa = (ac*aq)/(ac + aq)
c      pref = sqrt(PI/(ac+aq))
      
      if(.not. geom_ncent(geom, nNuc))
     $  call errquit("EqmE: geom_ncent", 1, GEOM_ERR)
      if(.not. ma_push_get(mt_dbl,3*nNuc,'QM xyz',l_qmxyz,k_qmxyz))
     $  call errquit("EqmE: malloc k_qmxyz", 1, MA_ERR)
      if(.not. geom_cart_coords_get(geom, dbl_mb(k_qmxyz)))
     $  call errquit("EqmE: qm coords get failed", 1, GEOM_ERR)
c      dimfactors = 0.0d0
c      do i = 1, nDIM
c        temp = 0.0d0
c        write(luout,*) "DIM"
c        write(luout,*) xyz(:, i)
c        do j = 1, nNuc
c          write(luout,*) "QM"
c          write(luout,*) dbl_mb(k_qmxyz+3*(j-1)),
c     $                   dbl_mb(k_qmxyz+3*(j-1)+1),
c     $                   dbl_mb(k_qmxyz+3*(j-1)+2)
c          ox = pref*exp(-aa*(dbl_mb(k_qmxyz+3*(j-1))-xyz(1,i))**2)
c          oy = pref*exp(-aa*(dbl_mb(k_qmxyz+3*(j-1)+1)-xyz(2,i))**2)
c          oz = pref*exp(-aa*(dbl_mb(k_qmxyz+3*(j-1)+2)-xyz(3,i))**2)
c          temp = temp + ox*oy*oz
c        end do
c        dimfactors(i) = temp
c      end do
c      write(luout,*) "factors:"
cc      write(luout,*) dimfactors
      call hnd_elfcon(basis, geom, g_dens, xyz, nDIM, fld, 1)
c      call dim_elfcon(basis, geom, g_dens, xyz, nDIM, fld, 1)
      call ga_sync()
c      call ga_dgop(2112, fld, 3*nDIM, '+')
      if(.not.ma_pop_stack(l_qmxyz))
     %  call errquit('problem destroying qmxyz in EqmE', 1, MA_ERR)
      end subroutine dimqm_EqmE
