      subroutine tce_aod1(rtdb,g_aod1)
c
c $Id: tce_aod1.F,v 1.1 2003-07-25 17:56:56 sohirata Exp $
c
c     Returns a GA handle for AO Fock matrix.
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "tce.fh"
      integer rtdb               ! Run-time database
      integer g_aod1(3)          ! GA handle for AO dipole matrices
c
c     ===============
c     Allocate g_aod1
c     ===============
c
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole X',
     1  -1,-1,g_aod1(1))) call errquit('tce_aod1: GA problem',0)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole Y',
     1  -1,-1,g_aod1(2))) call errquit('tce_aod1: GA problem',1)
      if (.not.ga_create(mt_dbl,nbf,nbf,'AO dipole Z',
     1  -1,-1,g_aod1(3))) call errquit('tce_aod1: GA problem',2)
c
c     ====================
c     Make dipole matrices
c     ====================
c
      call int_dip_ga(ao_bas_han,ao_bas_han,
     1  g_aod1(1),g_aod1(2),g_aod1(3))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(1))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(2))
      if (util_print('aod1',print_debug)) call ga_print(g_aod1(3))
      return
      end
