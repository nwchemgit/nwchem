      subroutine cphf_int_debug(geom, basis, g_dens)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer g_dens, geom, basis
      integer g_x, g_y, g_z
      logical ga_create_atom_blocked
      external ga_create_atom_blocked
c
      g_x = ga_create_atom_blocked(geom, basis, 'cphf: dip x')
      g_y = ga_create_atom_blocked(geom, basis, 'cphf: dip y')
      g_z = ga_create_atom_blocked(geom, basis, 'cphf: dip z')
c
      call int_dip_ga(basis, basis, g_x, g_y, g_z)
c
      write(6,*) ga_ddot(g_dens, g_x)
      write(6,*) ga_ddot(g_dens, g_y)
      write(6,*) ga_ddot(g_dens, g_z)
c
      if (.not. ga_destroy(g_x)) call errquit('cphf',0, GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit('cphf',0, GA_ERR)
      if (.not. ga_destroy(g_z)) call errquit('cphf',0, GA_ERR)
c
      end

