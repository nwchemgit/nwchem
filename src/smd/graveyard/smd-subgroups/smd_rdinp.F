      subroutine smd_rdinp(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer rtdb
c
c      character*32 pname
c      character*255 filename
cc
c      pname = "smd_rdinp"
cc
c      write(*,*) "in "//pname
cc
cc     read and allocate atom structures
cc     ---------------------------------
c      call smd_atom_init(namespace,rtdb)
cc
cc     read and allocate atom-specific parameters
cc     ------------------------------------------
c      call smd_param_init(namespace,rtdb)
cc
cc     read and allocate lattice parameters
cc     ------------------------------------
c      call smd_lat_init(namespace,rtdb)
cc
cc     read and allocate vdw parameters
cc     --------------------------------
c      call smd_pot_init(namespace,namespace,rtdb)
      return
      end

c $Id$
