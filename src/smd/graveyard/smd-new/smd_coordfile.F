      subroutine smd_coordfile_natoms(nt)
      implicit none
#include "errquit.fh"
#include "util.fh"
      integer nt
c
      character*(nw_max_path_len) filename
      logical result
      character*30 pname
c
      pname = "smd_coordfile_natoms"
c
      call smd_srtdb_get_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)

c
      call smd_pdb_natoms(filename,nt)

      end

      subroutine smd_coordfile_read_coords(nt,c)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(nt,3)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname
      character*(nw_max_path_len) filename
      logical result


      pname = "smd_coordfile_read_coords"
c
c
      call smd_srtdb_get_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)

c
      call smd_pdb_read_coords(filename,nt,c)

      return

      end

      subroutine smd_coordfile_read_atomres(nt,ta,tr,ir)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "global.fh"
      integer nt
      character*16 ta(nt)
      character*16 tr(nt)
      integer ir(nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname
      character*(nw_max_path_len) filename
      logical result
c
      pname = "smd_coordfile_read_atomres"
c
c
      call smd_srtdb_get_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)
c
      call smd_pdb_read_atomres(filename,nt,ta,tr,ir)

      return

      end

c $Id$
