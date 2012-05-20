      subroutine smd_coordfile_natoms_par(nt)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer nt
c
      character*(33) filename
      logical result
      character*30 pname
      logical rtdb_mode
      integer master
c
      pname = "smd_coordfile_natoms"
c
      master = 0
      rtdb_mode = rtdb_parallel(.false.)
      write(*,*) "in",pname
      call util_flush(6)
      if (ga_nodeid().eq.master) then
        call smd_coordfile(filename,result)
        if(.not.result)
     >    call errquit(
     >         pname//'no coord file found',0, 0)

c
        call smd_pdb_natoms(filename,nt)
      end if
      call ga_brdcst(msg_smd,nt,
     >     ma_sizeof(mt_int,1,mt_byte),master)

      call ga_sync()
      rtdb_mode = rtdb_parallel(rtdb_mode)
      return

      end

      subroutine smd_coordfile_read_coords(nt,c)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
      integer nt
      double precision c(nt,3)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname
      character*(33) filename
      logical result
c
      pname = "smd_coordfile_read_coords"
c
c
      call smd_coordfile_read_coords_seq(nt,c)

      return

      end
      subroutine smd_coordfile_read_coords_par(nt,c)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
      integer nt
      double precision c(nt,3)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname
      character*(33) filename
      logical result
      logical rtdb_mode
      integer master
c
      pname = "smd_coordfile_read_coords"
c
c
      master = 0
      rtdb_mode = rtdb_parallel(.false.)
      if (ga_nodeid().eq.master) then
        call smd_coordfile_read_coords_seq(nt,c)
      end if

      call ga_brdcst(msg_smd,c,
     >     3*nt*ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_sync()
      rtdb_mode = rtdb_parallel(rtdb_mode)

      return

      end

      subroutine smd_coordfile_read_atomres_seq(nt,ta,tr,ir)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "global.fh"
#include "smd_const_data.fh"
      integer nt
      character*(*) ta(nt)
      character*(*) tr(nt)
      integer ir(nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname
      character*(33) filename
      logical result
c
      pname = "smd_coordfile_read_atomres"
c
c
      call smd_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)
c
      call smd_pdb_read_atomres(filename,nt,ta,tr,ir)

      return

      end
      subroutine smd_coordfile_read_atomres(nt,ta,tr,ir)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "global.fh"
#include "smd_const_data.fh"
      integer nt
      character*(smd_string_size) ta(nt)
      character*(smd_string_size) tr(nt)
      integer ir(nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname
      character*(33) filename
      logical result
c
      pname = "smd_coordfile_read_atomres"
c
c
      call smd_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)
c
      call smd_pdb_read_atomres(filename,nt,ta,tr,ir)

      return

      end

      subroutine smd_coordfile_natoms_seq(nt)
      implicit none
#include "errquit.fh"
#include "util.fh"
      integer nt
c
      character*(33) filename
      logical result
      character*30 pname
c
      pname = "smd_coordfile_natoms"
c
      call smd_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)

c
      call smd_pdb_natoms(filename,nt)

      end

      subroutine smd_coordfile_natoms(nt)
      implicit none
#include "errquit.fh"
#include "util.fh"
      integer nt
c
      character*(33) filename
      logical result
      character*30 pname
c
      pname = "smd_coordfile_natoms"
c
      call smd_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)

c
      call smd_pdb_natoms(filename,nt)

      end

      subroutine smd_coordfile_read_coords_seq(nt,c)
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
      character*(33) filename
      logical result


      pname = "smd_coordfile_read_coords"
c
c
      call smd_coordfile(filename,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coord file found',0, 0)

c
      call smd_pdb_read_coords(filename,nt,c)

      return

      end

      subroutine smd_coordfile(filename,result)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "global.fh"

c
      character*(*) filename
      logical result
c
      character*30 pname

      pname = "smd_coordfile"

      result = .true.
      call smd_rtdb_get_string("smd:coordfile",1,
     >                           filename,result)
      write(*,*) "filename",filename

      end

c $Id$
