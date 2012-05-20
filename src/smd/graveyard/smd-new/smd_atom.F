      subroutine smd_atom_init()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_data.fh"
#include "stdio.fh"
#include "msgids.fh"
c     
      character*32 pname
      character*80 tag
      integer nt
      integer i_t,i_tr,i_ir,i_qf
      integer i_iconst
      logical result
      integer srtdb
      integer master
      logical mode
c
      master = 0
      pname = "smd_atom_init"
#ifdef SMD_DEBUG
      write(*,*) "in "//pname
#endif
c
c
c     get total number of atoms from pdb file
c     ---------------------------------------
      call smd_coordfile_natoms(nt)
#ifdef SMD_DEBUG
      write(*,*) "number of atoms",nt
#endif

      if(nt.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
       call smd_srtdb_get_handle(srtdb)
       if(.not.smd_data_create(srtdb,"atom:resid",nt,MT_INT)) then
         call smd_message_flush(luout)
         call errquit("",0,0)
       end if

       if(.not.smd_data_create(srtdb,"atom:name",16*nt,MT_BYTE)) then
         call smd_message_flush(luout)
         call errquit("",0,0)
       end if

       if(.not.smd_data_create(srtdb,"atom:resname",16*nt,MT_BYTE)) then
         call smd_message_flush(luout)
         call errquit("",0,0)
       end if

c
c     fill in names from pdb file
c     ---------------------------
      tag = "atom:name"
      if(.not.smd_data_get_index(srtdb,tag,i_t)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if

      tag = "atom:resname"
      if(.not.smd_data_get_index(srtdb,tag,i_tr)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if

      tag = "atom:resid"
      if(.not.smd_data_get_index(srtdb,tag,i_ir)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if
c
#ifdef SMD_DEBUG
       write(*,*) "0 smd_coordfile_read_atomres "//pname
#endif
      if(ga_nodeid().eq.master) then 
        call smd_coordfile_read_atomres(nt,
     +                           byte_mb(i_t),
     +                           byte_mb(i_tr),
     +                           int_mb(i_ir))

#ifdef SMD_DEBUG
       write(*,*) "1 smd_coordfile_read_atomres "//pname
#endif
      end if

      call ga_brdcst(msg_smd,byte_mb(i_t),
     >      16*nt,master)
      call ga_brdcst(msg_smd,byte_mb(i_tr),
     >      16*nt,master)
      call ga_brdcst(msg_smd,int_mb(i_ir),
     >      nt*ma_sizeof(mt_int,1,mt_byte),master)

      call ga_sync()

c      tag = "atom:name"
c      call smd_data_print(srtdb,tag,luout)
#ifdef SMD_DEBUG
       write(*,*) "out "//pname
#endif

      return
      end

      subroutine smd_atom_end()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "smd_data.fh"
#include "stdio.fh"
c     
      integer srtdb
c
       call smd_srtdb_get_handle(srtdb)
       if(.not.smd_data_destroy(srtdb,"atom:resname")) then
         call smd_message_flush(luout)
         call errquit("",0,0)
       end if
       if(.not.smd_data_destroy(srtdb,"atom:name")) then
         call smd_message_flush(luout)
         call errquit("",0,0)
       end if
      if(.not.smd_data_destroy(srtdb,"atom:resid")) then
         call smd_message_flush(luout)
         call errquit("",0,0)
       end if
c
      return
      end

c $Id$
