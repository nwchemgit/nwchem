      subroutine smd_param_init()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_data.fh"
#include "stdio.fh"
#include "msgids.fh"
c     
      character*32 pname
      character*80 tag
      character*255 filename
      integer na
      integer i_ta,i_tt,i_tr,i_q,i_m
      logical result
      integer db, master
c
c
      master = 0
      pname = "smd_param_init"
#ifdef SMD_DEBUG
       write(*,*) "in "//pname
#endif

c
      call smd_srtdb_get_handle(db)
c
c     get total number of atoms in parameter file
c     -------------------------------------------
      call smd_paramfile_natoms(na)
#ifdef SMD_DEBUG
      write(*,*) "number of atoms in parameter space",na
#endif
      if(na.le.0)
     >  call errquit(
     >       pname//'no parameters fond',0, 0)
c
c     create parameter data structures
c     ---------------------------
      tag = "param:atomname"
      if(.not.smd_data_create_get(db,tag,16*na,MT_BYTE,i_ta)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if
      tag = "param:typename"
      if(.not.smd_data_create_get(db,tag,16*na,MT_BYTE,i_tt)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if
      tag = "param:charge"
      if(.not.smd_data_create_get(db,tag,na,MT_DBL,i_q)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if
      tag = "param:mass"
      if(.not.smd_data_create_get(db,tag,na,MT_DBL,i_m)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if
      tag = "param:resname"
      if(.not.smd_data_create_get(db,tag,16*na,MT_BYTE,i_tr)) then
        call smd_message_flush(luout)
        call errquit(" ",0,0)
      end if

      return
      if(ga_nodeid().eq.master) then
      call smd_paramfile_read_atom(na,
     +                         byte_mb(i_ta),
     +                         byte_mb(i_tt),
     +                         byte_mb(i_tr),
     +                         dbl_mb(i_m),
     +                         dbl_mb(i_q))

      end if

      call ga_brdcst(msg_smd,byte_mb(i_ta),
     >      16*na,master)
      call ga_brdcst(msg_smd,byte_mb(i_tt),
     >      16*na,master)
      call ga_brdcst(msg_smd,byte_mb(i_tr),
     >      16*na,master)
      call ga_brdcst(msg_smd,dbl_mb(i_m),
     >      na*ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_brdcst(msg_smd,dbl_mb(i_q),
     >      na*ma_sizeof(mt_dbl,1,mt_byte),master)

      call ga_sync()

#ifdef SMD_DEBUG
       write(*,*) "out "//pname
#endif

      return
      end

c      subroutine smd_param_ntot(namespace,ntot)
c      implicit none
c#include "errquit.fh"
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "util.fh"
c#include "global.fh"
cc     
c      character*(*) namespace
c      integer ntot
cc
c      character*72 tag
c      character*30 pname
c      logical result
c
c      pname = "smd_param_ntot"
c      tag = "param:charge"
c      call smd_data_get_size(namespace,tag,ntot,result)
c      if(.not. result)  
c     >  call errquit(
c     >       pname//'error getting size for '//tag,0, 0)
c      ntot = ntot
c
c      return
c      end
c
c $Id$
