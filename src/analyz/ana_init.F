      subroutine ana_init(nsai,msai,nwmi,mwmi,nwai,mwai,ltask)
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer nsai,msai,nwmi,mwmi,nwai,mwai
      logical ltask
c
c     determine me and np
c     this is needed here for on-the-fly analysis since then there
c     is no call to the analysis main task routine
c
c     input: nsai : total number of solute atoms
c            msai : dimension arrays for solute atoms on this node
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      nsa=nsai
      msa=msai
      nwm=nwmi
      mwm=mwmi
      nwa=nwai
      mwa=mwai
c
c     read rtdb if on-the-fly
c     for analysis task calls this has already been done
c
      if(.not.ltask) call ana_rtdbin()
c
c     get filename from rtdb: only node 0 is doing i/o
c
      if(me.eq.0) then
      if(.not.rtdb_cget(irtdb,'ana:project',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'arg:project',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'prep:sysnam',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'file_prefix',1,filnam))
     + call md_abort('arg_rdname: error rtdb_cget project',0)
      endif
      endif
      endif
      filana=filnam(1:index(filnam,' ')-1)//'.ana'
      open(unit=lfnana,file=filana,form='formatted',status='unknown')
      endif
c
c     allocate memory
c
      if(.not.ma_push_get(mt_int,nsa,'hbnd',l_hbnd,i_hbnd))
     + call md_abort('Failed to allocate memory for hbnd',me)
      if(.not.ma_push_get(mt_int,nsa,'isel',l_isel,i_isel))
     + call md_abort('Failed to allocate memory for isel',me)
      if(.not.ma_push_get(mt_int,nsa,'ssel',l_ssel,i_ssel))
     + call md_abort('Failed to allocate memory for ssel',me)
      if(.not.ma_push_get(mt_int,nsa,'osel',l_osel,i_osel))
     + call md_abort('Failed to allocate memory for osel',me)
c
      if(.not.ma_push_get(mt_byte,48*nsa,'tag',l_tag,i_tag))
     + call md_abort('Failed to allocate memory for tag',me)
      if(.not.ma_push_get(mt_dbl,2*nsa,'val',l_val,i_val))
     + call md_abort('Failed to allocate memory for val',me)
c
      if(.not.ma_push_get(mt_dbl,3*nsa,'xref',l_xref,i_xref))
     + call md_abort('Failed to allocate memory for xref',me)
      if(.not.ma_push_get(mt_dbl,nsa,'xrms',l_xrms,i_xrms))
     + call md_abort('Failed to allocate memory for xrms',me)
      if(.not.ma_push_get(mt_dbl,nsa,'wt',l_wt,i_wt))
     + call md_abort('Failed to allocate memory for wt',me)
      if(.not.ma_push_get(mt_dbl,nsa,'swt',l_swt,i_swt))
     + call md_abort('Failed to allocate memory for swt',me)
      if(.not.ma_push_get(mt_dbl,nsa,'owt',l_owt,i_owt))
     + call md_abort('Failed to allocate memory for owt',me)
c
c     initialize rms and wt arrays
c
      call ana_inirms(dbl_mb(i_xref),dbl_mb(i_xrms),dbl_mb(i_wt))
c
c     initialize tags
c
      call ana_initag(byte_mb(i_tag),dbl_mb(i_val))
c
c     for task calls memory needs to be set aside for atom names,
c     coordinates and global atom numbers
c     for on-the-fly analyses these come directly from the md
c     module where these arrays are distributed.
c     for task analyses these arrays are replicated over all
c     nodes
c     currently only solute coordinates are analyzed
c
      if(ltask) then
      if(.not.ma_push_get(mt_byte,16*nsa,'sgmnam',l_snam,i_snam))
     + call md_abort('Failed to allocate memory for sgmnam',me)
      if(.not.ma_push_get(mt_dbl,3*nsa,'xdat',l_xdat,i_xdat))
     + call md_abort('Failed to allocate memory for xdat',me)
      if(.not.ma_push_get(mt_dbl,nsa,'qdat',l_qdat,i_qdat))
     + call md_abort('Failed to allocate memory for qdat',me)
      if(.not.ma_push_get(mt_dbl,3*nsa,'xadat',l_xadat,i_xadat))
     + call md_abort('Failed to allocate memory for xadat',me)
      if(.not.ma_push_get(mt_dbl,3*nsa,'xp',l_xp,i_xp))
     + call md_abort('Failed to allocate memory for xp',me)
      if(.not.ma_push_get(mt_int,nsa,'idat',l_idat,i_idat))
     + call md_abort('Failed to allocate memory for idat',me)
      if(.not.ma_push_get(mt_dbl,3*mwm*mwa,'wdat',l_wdat,i_wdat))
     + call md_abort('Failed to allocate memory for wdat',me)
      if(.not.ma_push_get(mt_dbl,nwa,'qwdat',l_qwdat,i_qwdat))
     + call md_abort('Failed to allocate memory for qwdat',me)
      if(.not.ma_push_get(mt_int,mwm,'wdat',l_wsel,i_wsel))
     + call md_abort('Failed to allocate memory for wsel',me)
      if(.not.ma_push_get(mt_int,mwm,'ndxw',l_ndxw,i_ndxw))
     + call md_abort('Failed to allocate memory for ndxw',me)
      endif
c
c     set the number of analyzed frames to zero
c
      ndata=0
c
      return
      end
