      subroutine ana_init(nsai,msai)
c
c $Id: ana_init.F,v 1.1 1998-12-29 18:01:02 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer nsai,msai,indx
c
      me=ga_nodeid()
      np=ga_nnodes()
c
c     get filename from rtdb
c
      if(me.eq.0) then
      if(.not.rtdb_cget(irtdb,'ana:project',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'arg:project',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'prep:sysnam',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'file_prefix',1,filnam))
     + call errquit('arg_rdname: error rtdb_cget project',0)
      endif
      endif
      endif
      lfnana=75
      filana=filnam(1:index(filnam,' ')-1)//'.ana'
      open(unit=lfnana,file=filana,form='formatted',status='unknown')
      endif
c
      nsa=nsai
      msa=msai
c
      if(me.eq.0) then
      write(*,'(a,i5)') 'Number of solute atoms',nsa
      write(*,'(a,i5)') 'Dimension solyte array',msa
      endif
c
      call ana_rtdbin()
c
      if(me.eq.0) then
      write(*,'(a,i5)') 'Input option irmsd    ',irmsd
      endif
c
      if(irmsd.gt.0) then
      if(.not.ma_push_get(mt_dbl,3*nsa,'xref',l_xref,i_xref))
     + call errquit('Failed to allocate memory for xref',me)
      if(.not.ma_push_get(mt_dbl,nsa,'xrms',l_xrms,i_xrms))
     + call errquit('Failed to allocate memory for xrms',me)
      if(.not.ma_push_get(mt_dbl,nsa,'wt',l_wt,i_wt))
     + call errquit('Failed to allocate memory for wt',me)
      call ana_inirms(dbl_mb(i_xref),dbl_mb(i_xrms),dbl_mb(i_wt))
      endif
c
      ndata=0
c
      return
      end
