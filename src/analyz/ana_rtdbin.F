      subroutine ana_rtdbin
c
c $Id: ana_rtdbin.F,v 1.3 1999-09-20 20:15:39 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer ibroad(10)
c
c     input options from rtdb
c
      if(me.eq.0) then
c
c     reference input file
c
      if(.not.rtdb_cget(irtdb,'ana:filref',1,filref))
     + call errquit('ana_input: rtdb_get failed',0)
c
c     trajectory input file
c
      if(.not.rtdb_cget(irtdb,'ana:filtrj',1,filtrj))
     + call errquit('ana_input: rtdb_get failed',0)
c
c     copy trajectory
c
      if(.not.rtdb_cget(irtdb,'ana:filcop',1,filcop))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_cget(irtdb,'ana:fmtcop',1,fmtcop))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'ana:lcop',mt_int,1,lcop))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'ana:mcopf',mt_int,1,mcopf))
     + call errquit('ana_input: rtdb_get failed',0)
c
c     superimpose trajectory
c
      if(.not.rtdb_cget(irtdb,'ana:filsup',1,filsup))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_cget(irtdb,'ana:fmtsup',1,fmtsup))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'ana:lsup',mt_int,1,lsup))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'ana:msupf',mt_int,1,msupf))
     + call errquit('ana_input: rtdb_get failed',0)
c
      if(.not.rtdb_get(irtdb,'ana:timoff',mt_dbl,1,timoff))
     + call errquit('ana_input: rtdb_get failed',0)
c
      if(.not.rtdb_get(irtdb,'ana:ifrfr',mt_int,1,ifrfr))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'ana:ifrto',mt_int,1,ifrto))
     + call errquit('ana_input: rtdb_get failed',0)
c
      if(.not.rtdb_get(irtdb,'ana:ifrst',mt_int,1,ifrst))
     + call errquit('ana_input: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'ana:ilast',mt_int,1,ilast))
     + call errquit('ana_input: rtdb_get failed',0)
c
c     root mean square deviation
c
      if(.not.rtdb_get(irtdb,'ana:irmsd',mt_int,1,irmsd))
     + call errquit('ana_rtdbin: rtdb_get failed',1)
c
      if(.not.rtdb_get(irtdb,'ana:lplt',mt_int,1,lplt))
     + call errquit('ana_rtdbin: rtdb_get failed',1)
      if(.not.rtdb_get(irtdb,'ana:spac',mt_dbl,1,spac))
     + call errquit('ana_rtdbin: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'ana:rcut',mt_dbl,1,rcut))
     + call errquit('ana_rtdbin: rtdb_get failed',0)
      if(.not.rtdb_cget(irtdb,'ana:filchg',1,filchg))
     + call errquit('ana_rtdbin: rtdb_get failed',0)
      if(.not.rtdb_cget(irtdb,'ana:filplt',1,filplt))
     + call errquit('ana_rtdbin: rtdb_get failed',0)
      if(.not.rtdb_cget(irtdb,'ana:filxyz',1,filxyz))
     + call errquit('ana_rtdbin: rtdb_get failed',0)
c
      endif
c
c     broadcast options to all nodes
c
      if(np.gt.1) then
      ibroad(1)=irmsd
      ibroad(2)=lcop
      ibroad(3)=lsup
      ibroad(4)=ifrfr
      ibroad(5)=ifrto
      ibroad(6)=lplt
      call ga_brdcst(mag_d00,ibroad,6*ma_sizeof(mt_int,1,mt_byte),0)
      irmsd=ibroad(1)
      lcop=ibroad(2)
      lsup=ibroad(3)
      ifrfr=ibroad(4)
      ifrto=ibroad(5)
      lplt=ibroad(6)
      endif
c
c     logical file numbers
c
      lfnana=75
      lfnref=76
      lfntrj=77
      lfncop=78
      lfnsup=79
      lfnchg=80
      lfnplt=81
      lfnxyz=82
c
c     set flags on what to do
c
      ltrj=irmsd
c
      if(lcop.gt.0) ltrj=ltrj+1
      if(lsup.gt.0) ltrj=ltrj+1
c
      return
      end
