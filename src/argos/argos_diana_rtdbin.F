      subroutine argos_diana_rtdbin
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
c      if(me.eq.0) then
c      if(.not.rtdb_cget(irtdb,'ana:filcmd',1,filcmd))
c     + call md_abort('Error rtdb_get filcmd',0)
c      endif
c
c      integer ibroad(10)
cc
cc     input options from rtdb
cc
c      if(me.eq.0) then
cc
cc     reference input file
cc
c      if(.not.rtdb_cget(irtdb,'ana:filref',1,filref))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
cc
cc     trajectory input file
cc
c      if(.not.rtdb_cget(irtdb,'ana:filtrj',1,filtrj))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
cc
cc     copy trajectory
cc
c      if(.not.rtdb_cget(irtdb,'ana:filcop',1,filcop))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_cget(irtdb,'ana:fmtcop',1,fmtcop))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_get(irtdb,'ana:lcop',mt_int,1,lcop))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_get(irtdb,'ana:mcopf',mt_int,1,mcopf))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
cc
cc     superimpose trajectory
cc
c      if(.not.rtdb_cget(irtdb,'ana:filsup',1,filsup))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_cget(irtdb,'ana:fmtsup',1,fmtsup))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_get(irtdb,'ana:lsup',mt_int,1,lsup))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_get(irtdb,'ana:msupf',mt_int,1,msupf))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
cc
c      if(.not.rtdb_get(irtdb,'ana:timoff',mt_dbl,1,timoff))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
cc
c      if(.not.rtdb_get(irtdb,'ana:ifrfr',mt_int,1,ifrfr))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_get(irtdb,'ana:ifrto',mt_int,1,ifrto))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
cc
c      if(.not.rtdb_get(irtdb,'ana:ifrst',mt_int,1,ifrst))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
c      if(.not.rtdb_get(irtdb,'ana:ilast',mt_int,1,ilast))
c     + call md_abort('argos_diana_input: rtdb_get failed',0)
cc
cc     root mean square deviation
cc
c      if(.not.rtdb_get(irtdb,'ana:irmsd',mt_int,1,irmsd))
c     + call md_abort('argos_diana_rtdbin: rtdb_get failed',1)
cc
c      if(.not.rtdb_get(irtdb,'ana:lplt',mt_int,1,lplt))
c     + call md_abort('argos_diana_rtdbin: rtdb_get failed',1)
c      if(.not.rtdb_get(irtdb,'ana:spac',mt_dbl,1,spac))
c     + call md_abort('argos_diana_rtdbin: rtdb_get failed',0)
c      if(.not.rtdb_get(irtdb,'ana:rcut',mt_dbl,1,rcut))
c     + call md_abort('argos_diana_rtdbin: rtdb_get failed',0)
c      if(.not.rtdb_cget(irtdb,'ana:filchg',1,filchg))
c     + call md_abort('argos_diana_rtdbin: rtdb_get failed',0)
c      if(.not.rtdb_cget(irtdb,'ana:filplt',1,filplt))
c     + call md_abort('argos_diana_rtdbin: rtdb_get failed',0)
c      if(.not.rtdb_cget(irtdb,'ana:filxyz',1,filxyz))
c     + call md_abort('argos_diana_rtdbin: rtdb_get failed',0)
cc
c      endif
cc
cc     broadcast options to all nodes
cc
c      if(np.gt.1) then
c      ibroad(1)=irmsd
c      ibroad(2)=lcop
c      ibroad(3)=lsup
c      ibroad(4)=ifrfr
c      ibroad(5)=ifrto
c      ibroad(6)=lplt
c      call ga_brdcst(mag_d00,ibroad,6*ma_sizeof(mt_int,1,mt_byte),0)
c      irmsd=ibroad(1)
c      lcop=ibroad(2)
c      lsup=ibroad(3)
c      ifrfr=ibroad(4)
c      ifrto=ibroad(5)
c      lplt=ibroad(6)
c      endif
cc
cc     logical file numbers
cc
c      lfnana=75
c      lfnref=76
c      lfntrj=77
c      lfncop=78
c      lfnsup=79
c      lfnchg=80
c      lfnplt=81
c      lfnxyz=82
cc
cc     set flags on what to do
cc
c      ltrj=irmsd
cc
c      if(lcop.gt.0) ltrj=ltrj+1
c      if(lsup.gt.0) ltrj=ltrj+1
c
      return
      end
