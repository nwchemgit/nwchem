      logical function pre_rtdbin(irtdb,ffield,
     + dir_s,dir_x,dir_u,dir_t,source,sysnam,calc,slvnam,slvmdl,
     + newtop,newseq,newrst,mcount,ncount,icount,mgrid,gdist,mnoe,
     + maxscf,qscale,
     + altloc,chain,icyren,iconst,model,nxlnk,mdold,ignore,scount,
     + mfract,nfract,ifract)
c
c $Id: pre_rtdbin.F,v 1.20 2000-04-26 20:32:20 d3j191 Exp $
c
c     function to read input for prepare module from rtdb
c
      implicit none
c
      logical util_nwchemrc_get
      external util_nwchemrc_get
c
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer irtdb,mcount,ncount
      character*3 slvnam
      character*80 source
      character*80 ffield,sysnam,calc
      character*10 slvmdl
      character*255 dir_s,dir_x,dir_u,dir_t
      integer newtop,newseq,newrst,icount(mcount),icyren,ignore
      integer len,ndx,mgrid,mnoe,maxscf,model,nxlnk,mdold
      integer mfract,nfract,iconst
      integer ifract(mfract)
      real*8 gdist,qscale
      character*1 altloc,chain
      character*4 scount(mcount)
c
      character*255 key,value
c
      ffield=' '
      dir_s=' '
      dir_x=' '
      dir_u=' '
      dir_t=' '
      source=' '
c
c     read from rtdb
c     --------------
c
      if(.not.rtdb_cget(irtdb,'prep:ffield',1,ffield)) continue
c
      if(.not.rtdb_cget(irtdb,'prep:dir_s',1,dir_s)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_x',1,dir_x)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_u',1,dir_u)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_t',1,dir_t)) continue
c
      if(.not.rtdb_cget(irtdb,'prep:sysnam',1,sysnam)) then
      if(.not.rtdb_cget(irtdb,'arg:project',1,sysnam)) then
      if(.not.rtdb_cget(irtdb,'file_prefix',1,sysnam))
     + call errquit('rtdb_get failed on sysnam',9999)
      endif
      endif
c
      if(.not.rtdb_cget(irtdb,'prep:source',1,source)) source=' '
c
      if(.not.rtdb_cget(irtdb,'prep:slvnam',1,slvnam)) slvnam='HOH'
      if(.not.rtdb_cget(irtdb,'prep:slvmdl',1,slvmdl)) slvmdl=' '
c
      mdold=0
      if(.not.rtdb_get(irtdb,'md:mdold',mt_int,1,mdold)) mdold=0
c
      if(.not.rtdb_get(irtdb,'prep:newtop',mt_int,1,newtop)) newtop=0
      if(.not.rtdb_get(irtdb,'prep:newseq',mt_int,1,newseq)) newseq=0
      if(.not.rtdb_get(irtdb,'prep:newrst',mt_int,1,newrst)) newrst=0
c
      if(.not.rtdb_get(irtdb,'prep:nfract',mt_int,1,nfract)) nfract=0
      if(nfract.gt.0) then
      if(.not.rtdb_get(irtdb,'prep:ifract',mt_int,nfract,ifract))
     + call errquit('Fraction input problem',0)
      endif
c
      if(.not.rtdb_get(irtdb,'prep:ncount',mt_int,1,ncount)) ncount=0
      if(ncount.gt.0) then
      if(.not.rtdb_get(irtdb,'prep:icount',mt_int,ncount,icount))
     + call errquit('Counter ion input problem',0)
      if(.not.rtdb_cget(irtdb,'prep:scount',ncount,scount))
     + call errquit('Counter ion input problem',0)
      endif
      if(.not.rtdb_get(irtdb,'prep:mgrid',mt_int,1,mgrid)) mgrid=24
      if(.not.rtdb_get(irtdb,'prep:rgrid',mt_dbl,1,gdist)) gdist=0.2d0
c
      if(.not.rtdb_get(irtdb,'prep:mnoe',mt_int,1,mnoe)) mnoe=1
c
      if(.not.rtdb_get(irtdb,'prep:ignore',mt_int,1,ignore)) ignore=0
c
      if(.not.rtdb_get(irtdb,'prep:maxscf',mt_int,1,maxscf)) maxscf=20
      if(.not.rtdb_get(irtdb,'prep:qscale',mt_dbl,1,qscale))
     + qscale=1.2d0
c
      if(.not.rtdb_get(irtdb,'prep:mnoe',mt_int,1,nxlnk)) nxlnk=1
c
      if(.not.rtdb_cget(irtdb,'prep:altloc',1,altloc)) altloc='A'
      if(.not.rtdb_cget(irtdb,'prep:chain',1,chain)) chain='A'
c
      if(.not.rtdb_get(irtdb,'prep:model',mt_int,1,model)) model=0
c
      if(.not.rtdb_get(irtdb,'prep:icyren',mt_int,1,icyren)) icyren=0
      if(.not.rtdb_get(irtdb,'prep:iconst',mt_int,1,iconst)) iconst=0
c
      if(slvnam.eq.'HOH'.and.slvmdl(1:1).eq.' ') slvmdl='spce      '
c
      ndx=index(sysnam,'_')
      len=index(sysnam,' ')
      if(ndx.gt.1.and.ndx.lt.len-1) then
      calc=sysnam(ndx+1:len)
      sysnam=sysnam(1:ndx-1)//' '
      else
      calc='md '
      endif
c
      if(newseq.eq.1) newtop=1
c
c     get defaults for unspecified variables from $HOME/.nwchemrc
c     -----------------------------------------------------------
c
      if(ffield(1:1).eq.' ') then
      key='ffield '
      if(util_nwchemrc_get(key,value)) then
      ffield=value(1:index(value,' ')-1)
      endif
      endif
      if(dir_s(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_s '
      if(util_nwchemrc_get(key,value))
     + dir_s=value(1:index(value,' '))
      endif
      if(dir_x(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_x '
      if(util_nwchemrc_get(key,value))
     + dir_x=value(1:index(value,' '))
      endif
      if(dir_u(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_u '
      if(util_nwchemrc_get(key,value))
     + dir_u=value(1:index(value,' '))
      endif
      if(dir_t(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_t '
      if(util_nwchemrc_get(key,value))
     + dir_t=value(1:index(value,' '))
      endif
      if(dir_t(1:1).eq.' ') dir_t='./ '
c
      pre_rtdbin=.true.
      return
      end
