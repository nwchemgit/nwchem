      logical function pre_rtdbin(irtdb,ffield,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + par_s,par_x,par_q,par_u,par_t,par_c,
     + source,sysnam,calc,slvnam,slvmdl,
     + newtop,newseq,newrst,mcount,ncount,icount,mgrid,gdist,mnoe,
     + maxscf,qscale,
     + altloc,chain,icyren,iconst,model,nxlnk,mdold,ignore,scount,
     + mfract,nfract,ifract,scale,cpk,his)
c
c $Id: pre_rtdbin.F,v 1.30 2002-02-13 00:53:36 d3j191 Exp $
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
      character*255 dir_s,dir_x,dir_q,dir_u,dir_t,dir_c
      character*255 par_s,par_x,par_q,par_u,par_t,par_c
      integer newtop,newseq,newrst,icount(mcount),icyren,ignore
      integer len,ndx,mgrid,mnoe,maxscf,model,nxlnk,mdold
      integer mfract,nfract,iconst
      integer ifract(mfract)
      real*8 gdist,qscale
      character*1 altloc,chain,his
      character*4 scount(mcount)
      real*8 scale,cpk
      integer i
c
      character*255 key,value
c
      ffield=' '
      dir_s=' '
      dir_x=' '
      dir_q=' '
      dir_u=' '
      dir_t=' '
      dir_c=' '
      par_s=' '
      par_x=' '
      par_q=' '
      par_u=' '
      par_t=' '
      par_c=' '
      source=' '
c
c     read from rtdb
c     --------------
c
      if(.not.rtdb_cget(irtdb,'prep:ffield',1,ffield)) continue
c
      if(.not.rtdb_cget(irtdb,'prep:dir_s',1,dir_s)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_x',1,dir_x)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_q',1,dir_q)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_u',1,dir_u)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_t',1,dir_t)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_c',1,dir_c)) continue
c
      if(.not.rtdb_cget(irtdb,'prep:par_s',1,par_s)) continue
      if(.not.rtdb_cget(irtdb,'prep:par_x',1,par_x)) continue
      if(.not.rtdb_cget(irtdb,'prep:par_x',1,par_q)) continue
      if(.not.rtdb_cget(irtdb,'prep:par_u',1,par_u)) continue
      if(.not.rtdb_cget(irtdb,'prep:par_t',1,par_t)) continue
      if(.not.rtdb_cget(irtdb,'prep:par_c',1,par_c)) continue
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
     + qscale=1.0d0
c
      if(.not.rtdb_get(irtdb,'prep:nxlnk',mt_int,1,nxlnk)) nxlnk=1
c
      if(.not.rtdb_cget(irtdb,'prep:altloc',1,altloc)) altloc='A'
      if(.not.rtdb_cget(irtdb,'prep:chain',1,chain)) chain='A'
      if(.not.rtdb_cget(irtdb,'prep:his',1,his)) his=' '
c
      if(.not.rtdb_get(irtdb,'prep:model',mt_int,1,model)) model=0
c
      if(.not.rtdb_get(irtdb,'prep:icyren',mt_int,1,icyren)) icyren=0
      if(.not.rtdb_get(irtdb,'prep:iconst',mt_int,1,iconst)) iconst=0
c
      if(.not.rtdb_get(irtdb,'prep:scale',mt_dbl,1,scale)) scale=1.0d0
      if(.not.rtdb_get(irtdb,'prep:cpk',mt_dbl,1,cpk)) cpk=0.0d0
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
      if(dir_q(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_q '
      if(util_nwchemrc_get(key,value))
     + dir_q=value(1:index(value,' '))
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
      if(dir_c(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_c '
      if(util_nwchemrc_get(key,value))
     + dir_c=value(1:index(value,' '))
      endif
      if(dir_c(1:1).eq.' ') dir_c='./ '
c
      len=index(dir_s,' ')-1
      if(len.gt.0) then
      if(dir_s(len-3:len).eq.'.par') then
      ndx=0
      do 1 i=1,len
      if(dir_s(i:i).eq.'/') ndx=i
    1 continue
      if(ndx.eq.0) then
      par_s=dir_s
      dir_s='./'
      elseif(ndx.eq.len) then
      par_s=' '
      else
      par_s=dir_s(ndx+1:len)
      dir_s=dir_s(1:ndx)
      endif
      endif
      endif
c
      len=index(dir_x,' ')-1
      if(len.gt.0) then
      if(dir_x(len-3:len).eq.'.par') then
      ndx=0
      do 2 i=1,len
      if(dir_x(i:i).eq.'/') ndx=i
    2 continue
      if(ndx.eq.0) then
      par_x=dir_x
      dir_x='./'
      elseif(ndx.eq.len) then
      par_x=' '
      else
      par_x=dir_x(ndx+1:len)
      dir_x=dir_x(1:ndx)
      endif
      endif
      endif
c
      len=index(dir_q,' ')-1
      if(len.gt.0) then
      if(dir_q(len-3:len).eq.'.par') then
      ndx=0
      do 3 i=1,len
      if(dir_q(i:i).eq.'/') ndx=i
    3 continue
      if(ndx.eq.0) then
      par_q=dir_q
      dir_q='./'
      elseif(ndx.eq.len) then
      par_q=' '
      else
      par_q=dir_q(ndx+1:len)
      dir_q=dir_q(1:ndx)
      endif
      endif
      endif
c
      len=index(dir_u,' ')-1
      if(len.gt.0) then
      if(dir_u(len-3:len).eq.'.par') then
      ndx=0
      do 4 i=1,len
      if(dir_u(i:i).eq.'/') ndx=i
    4 continue
      if(ndx.eq.0) then
      par_u=dir_u
      dir_u='./'
      elseif(ndx.eq.len) then
      par_u=' '
      else
      par_u=dir_u(ndx+1:len)
      dir_u=dir_u(1:ndx)
      endif
      endif
      endif
c
      len=index(dir_t,' ')-1
      if(len.gt.0) then
      if(dir_t(len-3:len).eq.'.par') then
      ndx=0
      do 5 i=1,len
      if(dir_t(i:i).eq.'/') ndx=i
    5 continue
      if(ndx.eq.0) then
      par_t=dir_t
      dir_t='./'
      elseif(ndx.eq.len) then
      par_t=' '
      else
      par_t=dir_t(ndx+1:len)
      dir_t=dir_t(1:ndx)
      endif
      endif
      endif
c
      len=index(dir_c,' ')-1
      if(len.gt.0) then
      if(dir_c(len-3:len).eq.'.par') then
      ndx=0
      do 6 i=1,len
      if(dir_c(i:i).eq.'/') ndx=i
    6 continue
      if(ndx.eq.0) then
      par_c=dir_c
      dir_c='./'
      elseif(ndx.eq.len) then
      par_c=' '
      else
      par_c=dir_c(ndx+1:len)
      dir_c=dir_c(1:ndx)
      endif
      endif
      endif
c
      key='histidine'
      if(his.eq.' ') then
      if(util_nwchemrc_get(key,value)) then
      if(value(1:3).eq.'HID') his='D'
      if(value(1:3).eq.'HIE') his='E'
      if(value(1:3).eq.'HIP') his='P'
      else
      his='E'
      endif
      endif
c
      pre_rtdbin=.true.
      return
      end
