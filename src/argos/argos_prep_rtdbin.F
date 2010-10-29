      logical function argos_prep_rtdbin(irtdb,title,
     + source,sysnam,calc,slvnam,slvmdl,
     + newtop,newseq,newrst,mcount,ncount,icount,
     + mgrid,nrgrid,iogrid,rogrid,gdist,mnoe,mpmf,
     + maxscf,qscale,
     + altloc,chain,icyren,iconst,model,nxlnk,mdold,ignore,scount,
     + fcount,mfract,nfract,ifract,scale,cpk,his,itopol,icomb)
c
c $Id$
c
c     function to read input for prepare module from rtdb
c
      implicit none
c
      logical util_nwchemrc_get,argos_prep_paths,argos_prep_paths2
      external util_nwchemrc_get,argos_prep_paths,argos_prep_paths2
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "argos_prep_common.fh"
c
      integer irtdb,mcount,ncount
      character*3 slvnam
      character*80 source
      character*80 sysnam,calc,title(2,3)
      character*10 slvmdl
      integer newtop,newseq,newrst,icount(mcount),icyren,ignore
      integer len,ndx,mgrid,mnoe,mpmf,maxscf,model,nxlnk,mdold
      integer mfract,nfract,iconst,itopol,icomb
      integer ifract(mfract)
      real*8 gdist,qscale
      character*1 altloc,chain,his
      character*4 scount(mcount)
      real*8 fcount
      real*8 scale,cpk
      integer i,j,nrgrid,iogrid(5),rogrid(2,5)
c
      character*255 key,value
c
      ffield=' '
c
      key='ffield '
      if(util_nwchemrc_get(key,value)) then
      ffield=value(1:index(value,' ')-1)
      if(.not.argos_prep_paths())
     + call md_abort('Error reading path info from .nwchemrc',0)
      endif
c
      source=' '
c
      do 1 i=1,2
      do 2 j=1,3
      title(i,j)=' '
    2 continue
    1 continue
c
c     read from rtdb
c     --------------
c
      if(.not.rtdb_cget(irtdb,'prep:title',6,title)) continue
 
      if(rtdb_cget(irtdb,'prep:ffield',1,ffield)) then
      if(.not.argos_prep_paths())
     + call md_abort('Error reading path info from .nwchemrc',0)
      endif
c
      if(.not.argos_prep_paths2(irtdb)) call md_abort('prepaths2 ',0)
c
      if(.not.rtdb_cget(irtdb,'prep:sysnam',1,sysnam)) then
      if(.not.rtdb_cget(irtdb,'arg:project',1,sysnam)) then
      if(.not.rtdb_cget(irtdb,'file_prefix',1,sysnam))
     + call md_abort('rtdb_get failed on sysnam',9999)
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
     + call md_abort('Fraction input problem',0)
      endif
c
      if(.not.rtdb_get(irtdb,'prep:ncount',mt_int,1,ncount)) ncount=0
      if(ncount.gt.0) then
      if(.not.rtdb_get(irtdb,'prep:icount',mt_int,ncount,icount))
     + call md_abort('Counter ion input problem',0)
      if(.not.rtdb_cget(irtdb,'prep:scount',ncount,scount))
     + call md_abort('Counter ion input problem',0)
      endif
      if(.not.rtdb_get(irtdb,'prep:fcount',mt_dbl,1,fcount))
     + call md_abort('Counter ion input problem',0)
      if(.not.rtdb_get(irtdb,'prep:mgrid',mt_int,1,mgrid)) mgrid=24
      if(.not.rtdb_get(irtdb,'prep:nrgrid',mt_int,1,nrgrid)) nrgrid=0
      if(nrgrid.gt.0) then
      if(.not.rtdb_get(irtdb,'prep:iogrid',mt_int,nrgrid,iogrid))
     + call md_abort('argos_prep_rtdb: rtdb_get failed',0)
      if(.not.rtdb_get(irtdb,'prep:rogrid',mt_int,2*nrgrid,rogrid))
     + call md_abort('argos_prep_rtdbin: rtdb_get failed',0)
      endif
      if(.not.rtdb_get(irtdb,'prep:rgrid',mt_dbl,1,gdist)) gdist=0.2d0
c
      if(.not.rtdb_get(irtdb,'prep:mnoe',mt_int,1,mnoe)) mnoe=0
c
      if(.not.rtdb_get(irtdb,'prep:mpmf',mt_int,1,mpmf)) mpmf=0
c
      if(.not.rtdb_get(irtdb,'prep:itopol',mt_int,1,itopol)) itopol=0
c
      if(.not.rtdb_get(irtdb,'prep:icomb',mt_int,1,icomb)) icomb=0
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
      lqhopmem=.false.
c
      call qhop_alloc()
      argos_prep_rtdbin=.true.
      return
      end
