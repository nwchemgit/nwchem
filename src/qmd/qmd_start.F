      subroutine qmd_start()
c
c $Id$
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      call qmd_dimens()
c
      return
      end
      subroutine qmd_dimens()
c
      implicit none
#include "errquit.fh"
c
#include "md_common.fh"
#include "qmd_common.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "nwmd_constants.fh"
c
      real*8 c(3),q,b
      character*16 tag
      integer i
c
      lgeom=0
      jdebug=idebug
      npbtyp=0
      nwm=0
      nwa=0
      msf=1
      nsf=1
      nsm=1
      msm=1
      nwa=0
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('arg_qinput: geom_create failed',0, GEOM_ERR)
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('arg_qinput: geom_rtdb_load failed',0, RTDB_ERR)
      if(.not.geom_ncent(igeom,nsa))
     + call errquit('arg_qinput: geom_ncent failed',0, GEOM_ERR)
      b=zero
      do 1 i=1,nsa
      if(.not.geom_cent_get(igeom,i,tag,c,q))
     + call errquit('arg_rdqdim: geom_cent_get failed',i, GEOM_ERR)
      if(abs(c(1)).gt.b) b=abs(c(1))
      if(abs(c(2)).gt.b) b=abs(c(2))
      if(abs(c(3)).gt.b) b=abs(c(3))
    1 continue
      box(1)=cau2nm*b
      box(2)=cau2nm*b
      box(3)=cau2nm*b
      if(.not.geom_destroy(igeom))
     + call errquit('arg_rdqdim: geom_destroy failed',0, GEOM_ERR)
      mwm=1
      mwa=3
      msa=nsa
c
c      print*,'Dimensions read from rtdb'
      return
      end
      subroutine qmd_setup(numwm,isl,xs,vs,gs,numsa)
c
      implicit none
#include "errquit.fh"
c
#include "qmd_common.fh"
#include "md_common.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "nwmd_constants.fh"
c
      logical task_energy, int_normalize
      external task_energy, int_normalize
c
      integer numwm,numsa
      integer isl(msa,mis2)
      real*8 xs(msa,3),vs(msa,3),gs(msa,3)
c
      real*8 c(3),v(3),q,b,u
      character*16 tag
      integer i,ibasis,junk
c
      numwm=0
      numsa=0
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('arg_qinput: geom_create failed',0, GEOM_ERR)
c
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('arg_qinput: geom_rtdb_load failed',0, RTDB_ERR)
      if(.not.geom_ncent(igeom,nsa))
     + call errquit('arg_qinput: geom_ncent failed',0, GEOM_ERR)
      b=zero
      do 1 i=1,nsa
      if(.not.geom_centv_get(igeom,i,tag,c,v,q))
     + call errquit('arg_rdqdim: geom_cent_get failed',i, GEOM_ERR)
      if(me.eq.0) then
      xs(i,1)=c(1)*cau2nm
      xs(i,2)=c(2)*cau2nm
      xs(i,3)=c(3)*cau2nm
      vs(i,1)=v(1)
      vs(i,2)=v(2)
      vs(i,3)=v(3)
      if(iguide.gt.0) then
      gs(i,1)=0.0d0
      gs(i,2)=0.0d0
      gs(i,3)=0.0d0
      endif
      isl(i,lsgan)=i
      isl(i,lsbox)=0
      isl(i,lsnod)=0
      isl(i,lsatt)=i
      isl(i,lsmol)=1
      isl(i,lssgm)=1
      isl(i,lsgrp)=1
      isl(i,lspgr)=1
      isl(i,lsdyn)=0
      isl(i,lsfrc)=1
      isl(i,lssss)=0
      isl(i,lsct1)=1
      isl(i,lsct2)=1
      isl(i,lsct3)=1
      endif
    1 continue
      if(me.eq.0) numsa=nsa
c
      if(.not.bas_create(ibasis,'ao basis'))
     + call errquit('arg_qenrgy: bas_create failed',0, BASIS_ERR)
c
      if(.not.bas_rtdb_load(irtdb,igeom,ibasis,'ao basis'))
     + call errquit('arg_qenrgy: bas_rtdb_load failed',0, RTDB_ERR)
c
      if (.not.int_normalize(irtdb,ibasis))
     + call errquit('arg_qenrgy: int_normalize',0, INT_ERR)
c
      call ga_sync()
      if(jdebug.gt.0) then
      call ma_summarize_allocated_blocks()
      endif
      call ga_sync()
c
      call int_init(irtdb,1,ibasis)
      call guess_dens_gen(irtdb,igeom,ibasis,.false.,junk,uqmatm)
      call int_terminate()
c
      uqmatm=uqmatm*cau2kj
c
      if(.not.geom_destroy(igeom))
     + call errquit('arg_qenrgy: geom_destroy failed',0, GEOM_ERR)
c
      if(.not.bas_destroy(ibasis))
     + call errquit('arg_qenrgy: bas_destroy failed',0, BASIS_ERR)
c
      call ga_sync()
      if(.not.task_energy(irtdb))
     + call errquit('arg_qenrgy: task_energy failed',0, CALC_ERR)
      call ga_sync()
c
      if(.not.rtdb_get(irtdb,'task:energy',mt_dbl,1,u))
     + call errquit('arg_qenrgy: rtdb_get energy failed',0, RTDB_ERR)
c
c      print*,'Coordinates read from rtdb'
      return
      end
