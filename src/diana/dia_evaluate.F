      subroutine dia_evaluate()
c
c $Id: dia_evaluate.F,v 1.1 2007-02-05 21:37:42 d3j191 Exp $
c
c     evaluate the defined properties
c     -------------------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer itval
      integer ilt,iht,jlt,jht,ndxt,ldt
c
c     allocate global array for the properties
c     ----------------------------------------
c
      if(ntval.gt.0) then
      if(.not.ga_destroy(ga_tval))
     + call md_abort('Failed to destroy ga_tval',0)
      endif
      ntval=ntors
c
      if(.not.ga_create(mt_dbl,ntval,nfrtot,'trj',ntval,nfrdim,
     + ga_tval)) call md_abort('Failed to create ga_tval',0)
c
      itval=1
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
      call ga_distribution(ga_tval,me,ilt,iht,jlt,jht)
      write(*,'(5i5)') me,ilt,iht,jlt,jht
      call ga_access(ga_tval,ilt,iht,jlt,jht,ndxt,ldt)
c
c     evaluate distances
c     evaluate angles
c
c     evaluate torsions
c     -----------------
c
      call dia_evaltors(dbl_mb(ndx),nsa,dbl_mb(ndxt),nfrme,itval)
c
      return
      end
      subroutine dia_evaltors(xtrj,ntrj,tval,nfr,itval)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
c
      real*8 dia_torsion
      external dia_torsion
c
      integer ntrj,nfr,itval
      real*8 xtrj(ntrj,3,nfr),tval(ntval,nfrdim)
c
      real*8 wtrj(3)
c
      integer i,k
c
      do 1 i=1,ntors
      do 2 k=1,nfr
      tval(itval+i-1,k)=dia_torsion(xtrj(1,1,k),ntrj,wtrj,1,3,
     + itors(i,1,1),itors(i,1,2),itors(i,2,1),itors(i,2,2),
     + itors(i,3,1),itors(i,3,2),itors(i,4,1),itors(i,4,2))
    2 continue
    1 continue
c
      return
      end
      subroutine dia_export()
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
#include "mafdecls.fh"
c
      integer i,j,ilt,iht,jlt,jht,nfr
      integer id1,id2,idexp,jtors
c
      if(me.eq.0) then
c
      idexp=0
c
      if(card(8:13).eq.'torsio') then
      read(card(15:19),'(i5)') jtors
      filexp=card(20:80)
      idexp=3
      if(jtors.eq.0) then
      id1=1
      id2=ntors
      else
      id1=0
      do 2 i=1,ntors
      if(itors(i,1,3).eq.jtors) then
      if(id1.eq.0) id1=i
      id2=i
      endif
    2 continue
      endif
      endif
c
      open(unit=lfnexp,file=filexp(1:index(filexp,' ')-1),
     + form='formatted',status='unknown')
c
      if(.not.ma_push_get(mt_dbl,nfrdim*ntval,'tv',l_tval,i_tval))
     + call md_abort('Failed to allocate tv',0)
c
      do 1 i=0,np-1
      call ga_distribution(ga_tval,i,ilt,iht,jlt,jht)
      call ga_get(ga_tval,ilt,iht,jlt,jht,dbl_mb(i_tval),ntval)
      nfr=min(nfrdim,nfrtot-i*nfrdim)
      if(idexp.eq.3) call dia_exptors(dbl_mb(i_tval),nfr,id1,id2)
    1 continue
c
      if(.not.ma_pop_stack(l_tval))
     + call md_abort('Failed to deallocate tv',0)
c
      close(unit=lfnexp)
c
      endif
c
      return
      end
      subroutine dia_exptors(tval,nfr,id1,id2)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
#include "mafdecls.fh"
c
      integer nfr,id1,id2
      real*8 tval(ntval,nfrdim)
c
      integer i,j
c
      do 1 j=1,nfr
      write(lfnexp,1000) (tval(i,j),i=id1,id2)
 1000 format(10f12.6)
    1 continue
c
      return
      end
