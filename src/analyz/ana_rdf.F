      subroutine ana_rdf_init()
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      if(.not.ma_push_get(mt_int,nrdf*nsel*mwa,'irdf',l_rdf,i_rdf))
     + call md_abort('Could not allocate irdf',0)
      print*,'rdf allocated in rdfhdr ',nrdf*nsel*mwa
c
      return
      end
      subroutine ana_rdfhdr(irdf)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer irdf(nsel,mwa,nrdf)
c
      character*255 fname
      integer i,j,k,lq,l,m
c
      fname=filtrj
      lq=index(filtrj,'?')
      if(lq.gt.0) then
      fname=filtrj(1:lq-1)//cnum//filtrj(lq+1:index(filtrj,' ')-1)
      endif
      lq=index(fname,'.trj')
      fname(lq:lq+3)='.rdf'
c
      open(unit=lfnrdf,file=fname(1:index(fname,' ')-1),
     + status='unknown')
c
      write(*,3333) fname(1:index(fname,' ')-1)
 3333 format(/,' Opening rdf file ',a)
      rewind(lfnrdf)
c
      numrdf=0
c
      do 1 i=1,nrdf
      do 2 j=1,mwa
      do 3 k=1,nsel
      irdf(k,j,i)=0
    3 continue
    2 continue
    1 continue
c
      print*,'rdf set to zero in rdfhdr'
c
      return
      end
      subroutine ana_rdf(isel,xs,xw,irdf)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(msa),irdf(nsel,mwa,nrdf)
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer i,j,k,l,m
      real*8 d
c
      i=0
      do 1 j=1,nsa
      if(isel(j).gt.0) then
      i=i+1
      do 2 k=1,nwm
      do 3 l=1,nwa
      d=sqrt((xs(j,1)-xw(k,l,1))**2+(xs(j,2)-xw(k,l,2))**2+
     + (xs(j,3)-xw(k,l,3))**2)
      m=int(dble(nrdf*rrdf)/d)
      if(m.le.nrdf) irdf(i,l,m)=irdf(i,l,m)+1
    3 continue
    2 continue
      endif
    1 continue
c
      numrdf=numrdf+1
c
      return
      end
      subroutine ana_rdfwrt(isel,irdf)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(msa)
      integer irdf(nsel,mwa,nrdf)
c
      integer i,j,k,l,m
      real*8 d
c
      write(*,'(a)') 'ANA_RDFWRT'
c
      do 4 i=1,nsel
      do 5 m=1,nrdf
      write(*,'(2i5,10i10)') i,m,(irdf(i,l,m),l=1,nwa)
    5 continue
    4 continue
c
c      i=0
c      do 1 j=1,nsa
c      if(isel(j).gt.0) then
c      i=i+1
c      do 2 k=1,nrdf
c      d=dble(k)*rrdf/dble(nrdf)
c      write(*,1000) d,(irdf(i,l,k),l=1,nwa)
c 1000 format(f12.6,10i5)
c    2 continue
c      endif
c    1 continue
c
      close(unit=lfnrdf,status='keep')
      write(*,'(a)') ' Closing rdf file '
c
      return
      end
