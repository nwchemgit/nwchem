      subroutine argos_diana_rdfwrt(isel,irdf)
c
      implicit none
c
#include "argos_diana_common.fh"
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
      write(*,'(a)') 'ARGOS_DIANA_RDFWRT'
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
c $Id$
