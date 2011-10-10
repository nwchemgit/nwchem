      subroutine argos_diana_rdfhdr(irdf)
c
      implicit none
c
#include "argos_diana_common.fh"
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
c $Id$
