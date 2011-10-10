      subroutine argos_diana_hbond_header(snam,isel,ihbnd)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      character*16 snam(msa)
      integer isel(msa),ihbnd(msa)
c
      integer i
c
      write(lfnhba,1000) nsa
 1000 format('hydrogen bond analysis',/,i7)
      do 1 i=1,nsa
      write(lfnhba,1001) i,snam(i),isel(i),ihbnd(i)
 1001 format(i7,3x,a,2i7)
    1 continue
c
      return
      end
c $Id$
