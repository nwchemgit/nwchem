      subroutine argos_diana_espplt()
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "util.fh"
c
      character*255 string
      integer i
c
      string=card(8:80)
      read(string(1:4),'(i4)') ngrid
      read(string(5:8),'(i4)') iesppb
      read(string(9:20),'(f12.6)') rcut
      string=card(28:80)
      i=index(string,' ')
      filchg=string(1:i)
      filplt=string(i+1:80)
      filepz=filplt(1:index(filplt,'.'))//'epz '
c
      if(filchg(1:1).ne.' ') call argos_diana_pltgrd(.true.,
     + int_mb(i_isel))
c
      return
      end
