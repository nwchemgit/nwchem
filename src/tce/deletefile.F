      subroutine deletefile(d_sf)
c
c $Id: deletefile.F,v 1.2 2002-10-23 01:38:51 sohirata Exp $
c
c     Creates a shared file and returns the handle
c
      implicit none
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      integer d_sf
c
      if (sf_destroy(d_sf).ne.0)
     1  call errquit('deletefile: sf problem',0)
      return
      end
