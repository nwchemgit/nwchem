      subroutine deletefile(d_sf)
c
c $Id: deletefile.F,v 1.1 2002-10-18 23:42:06 sohirata Exp $
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
      integer dummy
c
      dummy = sf_destroy(d_sf)
      return
      end
