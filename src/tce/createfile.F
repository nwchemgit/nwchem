      subroutine createfile(filename,d_sf,size)
c
c $Id: createfile.F,v 1.3 2002-10-23 22:44:32 sohirata Exp $
c
c     Creates a shared file and returns the handle
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      character*(*) filename
      integer d_sf
      integer size
c
      if (sf_create(filename,dfloat(size),dfloat(size),
     1  1.0d4,d_sf).ne.0)
     2  call errquit('createfile: sf problem',0)
      call tce_zero(d_sf,size)
      return
      end
