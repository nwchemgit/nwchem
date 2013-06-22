      subroutine pgcreatefile(filename,d_sf,size,pg)
c
c     Creates a shared file and returns the handle
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "eaf.fh"
#include "dra.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_mrcc.fh"
      character*(*) filename
      integer d_sf
      integer d_temp
      integer size
      integer i,pg
      logical nodezero
      logical parallel
      logical used

      nodezero = (ga_nodeid().eq.0)
      if (size.eq.0) then
        if (nodezero.and.util_print('files',print_high)) then
        if (nodezero) then
          write(LuOut,*) '*************** Warning ***************'
          write(LuOut,*) ' Create file size is zero. Calculation'
          write(LuOut,*) ' will continue by increasing the size.'
          write(LuOut,*) ' Use of a larger basis set is advised.'
        endif
        endif
        size = 1
      endif

      if (ioalg.eq.2) then
        if (.not.nga_create_config(MT_DBL,1,size,filename,
     1    -1,pg,d_sf)) then
          write(LuOut,*) ' available GA memory ',
     1    ga_memory_avail(),' bytes. Need: ',size
          call errquit ('createfile: failed ga_create size=',size,
     1      GA_ERR)
        endif
        call ga_zero(d_sf)
      else
          call errquit ('createfile: non-GA alg unsupported',0,GA_ERR)
      endif

      return
      end













c $Id$
