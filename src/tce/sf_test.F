      subroutine sf_test
c
c $Id: sf_test.F,v 1.3 2002-10-30 22:05:30 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "stdio.fh"
#include "tce.fh"
      integer d_test
      integer size
      integer l_sf,k_sf
      integer i
      logical nodezero
      character*255 filename
c
      nodezero = (ga_nodeid() .eq. 0)
      if (nodezero.and.util_print('sf',print_default))
     1  write(LuOut,9000)
      size = ga_nnodes()
      if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1  call errquit('sf_print: MA problem',0)
      call util_file_name('sftest',.true.,.false.,filename)
      call createfile(filename,d_test,size)
      dbl_mb(k_sf) = dfloat(ga_nodeid()+1)
      call put_block(d_test,dbl_mb(k_sf),1,ga_nodeid())
      call ga_sync()
      call get_block(d_test,dbl_mb(k_sf),size,0)
      do i = 1,size
        if (dbl_mb(k_sf+i-1) .ne. dfloat(i))
     1    call errquit('sf_test failed',i)
      enddo
      if (.not.ma_pop_stack(l_sf))
     1  call errquit('sf_print: MA problem',1)
      if (nodezero.and.util_print('sf',print_default))
     1  write(LuOut,9010)
      call deletefile(d_test)
 9000 format(/,1x,
     1 'This module uses a global file system accessible by all nodes')
 9010 format(/,1x,'I/O self test is being performed ... passed')
      return
      end
