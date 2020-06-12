c
c     Coherency
c
      program ga_test
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "tce_main.fh"
      integer d_test
      integer size
      integer l_sf,k_sf
      integer i
      logical nodezero
      character*255 filename
      call pbeginf()
      call ga_initialize()
      if (.not.ma_init(mt_dbl, 1000000, 1000000))
     1  call errquit('ma_init failed',0,0)
      ioalg=2
      nodezero = (ga_nodeid() .eq. 0)
      size = ga_nnodes()
      if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1  call errquit('sf_test: MA problem',0,0)
      call tce_filename('filetest',filename)
      call createfile(filename,d_test,size)
      dbl_mb(k_sf) = dble(ga_nodeid()+1)
      call put_block(d_test,dbl_mb(k_sf),1,ga_nodeid())
      call reconcilefile(d_test,size)
      call get_block(d_test,dbl_mb(k_sf),size,0)
      do i = 1,size
        if (dbl_mb(k_sf+i-1) .ne. dble(i))
     1    call errquit('sf_test failed',i,0)
      enddo
      if (.not.ma_pop_stack(l_sf))
     1  call errquit('sf_test: MA problem',1,0)
      if (nodezero)
     1  write(6,9010)
      call deletefile(d_test)
      call ga_terminate()
      call pend()
 9010 format(/,1x,'Parallel file system coherency ......... OK')
      end
