      subroutine sf_test
c
c $Id: sf_test.F,v 1.10 2003-08-05 21:55:45 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "eaf.fh"
#include "dra.fh"
#include "stdio.fh"
#include "tce.fh"
      integer d_test
      integer size
      integer l_sf,k_sf
      integer i
      logical nodezero
      character*255 filename
      double precision cpu, wall
      integer sizes(4)
      data sizes/100000,1000000,10000000,100000000/
c
      nodezero = (ga_nodeid() .eq. 0)
      if (nodezero.and.util_print('sf',print_default)) then
        if (ioalg.eq.0) then
          write(LuOut,9000)
        else if (ioalg.eq.1) then
          write(LuOut,9000)
        else if (ioalg.eq.2) then
          write(LuOut,9100)
        else if (ioalg.eq.3) then
          write(LuOut,9000)
        else if (ioalg.eq.4) then
          write(LuOut,9200)
        else if (ioalg.eq.5) then
          write(LuOut,9000)
        else if (ioalg.eq.6) then
          write(LuOut,9200)
        endif
      endif
c
c     Coherency
c
      size = ga_nnodes()
      if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1  call errquit('sf_test: MA problem',0)
      call tce_filename('filetest',filename)
      call createfile(filename,d_test,size)
      dbl_mb(k_sf) = dfloat(ga_nodeid()+1)
      call dratoga(d_test)
      call put_block(d_test,dbl_mb(k_sf),1,ga_nodeid())
      call reconcilefile(d_test,size)
      call get_block(d_test,dbl_mb(k_sf),size,0)
      call gatodra(d_test)
      do i = 1,size
        if (dbl_mb(k_sf+i-1) .ne. dfloat(i))
     1    call errquit('sf_test failed',i)
      enddo
      if (.not.ma_pop_stack(l_sf))
     1  call errquit('sf_test: MA problem',1)
      if (nodezero.and.util_print('sf',print_default))
     1  write(LuOut,9010)
      call deletefile(d_test)
c
c     Performance
c
      if (util_print('performance',print_debug)) then
        if (nodezero) write(LuOut,9020)
        do i = 1,4
          size = sizes(i)
          cpu  = -util_cpusec()
          wall = -util_wallsec()
          call tce_filename('filetest',filename)
          call createfile(filename,d_test,size)
          call reconcilefile(d_test,size)
          call deletefile(d_test)
          cpu  = cpu  + util_cpusec()
          wall = wall + util_wallsec()
          if (nodezero) then
            write(LuOut,9030) cpu, wall, size
            call util_flush(LuOut)
          endif
        enddo
      endif
 9000 format(/,1x,
     1 'Global files accessible by all nodes assumed')
 9100 format(/,1x,
     1 'Global array virtual files algorithm will be used')
 9200 format(/,1x,
     1 'Replicated distributed files algorithm will be used')
 9010 format(/,1x,'Parallel file system coherency ......... OK')
 9020 format(/,1x,'Parallel file system performance',
     1       /,1x,'-------------------------------------------')
 9030 format(1x,f8.2,' (',f8.2,') secs / ',i10,' words')
      return
      end
