      subroutine sf_test
c
c $Id: sf_test.F,v 1.8 2003-08-01 22:59:51 sohirata Exp $
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
      call put_block(d_test,dbl_mb(k_sf),1,ga_nodeid())
      call reconcilefile(d_test,size)
      call dratoga(d_test,size)
      call get_block(d_test,dbl_mb(k_sf),size,0)
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
c     if (nodezero) write(LuOut,9020)
c     do i = 1,4
c       size = sizes(i)
c       cpu  = -util_cpusec()
c       wall = -util_wallsec()
c       call tce_filename('filetest',filename)
c       call createfile(filename,d_test,size)
c       call reconcilefile(d_test,size)
c       call deletefile(d_test)
c       cpu  = cpu  + util_cpusec()
c       wall = wall + util_wallsec()
c       if (nodezero) then
c         write(LuOut,9030) cpu, wall, size
c         call util_flush(LuOut)
c       endif
c     enddo
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
