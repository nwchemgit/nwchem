      subroutine int_mem(max1e, maxg, mscratch_1e, mscratch_2e)
c $Id: int_mem.F,v 1.5 1995-05-09 18:00:59 d3e129 Exp $
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer max1e, maxg, mscratch_1e, mscratch_2e
c
      if (int_chk_init('int_mem')) then
        max1e       = isz_1e
        maxg        = isz_2e4c
        mscratch_1e = mem_1e
        mscratch_2e = mem_2e4c
      else
        call errquit('int_mem: int_init was not called' ,0)
      endif
      end
      subroutine int_mem_1e(max1e, mscratch_1e)
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer max1e, mscratch_1e
c
      if (int_chk_init('int_mem_1e')) then
        max1e       = isz_1e
        mscratch_1e = mem_1e
      else
        call errquit('int_mem_1e: int_init was not called' ,0)
      endif
      end
      subroutine int_mem_2e4c(maxg, mscratch_2e)
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxg, mscratch_2e
c
      if (int_chk_init('int_mem_2e4c')) then
        maxg        = isz_2e4c
        mscratch_2e = mem_2e4c
      else
        call errquit('int_mem_2e4c: int_init was not called' ,0)
      endif
      end
      subroutine int_mem_h1(maxh1, mscratch_h1)
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxh1, mscratch_h1
c
      if (int_chk_init('int_mem_h1')) then
        maxh1       = isz_1e
        mscratch_h1 = mem_1e_h1
      else
        call errquit('int_mem_h1: int_init was not called' ,0)
      endif
      end
      subroutine int_mem_2e3c(maxg, mscratch_2e3c)
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxg, mscratch_2e3c
c
      if (int_chk_init('int_mem_2e3c')) then
        maxg          = isz_2e3c
        mscratch_2e3c = mem_2e3c
      else
        call errquit('int_mem_2e3c: int_init was not called' ,0)
      endif
      end
      subroutine int_mem_2e2c(maxg, mscratch_2e2c)
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxg, mscratch_2e2c
c
      if (int_chk_init('int_mem_2e2c')) then
        maxg          = isz_2e2c
        mscratch_2e2c = mem_2e2c
      else
        call errquit('int_mem_2e2c: int_init was not called' ,0)
      endif
      end
      subroutine int_mem_3ov(maxbuf, mscratch_3ov)
      implicit none
#include "apiP.fh"
c::functions
      logical int_chk_init
      external int_chk_init
c::passed
      integer maxbuf, mscratch_3ov
c
      if (int_chk_init('int_mem_3ov')) then
        maxbuf       = isz_3ov
        mscratch_3ov = mem_3ov
      else
        call errquit('int_mem_3ov: int_init was not called' ,0)
      endif
      end
      subroutine int_mem_print()
      implicit none
#include "stdio.fh"
#include "apiP.fh"
c
      write(luout,*)' memory parameters for integral api:'
      write(luout,10000)
c
      if (init_int.eq.1)
     &    write(luout,*)' int_init was called '
      if (init_intd.eq.1)
     &    write(luout,*)' intd_init was called '

      write(luout,*)' mem_2e4c : ',mem_2e4c
      write(luout,*)' mem_1e   : ',mem_1e
      write(luout,*)' isz_2e4c : ',isz_2e4c
      write(luout,*)' isz_1e   : ',isz_1e
      write(luout,*)' mem_1e_h1: ',mem_1e_h1
      write(luout,*)' isz_2e3c : ',isz_2e3c
      write(luout,*)' isz_2e2c : ',isz_2e2c
      write(luout,*)' mem_2e3c : ',mem_2e3c
      write(luout,*)' mem_2e2c : ',mem_2e2c
      write(luout,*)' isz_3ov  : ',isz_3ov
      write(luout,*)' mem_3ov  : ',mem_3ov
10000 format(//)
      end
