      program nwchem
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "rtdb.fh"
      integer rtdb
      integer stack
      integer heap
c
c     Create parallel processes and initialize IPC layer
c
      call pbeginf()
c
c     Look for memory directive in the input ... must eventually
c     open the rtdb first so that can get memory directive out of that
c     if it is not in the input
c
      call input_mem_size(stack, heap)
c
c     Initialize local memory allocator
c
      if (.not. ma_init(mt_dbl, stack, heap))
     $     call errquit('nwchem: ma_init failed', -1)
c
c     Intialize global memory allocator
c
      call ga_initialize()
c
c     Read the input and open the RTDB
c
      call start(rtdb)
c
c     Call driver routine to do something useful
c
      call control(rtdb)
c
c     Close the RTDB
c
      if (.not. rtdb_par_close(rtdb, 'keep'))
     $     call errquit('nwchem: rtdb_close failed', rtdb)
c
c     Tidy up the global memory
c
      call ga_terminate()
c
c     Tidy up the parallel environment
c
      call pend()
c
c     Done
c
      end
      subroutine start(rtdb)
      implicit none
#include "mafdecls.h"
#include "rtdb.fh"
      integer rtdb              ! [output]
      character*255 rtdb_name
      character*255 file_prefix
      logical startup
c
c     Figure out if this is a start/restart and what the database
c     and file prefix are
c
      call input_start_opt(rtdb_name, file_prefix, startup)
c
c     Actually open the database and insert the file_prefix
c
      if (startup) then
         if (.not. rtdb_par_open(rtdb_name, 'empty', rtdb))
     $        call errquit('start: rtdb_open empty failed ', 0)
      else
         if (.not. rtdb_par_open(rtdb_name, 'old', rtdb))
     $        call errquit('start: rtdb_open old failed ', 0)
      endif
      if (.not. rtdb_cput(rtdb, 'file_prefix', 1, file_prefix))
     $     call errquit('start: rtdb_put of file_prefix failed', 0)
c      
c     Make sure input is all ready to go
c
      rewind 5
      call inp_init(5,6)
c
c     Parse input data and shove into the database
c
      call input_parse(rtdb)
c
c     Done
c
      end
      subroutine control(rtdb)
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb              ! [input]
      integer geom
      integer basis
c
      write(6,*)
      write(6,*)
      if (.not. geom_create(geom, 'tmp'))
     $     call errquit('control: geom_create failed', 0)
      if  (.not. geom_rtdb_load(rtdb, geom, 'test geom'))
     $     call errquit('control: geom_load failed', 0)
      if (.not. geom_print(geom))
     $     call errquit('control: geom_print failed', 0)
      if (.not. bas_create(basis,' '))
     $     call errquit('control: basis create failed', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $      call errquit('control: basis load failed', 0)
      if (.not. bas_print(basis))
     $      call errquit('control: basis print failed', 0)
      if (.not. geom_destroy(geom))
     $     call errquit('control: geom_destroy failed', 0)
      if (.not.bas_destroy(basis))
     $     call errquit('control: basis destroy failed', 0)
      if (.not. rtdb_print(rtdb, .true.))
     $     call errquit('control: rtdb_print failed', 0)
c
      end
