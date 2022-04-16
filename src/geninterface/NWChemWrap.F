C $Id$
      subroutine NWChem_nwchemStart(directory)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "util.fh"
#include "inp.fh"


      character*(nw_max_path_len) input_filename, rtdb_name
      character*(*) directory
      integer stack, heap, global
      logical overify, ohardfail, ostartup, ocontinue
      integer rtdb
      common /sidlinfo/ rtdb
      character*32 temp_directory1, temp_directory2

C  Create parallel processes
      call pbeginf()
    
C  Worry about input file.  For now use a dummy one
      write(input_filename, '('''',a,''/junk.nw'')')
     $      directory(1:inp_strlen(directory))
      if (nodeid().eq.0) then
         open(LuIn, file=input_filename, status='unknown', 
     $        form='formatted',
     $        err=10000)
         goto 10001
10000    call errquit('nwchem: failed to open the input file', 0)
10001    continue
c
      else
#if !(defined(KSR) || defined(IBM) || defined(FUJITSU) )
         close(LuIn)
#endif
      endif

C  Write scratch directory path into dummy input file

C      write(temp_directory1, '(''scratch_dir '',a,'' ; '')')
C     $      directory(1:inp_strlen(directory))
C
C      if (nodeid().eq.0) then
C         write (LuIn,*) temp_directory1
C      end if
C
C      write(temp_directory2, '(''permanent_dir '',a,'''')')
C     $      directory(1:inp_strlen(directory))
C
C      if (nodeid().eq.0) then
C         write (LuIn,*) temp_directory2
C         write (LuIn,*)
C      end if

C  Set up memory environment
      call input_mem_size(stack, heap, global, overify, ohardfail)
      call ga_initialize()
      if ( ga_uses_ma() ) then
        if (.not.ma_init(mt_dbl, stack, heap+global))
     $      call errquit
     $      ('initialize_qm: ma_init failed (ga_uses_ma=T)',555)
      else
        if (.not.ma_init(mt_dbl,stack,heap))
     $      call errquit
     $      ('initialize_qm: ma_init failed (ga_uses_ma=F)',555)
        call ga_set_memory_limit(ma_sizeof(mt_dbl,global,mt_byte))
      endif

C  Open rtdb
C  Need to get the input file name and the rtdb_name taken care of
      call input_file_info(input_filename, rtdb_name,
     $                     ostartup, ocontinue)
      if (.not. rtdb_open(rtdb_name, 'empty', rtdb))
     $     call errquit('initialize_qm: rtdb_open empty failed ', 555,0)
      call util_set_rtdb_state(rtdb,.true.,.false.,.false.)
      call util_file_info_rtdb(rtdb)
      call geom_hnd_parallel(.true.)

      return
      end

      subroutine NWChem_nwchemEnd()
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "util.fh"
      integer rtdb
      common /sidlinfo/ rtdb

C  Close the rtdb
      if (.not. rtdb_close(rtdb, 'delete'))
     $     call errquit('finalize_qm: rtdb_close failed', rtdb,0)

C  Close the global memory
      call ga_sync()
      call ga_terminate()

C  Close the parallel environment
      call pend()

      return
      end

      subroutine NWChem_setCoordinates(coords)
      implicit none
#include "geom.fh"

      integer geom
      double precision coords(*)
      integer rtdb
      common /sidlinfo/ rtdb

      if (.not.geom_create(geom,'geometry'))
     $      call errquit('setCoordinates: geom_create failed',555,0)

      if (.not.geom_rtdb_load(rtdb, geom,'geometry'))
     $      call errquit('setCoordinates: geom_rtdb_load failed',555,0)

      if (.not.geom_cart_coords_set(geom,coords))
     $      call errquit
     $      ('setCoordinates: geom_cart_coords_set failed',555)

      if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     $      call errquit('setCoordinates: geom_rtdb_store failed',555,0)

      if (.not.geom_destroy(geom))
     $      call errquit('setCoordinates: geom_destroy failed',555,0)

      return
      end

      subroutine NWChem_setBasisSet(basisSet)
      implicit none
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"

      external nw_inp_from_character
      logical nw_inp_from_character

      character*(*) basisSet
      character*32 temp_basis

      integer rtdb
      common /sidlinfo/ rtdb

      write(temp_basis, '(''basis; * library '',a,'' ; end'')')
     $      basisSet(1:inp_strlen(basisSet))

      if (.not. nw_inp_from_character(rtdb, temp_basis))
     $ call errquit('setBasisSet: nw_inp_from_character failed',555,0)

      return
      end

      subroutine NWChem_TaskEnergy(energy)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "inp.fh"

      double precision energy
      external task_energy
      logical task_energy

      integer rtdb
      common /sidlinfo/ rtdb

      if (.not.task_energy(rtdb))
     $     call errquit('NWChemMod:task_energy failed',555,0)

      if (.not.rtdb_get(rtdb, "task:energy", mt_dbl, 1, energy))
     $     call errquit('NWChemMod:rtdb_get failed',555,0)

      return
      end

      subroutine NWChem_TaskGradient(gradient)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "inp.fh"

      double precision gradient(*)
      external task_gradient
      logical task_gradient
      external task_energy
      logical task_energy

      integer rtdb
      common /sidlinfo/ rtdb
      integer ncent

      call NWChem_getNumCoordinates(ncent)

      if (.not.task_energy(rtdb))
     $     call errquit('NWChemMod:task_energy failed',555,0)

      if (.not.task_gradient(rtdb))
     $     call errquit('NWChemMod:task_gradient failed',555,0)

      if (.not.rtdb_get(rtdb, "task:gradient", mt_dbl, ncent,
     $     gradient))
     $     call errquit('NWChemMod:rtdb_get failed',555,0)

      return
      end

      subroutine NWChem_setCoordinatesFromFile(filename)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"

      logical nw_inp_from_file
      external nw_inp_from_file

      character*(nw_max_path_len) filename

      integer rtdb
      common /sidlinfo/ rtdb

      if (.not.nw_inp_from_file(rtdb,filename))
     $ call errquit('setCoordinatesFromFile: read inp from file failed'
     $ ,555)

      return
      end

      subroutine NWChem_setTheory(theoryName)
      implicit none
#include "rtdb.fh"

      integer rtdb
      common /sidlinfo/ rtdb

      character*(30) theoryName

      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theoryName))
     $     call errquit('setTheory: failed setting theory',555,0)

      return
      end

      subroutine NWChem_getCoordinates(coords)
      implicit none
#include "geom.fh"

      integer geom
      double precision coords(*)
      integer rtdb
      common /sidlinfo/ rtdb

      if (.not.geom_create(geom,'geometry'))
     $      call errquit('getCoordinates: geom_create failed',555,0)

      if (.not.geom_rtdb_load(rtdb, geom,'geometry'))
     $      call errquit('getCoordinates: geom_rtdb_load failed',555,0)

      if (.not.geom_cart_coords_get(geom,coords))
     $      call errquit
     $      ('getCoordinates: geom_cart_coords_get failed',555)

      if (.not.geom_destroy(geom))
     $      call errquit('getCoordinates: geom_destroy failed',555,0)

      return
      end

      subroutine NWChem_getNumCoordinates(ncent)
      implicit none
#include "geom.fh"

      integer geom
      integer ncent
      integer rtdb
      common /sidlinfo/ rtdb

      if (.not.geom_create(geom,'geometry'))
     $      call errquit('getNumCoordinates: geom_create failed',555,0)

      if (.not.geom_rtdb_load(rtdb, geom,'geometry'))
     $      call errquit('getNumCoordinates: geom_rtdb_load failed',
     $   555,0)

      if (.not.geom_ncent(geom,ncent))
     $      call errquit
     $      ('getNumCoordinates: geom_ncent failed',555)

      if (.not.geom_destroy(geom))
     $      call errquit('getNumCoordinates: geom_destroy failed',555,0)

      ncent=ncent*3;

      return
      end

      subroutine NWChem_TaskHessian(hessian)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "inp.fh"

      double precision hessian(*)
      external hess_read
      external task_hessian
      logical task_hessian
C      external task_energy
C      logical task_energy

      integer rtdb
      common /sidlinfo/ rtdb
      integer ncent

      call NWChem_getNumCoordinates(ncent)

!      if (.not.task_energy(rtdb))
!     $     call errquit('NWChem_TaskHessian:task_energy failed',555)

      if (.not.task_hessian(rtdb))
     $     call errquit('NWChem_TaskHessian:task_hessian failed',555,0)

      call hess_read(ncent/3,hessian)

      return
      end
