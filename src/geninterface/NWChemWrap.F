C $Id: NWChemWrap.F,v 1.4 2006-06-13 06:31:26 windus Exp $
      subroutine NWChem_nwchemStart(scratch_dir, permanent_dir)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "util.fh"
#include "inp.fh"


      character*(nw_max_path_len) input_filename, rtdb_name
      character*(nw_max_path_len) input_filename0
      character*(*) scratch_dir, permanent_dir
      integer stack, heap, global
      logical overify, ohardfail, ostartup, ocontinue
      integer rtdb
      common /sidlinfo/ rtdb
      character*(nw_max_path_len) temp_directory1, temp_directory2
      integer init_counter, sgroupid
      common /global_counter/ init_counter
      logical ma_status

C  Create parallel processes
!      call pbeginf()

      sgroupid=ga_pgroup_get_default()
      write (6,*) 'sgroupid=',sgroupid
      write (6,*) 'LuIn=',LuIn

C  Create unique input file name
      write(input_filename0, '('''',a,''/junk.'')')
     $     permanent_dir(1:inp_strlen(permanent_dir))

C temporary hack to resolve file name issues      
      if(sgroupid.eq.-1) then
         write(input_filename, '('''',a,i2,''.nw'')')
     $        input_filename0(1:inp_strlen(input_filename0)),sgroupid
      else if(sgroupid.gt.9) then
         write(input_filename, '('''',a,i2,''.nw'')')
     $        input_filename0(1:inp_strlen(input_filename0)),sgroupid
      else if(sgroupid.gt.99) then
         write(input_filename, '('''',a,i3,''.nw'')')
     $        input_filename0(1:inp_strlen(input_filename0)),sgroupid
      else
         write(input_filename, '('''',a,i1,''.nw'')')
     $        input_filename0(1:inp_strlen(input_filename0)),sgroupid
      endif
      
      write(6,*) 'input file name =', input_filename
      write(6,*) 'my node id is', ga_nodeid()
    
C  Worry about input file.  For now use a dummy one
c      write(input_filename, '('''',a,''/junk.nw'')')
c     $      permanent_dir(1:inp_strlen(permanent_dir))
      if (ga_nodeid().eq.0) then
         open(LuIn, file=input_filename, status='unknown', 
     $        form='formatted',
     $        err=10000)
         goto 10001
10000    call errquit('nwchem: failed to open the input file', 0)
10001    continue
c
      else
#if !(defined(KSR) || defined(IBM) || defined(FUJITSU_SOLARIS) ||defined(FUJITSU_VPP) ||defined(FUJITSU_VPP64))
         close(LuIn)
#endif
      endif

C  Write scratch directory path into dummy input file

      write(temp_directory1, '(''scratch_dir '',a,'''')')
     $      scratch_dir(1:inp_strlen(scratch_dir))

      if (ga_nodeid().eq.0) then
         write(6,*) 'writing to file:', input_filename
         write (LuIn,*) temp_directory1
      end if

      write(temp_directory2, '(''permanent_dir '',a,'''')')
     $      permanent_dir(1:inp_strlen(permanent_dir))

      if (ga_nodeid().eq.0) then
         write(6,*) 'writing to file:', input_filename
         write (LuIn,*) temp_directory2
      end if

C     Set up memory environment
      write (6,*) 'Set up memory environment'
      call input_mem_size(stack, heap, global, overify, ohardfail)
c      if (.not.(init_counter.eq.1)) then
!     call ga_initialize()
!     write (6,*) 'Called ga_initialize inside NWChemWrap.F'
      if ( ga_uses_ma() ) then
         ma_status=ma_init(mt_dbl, stack, heap+global) 
      else
         ma_status=ma_init(mt_dbl,stack,heap)
         call ga_set_memory_limit(ma_sizeof(mt_dbl,global,mt_byte))
      endif
      
C     Open rtdb
      write(6,*) 'Opening rtdb'
C     Need to get the input file name and the rtdb_name taken care of
      
      call input_file_info(input_filename, rtdb_name,
     $     ostartup, ocontinue)
      write(6,*) 'Before rtdb_open'
      if (.not. rtdb_open(rtdb_name, 'empty', rtdb))
     $     call errquit('initialize_qm: rtdb_open empty failed ',555,0)
      write(6,*) 'After rtdb_open'
      call util_set_rtdb_state(rtdb,.true.,.false.,.false.)
      call util_file_info_rtdb(rtdb)
      write(6,*) 'Before geom_hnd_parallel'
      call geom_hnd_parallel(.true.)
      
      write(6,*) 'end NWChem_nwchemStart'
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
!      call ga_sync()
!      call ga_terminate()

C  Close the parallel environment
!      call pend()

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
      character*30 temp_basis

      integer rtdb
      common /sidlinfo/ rtdb

      write(6,*) 'In NWChem_setBasisSet(): basisSet=',basisSet
      write(temp_basis, '(''basis; * library '',a,'' ; end'')')
     $      basisSet(1:inp_strlen(basisSet))
      
      write(6,*) 'In NWChem_setBasisSet(): temp_basis = ', temp_basis
      write(6,*) 'In NWChem_setBasisSet(): rtdb = ', rtdb
      if (.not. nw_inp_from_character(rtdb, temp_basis))
     $ call errquit('setBasisSet: nw_inp_from_character failed',555,0)

      write(6,*) 'End NWChem_setBasisSet()'
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

C      if (.not.task_energy(rtdb))
C     $     call errquit('NWChemMod:task_energy failed',555,0)

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

      character*(*) filename
      character*(nw_max_path_len) input_filename
      integer rtdb
      common /sidlinfo/ rtdb

      write(6,*) 'In NWChem_setCoordinatesFromFile'
      write(6,*) 'filename = ', filename
      if (.not.nw_inp_from_file(rtdb, filename))
     $ call errquit('setCoordinatesFromFile: read inp from file failed'
     $     ,555)

      return
      end

      subroutine NWChem_setTheory(theoryName)
      implicit none
#include "rtdb.fh"

      integer rtdb
      common /sidlinfo/ rtdb

      character*(*) theoryName

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
