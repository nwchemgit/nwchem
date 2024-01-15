C
C $Id$
C
C This is a first attempt at an interface code for Venus written in F77.
C
      Subroutine initialize_qm(theory, basis, natom, 
     $    coordinates, labels, charge, printName, stack, heap,
     $    global,venusnw_file)

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "geninterface.fh"
#include "inp.fh"

      character*256 theory
      character*256 basis
      integer natom,I3N
      double precision charge
      double precision coordinates (3*natom)
      character*16 labels(natom)
      character*256 printName
      character*256 venusnw_file
      character*(nw_max_path_len) input_filename, rtdb_name

      integer stack, heap, global, i
      integer mdtob1,mitob1
      logical overify, ohardfail, ostartup, ocontinue

      logical setTheory, setBasisSet, setInitCoordinates, setCharge
      external setTheory, setBasisSet, setInitCoordinates, setCharge

      logical nw_inp_from_file
      external nw_inp_from_file
 
      logical setPrint
      external setPrint

c      write (6,*) 'theory=',theory,' basis=',basis,' natom=',natom,
c     $       ' charge=',charge,' printName=',printName,' stack=',stack,
c     $       ' heap=',heap,' global=',global,' venusnw_file=',
c     $       venusnw_file

C  Create parallel processes
c - raj
      stack = 29582399
      heap = 29582386
      global = 29582386
c----
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      mdtob1=MA_sizeof(MT_DBL,1,MT_BYTE)
      call ga_brdcst(MT_INT, natom, mitob1, 0)
      call ga_brdcst(MT_DBL, coordinates, mdtob1*3*natom, 0)
C change MT_INT-256/4 to MT_BYTE-256
      call ga_brdcst(MT_BYTE, theory, 256, 0)
      call ga_brdcst(MT_BYTE, basis, 256, 0)
      call ga_brdcst(MT_INT, labels, 16*natom*mitob1/4, 0)
C      write(6,*) labels(1), labels(2), labels(3),labels(4)
C     * ,labels(5),labels(6)
C      call flush(6)
C YA: Broadcast charge
      call ga_brdcst(MT_DBL, charge, mdtob1, 0)
      call ga_brdcst(MT_BYTE, printName, 256, 0)
C PA Broadcast Memory parameters
      call ga_brdcst(MT_INT, stack, mitob1, 0)
      call ga_brdcst(MT_INT, heap, mitob1, 0)
      call ga_brdcst(MT_INT, global, mitob1, 0)
c--raj
       write (6,*) 'stack=',stack,' heap=',heap,' global=',global
c--
      numatoms = natom
C  Worry about input file.  For now use a dummy one.  Do I need to open it?

      input_filename = venusnw_file
      if (ga_nodeid().eq.0) then
         open(LuIn, file=input_filename, status='old', form='formatted',
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

C  Set up memory environment

!      call input_mem_size(stack, heap, global, overify, ohardfail)

C YA: Hardwire memory limits assuming total equal 1800.0 Mbytes
c       stack = 58982399
c       heap = 58982386
c       global = 117964800

       overify = .true.
       ohardfail = .true.
c       write (6,*) 'stack=',stack,' heap=',heap,' global=',global,
c     *            ' overify=', overify,' ohardfail=',ohardfail
!      call ga_initialize()
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

c     set no. threads for threaded BLAS to 1
      call util_blas_set_num_threads(1)

C  Open rtdb
C  Need to get the input file name and the rtdb_name taken care of

!     sets permanent_ and scratch_dir
!      rtdb_name = "junk.rtdb"
      call input_file_info(input_filename, rtdb_name,
     $                     ostartup, ocontinue)
      ostartup = .true.
      ocontinue = .false.
      if (.not. rtdb_open(rtdb_name, 'empty', rtdb))
     $     call errquit('initialize_qm: rtdb_open empty failed ',
     $   555,0)
      call util_set_rtdb_state(rtdb,.true.,.false.,.false.)
      call util_file_info_rtdb(rtdb)
      call geom_hnd_parallel(.true.)

c      write (6,*) 'In geninterface.F: venusnw_file=',venusnw_file
!  Test for .nw file from VENUS
      if (venusnw_file.eq.' ') then
       write (6,*) 'No NWChem file name is detected in VENUS input'
       write (6,*) 'Using only VENUS input data for NWChem init'

C  Process number of scf iterations

      if (.not. rtdb_put(rtdb, 'scf:maxiter', mt_int, 1, 200))
     $     call errquit('scf-input: rtdb_put failed', 0,0)

C  Process number of dft iterations

      if (.not. rtdb_put(rtdb, 'dft:iterations',
     &     mt_int, 1, 200))
     &     call errquit('dft_rdinput: rtdb_put failed', 101,0)

C  Process theory type and put into rtdb

      if (.not. setTheory(theory))
     $     call errquit('initialize_qm: setTheory failed', 555,0)

C  Process basis type and put into rtdb

      if (.not. setBasisSet(basis))
     $     call errquit('initialize_qm: setBasisSet failed', 555,0)

C  Process charge and put into rtdb

      if (.not. setCharge(charge))
     $     call errquit('initialize_qm: setCharge failed', 555,0)

C  Process print type and put into rtdb

c      write (6,*) 'printName=',printName, ga_nodeid()
C      i = rtdb_parallel(0)
C      if (ga_nodeid().eq.0) then
       if (.not. setPrint(printName))
     $     call errquit('initialize_qm: setPrint failed', 555)
C      end if
C      i = rtdb_parallel(i)
c      write (6,*) 'after printName=',printName, ga_nodeid()

C  Process molecule input and put into rtdb

      if (.not. setInitCoordinates(coordinates,labels,natom))
     $     call errquit
     $     ('initialize_qm:setInitCoordinates failed',555)

c      initialize_qm = .true.
c       write(6,*)'Initalize called'

      else
       write (6,*) 'Reading .nw file'
       if (.not.nw_inp_from_file(rtdb,venusnw_file))
     $ call errquit ('nw_inp_from_file: Error opening file',555)

C  Process molecule input and put into rtdb

      if (.not. setInitCoordinates(coordinates,labels,natom))
     $     call errquit
     $     ('initialize_qm:setInitCoordinates failed',555)

      end if

      return
      end
C
C
      logical function finalize_qm()

      implicit none

#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "geninterface.fh"

      finalize_qm = .false.

C  Close the rtdb

      if (.not. rtdb_close(rtdb, 'delete'))
     $     call errquit('finalize_qm: rtdb_close failed', rtdb,0)

C  Close the global memory

C      call ga_sync()
C      call ga_terminate()

C  Close the parallel environment

C      call pend()

      finalize_qm = .true.

      return
      end

C
C
      logical function setBasisSet(basisSet)

      implicit none

#include "inp.fh"
#include "util.fh"
#include "geninterface.fh"
      character*(*) basisSet
      character*32 temp_basis

      external nw_inp_from_character
      logical nw_inp_from_character

      setBasisSet = .false.

      temp_basis = basisSet(1:inp_strlen(basisSet))
	write(6,*)'basis in setbasis set function',temp_basis 
      write(basisSet, '(''basis; * library '',a,'' ; end'')')
     $     temp_basis(1:inp_strlen(temp_basis))

      if (.not. nw_inp_from_character(rtdb, basisSet))
     $     call errquit('setBasisSet: setting rtdb failed', 555,0)

      setBasisSet = .true.

      return 
      end
C
C
      logical function setTheory(theoryName)

      implicit none
#include "geninterface.fh"
#include "rtdb.fh"
#include "util.fh"

      character*(*) theoryName
      character*32 dft_func

      external nw_inp_from_character
      logical nw_inp_from_character

      setTheory = .false.

      if (theoryName.eq.'b3lyp') then 
         write(dft_func, '(''dft; XC b3lyp; end'')')
         if (.not. nw_inp_from_character(rtdb, dft_func))
     $   call errquit('setTheoryB3LYP: setting rtdb failed', 555,0)
         if (.not. rtdb_cput(rtdb, 'task:theory', 1, 'dft'))
     $     call errquit('setTheory: failed setting theory',555,0)

      else

         if (.not. rtdb_cput(rtdb, 'task:theory', 1, theoryName))
     $     call errquit('setTheory: failed setting theory',555,0)
      end if

      setTheory = .true.

      return
      end
C
C
      logical function setCoordinates(coords)

      implicit none

#include "geom.fh"
#include "geninterface.fh"

      integer geom
      character*12 units
      double precision scale
      double precision coords(*)
      integer i

      setCoordinates = .false.

      if (.not.geom_create(geom,'geometry'))
     $      call errquit('setCoordinates: geom_create failed',555,0)

      if (.not.geom_rtdb_load(rtdb, geom,'geometry'))
     $      call errquit('setCoordinates: geom_rtdb_load failed',555,0)

      if (.not. geom_set_user_units(geom,'angstroms'))
     $     call errquit
     $     ('setInitCoordinates: failed setting user units',0)

      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit
     $     ('setInitCoordinates:: failed getting user scale',0)

c      write(6,*) 'numatoms is in setCoordinates', numatoms
      call flush(6)
      if (scale.ne.1.0d0) 
     $   call dscal(3*numatoms, scale, coords, 1)

      if (.not.geom_cart_coords_set(geom,coords))
     $      call errquit
     $      ('setCoordinates: geom_cart_coords_set failed',555)

      if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     $      call errquit('setCoordinates: geom_rtdb_store failed',555,0)

      if (.not.geom_destroy(geom))
     $      call errquit('setCoordinates: geom_destroy failed',555,0)

c      write(6,*) 'numatoms at end of setCoordinates', numatoms
      call flush(6)
      setCoordinates = .true.

      return
      end
C
C
      logical function setInitCoordinates(coordinates,labels,natoms)

      implicit none

#include "geninterface.fh"
#include "geom.fh"
#include "rtdb.fh"

      integer natoms, geom, atn, i, j
      double precision coordinates(3,natoms)
      character*16 labels(natoms), element
      character*2 symbol
      character*12 units
      double precision scale
      double precision masses(natoms), charges(natoms)
      logical is_atom

      setInitCoordinates = .false.

      if (.not.geom_create(geom,'geometry'))
     $      call errquit('setCoordinates: geom_create failed',555,0)
      if (.not. geom_set_user_units(geom,'angstroms'))
     $     call errquit
     $     ('setInitCoordinates: failed setting user units',0)
      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit
     $     ('setInitCoordinates:: failed getting user scale',0)

c      if (.not.geom_rtdb_load(rtdb, geom,'geometry'))
c     $      call errquit
c     $      ('setInitCoordinates: geom_rtdb_load failed',555)

C  Map labels to atoms and charges

      do i = 1, natoms
        is_atom = geom_tag_to_element(labels(i), symbol,
     $             element, atn)

        if (.not. is_atom) call errquit
     $        ('setInitCoordinates:center is not an atom',555)

        if (.not.
     $       geom_atn_to_default_mass(atn,masses(i)))
     $       call errquit('setInitCoordinates: default mass failed',
     $       555)

        charges(i) = atn
      enddo

      call geom_check_input_quants
     $     (labels,charges,natoms,1.0d-12,'charge')

      call geom_check_input_quants
     $     (labels,masses,natoms,1.0d-8,'mass')

      if (scale.ne.1.0d0) 
     $   call dscal(3*natoms, scale, coordinates, 1)

      if (.not.geom_cart_set(geom,natoms,labels,coordinates,charges))
     $     call errquit
     $     ('setInitCoordinates: geom_cart_coords_set failed',555)

      if (.not. geom_masses_set(geom, natoms, masses))
     $     call errquit
     $     ('setInitCoordinates: geom_masses_set failed',555)

      if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     $     call errquit
     $     ('setInitCoordinates: geom_rtdb_store failed',555)

      if (.not.geom_destroy(geom))
     $     call errquit('setInitCoordinates: geom_destroy failed',555,0)

      setInitCoordinates = .true.

      return
      end
C
      logical function evaluateObjective(coords, energy)

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geninterface.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"

      common/threads/myid
      integer myid 
      double precision coords(*), energy
      integer i
      external setCoordinates, task_energy
      logical setCoordinates, task_energy
      character*256 movecs_in

      evaluateObjective = .false.

      call ga_brdcst(MT_DBL, coords,
     $     MA_sizeof(MT_DBL,(3*numatoms),MT_BYTE), 0)
      call flush(6)

      if (.not.setCoordinates(coords))
     $     call errquit('evaluateObjective: setCoordinates failed',
     $    555,0)

      if (.not.task_energy(rtdb)) then
         write (6,*) 'SCF failed running rest.scf script'
         movecs_in = 'atomic'
        if (.not. rtdb_cput(rtdb, 'dft:input vectors', 1, movecs_in))
     $  call errquit('name_rout: failed setting movecs_in',0,0)
        if (.not.task_energy(rtdb))
     &      call errquit('task_energy failed',555,0)
      endif

      if (.not.rtdb_get(rtdb, "task:energy", mt_dbl, 1, energy))
     $     call errquit('evaluateObjective:rtdb_get failed',555,0)

      evaluateObjective = .true.

      return
      end
C
      logical function evaluateGradient(coords, gradient)

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geninterface.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "tcgmsg.fh"

      common/threads/myid
      integer myid
      double precision coords(3*numatoms), gradient(3*numatoms)
      integer i

      external setCoordinates, task_gradient
      logical setCoordinates, task_gradient
      external task_energy
      logical task_energy
      character*256 movecs_in 

      evaluateGradient = .false.

      call ga_brdcst(MT_DBL, coords, mdtob(3*numatoms), 0)

      if (.not.setCoordinates(coords))
     $     call errquit('evaluateGradient: setCoordinates failed',555,0)

cc
      if (.not.task_energy(rtdb)) then
        write (6,*) 'SCF failed running rest.scf script'
        movecs_in = 'atomic'
        if (.not. rtdb_cput(rtdb, 'dft:input vectors', 1, movecs_in))
     $      call errquit('name_rout: failed setting movecs_in',0,0)
        if (.not.task_energy(rtdb))
     &      call errquit('task_energy failed',555,0)
      endif
      if (.not.task_gradient(rtdb))
     $     call errquit('evaluateGradient:task_gradient failed',555,0)

C      if (.not.rtdb_print(rtdb, .true.))
C     $     call errquit('evaluateGradient:rtdb_print failed',555)
C
      if (.not.rtdb_get(rtdb, "task:gradient", mt_dbl, numatoms*3, 
     $     gradient))
     $     call errquit('evaluateGradient:rtdb_get failed',555,0)

c
c       call ga_sync()

      evaluateGradient = .true.

      return
      end
C
      logical function evaluateObjectiveAndGradient(coords, 
     $       energy, gradient)

      implicit none
      integer i

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geninterface.fh"
c

      double precision coords(*), energy, gradient(*)

      external evaluateGradient
      logical evaluateGradient

      evaluateObjectiveAndGradient = .false.
c
      if (.not.evaluateGradient(coords, gradient))
     $     call errquit
     $     ('evaluateObjectiveAndGradient:evaluateGradient failed',555)

      if (.not.rtdb_get(rtdb, "task:energy", mt_dbl, 1, energy))
     $     call errquit
     $     ('evaluateObjectiveAndGradient:rtdb_get failed',555)
c
      evaluateObjectiveAndGradient = .true.

      return
      end

      logical function evaluateHessian(coords,
     $       hessian)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geninterface.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "tcgmsg.fh"
c--- added 11/15/05
#include "msgids.fh"
#include "msgtypesf.h"
c---

      double precision coords(*), hessian(*)

      external hess_read
      external task_hessian
      external setCoordinates
      logical task_hessian
      logical setCoordinates
c
      call ga_brdcst(MT_DBL, coords, mdtob(3*numatoms), 0)
c
      evaluateHessian = .false.

      if (.not.setCoordinates(coords))
     $     call errquit('NWChem_TaskHessian: setCoordinates failed',
     $    555,0)

      if (.not.task_hessian(rtdb))
     $     call errquit('NWChem_TaskHessian:task_hessian failed',555,0)

c--- changed 11/15/05 
       if (ga_nodeid().eq.0) then
         call hess_read(numatoms, hessian)
       endif
c
c  Broadcast results to all of the processors 
c
       call ga_brdcst(msg_hess_exy+MSGDBL, hessian,
     &               mdtob(3*numatoms*3*numatoms), 0)
c---
      evaluateHessian = .true.

      return
      end

      subroutine net_init (myid, nnodes1)

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "util.fh"
#include "geninterface.fh"

      integer myid       
      integer nnodes1
      call pbeginf()
      call ga_initialize()
      myid=nodeid()
      nnodes1=nnodes()

      return
      end

      subroutine net_final() 
      
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "util.fh"
#include "geninterface.fh"

C      call ga_sync()
      call ga_terminate()

C  Close the parallel environment

      call pend()

      return
      end

CBERT subroutine net_broadcast(myid,idop)
      subroutine net_broadcast(idop)
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "geninterface.fh"

      integer myid,idop

      call ga_brdcst(MT_INT, idop, MA_sizeof(MT_INT,1,MT_BYTE), 0)

      return
      end

C
C
      logical function setCharge(charge)

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "util.fh"
#include "geninterface.fh"
#include "inp.fh"

      double precision charge

      setCharge = .false.

      if (.not. rtdb_put(rtdb, 'charge', MT_DBL, 1, charge))
     $     call errquit
     $     ('input_charge: failed to write to rtdb', 0)

      setCharge = .true.

      return
      end

      logical function setPrint(printName)
      implicit none
#include "geninterface.fh"
#include "rtdb.fh"
#include "util.fh"

      character*256 printName
      integer level

      setPrint = .false.

      if (printName.eq."none") then
         level=print_none
      else if (printName.eq."low") then
         level=print_low
      else if (printName.eq."high") then
         level=print_high
      else if (printName.eq."debug") then
         level=print_debug
      else
         level=print_medium
      end if

      call util_print_set_level(level)

      setPrint = .true.

      return
      end


      subroutine clearOld()

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geninterface.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "tcgmsg.fh"

      common/threads/myid
      integer myid
      character*256 movecs_in


      movecs_in = 'atomic'
      if (.not. rtdb_cput(rtdb, 'dft:input vectors', 1, movecs_in))
     $  call errquit('name_rout: failed setting movecs_in',0,0)

      return
      end

