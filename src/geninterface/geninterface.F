C
C $Id: geninterface.F,v 1.5 2003-10-09 23:57:05 yuri Exp $
C
C This is a first attempt at an interface code for Venus written in F77.
C
      Subroutine initialize_qm(theory, basis, natom, 
     $    coordinates, labels, charge, printName)

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "util.fh"
#include "geninterface.fh"

      character*256 theory
      character*256 basis
      integer natom,I3N
      double precision charge
      double precision coordinates (3*natom)
      character*16 labels(natom)
      character*256 printName
      character*(nw_max_path_len) input_filename, rtdb_name

      integer stack, heap, global, i 
      logical overify, ohardfail, ostartup, ocontinue

      logical setTheory, setBasisSet, setInitCoordinates, setCharge
      external setTheory, setBasisSet, setInitCoordinates, setCharge
 
      logical setPrint
      external setPrint

C  Create parallel processes

      call ga_brdcst(MT_INT, natom, mitob(1), 0)
      write (6,*) 'after natom ', natom, ga_nodeid()
      write (6,*) '3*natom ', 3*natom, ga_nodeid()
      call flush(6)
      call ga_brdcst(MT_DBL, coordinates, mdtob(3*natom), 0)
      write (6,*) 'after coordinates ', ga_nodeid()
      I3N = 3*natom
      do i = 1,I3N
         write(6,*) coordinates(i), i, ga_nodeid()
      enddo
      call flush(6)
      call ga_brdcst(MT_INT, theory, 256/4, 0)
      write (6,*) 'after theory ', theory, ga_nodeid()
      call flush(6)
      call ga_brdcst(MT_INT, basis, 256/4, 0)
      write (6,*) 'after basis ', basis, ga_nodeid()
      call flush(6)
      call ga_brdcst(MT_INT, labels, 16*mitob(natom)/4, 0)
      write (6,*) 'after labels ', ga_nodeid()
      write(6,*) labels(1), labels(2), labels(3),labels(4)
     * ,labels(5),labels(6)
      call flush(6)
C YA: Broadcast charge
      call ga_brdcst(MT_DBL, charge, mdtob(1), 0)
      call ga_brdcst(MT_INT,printName, 256/4, 0)

      numatoms = natom
C  Worry about input file.  For now use a dummy one.  Do I need to open it?

      input_filename = "junk.nw"
c      if (nodeid().eq.0) then
c         open(LuIn, file=input_filename, status='old', form='formatted',
c     $        err=10000)
c         goto 10001
c10000    call errquit('nwchem: failed to open the input file', 0)
c10001    continue
c
c      else
c#if !(defined(KSR) || defined(IBM) || defined(FUJITSU_SOLARIS) ||defined(FUJITSU_VPP) ||defined(FUJITSU_VPP64))
c         close(LuIn)
c#endif
c      endif

C  Set up memory environment

!      call input_mem_size(stack, heap, global, overify, ohardfail)

C YA: Hardwire memory limits assuming total equal 1800.0 Mbytes
       stack = 58982399
       heap = 58982386
       global = 117964800

       overify = .true.
       ohardfail = .true.
       write (6,*) 'stack=',stack,' heap=',heap,' global=',global,
     *            ' overify=', overify,' ohardfail=',ohardfail
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

C  Open rtdb
C  Need to get the input file name and the rtdb_name taken care of

!      call input_file_info(input_filename, rtdb_name, 
!     $                     ostartup, ocontinue)
      rtdb_name = "junk.rtdb"
      ostartup = .true.
      ocontinue = .false.
      if (.not. rtdb_open(rtdb_name, 'empty', rtdb))
     $     call errquit('initialize_qm: rtdb_open empty failed ', 555)
      call util_set_rtdb_state(rtdb,.true.,.false.,.false.)
      call util_file_info_rtdb(rtdb)
      call geom_hnd_parallel(.true.)

C  Process number of scf iterations

      if (.not. rtdb_put(rtdb, 'scf:maxiter', mt_int, 1, 200))
     $     call errquit('scf-input: rtdb_put failed', 0)

C  Process number of dft iterations

      if (.not. rtdb_put(rtdb, 'dft:iterations',
     &     mt_int, 1, 200))
     &     call errquit('dft_rdinput: rtdb_put failed', 101)

C  Process theory type and put into rtdb

      if (.not. setTheory(theory))
     $     call errquit('initialize_qm: setTheory failed', 555)

C  Process basis type and put into rtdb

      if (.not. setBasisSet(basis))
     $     call errquit('initialize_qm: setBasisSet failed', 555)

C  Process charge and put into rtdb

      if (.not. setCharge(charge))
     $     call errquit('initialize_qm: setCharge failed', 555)

C  Process molecule input and put into rtdb

      if (.not. setInitCoordinates(coordinates,labels,natom))
     $     call errquit
     $     ('initialize_qm:setInitCoordinates failed',555)

c      initialize_qm = .true.
       write(6,*)'Initalize called'
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
     $     call errquit('finalize_qm: rtdb_close failed', rtdb)

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
     $     call errquit('setBasisSet: setting rtdb failed', 555)

      setBasisSet = .true.

      return 
      end
C
C
      logical function setTheory(theoryName)

      implicit none
#include "geninterface.fh"
#include "rtdb.fh"

      character*(*) theoryName

      setTheory = .false.

      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theoryName))
     $     call errquit('setTheory: failed setting theory',555)

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
     $      call errquit('setCoordinates: geom_create failed',555)

      if (.not.geom_rtdb_load(rtdb, geom,'geometry'))
     $      call errquit('setCoordinates: geom_rtdb_load failed',555)

      if (.not. geom_set_user_units(geom,'angstroms'))
     $     call errquit
     $     ('setInitCoordinates: failed setting user units',0)

      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit
     $     ('setInitCoordinates:: failed getting user scale',0)

      write(6,*) 'numatoms is in setCoordinates', numatoms
      call flush(6)
      if (scale.ne.1.0d0) 
     $   call dscal(3*numatoms, scale, coords, 1)

      if (.not.geom_cart_coords_set(geom,coords))
     $      call errquit
     $      ('setCoordinates: geom_cart_coords_set failed',555)

      if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     $      call errquit('setCoordinates: geom_rtdb_store failed',555)

      if (.not.geom_destroy(geom))
     $      call errquit('setCoordinates: geom_destroy failed',555)

      write(6,*) 'numatoms at end of setCoordinates', numatoms
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
     $      call errquit('setCoordinates: geom_create failed',555)
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
     $     call errquit('setInitCoordinates: geom_destroy failed',555)

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
#include "tcgmsg.fh"

      double precision coords(*), energy
      integer i
      external setCoordinates, task_energy
      logical setCoordinates, task_energy
   
      evaluateObjective = .false.

      call ga_brdcst(MT_DBL, coords, mdtob(3*numatoms), 0)
      call flush(6)

      if (.not.setCoordinates(coords))
     $     call errquit('evaluateObjective: setCoordinates failed',555)

      if (.not.task_energy(rtdb))
     $     call errquit('evaluateObjective:task_energy failed',555)

      if (.not.rtdb_get(rtdb, "task:energy", mt_dbl, 1, energy))
     $     call errquit('evaluateObjective:rtdb_get failed',555)

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

      double precision coords(3*numatoms), gradient(3*numatoms)
      integer i

      external setCoordinates, task_gradient
      logical setCoordinates, task_gradient
      external task_energy
      logical task_energy

      evaluateGradient = .false.

      call ga_brdcst(MT_DBL, coords, mdtob(3*numatoms), 0)

      if (.not.setCoordinates(coords))
     $     call errquit('evaluateGradient: setCoordinates failed',555)

      if (.not.task_energy(rtdb))
     $     call errquit('evaluateGradient:task_gradient failed',555)
      if (.not.task_gradient(rtdb))
     $     call errquit('evaluateGradient:task_gradient failed',555)

C      if (.not.rtdb_print(rtdb, .true.))
C     $     call errquit('evaluateGradient:rtdb_print failed',555)
C
      if (.not.rtdb_get(rtdb, "task:gradient", mt_dbl, numatoms*3, 
     $     gradient))
     $     call errquit('evaluateGradient:rtdb_get failed',555)

      evaluateGradient = .true.

      return
      end
C
      logical function evaluateObjectiveAndGradient(coords, 
     $       energy, gradient)

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geninterface.fh"

      double precision coords(*), energy, gradient(*)

      external evaluateGradient
      logical evaluateGradient

      evaluateObjectiveAndGradient = .false.

      if (.not.evaluateGradient(coords, gradient))
     $     call errquit
     $     ('evaluateObjectiveAndGradient:evaluateGradient failed',555)

      if (.not.rtdb_get(rtdb, "task:energy", mt_dbl, 1, energy))
     $     call errquit
     $     ('evaluateObjectiveAndGradient:rtdb_get failed',555)

      evaluateObjectiveAndGradient = .true.

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

      subroutine net_broadcast(myid,idop)
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "util.fh"
#include "geninterface.fh"

      integer myid,idop

      call ga_brdcst(MT_INT, idop, mitob(1), 0)

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

