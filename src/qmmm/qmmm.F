c
c $Id: qmmm.F,v 1.11 2004-08-31 19:45:48 marat Exp $
c
c
      function qmmm_init(irtdb)
      implicit none
c
#include "inp.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
      logical qmmm_init
c
      integer nqm
      integer i_ind

      qmmm_init = .true.

c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)

c     get total number of quantum atoms
c     --------------------------------
      call mm_get_tot_nqm(nqm)

c     allocate and initialize qmmm data structure
c     -------------------------------------------
      call qmmm_data_alloc(nqm)
      call qmmm_data_init(irtdb)

c     establish mapping between mm and qm atoms
c     ----------------------------------------
      i_ind = qmmm_get_i_ai()
      call mm_get_ind(nqm,int_mb(i_ind))

c     create initial qm geometry
c     --------------------------
      call qmmm_create_qm_geom(irtdb)

      return
      end
c
      function qmmm_end(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"

      integer rtdb
      logical qmmm_end

      call mm_end()
      call qmmm_data_release()

      qmmm_end = .true.
      return
      end
c
      function qmmm_energy(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
c
      integer rtdb
      logical qmmm_energy
c     local variables
      logical status
      logical do_qm, do_mm
      logical mm_only, qm_only
      character*255 theory
      double precision energy
      double precision mm_energy
      double precision eatoms
c     external functions
      logical  task_energy_doit
      external task_energy_doit
      
      do_qm=.true.
      do_mm=.true.

      if (.not. rtdb_get(rtdb,'qmmm:mm_only',mt_log,1,mm_only)) then
        mm_only=.false.
      end if
      if (.not. rtdb_get(rtdb,'qmmm:qm_only',mt_log,1,qm_only)) then
        qm_only=.false.
      end if
      if(mm_only) then
       do_mm=.true.
       do_qm=.false.
      end if
      if(qm_only) then
       do_qm=.true.
       do_mm=.false.
      end if
c
c     sync qm coordinates with the geometry object
c     that was potentially modified outside qmmm (e.g.driver)     
      call qmmm_sync_qm_geom(rtdb)

      call qmmm_create_bq_geom(rtdb)

      qmmm_energy = .false.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

c     QM energy calculation
c     --------------------
      status = task_energy_doit(rtdb,theory,energy)
      
      call qmmm_restore_qm_geom(rtdb)
 
c     update QM coordinates for MM module
c     very important!!!
c     ----------------------------------
      call qmmm_update_mm_coord(rtdb)
 
c     single point MM calculation
c     ----------------------------------
       call mm_task_sp()

c     combine QM and MM energies
c     -------------------------
      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)

      if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy', 0, RTDB_ERR)

      mm_energy = mm_energy/cau2kj + 0.5d0*eatoms
      if(mm_only) then
        energy=0.0d0
      end if
      if(qm_only) then
        mm_energy=0.0d0
      end if

      energy = energy + mm_energy

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      qmmm_energy=status

      end
c
      function qmmm_gradient(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
      logical qmmm_gradient
      logical status
      logical task_gradient_doit
      logical do_qm, do_mm
      logical mm_only, qm_only
      external task_gradient_doit
      character*255 theory
      double precision energy
      double precision mm_energy
      integer i_g,h_g 
      integer i_gmm,h_gmm 
      integer i_ai
      double precision eatoms
      integer nqm
      integer ntot
c
      do_qm=.true.
      do_mm=.true.

      if (.not. rtdb_get(rtdb,'qmmm:mm_only',mt_log,1,mm_only)) then
        mm_only=.false.
      end if
      if (.not. rtdb_get(rtdb,'qmmm:qm_only',mt_log,1,qm_only)) then
        qm_only=.false.
      end if
      if(mm_only) then
       do_mm=.true.
       do_qm=.false.
      end if
      if(qm_only) then
       do_qm=.true.
       do_mm=.false.
      end if


      status = .true.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

c
c     sync qm coordinates with the geometry object
c     that was potentially modified outside qmmm (e.g.driver)     
      call qmmm_sync_qm_geom(rtdb)

c
c    first do quantum part
c    --------------------
      call qmmm_create_bq_geom(rtdb)

      nqm  = qmmm_get_nqm()
      ntot = qmmm_get_ntot()

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)

      if(.not.ma_alloc_get(MT_DBL, 3*nqm, 'qmmm mm grad array',
     &      h_gmm, i_gmm) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(3*nqm,0,dbl_mb(i_gmm),1)
 

      status = task_gradient_doit(rtdb,theory,energy,dbl_mb(i_g))

c
c     now do mm part
c     ----------------------------------
      call qmmm_update_mm_coord(rtdb)
      call mm_task_sp()
      i_ai = qmmm_get_i_ai()

      call dfill(3*nqm,0,dbl_mb(i_gmm),1)

      call mm_get_force(nqm,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_gmm))

      if(mm_only) then
        call dfill(3*nqm,0,dbl_mb(i_g),1)
      end if
      if(qm_only) then
        call dfill(3*nqm,0,dbl_mb(i_gmm),1)
      end if


      call daxpy(3*nqm,cau2nm/cau2kj,dbl_mb(i_gmm),1,dbl_mb(i_g),1)

      if(.not.rtdb_put(rtdb,'qmmm:gradient',mt_dbl,3*nqm,dbl_mb(i_g)))
     + call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)

      if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy', 0, RTDB_ERR)

      mm_energy = mm_energy/cau2kj + 0.5d0*eatoms
      if(mm_only) then
        energy=0.0d0
      end if
      if(qm_only) then
        mm_energy=0.0d0
      end if

      energy = energy + mm_energy

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      call qmmm_restore_qm_geom(rtdb)

      if(.not.ma_free_heap(h_gmm))
     & call errquit('qmmm g: Failed to deallocate heap',nqm,
     &       MA_ERR)
      if(.not.ma_free_heap(h_g))
     & call errquit('qmmm g: Failed to deallocate heap',nqm,
     &       MA_ERR)
 
      qmmm_gradient=status

      end
c
      function qmmm(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer rtdb
      logical qmmm
      character*30 operation
      logical status
      logical task_energy
      external task_energy
      logical qmmm_init
      external qmmm_init
      logical qmmm_end
      external qmmm_end
c
c
      status = .true.
      write(6,*) "initializing qmmm"
      if( .not. qmmm_init(rtdb)) then
       call errquit('qmmm_init',0,INPUT_ERR)
      end if
c
c     determine calculation type
c     --------------------------
      if(.not.rtdb_cget(rtdb,'task:operation',1,operation))
     + call errquit('qmmm: no calculation type',0,INPUT_ERR)
c
      if (inp_compare(.false., 'energy', operation)) then
         status = task_energy(rtdb)
      end if
      if( .not. qmmm_end(rtdb)) then
       call errquit('qmmm_end',0,INPUT_ERR)
      end if
c
      qmmm=status
      return
      end

      subroutine qmmm_update_mm_coord(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer nqm
      integer i_ai
      integer i_c

      nqm  = qmmm_get_nqm()
      i_ai = qmmm_get_i_ai()
      i_c = qmmm_get_i_c()


      call mm_set_geom(nqm,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_c))

      end

      subroutine qmmm_update_qm_forces(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
      integer irtdb
c     local variables
      integer nt
      integer i_ai
      integer i_gmm
      integer i_g

      nt   = qmmm_get_nqm()
      i_ai = qmmm_get_i_ai()
      i_g = qmmm_get_i_g()
      i_gmm = qmmm_get_i_gmm()

      call dfill(3*nt,0,dbl_mb(i_gmm),1)

      call mm_get_force(nt,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_gmm))


      if(.not.rtdb_get(irtdb,'qmmm:gradient',mt_dbl,3*nt,dbl_mb(i_g)))
     + call errquit('qmmm: rtdb_get gradient failed',0, RTDB_ERR)

      call daxpy(3*nt,cau2nm/cau2kj,dbl_mb(i_gmm),1,dbl_mb(i_g),1)

      if(.not.rtdb_put(irtdb,'qmmm:gradient',mt_dbl,3*nt,dbl_mb(i_g)))
     + call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

      end


