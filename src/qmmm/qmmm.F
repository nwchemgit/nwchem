c
c $Id: qmmm.F,v 1.17 2004-09-13 19:57:14 marat Exp $
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
      qmmm_init = .true.

c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)

c     allocate and initialize qmmm data structure
c     -------------------------------------------
      call qmmm_data_init(irtdb)

c     initialize bq settings 
      call qmmm_bq_data_init(irtdb)

c     create initial qm geometry
c     --------------------------
      call qmmm_geom_init(irtdb)

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
c     ------------------------------------------------------
      call qmmm_geom_prepare(rtdb)

      qmmm_energy = .false.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

c     QM energy calculation
c     --------------------
      status = task_energy_doit(rtdb,theory,energy)
      call qmmm_geom_finish(rtdb)
 
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
      integer nqm,nbq
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
c     ------------------------------------------------------
      call qmmm_geom_prepare(rtdb)

      nqm  = qmmm_get_nqm()
      nbq  = qmmm_get_nbq()
      ntot = nqm+nbq

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
      call qmmm_geom_finish(rtdb)
      call qmmm_update_mm_coord(rtdb)
      call mm_task_sp()
      i_ai = qmmm_get_i_ai()

      call mm_get_solute_force_quant(nqm,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_gmm))
 
c
      if(mm_only) then
        call dfill(3*nqm,0,dbl_mb(i_g),1)
      end if
      if(qm_only) then
        call dfill(3*nqm,0,dbl_mb(i_gmm),1)
      end if


      call daxpy(3*nqm,1.0d0,dbl_mb(i_gmm),1,dbl_mb(i_g),1)


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


      if(.not.ma_free_heap(h_gmm))
     & call errquit('qmmm g: Failed to deallocate heap',nqm,
     &       MA_ERR)
      if(.not.ma_free_heap(h_g))
     & call errquit('qmmm g: Failed to deallocate heap',nqm,
     &       MA_ERR)
 
      qmmm_gradient=status

      end
c
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

      call mm_set_solute_coord_quant(nqm,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_c))

      end



