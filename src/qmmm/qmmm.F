c
c $Id: qmmm.F,v 1.29 2004-10-18 23:10:19 marat Exp $
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
#include "global.fh"
      integer irtdb
      logical qmmm_init
c
      character*32 pname
c
      pname = "qmmm_init"


      write(*,*) ga_nodeid(),"in "//pname
    
      qmmm_init = .true.

c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)
      call qmmm_init_qm(irtdb)

      if (.not. rtdb_put(irtdb,'qmmm:init',mt_log,1,.true.))
     $     call errquit('qmmm: failed qmmm:init', 0, RTDB_ERR)

      if(qmmm_master()) then
        call qmmm_print_header(irtdb)
        call qmmm_print_params(irtdb)
      end if

      write(*,*) ga_nodeid(),"out "//pname

      return
      end
c
      subroutine qmmm_init_qm(irtdb)
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
c
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
#include "errquit.fh"

      integer rtdb
      logical qmmm_end
      qmmm_end = .true.
c
      call mm_end()
      call qmmm_data_release()
      call qmmm_bq_data_dealloc()

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
#include "qmmm.fh"
c
      integer rtdb
      logical qmmm_energy
c     local variables
      logical status
      character*255 theory
      double precision energy
      double precision qm_energy
      double precision mm_energy
      double precision eatoms
c     external functions
      logical  task_energy_doit
      external task_energy_doit
      
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      call qmmm_get_qm_coord_qm(rtdb)
      call qmmm_push_qm_coord_mm(rtdb)
      call mm_task_sp()

      call qmmm_create_bq_geom(rtdb)

c     QM energy calculation
c     --------------------
      status = task_energy_doit(rtdb,theory,qm_energy)
c
      call qmmm_geom_finish(rtdb)

c     combine QM and MM energies
c     -------------------------
      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)
      mm_energy = mm_energy/cau2kj
      energy = qm_energy + mm_energy-eatoms

      if (.not. rtdb_put(rtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)


      call qmmm_print_energy(rtdb)

      qmmm_energy=status

      end

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
      external task_gradient_doit
      character*255 theory
      double precision energy
      double precision mm_energy
      double precision qm_energy
      integer i_g,h_g 
      integer i_gmm,h_gmm 
      integer i_ai
      double precision eatoms
      integer nqm,nbq
      logical link
      integer ntot
      integer i


      link = .not.qmmm_h_link()


      status = .true.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      call qmmm_get_qm_coord_qm(rtdb)
      call qmmm_push_qm_coord_mm(rtdb)
      call mm_task_sp()

      call qmmm_create_bq_geom(rtdb)

      nqm  = qmmm_get_nqm()
      nbq  = qmmm_get_nbq()
      ntot = nqm+nbq
      i_ai = qmmm_get_i_ai()

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

      status = task_gradient_doit(rtdb,theory,qm_energy,dbl_mb(i_g))

      call mm_get_solute_force_quant(nqm,
     &                     link,
     &                     int_mb(i_ai),
     &                     dbl_mb(i_gmm))
 
c
      call daxpy(3*nqm,1.0d0,dbl_mb(i_gmm),1,dbl_mb(i_g),1)

      if(.not.rtdb_put(rtdb,'qmmm:gradient',mt_dbl,3*nqm,dbl_mb(i_g)))
     +   call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy',0,RTDB_ERR)

      mm_energy = mm_energy/cau2kj
      energy = qm_energy + mm_energy-eatoms

      if (.not. rtdb_put(rtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

       if(.not.ma_free_heap(h_gmm))
     & call errquit('qmmm g: Failed to deallocate heap',nqm,
     &       MA_ERR)
      if(.not.ma_free_heap(h_g))
     & call errquit('qmmm g: Failed to deallocate heap',nqm,
     &       MA_ERR)

      call qmmm_print_energy(rtdb)
 
      qmmm_gradient=status
      call qmmm_geom_finish(rtdb)

      end

