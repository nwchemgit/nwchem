c
c $Id: task_qmmm_gradient.F,v 1.2 2004-11-24 22:57:51 marat Exp $
c
c
      function task_qmmm_gradient(rtdb)
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
      logical  task_qmmm_gradient
c
      logical  qmmm_gradient
      external qmmm_gradient
c     local variables
      character*32 interface
      double precision energy

      if (.not. rtdb_cget(rtdb, 'qmmm:interface', 1, interface))
     $     call errquit('qmmm_interface: failed ',0,0)
c
      
      task_qmmm_gradient = .true.
      if (interface.eq.'mm') then
         call qmmm_get_qm_coord_qm(rtdb)
         call qmmm_push_qm_coord_mm(rtdb)
         call md_sp()
         if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)

         if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy/cau2kj))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      else if (interface.eq.'qm') then
         task_qmmm_gradient = qmmm_gradient(rtdb)
      endif
c
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

      call util_print_push() 
      call util_print_rtdb_load(rtdb, 'qmmm')

      link = .not.qmmm_h_link()


      status = .true.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

c     get qm coordinates from geometry object
      call qmmm_get_qm_coord_qm(rtdb)
c     push qm coordinates to mm module
      call qmmm_push_qm_coord_mm(rtdb)
c     calculate mm gradients
      call mm_task_sp()

c     h_qmm will hold mm contribution to forces on quantum atoms
      nqm  = qmmm_get_nqm()
      if(.not.ma_alloc_get(MT_DBL, 3*nqm, 'qmmm mm grad array',
     &      h_gmm, i_gmm) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(3*nqm,0,dbl_mb(i_gmm),1)

      i_ai = qmmm_get_i_ai()
      call mm_get_solute_force_quant(nqm,
     &                     link,
     &                     int_mb(i_ai),
     &                     dbl_mb(i_gmm))
 

c     create full geometry object with Bq's
      call qmmm_create_bq_geom(rtdb)

      nbq  = qmmm_get_nbq()
      ntot = nqm+nbq

c     h_g will hold QM gradient for the full system
      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)

c     calculate gradients from qm
      status = task_gradient_doit(rtdb,theory,qm_energy,dbl_mb(i_g))

c     add on classical contributions
      call daxpy(3*nqm,1.0d0,dbl_mb(i_gmm),1,dbl_mb(i_g),1)

c     push final gradients onto qm module
      call qmmm_push_qm_grad_qm(rtdb,3*nqm,dbl_mb(i_g))

c     take care of the energies
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

      call util_print_pop()
      call mm_write_restart()

      end

      subroutine qmmm_push_qm_grad_qm(irtdb,n,g)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
      integer irtdb
      integer n
      double precision g(n)
c     local variables
      integer nqm
      integer ncent
      integer icent
      integer i
      integer k

      integer igeom
      integer i_c
      integer i_act
      integer i_gtmp,h_gtmp
      character*32 pname

      pname = "qmmm_push_qm_coord_qm"

      nqm   = qmmm_get_nqm()
      i_c   = qmmm_get_i_c()
      i_act = qmmm_get_i_act()

      ncent = 0
      do i=1,nqm
        if(log_mb(i_act+i-1)) ncent = ncent + 1
      end do

      if(.not.ma_push_get(mt_dbl,3*ncent,'gtmp',h_gtmp,i_gtmp))
     + call errquit( pname//'Failed to allocate memory for gtmp',
     + 3*ncent, MA_ERR)

      icent = 0
      do i=1,n/3
        if(log_mb(i_act+i-1)) then
          icent = icent + 1
          do k=1,3
            dbl_mb(i_gtmp+3*(icent-1)+k-1) = g(3*(i-1)+k)
            write(*,*) "qmmm_push_qm_grad_qm", g(3*(i-1)+k)
          end do
        end if
      end do

      if(.not.rtdb_put(irtdb,'qmmm:gradient',mt_dbl,
     +                 3*ncent,dbl_mb(i_gtmp)))
     +   call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

      if(.not.ma_pop_stack(h_gtmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack gtmp',ncent,
     &       MA_ERR)

      end

