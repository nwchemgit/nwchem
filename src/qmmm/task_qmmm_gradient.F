c
c $Id: task_qmmm_gradient.F,v 1.16 2005-02-15 06:39:44 marat Exp $
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
      logical ignore_mm
      logical ignore_qm
      character*30 operation
      logical do_actlist
      logical save_rst

      call util_print_push() 
      call util_print_rtdb_load(rtdb, 'qmmm')

      link = .false.

      
      if (.not. rtdb_get(rtdb,'qmmm:ignore_mm',mt_log,1,ignore_mm)) then
         ignore_mm = .false.
      end if

      if (.not. rtdb_get(rtdb,'qmmm:ignore_qm',mt_log,1,ignore_qm)) then
         ignore_qm = .false.
      end if


      status = .true.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = 'energy'

      do_actlist = .false.
      if (operation.eq.'hessian' .or.
     >    operation.eq.'frequencies' .or.
     >    operation.eq.'freq') do_actlist=.true.

      save_rst = .false.
      if (operation.eq.'optimize' .or.
     >    operation.eq.'saddle') save_rst = .true.


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
      call qmmm_create_full_geom(rtdb)
      call qmmm_print_current_geom(rtdb)

      nbq  = qmmm_get_nbq()
      ntot = nqm+nbq

c     h_g will hold QM gradient for the full system
      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)
c
      if(do_actlist) 
     >  call qmmm_load_fixed_atoms(rtdb)


c     calculate gradients from qm
      status = task_gradient_doit(rtdb,theory,qm_energy,dbl_mb(i_g))

      if(ignore_qm) call dfill(3*ntot,0,dbl_mb(i_g),1)
      if(ignore_mm) call dfill(3*nqm,0,dbl_mb(i_gmm),1)

c     add on classical contributions
      call daxpy(3*nqm,1.0d0,dbl_mb(i_gmm),1,dbl_mb(i_g),1)
      
      call qmmm_link_forces_qm(ntot,dbl_mb(i_g))

c     push final gradients onto qm module
      call qmmm_push_qm_grad_qm(rtdb,3*nqm,dbl_mb(i_g))

c     take care of the energies
      mm_energy = 0.0d0
      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy',0,RTDB_ERR)

      mm_energy = mm_energy/cau2kj
      if(ignore_qm) qm_energy = 0.0d0
      if(ignore_mm) mm_energy = 0.0d0
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

c
      call qmmm_restore_qm_geom(rtdb)

      call util_print_pop()
      if(save_rst)
     > call mm_write_restart()

      if(do_actlist) 
     >  call qmmm_unload_fixed_atoms(rtdb)

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

      subroutine qmmm_load_fixed_atoms1(rtdb)
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
c
      logical prior_cons
      character*32 namespace
      character*32 pname
      integer nqm,nact,nfix
      integer i,j
      integer i_cons,h_cons
      integer i_act
      integer icent
      integer i_lfix,h_lfix
      integer ma_type
     
      pname = "qmmm_create_fixed_atoms"
     
      call cons_inquire_fixed_atoms
     >   (
     >   prior_cons,
     >   rtdb
     >   )

      if(.not.prior_cons) return

      write(*,*) "found prior constraints"
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact, h_cons)) then
            if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)
      end if

      call cons_sort_i(nact,int_mb(i_cons))

      namespace = "qmmm_old_cons"
      call cons_save_fixed_atoms(namespace,rtdb)
      call cons_unload_fixed_atoms(rtdb)

      nqm  = qmmm_get_nqm()
      i_act = qmmm_get_i_act()


      if(.not.ma_push_get(MT_LOG, nqm, 'qmmm lfix',
     &      h_lfix, i_lfix) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nqm, MA_ERR)


      do i=1,nqm
        log_mb(i_lfix+i-1)=.true.
      end do

      do i=1,nact
        j = int_mb(i_cons+i-1)
        log_mb(i_lfix+j-1)=.false.
        write(*,*) "j=",j
      end do

      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)

      if(.not.ma_push_get(MT_INT, nqm-nact, 'qmmm cons',
     &      h_cons, i_cons) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nqm, MA_ERR)



      do i=1,nact
        int_mb(i_cons+i-1) = -1
      end do


      icent = 0
      nfix = 0
      do i=1,nqm
        if(log_mb(i_act+i-1)) then
          icent = icent + 1
          if(log_mb(i_lfix+icent-1)) then
            nfix = nfix+1
            int_mb(i_cons+nfix-1) = i
            write(*,*) "i_cons,nfix",int_mb(i_cons+nfix-1),nfix
          end if
        end if
      end do


      namespace = "qmmm_new_cons"
      write(*,*) "nfix=",nfix
      write(*,*) (int_mb(i_cons+i-1),i=1,nfix)

      call cons_add_fixed_atoms(namespace,nfix,int_mb(i_cons),rtdb)
      call cons_load_fixed_atoms(namespace,rtdb)

      if (.not. ma_pop_stack(h_cons) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)


      if (.not. ma_pop_stack(h_lfix) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)

      end

      subroutine qmmm_load_fixed_atoms(rtdb)
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
c
      logical prior_cons
      character*32 pname
      integer nqm,nact,nfix
      integer nact1
      integer i,j
      integer i_cons,h_cons
      integer i_act
      integer icent
      integer i_lfix,h_lfix
      integer ma_type
     
      pname = "qmmm_create_fixed_atoms"
     
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact, h_cons)) then
            if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)
      else
       return
      end if

      if (.not.rtdb_delete(rtdb, 'geometry:actlist'))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not.rtdb_put(rtdb, 'qmmm:actlist',
     >       mt_int,nact,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      call cons_sort_i(nact,int_mb(i_cons))

      nqm  = qmmm_get_nqm()
      i_act = qmmm_get_i_act()


      if(.not.ma_push_get(MT_LOG, nqm, 'qmmm lfix',
     &      h_lfix, i_lfix) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nqm, MA_ERR)


      do i=1,nqm
        log_mb(i_lfix+i-1)=.true.
      end do

      do i=1,nact
        j = int_mb(i_cons+i-1)
        log_mb(i_lfix+j-1)=.false.
      end do


      icent = 0
      nact1 = 0
      do i=1,nqm
        if(log_mb(i_act+i-1)) then
          icent = icent + 1
          if(.not.log_mb(i_lfix+icent-1)) then
            nact1=nact1+1
            int_mb(i_cons+nact1-1) = i
          end if
        end if
      end do

      if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >        mt_int,nact1,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)


      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)


      if (.not. ma_pop_stack(h_lfix) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)

      end

      subroutine qmmm_unload_fixed_atoms(rtdb)
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
c
      character*32 pname
      integer nact
      integer i_cons,h_cons
      integer ma_type
     
      pname = "qmmm_unload_fixed_atoms"
     
      if (rtdb_ma_get(rtdb, 'qmmm:actlist', ma_type,
     $        nact, h_cons)) then
            if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)
      else
        return
      end if

      if (.not.rtdb_delete(rtdb, 'qmmm:actlist'))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >       mt_int,nact,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)

      end

      subroutine qmmm_create_fixed_atoms(rtdb)
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
c
      logical prior_cons
      character*32 namespace
      character*32 pname
      integer nqm,nbq
      integer i
      integer i_cons,h_cons
     
      pname = "qmmm_create_fixed_atoms"
     
      call cons_inquire_fixed_atoms
     >   (
     >   prior_cons,
     >   rtdb
     >   )

      if(prior_cons) then
        write(*,*) "found prior constraints"
        namespace = "qmmm_old_cons"
        call cons_save_fixed_atoms(namespace,rtdb)
        namespace = "qmmm_new_cons"
        call cons_save_fixed_atoms(namespace,rtdb)
        call cons_unload_fixed_atoms(rtdb)
      end if

      nqm  = qmmm_get_nqm()
      nbq  = qmmm_get_nbq()

      if(.not.ma_push_get(MT_INT, nbq, 'qmmm cons list',
     &      h_cons, i_cons) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nbq, MA_ERR)


      do i=1,nbq
        int_mb(i_cons+i-1)=nqm+i
      end do

      namespace = "qmmm_new_cons"
      call cons_add_fixed_atoms(namespace,nbq,int_mb(i_cons),rtdb)
      call cons_load_fixed_atoms(namespace,rtdb)

      if (.not. ma_pop_stack(h_cons) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)



      end

      subroutine qmmm_link_forces_qm(nt,g)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
      integer nt
      double precision g(3,nt)
c     local variables
      integer nqm
      integer nlink
      integer i,k
      integer i_lb
      integer i_c
      integer h_cl,i_cl
      integer i_ai
      integer il,iq,ih
      integer i_mqm
      integer h_mqm
      double precision l_qh,l_ql
      double precision rh(3)
      double precision rq(3)
      double precision rl(3)
      double precision rql(3)
      double precision fh(3)
      double precision fq(3)
      double precision fh_rql
      double precision a,b
      character*32 pname

      double precision qmmm_geom_distance
      external qmmm_geom_distance

      double precision ddot
      external ddot


      pname = "qmmm_link_forces"

      if(qmmm_print_debug())
     >   write(*,*) "in "//pname

      nlink = qmmm_get_nlink()

c
c     nothing to do if no links
c     -------------------------
      if(nlink.eq.0) return

      nqm = qmmm_get_nqm()

c
c     nothing to do if links are not hydrogens
c     ----------------------------------------
      if(.not.qmmm_h_link()) return


      if(.not.ma_push_get(mt_dbl,3*nlink,'cl',h_cl,i_cl))
     + call errquit( pname//'Failed to allocate memory for cl',
     + 3*nlink, MA_ERR)

      if(.not.ma_push_get(mt_int,nlink,'mqm',h_mqm,i_mqm))
     + call errquit( pname//'Failed to allocate memory for mqm',
     + nlink, MA_ERR)

c
c     get pointer to index of link pairs
c     ----------------------------------
      i_lb = qmmm_get_i_lb()
      i_ai = qmmm_get_i_ai()
c
c     map to a local index

      if(.not.qmmm_map(nqm,int_mb(i_ai),
     >                 nlink,int_mb(i_lb),
     >                 int_mb(i_mqm)))  
     >     call errquit( pname//'Failed to map',
     > 0, 0)


c     get pointer to coord array
c     --------------------------
      i_c   = qmmm_get_i_c()
c
c     get coordinates of mm link atoms
c     -------------------------------
      call mm_get_solute_coord_quant(nlink,
     &                   .true.,
     &                   int_mb(i_lb+nlink),
     &                   dbl_mb(i_cl))

c
c     loop over link atom pairs
c     -------------------------
      do i=1,nlink
        ih = nqm-nlink+i
        iq = int_mb(i_mqm+i-1)
        il = i
        do k=1,3
           fh(k)  = g(k,ih)
           rh(k)  = dbl_mb(i_c+(ih-1)*3+k-1)
           rq(k)  = dbl_mb(i_c+(iq-1)*3+k-1)
           rl(k)  = dbl_mb(i_cl+(il-1)*3+k-1)
           rql(k) = rq(k)-rl(k)
        end do
        l_qh=qmmm_geom_distance(rh,rq)
        l_ql=qmmm_geom_distance(rl,rq)
        fh_rql = ddot(3,fh,1,rql,1)
c        a = 1.0d0 - l_qh/l_ql
c        b = fh_rql*l_qh/l_ql**3
        a = 0.709d0
        do k=1,3
c          fq(k) = fh(k)*a + rql(k)*b
          fq(k) = fh(k)*(1.0d0-a)
          g(k,iq) = g(k,iq) + fq(k)
        end do
      end do
    
      if(.not.ma_pop_stack(h_mqm))
     + call errquit( pname//'Failed to deallocate memory for mqm',
     + nlink, MA_ERR)

      if(.not.ma_pop_stack(h_cl))
     + call errquit( pname//'Failed to deallocate memory for cl',
     + 3*nlink, MA_ERR)

      end
c
      subroutine qmmm_link_forces_mm(nt,nb,g,ibq,gbq)
c     nt      - number of quantum atoms
c     nb      - number of bq atoms
c     g(nt)    - forces on quantum atoms
c     ibq(nb) - global index for bq atoms
c     gbq(nb) - forces on quantum atoms 
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
      integer nt
      integer nb
      double precision g(3,nt)
      double precision gbq(3,nb)
      integer ibq(nb)
c     local variables
      integer nqm
      integer nlink
      integer i,k
      integer i_lb
      integer i_c
      integer h_cl,i_cl
      integer i_ai
      integer ib
      integer il,iq,ih
      integer i_mqm,h_mqm
      integer i_mbq,h_mbq
      integer k1
      double precision l_qh,l_ql
      double precision rhq(3)
      double precision rh(3)
      double precision rq(3)
      double precision rl(3)
      double precision rlq(3)
      double precision fh(3)
      double precision fl(3)
      double precision fh_rql
      double precision a,b
      character*32 pname

      double precision qmmm_geom_distance
      external qmmm_geom_distance

      double precision ddot
      external ddot


      pname = "qmmm_link_forces_mm"

      if(qmmm_print_debug())
     >   write(*,*) "in "//pname
c
c     number of link atoms
c     -------------------
      nlink = qmmm_get_nlink()
c
c     nothing to do if no links
c     -------------------------
      if(nlink.eq.0) return
c
c     number of quantum atoms
c     -----------------------
      nqm = qmmm_get_nqm()
c
c     check for consistency
c     ---------------------
      if(nqm.ne.nt)
     > call errquit( pname//'inconsistent number of atoms',
     > 3*nlink, MA_ERR)
c
c     nothing to do if links are not hydrogens
c     ----------------------------------------
      if(.not.qmmm_h_link()) return
c
c     coordinates for link atoms
c     --------------------------
      if(.not.ma_push_get(mt_dbl,3*nlink,'cl',h_cl,i_cl))
     + call errquit( pname//'Failed to allocate memory for cl',
     + 3*nlink, MA_ERR)
c
c     mapping from global to local index for quantum atoms
c     ----------------------------------------------------
      if(.not.ma_push_get(mt_int,nlink,'mqm',h_mqm,i_mqm))
     + call errquit( pname//'Failed to allocate memory for mqm',
     + nlink, MA_ERR)
c
c     mapping from global to local index for bq atoms
c     ----------------------------------------------------
      if(.not.ma_push_get(mt_int,nlink,'mbq',h_mbq,i_mbq))
     + call errquit( pname//'Failed to allocate memory for mbq',
     + nlink, MA_ERR)
c
c     get pointer to index of link pairs
c     ----------------------------------
      i_lb = qmmm_get_i_lb()
c
c     get pointer to global atom index for quantum atoms
c     --------------------------------------------------
      i_ai = qmmm_get_i_ai()
c
c     map to a local index for quantum atoms
c     --------------------------------------
      if(.not.qmmm_map(nqm,int_mb(i_ai),
     >                 nlink,int_mb(i_lb),
     >                 int_mb(i_mqm)))  
     >     call errquit( pname//'Failed to map quantum atoms',0,0)
c
c     map to a local index for bq atoms
c     ---------------------------------
      if(.not.qmmm_map(nb,ibq,
     >                 nlink,int_mb(i_lb+nlink),
     >                 int_mb(i_mbq)))  
     >     call errquit( pname//'Failed to map bq atoms',0,0)


c     get pointer to quantum coords
c     -----------------------------
      i_c   = qmmm_get_i_c()
c
c     get coordinates of mm link atoms
c     -------------------------------
      call mm_get_solute_coord_quant(nlink,
     &                   .true.,
     &                   int_mb(i_lb+nlink),
     &                   dbl_mb(i_cl))

c
c     actually calculate link forces
c     ------------------------------
      do i=1,nlink
        ih = nqm-nlink+i
        iq = int_mb(i_mqm+i-1)
        ib = int_mb(i_mbq+i-1)
        il = i
        do k=1,3
           fh(k)  = g(k,ih)
           rh(k)  = dbl_mb(i_c+(ih-1)*3+k-1)
           rq(k)  = dbl_mb(i_c+(iq-1)*3+k-1)
           rl(k)  = dbl_mb(i_cl+(il-1)*3+k-1)
           rlq(k) = rl(k)-rq(k)
           rhq(k) = rq(k)-rh(k)
        end do
        l_qh=qmmm_geom_distance(rh,rq)
        l_ql=qmmm_geom_distance(rl,rq)
        fh_rql = ddot(3,fh,1,rlq,1)
         a = 0.709d0
c        a = l_qh/l_ql
c        b = fh_rql*l_qh/l_ql**3
        do k=1,3
c          fl(k) = fh(k)*a - rlq(k)*b
          fl(k) = fh(k)*a 
          gbq(k,ib) = gbq(k,ib) + fl(k)
        end do
c        write(*,*) "debug coord link",(rl(k1),k1=1,3)
c        write(*,*) "debug coord hydr",(rh(k1),k1=1,3)
c        write(*,*) "debug force hydr",(fh(k1),k1=1,3)
c        write(*,*) "forces on link bond",iq,ih
c        write(*,*) (fl(k),k=1,3)
c        write(*,*) "debug l_qh",l_qh
c        write(*,*) "debug angle",ddot(3,rhq,1,rlq,1)/(l_ql*l_qh)
        
      end do

      if(.not.ma_pop_stack(h_mbq))
     + call errquit( pname//'Failed to deallocate memory for il',
     + nlink, MA_ERR)

       if(.not.ma_pop_stack(h_mqm))
     + call errquit( pname//'Failed to deallocate memory for il',
     + nlink, MA_ERR)
    
      if(.not.ma_pop_stack(h_cl))
     + call errquit( pname//'Failed to deallocate memory for cl',
     + 3*nlink, MA_ERR)

      end
 
