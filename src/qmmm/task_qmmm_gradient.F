c
c $Id: task_qmmm_gradient.F,v 1.20 2005-05-29 18:53:07 marat Exp $
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
#include "qmmm.fh"
c
      integer rtdb
      logical  task_qmmm_gradient
c
      logical  qmmm_gradient
      external qmmm_gradient
c     local variables
      character*32 interface
      double precision energy

      interface = qmmm_get_interface()
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
c         call qmmm_push_grad_qm(rtdb)
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
      logical task_energy_doit
      external task_energy_doit
      logical task_gradient_doit
      external task_gradient_doit
      character*255 theory
      double precision energy
      double precision mm_energy
      double precision qm_energy
      integer i_g,h_g 
      integer i_gt,h_gt
      integer i_ai
      double precision eatoms
      logical link
      integer ntot
      character*30 operation
      logical do_actlist
      logical do_print
      logical save_rst
      integer nfg,ng
      integer i_ifg,i_ig
      integer nact
      integer ngs,ngw
      integer i

      status = .true.

      call util_print_push() 
      call util_print_rtdb_load(rtdb, 'qmmm')
      
      ng  = qmmm_geom_get_ng()
      ngs  = qmmm_geom_get_ngs()
      ngw  = qmmm_geom_get_ngw()
      i_ig = qmmm_geom_get_i_ig()
      nfg = qmmm_geom_get_nfg()
      i_ifg = qmmm_geom_get_i_ifg()

      if(.not.ma_alloc_get(MT_DBL, 3*nfg, 'qmmm grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*nfg,0,dbl_mb(i_g),1)

      if(.not.ma_alloc_get(MT_DBL, 3*ng, 'qmmm grad array',
     &      h_gt, i_gt) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ng,0,dbl_mb(i_gt),1)

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

      do_print = .false.
      if (operation.eq.'gradient') do_print = .true.

       if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      call qmmm_push_coord_mm(rtdb)

c     calculate mm gradients
      call mm_task_sp()
c
c     create full geometry
c     -------------------
      call qmmm_create_full_geom(rtdb)
      call qmmm_bq_coord_update()
c
      call qmmm_map_fixed_atoms(rtdb,nact)
c
c     QM energy calculation
c     --------------------
      if(nact.ne.0) then
        status = task_gradient_doit(rtdb,theory,qm_energy,dbl_mb(i_g))
        call qmmm_links_adjust_forces(nfg,
     &                       int_mb(i_ifg),
     &                       dbl_mb(i_g))


        call mm_add_solute_force(nfg,
     &                     int_mb(i_ifg),
     &                     dbl_mb(i_g))

        call qmmm_bq_gradient(rtdb)
       else
        status = task_energy_doit(rtdb,theory,qm_energy)
      end if


      call mm_get_solute_force(ngs,
     &                     int_mb(i_ig),
     &                     dbl_mb(i_gt))

      call mm_get_solvent_force(ngw,
     &                     int_mb(i_ig+ngs),
     &                     dbl_mb(i_gt+3*ngs))

      
      write(*,*) "final forces"
      do i=1,3*ng
      write(*,*) dbl_mb(i_gt+i-1)
      end do

      if(.not.rtdb_put(rtdb,'qmmm:gradient',mt_dbl,
     +                 3*ng,dbl_mb(i_gt)))
     +   call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

c     take care of the energies
c
      mm_energy = 0.0d0
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

      if(.not.ma_free_heap(h_gt))
     & call errquit('qmmm g: Failed to deallocate heap',nfg,
     &       MA_ERR)

      if(.not.ma_free_heap(h_g))
     & call errquit('qmmm g: Failed to deallocate heap',nfg,
     &       MA_ERR)

      call qmmm_print_energy(rtdb)
 
      qmmm_gradient=status

c
      call qmmm_restore_geom(rtdb)

      call util_print_pop()
      if(save_rst)
     > call mm_write_restart()

      call qmmm_restore_fixed_atoms(rtdb)

      if(do_print)
     >  call mm_print_info()

      end

      subroutine qmmm_push_grad_qm(rtdb)
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
      integer i_ig
      integer i_gt,h_gt
      integer ng

      
      ng   = qmmm_geom_get_ng()
      i_ig = qmmm_geom_get_i_ig()

      if(.not.ma_alloc_get(MT_DBL, 3*ng, 'qmmm grad array',
     &      h_gt, i_gt) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ng, MA_ERR)
      call dfill(3*ng,0,dbl_mb(i_gt),1)

      call mm_get_solute_force(ng,
     &                     int_mb(i_ig),
     &                     dbl_mb(i_gt))

      if(.not.rtdb_put(rtdb,'qmmm:gradient',mt_dbl,
     +                 3*ng,dbl_mb(i_gt)))
     +   call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

      if(.not.ma_free_heap(h_gt))
     & call errquit('qmmm g: Failed to deallocate heap',ng,
     &       MA_ERR)

      end

      subroutine qmmm_bq_gradient(rtdb)
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
      integer i_c,h_c
      integer i_q,h_q
      integer i_i,h_i
      integer i_g,h_g
  
      integer i,k,na
      character*30 pname


      pname = "qmmm_bq_gradient"

      na  = qmmm_get_nbqa() 
     
      if(na.eq.0) return

      if(.not.ma_push_get(MT_DBL, na, 'active bq charges',
     &      h_q, i_q) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      na, MA_ERR)
      call dfill(na,0.0d0,dbl_mb(i_q),1)

      if(.not.ma_push_get(MT_DBL,3*na, 'active bq coord',
     &      h_c, i_c) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*na, MA_ERR)
      call dfill(3*na,0.0d0,dbl_mb(i_c),1)

      if(.not.ma_push_get(MT_DBL,3*na, 'active bq forces',
     &      h_g, i_g) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*na, MA_ERR)
      call dfill(3*na,0.0d0,dbl_mb(i_g),1)


      call qmmm_get_cqbqa(na,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q))


      call qmmm_bq_forces(rtdb,
     >                    na,
     >                    dbl_mb(i_c),
     >                    dbl_mb(i_q),
     >                    dbl_mb(i_g))

      write(*,*) "bq forces"
      do i=1,3*na
      write(*,*) dbl_mb(i_g+i-1)
      end do

c      call mm_add_solute_force(na,
c     &                     int_mb(i_i),
c     &                     dbl_mb(i_g))
 
       call qmmm_bq_push_grad(na,dbl_mb(i_g))

      if(.not.ma_pop_stack(h_g))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      if(.not.ma_pop_stack(h_c))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      if(.not.ma_pop_stack(h_q))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      call qmmm_restore_geom(rtdb)

      
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

      subroutine qmmm_map_fixed_atoms(rtdb,nact1)
      implicit none
c
c this will map active atom list from auxiliary geometry in qm modules
c to a  geometry used in task gradient. 
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
      integer nact1
c
      logical prior_cons
      character*32 pname
      integer nact
      integer i,j
      integer i_cons,h_cons
      integer i_act,h_act
      integer ma_type
      integer ng,nfg,nbqa
      integer i_ig,i_ifg,i_ibqa
      integer jj0,jj,ii
      integer nlink
      integer i_lb
     
      pname = "qmmm_map_fixed_atoms"
c
c     pointers to global index for auxiliary (i_ig)
c     and gradient geometry (i_ifg)
c     --------------------------------------------
      ng   = qmmm_geom_get_ng()
      i_ig = qmmm_geom_get_i_ig()

      nfg   = qmmm_geom_get_nfg()
      i_ifg = qmmm_geom_get_i_ifg()

      nlink  = qmmm_get_nlink()
      i_lb   = qmmm_get_i_lb()
     
c
c     get active list
c     --------------
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact, h_cons)) then

        if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)

        if (.not.rtdb_delete(rtdb, 'geometry:actlist'))
     $       call errquit(pname,0,
     &       RTDB_ERR)

        if (.not.rtdb_put(rtdb, 'qmmm:actlist',
     >       mt_int,nact,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      else
        nact = ng
        if(.not.ma_alloc_get(mt_int, nact, 'qmmm actlist',
     &      h_cons, i_cons) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nact, MA_ERR)
        do i=1,nact
           int_mb(i_cons+i-1) = i
        end do
      end if

c
c     create new active list
c     ----------------------
      if(.not.ma_push_get(MT_INT, nfg, 'qmmm fixed atom list',
     &      h_act, i_act) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nact, MA_ERR)

      nact1=0
      do i=1,nfg
        ii = int_mb(i_ifg+i-1)
        do j=1,ng
         jj0 = int_mb(i_cons+j-1)
         jj  = int_mb(i_ig+jj0-1)
         if(ii.eq.jj) then
            nact1 = nact1 + 1
            int_mb(i_act+nact1-1) = i
            goto 1
          end if
        end do
        do j=1,nlink
         jj  = int_mb(i_lb+nlink+j-1)
         if(ii.eq.jj) then
            nact1 = nact1 + 1
            int_mb(i_act+nact1-1) = i
            goto 1
          end if
        end do
1       continue
      end do


      if(nact1.ne.0) then
        if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >          mt_int,nact1,int_mb(i_act)))
     >         call errquit(pname,0,
     >         RTDB_ERR)
      end if

      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)


      if (.not. ma_pop_stack(h_act) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)

      end

      subroutine qmmm_restore_fixed_atoms(rtdb)
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
      logical ignore
     
      pname = "qmmm_restore_fixed_atoms"
     
      ignore = rtdb_delete(rtdb, 'geometry:actlist')
     
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
 
