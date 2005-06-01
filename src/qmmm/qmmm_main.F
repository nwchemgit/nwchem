      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy)
c
c $Id: qmmm_main.F,v 1.41 2005-06-01 16:27:03 marat Exp $
c
      implicit none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
c
      external task_gradient,esp
      logical task_gradient,esp
c
      integer cf_element
      external cf_element
c
      integer inrtdb
      integer mwm,nwm,mwa,nwa,msa,nsa
      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
      real*8 energy,xw(mwm,3,mwa),fw(mwm,3,mwa),xs(msa,3),fs(msa,3)
c
      character*32 tag
      character*32 interface
      character*32 optimization
      character*255 theory
      logical qmmm_init
      character*32 pname
      integer nbq,nqm,ntot,nlink
      integer i_ai,i_ibq
      integer i_g,h_g 
      integer i_gqm,h_gqm 
      integer i_gbq,h_gbq
      double precision eatoms
      logical status
      integer i
      logical fixqm
      logical ignore_mm
      logical ignore_qm
      logical static_qm
c
      logical task_energy_doit
      external task_energy_doit
      logical task_gradient_doit
      external task_gradient_doit
      double precision norm
      integer nbqs, nbqw
      logical link
c 
      integer i_cbq,i_qbq
c
      integer igeom

      pname = 'qmmm_main'
      if(qmmm_print_debug()) 
     >  write(*,*) "in "//pname
cc
cc     make sure qmmm is initialized
cc     -----------------------------
c      if (.not. rtdb_get(inrtdb,'qmmm:init',mt_log,1,qmmm_init))
c     $     qmmm_init=.false.
c      if(.not.qmmm_init) return
cc
cc     nothing should be calculated if
cc     qm module (e.g. driver is in charge)
cc     -----------------------------------
c      interface = qmmm_get_interface()
c      if(interface.eq.'qm') then
c      if(qmmm_print_debug())
c     > write(*,*) "exiting out of qmmm_forces since interface is set qm"
c        energy = 0
c        return
c      end if
c
c      if (.not. rtdb_cget(inrtdb, 'qmmm:optimization', 1, optimization))
c     $     call errquit(pname//'failed to get optimization flag',0,0)
c     
cc 
cc     get other options from rtdb
cc     --------------------------
c      tag ="qmmm:fixqm" 
c      if (.not.rtdb_get(inrtdb,tag,mt_log,1,fixqm)) 
c     >      call errquit(pname//'failed to get fixqm',0,
c     >       RTDB_ERR)
c      
c      if (.not.rtdb_get(inrtdb,'qmmm:ignore_mm',mt_log,1,ignore_mm))
c     >   ignore_mm = .false.
c
c      if (.not.rtdb_get(inrtdb,'qmmm:static_qm',mt_log,1,static_qm))
c     >   static_qm = .false.
c
c      if (.not.rtdb_get(inrtdb,'qmmm:ignore_qm',mt_log,1,ignore_qm)) 
c     >   ignore_qm = .false.
c
c       if (.not.rtdb_get(inrtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
c     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)
cc
cc     link are now treated as classical atoms
cc     ----------------------------------------
c      link = .false.
cc
cc     get coordinates for qm region from mm module
cc     -------------------------------------------- 
c      call qmmm_get_qm_coord_mm(inrtdb)
cc
cc     propagate coordinates for qm region into qm modules
cc     ---------------------------------------------------
c      call qmmm_push_qm_coord_qm(inrtdb)
cc
cc     create full geometry with bq charges for qm module
cc     -------------------------------------------------- 
c      call qmmm_create_full_geom(inrtdb)
cc
cc     get qmmm data
cc     ------------
c      nqm   = qmmm_get_nqm()
c      nbq   = qmmm_get_nbq ()
c      nbqs  = qmmm_get_nbqs()
c      nbqw  = qmmm_get_nbqw()
c      nlink = qmmm_get_nlink()
c      ntot = nqm+nbq
cc
c      i_ai  = qmmm_get_i_ai()
c      i_ibq = qmmm_get_i_ibq()
c      i_cbq   = qmmm_get_i_cbq()
c      i_qbq   = qmmm_get_i_qbq()
c
cc
cc     allocate gradient arrays
cc     ------------------------
c      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm grad array',
c     &      h_g, i_g) ) call errquit(
c     &      pname//'unable to allocate heap space',
c     &      3*ntot, MA_ERR)
c      call dfill(3*ntot,0,dbl_mb(i_g),1)
c
c      if(.not.ma_alloc_get(MT_DBL, 3*nqm, 'qmmm grad array',
c     &      h_gqm, i_gqm) ) call errquit(
c     &      pname//'unable to allocate heap space',
c     &      3*nqm, MA_ERR)
c      call dfill(3*nqm,0,dbl_mb(i_gqm),1)
c
c      if(.not.ma_alloc_get(MT_DBL, 3*nbq, 'qmmm grad array',
c     &      h_gbq, i_gbq) ) call errquit(
c     &      pname//'unable to allocate heap space',
c     &      3*nqm, MA_ERR)
c      call dfill(3*nbq,0,dbl_mb(i_gbq),1)
cc
cc     calculate qm contributions to gradients and energy
cc     --------------------------------------------------
c      if (.not. rtdb_cget(inrtdb, 'task:theory', 1, theory))
c     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
c      if(ignore_qm) then
c        energy = 0.0d0
c        call dfill(3*ntot,0,dbl_mb(i_g),1)
c      else
c        if(optimization.eq."mm") then 
c          if(.not.static_qm)  
c     >    status = task_energy_doit(inrtdb,theory,energy)
c        else
c          status = task_gradient_doit(inrtdb,theory,energy,dbl_mb(i_g))
c        end if
c      end if
cc
cc
cc     do forces on bq's
cc     ----------------
c      if(.not.ignore_qm) then
c      call qmmm_bq_forces(inrtdb,nbq,
c     >                    dbl_mb(i_cbq),
c     >                    dbl_mb(i_qbq),
c     >                    dbl_mb(i_g+3*nqm))       
c      end if
cc
cc
c      if(ga_nodeid().eq.0) 
c     >  write(*,*) "qm portion of qmmm energy",energy
cc
c      energy = energy-eatoms
c      energy=energy*cau2kj
cc
cc     get classical forces
cc     -------------------
c      call mm_get_solute_force_quant(nqm,
c     >                               link,
c     >                               int_mb(i_ai),
c     >                               dbl_mb(i_gqm))
c
c      call mm_get_solute_force_bq(nbqs,
c     >                               int_mb(i_ibq),
c     >                               dbl_mb(i_gbq))
c
c
c      call mm_get_solvent_force(nbqw,
c     >                          int_mb(i_ibq+nbqs),
c     >                          dbl_mb(i_gbq+3*nbqs))
c
c
cc
cc     add quantum contributions
cc     or zero out all forces on quantum region
cc     ----------------------------------------
c      norm = 0.0d0
c      if(fixqm) then
c        call dfill(3*nqm,0,dbl_mb(i_gqm),1)
c      else
c        do i=1,3*nqm
c           dbl_mb(i_gqm+i-1)=dbl_mb(i_gqm+i-1)+dbl_mb(i_g+i-1)
c           norm = norm + dbl_mb(i_gqm+i-1)*dbl_mb(i_gqm+i-1)
c        end do
c      end if
c      norm = sqrt(norm)
c      if(qmmm_print_debug())
c     > write(*,*) "qnorm",norm
cc
c      call qmmm_push_qm_grad_qm(inrtdb,3*nqm,dbl_mb(i_gqm))
cc
cc     add quantum contributions on bq region
cc     --------------------------------------
c      norm = 0.0d0
c      do i=1,3*nbq
c         dbl_mb(i_gbq+i-1)=dbl_mb(i_gbq+i-1)+dbl_mb(i_g+3*nqm+i-1)
c         norm = norm + dbl_mb(i_gbq+i-1)*dbl_mb(i_gbq+i-1)
c      end do
c      norm = sqrt(norm/nbq)
c      if(qmmm_print_debug())
c     > write(*,*) "bqnorm",norm
c
cc
cc    add link forces
cc    ---------------
c      call qmmm_link_forces_qm(nqm,dbl_mb(i_gqm))
c      call qmmm_link_forces_mm(nqm,nbq,dbl_mb(i_gqm),
c     >                         int_mb(i_ibq),
c     >                         dbl_mb(i_gbq))
c
cc
cc    propagate final forces into mm module
cc    -------------------------------------
c      call mm_set_solute_force_quant(nqm,
c     >                               .true.,
c     >                               int_mb(i_ai),
c     >                               dbl_mb(i_gqm))
c
c      call mm_set_solute_force_bq(nbqs,
c     >                               int_mb(i_ibq),
c     >                               dbl_mb(i_gbq))
c
c
c      call mm_set_solvent_force(nbqw,
c     >                          int_mb(i_ibq+nbqs),
c     >                          dbl_mb(i_gbq+3*nbqs))
c
cc
cc    restore normal qm geometry (i.e. without bq's)
cc    ---------------------------------------------
c      call qmmm_restore_geom(inrtdb)
c      call qmmm_unset_actlist(inrtdb)
c
c      if(qmmm_print_debug())
c     > write(*,*) "out "//pname
c
c      if(.not.ma_free_heap(h_gbq))
c     & call errquit(pname//'gbq: Failed to deallocate heap',nqm,
c     &       MA_ERR)
c      if(.not.ma_free_heap(h_gqm))
c     & call errquit(pname//' gqm: Failed to deallocate heap',nqm,
c     &       MA_ERR)
c      if(.not.ma_free_heap(h_g))
c     & call errquit(pname//' g: Failed to deallocate heap',nqm,
c     &       MA_ERR)
c 
      return
      end

      subroutine qmmm_set_actlist(n1,n2,rtdb)
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
      integer n1,n2
      integer rtdb
      logical ignore
c
      character*32 pname
      integer nact
      integer i
      integer j
      integer i_cons,h_cons
     
      pname = "qmmm_set_actlist"
      nact = n2-n1+1

      if(.not.ma_push_get(MT_INT, nact, 'qmmm hcons',
     &      h_cons, i_cons) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nact, MA_ERR)

      j = n1
      do i=1,nact
       int_mb(i_cons+i-1) = j
       j = j + 1
      end do

      ignore = rtdb_delete(rtdb, 'geometry:actlist')


      if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >        mt_int,nact,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not. ma_pop_stack(h_cons) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)


      end

      subroutine qmmm_unset_actlist(rtdb)
      implicit none
c
#include "rtdb.fh"
c
      integer rtdb
      logical ignore

      ignore = rtdb_delete(rtdb, 'geometry:actlist')

      end

