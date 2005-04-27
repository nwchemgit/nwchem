      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy)
c
c $Id: qmmm_main.F,v 1.38 2005-04-27 00:54:07 marat Exp $
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
c
c     make sure qmmm is initialized
c     -----------------------------
      if (.not. rtdb_get(inrtdb,'qmmm:init',mt_log,1,qmmm_init))
     $     qmmm_init=.false.
      if(.not.qmmm_init) return
c
c     nothing should be calculated if
c     qm module (e.g. driver is in charge)
c     -----------------------------------
      interface = qmmm_get_interface()
      if(interface.eq.'qm') then
      if(qmmm_print_debug())
     > write(*,*) "exiting out of qmmm_forces since interface is set qm"
        energy = 0
        return
      end if

      if (.not. rtdb_cget(inrtdb, 'qmmm:optimization', 1, optimization))
     $     call errquit(pname//'failed to get optimization flag',0,0)
     
c 
c     get other options from rtdb
c     --------------------------
      tag ="qmmm:fixqm" 
      if (.not.rtdb_get(inrtdb,tag,mt_log,1,fixqm)) 
     >      call errquit(pname//'failed to get fixqm',0,
     >       RTDB_ERR)
      
      if (.not.rtdb_get(inrtdb,'qmmm:ignore_mm',mt_log,1,ignore_mm))
     >   ignore_mm = .false.

      if (.not.rtdb_get(inrtdb,'qmmm:static_qm',mt_log,1,static_qm))
     >   static_qm = .false.

      if (.not.rtdb_get(inrtdb,'qmmm:ignore_qm',mt_log,1,ignore_qm)) 
     >   ignore_qm = .false.

       if (.not.rtdb_get(inrtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)
c
c     link are now treated as classical atoms
c     ----------------------------------------
      link = .false.
c
c     get coordinates for qm region from mm module
c     -------------------------------------------- 
      call qmmm_get_qm_coord_mm(inrtdb)
c
c     propagate coordinates for qm region into qm modules
c     ---------------------------------------------------
      call qmmm_push_qm_coord_qm(inrtdb)
c
c     create full geometry with bq charges for qm module
c     -------------------------------------------------- 
      call qmmm_create_full_geom(inrtdb)
c
c     get qmmm data
c     ------------
      nqm   = qmmm_get_nqm()
      nbq   = qmmm_get_nbq ()
      nbqs  = qmmm_get_nbqs()
      nbqw  = qmmm_get_nbqw()
      nlink = qmmm_get_nlink()
      ntot = nqm+nbq
c
      i_ai  = qmmm_get_i_ai()
      i_ibq = qmmm_get_i_ibq()
      i_cbq   = qmmm_get_i_cbq()
      i_qbq   = qmmm_get_i_qbq()

c
c     allocate gradient arrays
c     ------------------------
      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm grad array',
     &      h_g, i_g) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)

      if(.not.ma_alloc_get(MT_DBL, 3*nqm, 'qmmm grad array',
     &      h_gqm, i_gqm) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*nqm, MA_ERR)
      call dfill(3*nqm,0,dbl_mb(i_gqm),1)

      if(.not.ma_alloc_get(MT_DBL, 3*nbq, 'qmmm grad array',
     &      h_gbq, i_gbq) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*nqm, MA_ERR)
      call dfill(3*nbq,0,dbl_mb(i_gbq),1)
c
c     calculate qm contributions to gradients and energy
c     --------------------------------------------------
      if (.not. rtdb_cget(inrtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
      if(ignore_qm) then
        energy = 0.0d0
        call dfill(3*ntot,0,dbl_mb(i_g),1)
      else
        if(optimization.eq."mm") then 
          if(.not.static_qm)  
     >    status = task_energy_doit(inrtdb,theory,energy)
        else
          status = task_gradient_doit(inrtdb,theory,energy,dbl_mb(i_g))
        end if
      end if
c
c
c     do forces on bq's
c     ----------------
      call qmmm_bq_forces(inrtdb,nbq,
     >                    dbl_mb(i_cbq),
     >                    dbl_mb(i_qbq),
     >                    dbl_mb(i_g+3*nqm))       
c
c
      if(ga_nodeid().eq.0) 
     >  write(*,*) "qm portion of qmmm energy",energy
c
      energy = energy-eatoms
      energy=energy*cau2kj
c
c     get classical forces
c     -------------------
      call mm_get_solute_force_quant(nqm,
     >                               link,
     >                               int_mb(i_ai),
     >                               dbl_mb(i_gqm))

      call mm_get_solute_force_bq(nbqs,
     >                               int_mb(i_ibq),
     >                               dbl_mb(i_gbq))


      call mm_get_solvent_force(nbqw,
     >                          int_mb(i_ibq+nbqs),
     >                          dbl_mb(i_gbq+3*nbqs))


c
c     add quantum contributions
c     or zero out all forces on quantum region
c     ----------------------------------------
      norm = 0.0d0
      if(fixqm) then
        call dfill(3*nqm,0,dbl_mb(i_gqm),1)
      else
        do i=1,3*nqm
           dbl_mb(i_gqm+i-1)=dbl_mb(i_gqm+i-1)+dbl_mb(i_g+i-1)
           norm = norm + dbl_mb(i_gqm+i-1)*dbl_mb(i_gqm+i-1)
        end do
      end if
      norm = sqrt(norm)
      if(qmmm_print_debug())
     > write(*,*) "qnorm",norm
c
      call qmmm_push_qm_grad_qm(inrtdb,3*nqm,dbl_mb(i_gqm))
c
c     add quantum contributions on bq region
c     --------------------------------------
      norm = 0.0d0
      do i=1,3*nbq
         dbl_mb(i_gbq+i-1)=dbl_mb(i_gbq+i-1)+dbl_mb(i_g+3*nqm+i-1)
         norm = norm + dbl_mb(i_gbq+i-1)*dbl_mb(i_gbq+i-1)
      end do
      norm = sqrt(norm/nbq)
      if(qmmm_print_debug())
     > write(*,*) "bqnorm",norm

c
c    add link forces
c    ---------------
      call qmmm_link_forces_qm(nqm,dbl_mb(i_gqm))
      call qmmm_link_forces_mm(nqm,nbq,dbl_mb(i_gqm),
     >                         int_mb(i_ibq),
     >                         dbl_mb(i_gbq))

c
c    propagate final forces into mm module
c    -------------------------------------
      call mm_set_solute_force_quant(nqm,
     >                               .true.,
     >                               int_mb(i_ai),
     >                               dbl_mb(i_gqm))

      call mm_set_solute_force_bq(nbqs,
     >                               int_mb(i_ibq),
     >                               dbl_mb(i_gbq))


      call mm_set_solvent_force(nbqw,
     >                          int_mb(i_ibq+nbqs),
     >                          dbl_mb(i_gbq+3*nbqs))

c
c    restore normal qm geometry (i.e. without bq's)
c    ---------------------------------------------
      call qmmm_restore_qm_geom(inrtdb)
      call qmmm_unset_actlist(inrtdb)

      if(qmmm_print_debug())
     > write(*,*) "out "//pname

      if(.not.ma_free_heap(h_gbq))
     & call errquit(pname//'gbq: Failed to deallocate heap',nqm,
     &       MA_ERR)
      if(.not.ma_free_heap(h_gqm))
     & call errquit(pname//' gqm: Failed to deallocate heap',nqm,
     &       MA_ERR)
      if(.not.ma_free_heap(h_g))
     & call errquit(pname//' g: Failed to deallocate heap',nqm,
     &       MA_ERR)
 
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

