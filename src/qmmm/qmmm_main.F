      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy)
c
c $Id: qmmm_main.F,v 1.26 2005-01-06 23:55:46 edo Exp $
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
      character*255 theory
      logical qmmm_init
      character*32 pname
      integer nbq,nqm,ntot
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
c
      logical task_gradient_doit
      external task_gradient_doit
      double precision norm
      integer nbqs, nbqw
      logical link

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
      tag ="qmmm:interface"
      if (.not.rtdb_cget(inrtdb,tag,1,interface)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)
      if(interface.eq.'qm') then
      if(qmmm_print_debug())
     > write(*,*) "exiting out of qmmm_forces since interface is set qm"
        energy = 0
        return
      end if
     
c 
c     get other options from rtdb
c     --------------------------
      tag ="qmmm:fixqm" 
      if (.not.rtdb_get(inrtdb,tag,mt_log,1,fixqm)) 
     >      call errquit(pname//'failed to get fixqm',0,
     >       RTDB_ERR)
      
      if (.not.rtdb_get(inrtdb,'qmmm:ignore_mm',mt_log,1,ignore_mm))
     >   ignore_mm = .false.

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
      ntot = nqm+nbq
c
      i_ai  = qmmm_get_i_ai()
      i_ibq = qmmm_get_i_ibq()

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
        status = task_gradient_doit(inrtdb,theory,energy,dbl_mb(i_g))
      end if
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
      if(.not.rtdb_put(inrtdb,'qmmm:gradient',mt_dbl,
     >                 3*nqm,dbl_mb(i_gqm)))
     >   call errquit(pname//'rtdb_put gradient failed',0, RTDB_ERR)

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

