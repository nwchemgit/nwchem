      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy)
c
c $Id: qmmm_main.F,v 1.21 2004-10-19 20:01:49 marat Exp $
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
c      character*32 optimization
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
c
      logical task_gradient_doit
      external task_gradient_doit
      double precision norm

      pname = 'qmmm_forces'
      if(qmmm_print_debug()) 
     >  write(*,*) "in "//pname
      if (.not. rtdb_get(inrtdb,'qmmm:init',mt_log,1,qmmm_init))
     $     qmmm_init=.false.

      if(.not.qmmm_init) return

      tag ="qmmm:fixqm" 
      if (.not.rtdb_get(inrtdb,tag,mt_log,1,fixqm)) 
     >      call errquit(pname//'failed to get fixqm',0,
     >       RTDB_ERR)


      tag ="qmmm:interface"
      if (.not.rtdb_cget(inrtdb,tag,1,interface)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)


c      if (.not. rtdb_cget(inrtdb,'qmmm:optimization',1,optimization))
c     $     call errquit('qmmm_forces: failed to get qmmm:optimization',
c     &       0,RTDB_ERR)
c
c      fixqm = .false.
c      interface = "qm"
c      if (optimization.eq.'all') then
c         fixqm = .false.
c         interface = "mm"
c      else if (optimization.eq.'mm') then
c         fixqm = .true.
c         interface = "mm"
c      else if (optimization.eq.'qm') then
c         fixqm = .false.
c         interface = "qm"
c      endif
c
      if(interface.eq.'qm') then
      if(qmmm_print_debug())
     > write(*,*) "exiting out of qmmm_forces since interface is set qm"
        energy = 0
        return
      end if
      
      call qmmm_get_qm_coord_mm(inrtdb)
      call qmmm_push_qm_coord_qm(inrtdb)
 
      call qmmm_create_bq_geom(inrtdb)
 
      nqm  = qmmm_get_nqm()
      nbq  = qmmm_get_nbqs()
      ntot = nqm+nbq
      i_ai = qmmm_get_i_ai()
      i_ibq = qmmm_get_i_ibq()

      if(.not.ma_alloc_get(MT_DBL, 3*ntot, 'qmmm grad array',
     &      h_g, i_g) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*ntot, MA_ERR)
      call dfill(3*ntot,0,dbl_mb(i_g),1)

      if(.not.ma_alloc_get(MT_DBL, 3*nqm, 'qmmm grad array',
     &      h_gqm, i_gqm) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*nqm, MA_ERR)
      call dfill(3*nqm,0,dbl_mb(i_gqm),1)

      if(.not.ma_alloc_get(MT_DBL, 3*nbq, 'qmmm grad array',
     &      h_gbq, i_gbq) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      3*nqm, MA_ERR)
      call dfill(3*nbq,0,dbl_mb(i_gbq),1)

      if (.not. rtdb_cget(inrtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
       if (.not.rtdb_get(inrtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      status = task_gradient_doit(inrtdb,theory,energy,dbl_mb(i_g))
      energy = energy-eatoms
      energy=energy*cau2kj

      call mm_get_solute_force_quant(nqm,
     >                               .true.,
     >                               int_mb(i_ai),
     >                               dbl_mb(i_gqm))

      call mm_get_solute_force_bq(nbq,
     >                               int_mb(i_ibq),
     >                               dbl_mb(i_gbq))


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
      if(.not.rtdb_put(inrtdb,'qmmm:gradient',mt_dbl,
     >                 3*nqm,dbl_mb(i_gqm)))
     +   call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)


      norm = 0.0
      do i=1,3*nbq
         dbl_mb(i_gbq+i-1)=dbl_mb(i_gbq+i-1)+dbl_mb(i_g+3*nqm+i-1)
         norm = norm + dbl_mb(i_gbq+i-1)*dbl_mb(i_gbq+i-1)
      end do
      norm = sqrt(norm/nbq)
      if(qmmm_print_debug())
     > write(*,*) "bqnorm",norm


      call mm_set_solute_force_quant(nqm,
     >                               .true.,
     >                               int_mb(i_ai),
     >                               dbl_mb(i_gqm))

      call mm_set_solute_force_bq(nbq,
     >                               int_mb(i_ibq),
     >                               dbl_mb(i_gbq))


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

