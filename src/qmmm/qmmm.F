c
c $Id: qmmm.F,v 1.1 2004-02-10 01:07:32 marat Exp $
c
      logical function qmmm(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer rtdb
      character*30 operation
      logical status
      logical task_energy
      external task_energy
c
c
      write(6,*) "initializing qmmm"
      call qmmm_init(rtdb)
c
c     determine calculation type
c     --------------------------
      if(.not.rtdb_cget(rtdb,'task:operation',1,operation))
     + call errquit('qmmm: no calculation type',0,INPUT_ERR)
c
      if (inp_compare(.false., 'energy', operation)) then
         status = task_energy(rtdb)
      end if
      call qmmm_end(rtdb)
c
      qmmm=status
      return
      end
c
      subroutine qmmm_init(rtdb)
      implicit none
c
#include "inp.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
c
#include "qmmm_data.fh"
      integer rtdb
c
      character*(nw_max_path_len) filename
      integer i
      integer i_m,l_m
      integer i_c,l_c
      integer i_q,l_q
      integer i_t,l_t
      integer i_i,l_i
      integer igeom
c
c     initialize mm part
c     (note it will also want operation, see below)
c     --------------------------------------------
      call mm_init(rtdb)

c     initialize qm geometry
c     ---------------------
      call mm_get_tot_nqm(tot_nqm)
c
c     set up indexing array
c     ---------------------
      if(.not.ma_alloc_get(MT_INT, tot_nqm, 'qmmm index array',
     &      h_ai, i_ai) ) call errquit(
     &      'qmmm_init: unable to allocate heap space',
     &      tot_nqm, MA_ERR)
      call ifill(tot_nqm,0,int_mb(i_ai),1)
c
      call mm_qm_mapping(tot_nqm,int_mb(i_ai))
c
c     set up qm geometry
c     ------------------
      if(.not.ma_push_get(mt_dbl,3*tot_nqm,'m',l_m,i_m))
     + call errquit('qmmm: Failed to allocate for c',tot_nqm, MA_ERR)
      if(.not.ma_push_get(mt_dbl,3*tot_nqm,'c',l_c,i_c))
     + call errquit('qmmm: Failed to allocate for c',tot_nqm, MA_ERR)
      if(.not.ma_push_get(mt_dbl,tot_nqm,'q',l_q,i_q))
     + call errquit('qmmm: Failed to allocate for q',tot_nqm, MA_ERR)
      if(.not.ma_push_get(mt_byte,16*tot_nqm,'t',l_t,i_t))
     + call errquit('qmmm: Failed to allocate for t',tot_nqm, MA_ERR)
      if(.not.ma_push_get(mt_int,tot_nqm,'i',l_i,i_i))
     + call errquit('qmmm: Failed to allocate for i',tot_nqm, MA_ERR)
  
      call mm_get_qm_geom(tot_nqm,
     &                    int_mb(i_ai),
     &                    int_mb(i_i),
     &                    dbl_mb(i_c),
     &                    dbl_mb(i_q),
     &                    dbl_mb(i_m),
     &                    byte_mb(i_t))
 
       if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
      if(.not.geom_cart_set(igeom,tot_nqm,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)
      if(.not.geom_masses_set(igeom,tot_nqm,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

c      if (.not.rtdb_get(inrtdb,'qmmm:xyz',mt_log,1,xyz)) xyz=.false.
c      if (.not.rtdb_get(inrtdb,'qmmm:xyzq',mt_log,1,xyzq)) xyzq=.false.
c
cc      if(xyz) then
        call util_file_name('xyz', .false., .false., filename)
        open(unit=43,status="unknown",form="formatted",file=filename)
        if(.not.geom_print_xyz(igeom,43))
     +    call errquit('qmmm: Failed to xyz file',0)
        close(43)
cc      end if

      if(.not.geom_rtdb_store(rtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)


      if(.not.ma_pop_stack(l_i))
     & call errquit('qmmm: Failed to deallocate stack',tot_nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(l_t))
     & call errquit('qmmm: Failed to deallocate stack',tot_nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(l_q))
     & call errquit('qmmm: Failed to deallocate stack',tot_nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(l_c))
     & call errquit('qmmm: Failed to deallocate stack',tot_nqm,
     &       MA_ERR)
      if(.not.ma_pop_stack(l_m))
     & call errquit('qmmm: Failed to deallocate stack',tot_nqm,
     &       MA_ERR)
      
      end
c
      subroutine qmmm_end(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
c
      call mm_end()
c
      end

