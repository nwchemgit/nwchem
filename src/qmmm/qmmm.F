c
c $Id: qmmm.F,v 1.5 2004-06-15 18:17:34 marat Exp $
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
      integer nqm
      integer i_ind

      qmmm_init = .true.

c     initialize mm part
c     (note it will also want operation)
c     --------------------------------------------
      call mm_init(irtdb)

c     get total number of quantum atoms
c     --------------------------------
      call mm_get_tot_nqm(nqm)

c     allocate qmmm data structure
c     ---------------------------
      call qmmm_data_alloc(nqm)

c     establish mapping between mm and qm atoms
c     ----------------------------------------
      i_ind = qmmm_get_i_ai()
      call mm_get_ind(nqm,int_mb(i_ind))

c     create qm geometry
c     ------------------
      call qmmm_create_geom(irtdb)
    
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
      character*255 theory
      double precision energy
      double precision mm_energy
      double precision eatoms
      integer tot_nbq
c     external functions
      logical  task_energy_doit
      external task_energy_doit
      logical  qmmm_init
      external qmmm_init
      logical  qmmm_end
      external qmmm_end

      qmmm_energy = .false.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

c     QM energy calculation
c     --------------------
      status = task_energy_doit(rtdb,theory,energy)

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed storing energy', 0, RTDB_ERR)

 
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
      write(*,*) "qmmm energy",energy,mm_energy,energy+mm_energy,eatoms
      energy = energy + mm_energy

      write(*,*) "qmmm energy",energy,mm_energy

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
      external task_gradient_doit
      character*255 theory
      double precision energy
      double precision mm_energy
      integer i_g 
      double precision eatoms
      integer nat
c
      status = .true.
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

c
c    first do quantum part
c    --------------------
      nat = qmmm_get_nat()
      i_g = qmmm_get_i_g()
      status = task_gradient_doit(rtdb,theory,energy,dbl_mb(i_g))
      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed storing energy', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'qmmm:gradient',mt_dbl,3*nat,dbl_mb(i_g)))
     $        call errquit('oniom: failed storing gradient', 0,
     &       RTDB_ERR)
 
c
c     now do mm part
c     ----------------------------------
      call qmmm_update_mm_coord(rtdb)
      call mm_task_sp()
      call qmmm_update_qm_forces(rtdb)

      if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)

      if (.not.rtdb_get(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy', 0, RTDB_ERR)

      mm_energy = mm_energy/cau2kj + 0.5d0*eatoms
      energy = energy + mm_energy

      write(*,*) "qmmm energy",energy,mm_energy

      if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)

      qmmm_gradient=status

      end
c
      function qmmm(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer rtdb
      logical qmmm
      character*30 operation
      logical status
      logical task_energy
      external task_energy
      logical qmmm_init
      external qmmm_init
      logical qmmm_end
      external qmmm_end
c
c
      status = .true.
      write(6,*) "initializing qmmm"
      if( .not. qmmm_init(rtdb)) then
       call errquit('qmmm_init',0,INPUT_ERR)
      end if
c
c     determine calculation type
c     --------------------------
      if(.not.rtdb_cget(rtdb,'task:operation',1,operation))
     + call errquit('qmmm: no calculation type',0,INPUT_ERR)
c
      if (inp_compare(.false., 'energy', operation)) then
         status = task_energy(rtdb)
      end if
      if( .not. qmmm_end(rtdb)) then
       call errquit('qmmm_end',0,INPUT_ERR)
      end if
c
      qmmm=status
      return
      end

      subroutine qmmm_update_mm_coord(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer nt
      integer igeom
      integer i_ai
      integer i_c

      nt   = qmmm_get_nat()
      i_ai = qmmm_get_i_ai()
      i_c = qmmm_get_i_c()


      if(.not.geom_create(igeom,'geometry'))
     & call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     & call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

      if(.not.geom_cart_coords_get(igeom, dbl_mb(i_c)))
     & call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_destroy(igeom))
     & call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      call mm_set_geom(nt,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_c))

      end

      subroutine qmmm_update_qm_forces(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
      integer irtdb
c     local variables
      integer nt
      integer i_ai
      integer i_gmm
      integer i_g

      nt   = qmmm_get_nat()
      i_ai = qmmm_get_i_ai()
      i_g = qmmm_get_i_g()
      i_gmm = qmmm_get_i_gmm()

      call dfill(3*nt,0,dbl_mb(i_gmm),1)

      call mm_get_force(nt,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_gmm))


      if(.not.rtdb_get(irtdb,'qmmm:gradient',mt_dbl,3*nt,dbl_mb(i_g)))
     + call errquit('qmmm: rtdb_get gradient failed',0, RTDB_ERR)

      call daxpy(3*nt,cau2nm/cau2kj,dbl_mb(i_gmm),1,dbl_mb(i_g),1)

      if(.not.rtdb_put(irtdb,'qmmm:gradient',mt_dbl,3*nt,dbl_mb(i_g)))
     + call errquit('qmmm: rtdb_out gradient failed',0, RTDB_ERR)

      end

c      subroutine qmmm_create_geom(irtdb)
c      implicit none
c#include "mafdecls.fh"
c#include "errquit.fh"
c#include "geom.fh"
c#include "qmmm.fh"
c      integer irtdb
cc     local variables
c      integer i
c      integer ntot
c      integer nqm
c      integer nbqs
c      integer nbqw
c      integer igeom
c      integer i_c
c      integer i_q
c      integer i_t
c      integer i_m
c      integer i_ind
c      integer i_num
c      integer i_cbq,h_cbq
c      integer i_qbq,h_qbq
c      integer i_cbqw,h_cbqw
c      integer i_qbqw,h_qbqw
c      integer i_qf,h_qf
c      integer i_c_all, h_c_all
c      integer i_q_all, h_q_all
c      integer i_t_all, h_t_all
c      character*16 myname
c      character*4 bqtag
c
cc****************************************************
cc     Stage 1 get cart info for quantum atoms
cc     Stage 2 get cart info for solute point charges
cc     Stage 3 get cart info for solvent point charges
cc     Stage 4 combine all three into geometry object
cc***************************************************
c
cc     -----------------------------------------
cc     Stage 1: collecting info on quantum atoms
cc     -----------------------------------------
c      nqm   = qmmm_get_nat()
c      i_ind = qmmm_get_i_ai()
c      i_num = qmmm_get_i_num()
c      i_c   = qmmm_get_i_c()
c      i_q   = qmmm_get_i_q()
c      i_t   = qmmm_get_i_t()
c      i_m   = qmmm_get_i_m()
c
c
c      call mm_get_geom(nqm,
c     >                 int_mb(i_ind),
c     >                 int_mb(i_num),
c     >                 byte_mb(i_t),
c     >                 dbl_mb(i_c),
c     >                 dbl_mb(i_q),
c     >                 dbl_mb(i_m))
c
c      myname="test.pdb"
c      call qmmm_print_pdb(nqm,myname,dbl_mb(i_c),
c     >        dbl_mb(i_q),byte_mb(i_t))
c
cc     ------------------------------------------------
cc     Stage 2: collecting info on solute point charges
cc     ------------------------------------------------
c
cc      call mm_create_bq_geom(rtdb)
c       call mm_get_tot_nbqs(nbqs)
c       write(*,*) "tot_nbq=",nbqs
c      if(.not.ma_push_get(mt_dbl,3*nbqs,'cbq',h_cbq,i_cbq))
c     + call errquit('qmmm: Failed to allocate memory for cbq',3*nbqs,
c     &       MA_ERR)
c      if(.not.ma_push_get(mt_dbl,nbqs,'qbq',h_qbq,i_qbq))
c     + call errquit('qmmm: Failed to allocate memory for qbq',nbqs,
c     &       MA_ERR)
c      
c
c      call mm_get_geom_bqs(nbqs,dbl_mb(i_cbq),dbl_mb(i_qbq))
c      write(*,*) "after mm_get_geom_bqs"
c      myname = "test-bq.pdb"
c      bqtag="Bq"
c      call qmmm_print_bq_pdb(nbqs,myname,dbl_mb(i_cbq),
c     >      dbl_mb(i_qbq),bqtag)
c
c
c      if(.not.ma_pop_stack(h_qbq))
c     & call errquit('qmmm: Failed to deallocate stack',nbqs,
c     &       MA_ERR)
c
c      if(.not.ma_pop_stack(h_cbq))
c     & call errquit('qmmm: Failed to deallocate stack',nbqs,
c     &       MA_ERR)
c
cc     -------------------------------------------------
cc     Stage 3: collecting info on solvent point charges
cc     -------------------------------------------------
c       call mm_get_tot_nbqw(nbqw)
c       write(*,*) "tot_nbqw=",nbqw
c      if(.not.ma_push_get(mt_dbl,3*nbqw,'cbqw',h_cbqw,i_cbqw))
c     + call errquit('qmmm: Failed to allocate memory for cbqw',3*nbqw,
c     &       MA_ERR)
c      if(.not.ma_push_get(mt_dbl,nbqw,'qbq',h_qbqw,i_qbqw))
c     + call errquit('qmmm: Failed to allocate memory for qbq',nbqw,
c     &       MA_ERR)
c      
c
c      write(*,*) "before mm_get_geom_bqw"
c      call mm_get_geom_bqw(nbqw,dbl_mb(i_cbqw),dbl_mb(i_qbqw))
c      write(*,*) "after mm_get_geom_bqw"
c      myname = "test-bqw.pdb"
c      bqtag="Bqw"
c      call qmmm_print_bq_pdb(nbqw,myname,dbl_mb(i_cbqw),
c     >      dbl_mb(i_qbqw),bqtag)
c
c      if(.not.ma_pop_stack(h_qbqw))
c     & call errquit('qmmm: Failed to deallocate stack',nbqs,
c     &       MA_ERR)
c
c      if(.not.ma_pop_stack(h_cbqw))
c     & call errquit('qmmm: Failed to deallocate stack',nbqs,
c     &       MA_ERR)
c
cc
cc     Stage 4 creating geometry
cc     ------------------------
c      ntot = nqm+nbqs+nbqw
c      if(.not.ma_push_get(mt_dbl,3*nbqw,'c_all',h_c_all,i_c_all))
c     + call errquit('qmmm: Failed to allocate memory for c_all',
c     + 3*ntot, MA_ERR)
c      if(.not.ma_push_get(mt_dbl,nbqw,'q_all',h_q_all,i_q_all))
c     + call errquit('qmmm: Failed to allocate memory for q_all',ntot,
c     &       MA_ERR)
c
c       do i=1,ntot
c
c
c       end do
c
c
c      if(.not.ma_pop_stack(h_c_all))
c     & call errquit('qmmm: Failed to deallocate stack',nbqs,
c     &       MA_ERR)
c
c      if(.not.ma_pop_stack(h_q_all))
c     & call errquit('qmmm: Failed to deallocate stack',nbqs,
c     &       MA_ERR)
c
c
c   
c      if(.not.geom_create(igeom,'geometry'))
c     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
c
c      if(.not.ma_push_get(mt_dbl,nqm,'qf',h_qf,i_qf))
c     + call errquit('qmmm: Failed to allocate memory for qf',3*nqm,
c     &       MA_ERR)
c
c      do i=1,nqm
c       dbl_mb(i_qf+i-1)=dble(int_mb(i_num+i-1))
c      end do
c 
c      if(.not.geom_cart_set(igeom,nqm,byte_mb(i_t),
c     + dbl_mb(i_c),dbl_mb(i_qf)))
c     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)
c
c      if(.not.ma_pop_stack(h_qf))
c     & call errquit('qmmm: Failed to deallocate stack',nqm,
c     &       MA_ERR)
c 
c      if(.not.geom_masses_set(igeom,nqm,dbl_mb(i_m)))
c     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
c      call geom_compute_values(igeom)
c
c      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
c     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)
c
c       if(.not.geom_destroy(igeom))
c     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)
c
c  
c      end
c
c      subroutine qmmm_print_pdb(nt,myname,c,q,t)
c      implicit none
c#include "mafdecls.fh"
c#include "errquit.fh"
c#include "qmmm.fh"
c#include "util.fh"
c#include "inp.fh"
c      integer nt
c      double precision c(3,nt)
c      double precision q(nt)
c      character*16 t(nt)
c      character*16 myname
c
c      integer i
c      integer n
c      character*(nw_max_path_len) filename
c
c      n=inp_strlen(myname)
c      call util_file_name(myname(1:n), .false.,
c     >               .false., filename)
c      open(unit=46,status="unknown",form="formatted",file=myname)
c
c
c      do i=1,nt
c       write(46,FMT=9000)
c     >           i,t(i),c(1,i)*0.529177249d00,
c     >            c(2,i)*0.529177249d00,
c     >            c(3,i)*0.529177249d00, q(i)
c
c      end do
c9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)
c
c      call util_flush(46)
c      close(46)
c
c
c      end
c
c      subroutine qmmm_print_bq_pdb(nt,myname,c,q,t)
c      implicit none
c#include "mafdecls.fh"
c#include "errquit.fh"
c#include "qmmm.fh"
c#include "util.fh"
c#include "inp.fh"
c      integer nt
c      double precision c(3,nt)
c      double precision q(nt)
c      character*4 t
c      character*16 myname
c
c      integer i
c      integer n
c      character*(nw_max_path_len) filename
c
c      n=inp_strlen(myname)
c      call util_file_name(myname(1:n), .false.,
c     >               .false., filename)
c      open(unit=46,status="unknown",form="formatted",file=myname)
c
c
c      do i=1,nt
c       write(46,FMT=9000)
c     >           i,t,c(1,i)*0.529177249d00,
c     >            c(2,i)*0.529177249d00,
c     >            c(3,i)*0.529177249d00, q(i)
c
c      end do
c9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)
c
c      call util_flush(46)
c      close(46)
c
c
c      end
c
