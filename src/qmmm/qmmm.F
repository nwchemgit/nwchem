c
c $Id: qmmm.F,v 1.3 2004-03-12 23:37:05 marat Exp $
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

      subroutine qmmm_create_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer i
      integer nqm
      integer nbqs
      integer nbqw
      integer igeom
      integer i_c
      integer i_q
      integer i_t
      integer i_m
      integer i_ind
      integer i_num
      integer i_cbq,h_cbq
      integer i_qbq,h_qbq
      integer i_cbqw,h_cbqw
      integer i_qbqw,h_qbqw
      integer i_qf,h_qf
      character*16 myname

      nqm   = qmmm_get_nat()

      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c = qmmm_get_i_c()
      i_q = qmmm_get_i_q()
      i_t = qmmm_get_i_t()
      i_m = qmmm_get_i_m()


      call mm_get_geom(nqm,
     >                 int_mb(i_ind),
     >                 int_mb(i_num),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q),
     >                 dbl_mb(i_m))

      myname="test.pdb"
      call qmmm_print_pdb(nqm,myname,dbl_mb(i_c),
     >        dbl_mb(i_q),byte_mb(i_t))
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.ma_push_get(mt_dbl,nqm,'qf',h_qf,i_qf))
     + call errquit('qmmm: Failed to allocate memory for qf',3*nqm,
     &       MA_ERR)

      do i=1,nqm
       dbl_mb(i_qf+i-1)=dble(int_mb(i_num+i-1))
      end do
 
      if(.not.geom_cart_set(igeom,nqm,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_qf)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.ma_pop_stack(h_qf))
     & call errquit('qmmm: Failed to deallocate stack',nqm,
     &       MA_ERR)
 
      if(.not.geom_masses_set(igeom,nqm,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)


c     create bq geometry
c     ------------------
c      call mm_create_bq_geom(rtdb)
       call mm_get_tot_nbqs(nbqs)
       write(*,*) "tot_nbq=",nbqs
      if(.not.ma_push_get(mt_dbl,3*nbqs,'cbq',h_cbq,i_cbq))
     + call errquit('qmmm: Failed to allocate memory for cbq',3*nbqs,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbqs,'qbq',h_qbq,i_qbq))
     + call errquit('qmmm: Failed to allocate memory for qbq',nbqs,
     &       MA_ERR)
      

      write(*,*) "before mm_get_geom_bqs"
      call mm_get_geom_bqs(nbqs,dbl_mb(i_cbq),dbl_mb(i_qbq))
      write(*,*) "after mm_get_geom_bqs"
      myname = "test-bq.pdb"
      call qmmm_print_bq_pdb(nbqs,myname,dbl_mb(i_cbq),
     >      dbl_mb(i_qbq),"Bq")


      if(.not.ma_pop_stack(h_qbq))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_cbq))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)

       call mm_get_tot_nbqw(nbqw)
       write(*,*) "tot_nbqw=",nbqw
      if(.not.ma_push_get(mt_dbl,3*nbqw,'cbqw',h_cbqw,i_cbqw))
     + call errquit('qmmm: Failed to allocate memory for cbqw',3*nbqw,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbqw,'qbq',h_qbqw,i_qbqw))
     + call errquit('qmmm: Failed to allocate memory for qbq',nbqw,
     &       MA_ERR)
      

      write(*,*) "before mm_get_geom_bqw"
      call mm_get_geom_bqw(nbqw,dbl_mb(i_cbqw),dbl_mb(i_qbqw))
      write(*,*) "after mm_get_geom_bqw"
      myname = "test-bqw.pdb"
      call qmmm_print_bq_pdb(nbqw,myname,dbl_mb(i_cbqw),
     >      dbl_mb(i_qbqw),"Bqw")

      if(.not.ma_pop_stack(h_qbqw))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_cbqw))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)
  
      end

      subroutine qmmm_print_pdb(nt,myname,c,q,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(3,nt)
      double precision q(nt)
      character*16 t(nt)
      character*16 myname

      integer i
      integer n
      character*(nw_max_path_len) filename

      n=inp_strlen(myname)
      call util_file_name(myname(1:n), .false.,
     >               .false., filename)
      open(unit=46,status="unknown",form="formatted",file=myname)


      do i=1,nt
       write(46,FMT=9000)
     >           i,t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      call util_flush(46)
      close(46)


      end

      subroutine qmmm_print_bq_pdb(nt,myname,c,q,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(3,nt)
      double precision q(nt)
      character*4 t
      character*16 myname

      integer i
      integer n
      character*(nw_max_path_len) filename

      n=inp_strlen(myname)
      call util_file_name(myname(1:n), .false.,
     >               .false., filename)
      open(unit=46,status="unknown",form="formatted",file=myname)


      do i=1,nt
       write(46,FMT=9000)
     >           i,t,c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      call util_flush(46)
      close(46)


      end

