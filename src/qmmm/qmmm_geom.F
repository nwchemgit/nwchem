c
c $Id: qmmm_geom.F,v 1.1 2004-06-15 18:17:34 marat Exp $
c
c
      subroutine qmmm_create_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer i
      integer ntot
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
      integer i_c_all, h_c_all
      integer i_q_all, h_q_all
      integer i_t_all, h_t_all
      character*16 myname
      character*4 bqtag

c****************************************************
c     Stage 1 get cart info for quantum atoms
c     Stage 2 get cart info for solute point charges
c     Stage 3 get cart info for solvent point charges
c     Stage 4 combine all three into geometry object
c***************************************************

c     -----------------------------------------
c     Stage 1: collecting info on quantum atoms
c     -----------------------------------------
      nqm   = qmmm_get_nat()
      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()


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

c     ------------------------------------------------
c     Stage 2: collecting info on solute point charges
c     ------------------------------------------------

c      call mm_create_bq_geom(rtdb)
       call mm_get_tot_nbqs(nbqs)
       write(*,*) "tot_nbq=",nbqs
      if(.not.ma_push_get(mt_dbl,3*nbqs,'cbq',h_cbq,i_cbq))
     + call errquit('qmmm: Failed to allocate memory for cbq',3*nbqs,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbqs,'qbq',h_qbq,i_qbq))
     + call errquit('qmmm: Failed to allocate memory for qbq',nbqs,
     &       MA_ERR)
      

      call mm_get_geom_bqs(nbqs,dbl_mb(i_cbq),dbl_mb(i_qbq))
      write(*,*) "after mm_get_geom_bqs"
      myname = "test-bq.pdb"
      bqtag="Bq"
      call qmmm_print_bq_pdb(nbqs,myname,dbl_mb(i_cbq),
     >      dbl_mb(i_qbq),bqtag)


      if(.not.ma_pop_stack(h_qbq))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_cbq))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)

c     -------------------------------------------------
c     Stage 3: collecting info on solvent point charges
c     -------------------------------------------------
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
      bqtag="Bqw"
      call qmmm_print_bq_pdb(nbqw,myname,dbl_mb(i_cbqw),
     >      dbl_mb(i_qbqw),bqtag)

      if(.not.ma_pop_stack(h_qbqw))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_cbqw))
     & call errquit('qmmm: Failed to deallocate stack',nbqs,
     &       MA_ERR)

c
c     Stage 4 creating geometry
c     ------------------------
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

