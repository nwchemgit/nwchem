c
c $Id: qmmm_geom.F,v 1.7 2004-08-31 19:45:48 marat Exp $
c
c
      subroutine qmmm_create_qm_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer i
      integer nqm
      integer igeom
      integer i_c
      integer i_t
      integer i_m
      integer i_ind
      integer i_num
      integer i_q
      integer i_qf
      character(255) basis_string

      nqm   = qmmm_get_nqm()
      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
      i_qf  = qmmm_get_i_qf()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()


      call mm_get_geom(nqm,
     >                 int_mb(i_ind),
     >                 int_mb(i_num),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q),
     >                 dbl_mb(i_m))

      do i=1,nqm
       dbl_mb(i_qf+i-1)=dble(int_mb(i_num+i-1))
      end do

      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set(igeom,nqm,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_qf)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,nqm,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)
 
c     print geometry
c     -------------
      call util_print_centered(6,"QM/MM geometry", 36, .true.)
      if(.not.geom_print(igeom))
     +    call errquit('qmmm: Failed to print geometry',0)

      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)
       
      basis_string = 'ecp; F library "qmmm_zhang_6-31Gs_ecp"; end'
      call nw_inp_from_character(irtdb, basis_string)

      end

      subroutine qmmm_sync_qm_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer nqm

      integer igeom
      integer i_c

      nqm   = qmmm_get_nqm()
      i_c   = qmmm_get_i_c()

c    --------------------------------------
c    get qm coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_cart_coords_get(igeom,dbl_mb(i_c))) 
     >    call errquit("qmmm:geom_cart_coords_get",0,0)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      end
     
      subroutine qmmm_create_bq_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer i
      integer nqm
      integer nbq
      integer ntot

      integer igeom
      integer i_c
      integer i_t
      integer i_m
      integer i_ind
      integer i_num
      integer i_qf
      integer i_cbq
      integer i_qbq
      integer i_q_all,h_q_all
      integer i_c_all,h_c_all
      integer i_t_all,h_t_all
      integer i_m_all,h_m_all
      logical nobq
      double precision charge

      if (.not. rtdb_get(irtdb,'qmmm:remove_bq',mt_log,1,nobq)) then
        nobq=.false.
      end if
 
c
c     load bq information (coord and charges)
c     into qmmm common arrays from md module
c     --------------------------------------
      if(.not.nobq) then
         call qmmm_load_bq()
      end if

      nqm   = qmmm_get_nqm()
      nbq   = qmmm_get_nbq()
      ntot  = qmmm_get_ntot()

c     nothing to do if no bq's
c     -----------------------
      if(nbq.eq.0) then
        write(*,*) "qmmm: no bq charges were found"
        return
      end if

      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_qf  = qmmm_get_i_qf()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()

      i_cbq = qmmm_get_i_cbq()
      i_qbq = qmmm_get_i_qbq()


c    --------------------------------------
c    get qm coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_store(irtdb,igeom,'qm_geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_cart_coords_get(igeom,dbl_mb(i_c))) 
     >    call errquit("qmmm:geom_cart_coords_get",0,0)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_delete(irtdb,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

c    --------------------------------------
c    create combined qm bq geometry
c    --------------------------------------
      if(.not.ma_push_get(mt_dbl,3*ntot,'c_all',h_c_all,i_c_all))
     + call errquit('qmmm: Failed to allocate memory for c_all',
     + 3*ntot, MA_ERR)
      if(.not.ma_push_get(mt_dbl,ntot,'q_all',h_q_all,i_q_all))
     + call errquit('qmmm: Failed to allocate memory for q_all',ntot,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,ntot,'m_all',h_m_all,i_m_all))
     + call errquit('qmmm: Failed to allocate memory for q_all',ntot,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*ntot,'t_all',h_t_all,i_t_all))
     + call errquit('qmmm: Failed to allocate memory for q_all',ntot,
     &       MA_ERR)

      call qmmm_comb_geom(ntot,
     >                    dbl_mb(i_c_all),
     >                    dbl_mb(i_q_all),
     >                    byte_mb(i_t_all),
     >                    dbl_mb(i_m_all),
     >                    nqm,
     >                    dbl_mb(i_c),
     >                    dbl_mb(i_qf),
     >                    byte_mb(i_t),
     >                    dbl_mb(i_m),
     >                    nbq,
     >                    dbl_mb(i_cbq),
     >                    dbl_mb(i_qbq))


      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set(igeom,ntot,byte_mb(i_t_all),
     + dbl_mb(i_c_all),dbl_mb(i_q_all)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,ntot,dbl_mb(i_m_all)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

      if(.not.geom_print_xyzq(igeom,67))
     +    call errquit('qmmm: Failed to print geometry',0)

      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
c     release temporary memory
c     ------------------------
      if(.not.ma_pop_stack(h_t_all))
     & call errquit('qmmm: Failed to deallocate stack',ntot,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m_all))
     & call errquit('qmmm: Failed to deallocate stack',ntot,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q_all))
     & call errquit('qmmm: Failed to deallocate stack',ntot,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c_all))
     & call errquit('qmmm: Failed to deallocate stack',ntot,
     &       MA_ERR)


c    
c     set total charge
c     -----------------------------
      charge =  qmmm_get_qcharge()
      do i=1,nbq
       charge = charge + dbl_mb(i_qbq+i-1)
      end do
      if(.not.rtdb_put(irtdb,'charge',mt_dbl,1,charge)) then
          call errquit('qmmm: Failed to set charge',0,
     &       RTDB_ERR)
      end if

      end

      subroutine qmmm_restore_qm_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer i
      integer ntot
      integer nqm
      integer nbq
      integer igeom
      integer i_c
      integer i_c_all, h_c_all
      double precision charge
      nbq  = qmmm_get_nbq()
      
c     nothing to do if no bq's
c     -----------------------
      if(nbq.eq.0) return

      nqm  = qmmm_get_nqm()
      ntot = qmmm_get_ntot()
      i_c  = qmmm_get_i_c()

c
c     get coord information from combined
c     qm and bq geometry
c     -----------------------------------
      if(.not.ma_push_get(mt_dbl,3*ntot,'c_all',h_c_all,i_c_all))
     + call errquit('qmmm: Failed to allocate memory for c_all',
     + 3*ntot, MA_ERR)

      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,"geometry"))
     >    call errquit("qmmm:load geometry",0,0)

      if(.not. geom_cart_coords_get(igeom,dbl_mb(i_c_all))) 
     >    call errquit("qmmm:geom_cart_coords_get",0,0)

      if(.not.geom_destroy(igeom))
     >    call errquit("qmmm:load geometry",0,0)

      if(.not.geom_rtdb_delete(irtdb,"geometry"))
     >    call errquit("qmmm:delete geometry",0,0)

c
c     copy qm coordinates
c     ------------------
      do i=1,3*nqm
       dbl_mb(i_c+i-1)= dbl_mb(i_c_all+i-1)
      end do

      if(.not.ma_pop_stack(h_c_all))
     & call errquit('qmmm: Failed to deallocate stack',ntot,
     &       MA_ERR)

c
c     create geometry out qm coordnates only
c     --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,"qm_geometry"))
     >    call errquit("qmmm:load geometry",0,0)

      if(.not. geom_cart_coords_set(igeom,dbl_mb(i_c))) 
     >    call errquit("qmmm:geom_cart_coords_set",0,0)

      if(.not.geom_rtdb_store(irtdb,igeom,"geometry"))
     >    call errquit("qmmm:load geometry",0,0)

      if(.not.geom_destroy(igeom))
     >    call errquit("qmmm:load geometry",0,0)

c    
c     return charge to normal
c     -----------------------------
      charge =  qmmm_get_qcharge()
      if(.not.rtdb_put(irtdb,'charge',mt_dbl,1,charge)) then
          call errquit('qmmm: Failed to set charge',0,
     &       RTDB_ERR)
      end if


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

      subroutine qmmm_comb_geom(ntot,c_all,q_all,t_all,m_all,
     >                          nqm,c,q,t,m,
     >                          nbq,cbq,qbq)

      integer ntot
      double precision c_all(3,ntot)
      double precision q_all(ntot)
      character*16 t_all(ntot)
      double precision m_all(ntot)
      integer nqm
      double precision c(3,nqm)
      double precision q(nqm)
      character*16 t(nqm)
      double precision m(nqm)
      integer nbq
      double precision cbq(3,nbq)
      double precision qbq(nbq)

c     local variables
c     ---------------
      integer i

c
c     coordinates 
c     ----------
      do i=1,nqm
        c_all(1,i) = c(1,i)  
        c_all(2,i) = c(2,i)  
        c_all(3,i) = c(3,i)  
        q_all(i) = q(i)  
        t_all(i) = t(i)  
        m_all(i) = m(i)  
      end do
      do i=1,nbq
        c_all(1,nqm+i)=cbq(1,i)
        c_all(2,nqm+i)=cbq(2,i)
        c_all(3,nqm+i)=cbq(3,i)
        q_all(nqm+i)=qbq(i)
        t_all(nqm+i)="Bq"
        m_all(nqm+i)=0.0d0
      end do

      end


