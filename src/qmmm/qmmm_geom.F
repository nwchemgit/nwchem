c
c $Id: qmmm_geom.F,v 1.16 2004-09-17 20:07:00 marat Exp $
c
c
      subroutine qmmm_geom_init(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
      integer irtdb
c     local variables
      integer i,i1
      integer nqm,nlink
      integer igeom
      integer i_c
      integer i_t
      integer i_m
      integer i_ind
      integer i_num
      integer i_q
      character*255 basis_string
      integer qmmm_geom_heavy_neighbor
      external qmmm_geom_heavy_neighbor


      nqm   = qmmm_get_nqm()
      nlink = qmmm_get_nlink()
      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()

      call mm_get_solute_geom_quant(nqm,
     >                 int_mb(i_ind),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q),
     >                 dbl_mb(i_m))


c
c     adjust link atom names and
c     basis set if needed
c     --------------------------
      do i=1,nqm
        int_mb(i_num+i-1)=int(dbl_mb(i_q+i-1))
      end do
      if(nlink.gt.0) then
c       call qmmm_link_atoms(nlink,
c     >                      int_mb(i_num+(nqm-nlink)),
c     >                      dbl_mb(i_q+(nqm-nlink)),
c     >                      byte_mb(i_t+16*(nqm-nlink)),
c     >                      dbl_mb(i_c+3*(nqm-nlink)),
c     >                      dbl_mb(i_m+(nqm-nlink)))
       call qmmm_link_atoms1(nlink,
     >                       nqm,
     >                      int_mb(i_num),
     >                      dbl_mb(i_q),
     >                      byte_mb(i_t),
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_m))

c       do i=1,nlink
c          i1=qmmm_geom_heavy_neighbor(nqm-nlink+i,nqm,int_mb(i_num),
c     >                                dbl_mb(i_c))
c          write(*,*) "qmmm link",i1,int_mb(i_num+i1-1)
c       end do
      end if

      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set(igeom,nqm,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,nqm,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)
 
c     print geometry
c     -------------
c      call util_print_centered(6,"QM/MM geometry", 36, .true.)
c      if(.not.geom_print(igeom))
c     +    call errquit('qmmm: Failed to print geometry',0)

      if(.not.geom_print_xyz(igeom,22))
     +    call errquit('qmmm: Failed to print geometry',0)
      call util_flush(22)

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

      subroutine qmmm_sync_geom(irtdb)
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
      integer i_cbq
      integer i_c_all, h_c_all
      double precision charge
      nbq  = qmmm_get_nbq()
      
      nqm  = qmmm_get_nqm()
      nbq  = qmmm_get_nbq()
      ntot = nbq+nqm
      i_c  = qmmm_get_i_c()
      i_cbq = qmmm_get_i_cbq()

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

c
c     copy coordinates
c     ------------------
      do i=1,3*nqm
       dbl_mb(i_c+i-1)= dbl_mb(i_c_all+i-1)
      end do

      do i=1,3*nbq
       dbl_mb(i_cbq+i-1)= dbl_mb(i_c_all+3*nqm+i-1)
      end do


      if(.not.ma_pop_stack(h_c_all))
     & call errquit('qmmm: Failed to deallocate stack c_all',ntot,
     &       MA_ERR)


      end
     
      subroutine qmmm_create_bq_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
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
      integer i_q
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
         call qmmm_bq_data_load()
      end if

      nqm   = qmmm_get_nqm()
      nbq   = qmmm_get_nbq()
      ntot  = nbq+nqm

c     nothing to do if no bq's
c     -----------------------
      if(nbq.eq.0) then
        write(*,*) "qmmm: no bq charges were found"
        return
      end if

      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
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
     >                    dbl_mb(i_q),
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

      if(.not.geom_print_pdb(igeom,23))
     +    call errquit('qmmm: Failed to print pdb geometry',0)

c     print geometry
c     -------------
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,"QM/MM geometry", 36, .true.)
      if(.not.geom_print(igeom))
     +    call errquit('qmmm: Failed to print geometry',0)
      end if


      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
c     release temporary memory
c     ------------------------
      if(.not.ma_pop_stack(h_t_all))
     & call errquit('qmmm: Failed to deallocate stack t_all',ntot,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m_all))
     & call errquit('qmmm: Failed to deallocate stack m_all',ntot,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q_all))
     & call errquit('qmmm: Failed to deallocate stack q_all',ntot,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c_all))
     & call errquit('qmmm: Failed to deallocate stack c_all',ntot,
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
      nbq  = qmmm_get_nbq()
      ntot = nbq+nqm
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
     & call errquit('qmmm: Failed to deallocate stack c_all',ntot,
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

      implicit none
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

      subroutine qmmm_link_atoms(n,inum,q,t,c,m)
      implicit none
#include "geom.fh"
#include "errquit.fh"
      integer n
      double precision inum(n)
      double precision q(n)
      double precision c(3,n)
      character*16 t(n)
      double precision m(n)

c     local variables
c     ---------------
      integer i
      integer i1

      do i=1,n
        q(i)=9.0
        t(i)="F_L"
        if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do

      end

      subroutine qmmm_link_atoms1(nl,n,inum,q,t,c,m)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "qmmm.fh"
      integer nl
      integer n
      double precision inum(n)
      double precision q(n)
      double precision c(3,n)
      character*16 t(n)
      double precision m(n)

c     local variables
c     ---------------
      integer i
      integer i1
      integer k
      logical h_link
      integer qmmm_geom_heavy_neighbor
      external qmmm_geom_heavy_neighbor
      double precision r1(3),r2(3)
      double precision qmmm_geom_distance
      external qmmm_geom_distance

      h_link=qmmm_h_link()
      if(h_link) then
      do i=n-nl+1,n
        i1=qmmm_geom_heavy_neighbor(i,n,inum,c)
        do k=1,3
         r1(k)=c(k,i)
         r2(k)=c(k,i1)
        end do
        write(*,*) "qmmm link bond",t(i),t(i1),qmmm_geom_distance(r1,r2)
        call qmmm_adjust_link(r1,r2)
        do k=1,3
         c(k,i) =r1(k)
         c(k,i1)=r2(k)
        end do
        t(i)="H_L"
        q(i)=1.0
        write(*,*) "qmmm link bond",t(i),t(i1),qmmm_geom_distance(r1,r2)
      end do
      else
      do i=n-nl+1,n
        q(i)=9.0
        t(i)="F_L"
        if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do
      end if

      end

      subroutine qmmm_geom_prepare(irtdb)
      implicit none
      integer irtdb
#include "qmmm.fh"

      call qmmm_sync_qm_geom(irtdb)
      call qmmm_create_bq_geom(irtdb)


      end

      subroutine qmmm_geom_finish(irtdb)
      implicit none
      integer irtdb
#include "qmmm.fh"

      call qmmm_restore_qm_geom(irtdb)

      end

      function qmmm_geom_heavy_neighbor(i,ntot,atn,c)

      implicit none
      integer i
      integer ntot
      double precision c(3,ntot)
      integer atn(ntot)
      integer qmmm_geom_heavy_neighbor
c 
      double precision qmmm_geom_distance
      external qmmm_geom_distance
      integer i1,k,j
      double precision r,rmin

      rmin=1000.0d0
      do i1=1,ntot
        if(atn(i1).ne.1) then
           r=qmmm_geom_distance(c(1,i),c(1,i1))
c          r=0.0d0
c          do k=1,3
c            r=r+(c(k,i)-c(k,i1))*(c(k,i)-c(k,i1))
c          end do
c          r=sqrt(r)
          write(*,*) "qmmm nei r",r
          if(r.lt.rmin) then
           if(r.gt. 1.0e-2) then
             write(*,*) "qmmm nei",r,atn(i1)
             j=i1
             rmin=r
           endif
          endif
        endif
      end do
      qmmm_geom_heavy_neighbor=j

      end

      function qmmm_geom_distance(r1,r2)

      implicit none
      double precision r1(3),r2(3)
      double precision qmmm_geom_distance
c 
      integer k
      double precision r

      r=0.0d0
      do k=1,3
        r=r+(r1(k)-r2(k))*(r1(k)-r2(k))
      end do
      r=sqrt(r)

      qmmm_geom_distance=r

      end

      subroutine qmmm_adjust_link(r1,r2)

      implicit none
      double precision r1(3),r2(3)
      double precision qmmm_geom_distance
      external  qmmm_geom_distance
c 
      integer k
      double precision dr(3)
      double precision r
      double precision rb

      r=qmmm_geom_distance(r1,r2)
      rb=1.0d0/0.529
      do k=1,3
        dr(k)=r2(k)-r1(k)
        dr(k)=dr(k)*(r-rb)/r
      end do

      do k=1,3
        r1(k)=r1(k)+dr(k)
      end do

      end

