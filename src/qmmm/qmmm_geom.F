c
c $Id: qmmm_geom.F,v 1.55 2005-04-09 01:42:01 marat Exp $
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
#include "qmmm_utils.fh"
#include "inp.fh"
      integer irtdb
c     local variables
      integer i,i1
      integer nqm,nlink,nps
      integer igeom
      integer i_c
      integer i_t
      integer i_m
      integer i_ind
      integer i_act
      integer i_num
      integer i_q
      integer i_ai
      integer i_ips
      integer i_lb
      character*255 basis_string,namespace
      character*255 basis_name
      character*255 link_ecp
      character*255 geom
      integer qmmm_geom_heavy_neighbor
      external qmmm_geom_heavy_neighbor
      integer tmp(2)
      integer nb
      integer k
      integer h_mqm,i_mqm
      character*32 pname

      pname = "qmmm_geom_init"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      nqm   = qmmm_get_nqm()
      nps   = qmmm_get_nps()
      nlink = qmmm_get_nlink()
      i_ind = qmmm_get_i_ai()
      i_act = qmmm_get_i_act()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()
      i_ai  = qmmm_get_i_ai()
      i_lb  = qmmm_get_i_lb()
      i_ips = qmmm_get_i_ips()

c
c     activate all atoms by default
c     -----------------------------
      do i=1,nqm
        log_mb(i_act+i-1)=.true.
      end do


      call mm_get_solute_geom_quant(nqm,
     >                 int_mb(i_ind),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q),
     >                 dbl_mb(i_m))

      if(qmmm_print_debug())
     >    write(*,*) "mm_get_solute_cons_quant",pname

      call mm_get_solute_cons_quant(nqm,
     >                       int_mb(i_ind),
     >                       log_mb(i_act))

      if(qmmm_print_debug())
     >    write(*,*) "out mm_get_solute_cons_quant",pname

      do i=1,nqm
        int_mb(i_num+i-1)=int(dbl_mb(i_q+i-1))
      end do

c     adjust link atom names and
c     basis set if needed
c     --------------------------
      if(nlink.gt.0) then

       call qmmm_adjust_link_params(nlink,
     >                      nqm,
     >                      int_mb(i_num),
     >                      dbl_mb(i_q),
     >                      byte_mb(i_t),
     >                      dbl_mb(i_m))

      end if

      if(nps.gt.0) then
       call qmmm_adjust_flinks(nqm,
     >                         int_mb(i_ai),
     >                         nps,
     >                         int_mb(i_ips),
     >                         int_mb(i_num),
     >                         dbl_mb(i_q),
     >                         byte_mb(i_t),
     >                         dbl_mb(i_m))

      
      end if

      call qmmm_adjust_link_coord()

      if(qmmm_print_debug()) then
        call util_print_centered(6,"QM system information",
     >              36, .true.)
        write(*,'(/,A,T20,A,T30,A,/)') "global index",
     >                                 "name","coordinates(au)"
        do i=1,nqm
         write(*,'(I5,T20,A,T30,3F12.6)') int_mb(i_ind+i-1),
     >              byte_mb(i_t+16*(i-1)),
     >              (dbl_mb(i_c+3*(i-1)+k-1),k=1,3)
        end do
      end if

c 
c     force inactivation of hydrogen links
c     ------------------------------------
      if(qmmm_h_link() .and. nlink.gt.0) then
        do i=1,nlink
          i1=nqm-nlink+i
          log_mb(i_act+i1-1)=.false.
        end do
      end if

      if(qmmm_fixlinks() .and. nlink.gt.0) then

        if(.not.ma_push_get(mt_int,nlink,'mqm',h_mqm,i_mqm))
     +     call errquit( pname//'Failed to allocate memory for mqm',
     +     nlink, MA_ERR)
         if(.not.qmmm_map(nqm,int_mb(i_ai),
     >                 nlink,int_mb(i_lb),
     >                 int_mb(i_mqm)))  
     >     call errquit( pname//'Failed to map', 0, 0)


        do i=1,nlink
          i1=int_mb(i_mqm+i-1)
          log_mb(i_act+i1-1)=.false.
        end do

        if(.not.ma_pop_stack(h_mqm))
     +   call errquit( pname//'Failed to deallocate memory for mqm',
     +   nlink, MA_ERR)

      end if


      if(qmmm_print_debug()) then
        write(*,*) "list of active atoms"
        do i=1,nqm
          write(*,*) i, log_mb(i_act+i-1)
        end do
      end if

      if (.not.rtdb_cget(irtdb,"qmmm:geom",1,geom)) then
        call  qmmm_create_qm_geom(irtdb)
      else 
        call  qmmm_load_qm_geom(geom,irtdb)
        call  qmmm_get_qm_coord_qm(irtdb)
        call  qmmm_push_qm_coord_mm(irtdb)
      end if
      
c
c     set link atoms ecp if any
c     ------------------------
      if (.not.rtdb_cget(irtdb,"qmmm:link_ecp",1,link_ecp)) 
     >    call errquit('qmmm_geom: failed to get qmmm:link_ecp',
     &       0,RTDB_ERR)
c
      if(link_ecp.eq."auto") then
        if(qmmm_halogen_link()) then
          if(qmmm_master()) then
            write(6,'(//)')
            call util_print_centered(6,
     >            "Setting up ecp basis set on psedo-carbon links",
     >            36, .true.)
            write(6,'(//)')
          end if
          basis_name = "qmmm_zhang_6-31Gs_ecp"
          nb = inp_strlen(basis_name)
          basis_string='ecp; F_L library "'//basis_name(1:nb)//'"; end'
          call nw_inp_from_character(irtdb, basis_string)
        end if
      end if


       call qmmm_print_pdb(nqm,"-qm-start.pdb",
     >                    dbl_mb(i_c),
     >                    dbl_mb(i_q),
     >                    byte_mb(i_t))

       call qmmm_bq_data_load()

      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_create_qm_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
c
      integer i
      integer j
      integer icent
      integer i_c
      integer i_t
      integer i_q
      integer i_m
      integer i_ind
      integer nact
      integer nqm
      integer i_act
      integer igeom
      character*32 pname
      integer i_q_all,h_q_all
      integer i_c_all,h_c_all
      integer i_t_all,h_t_all
      integer i_m_all,h_m_all
 
      pname = "qmmm_create_qm_geom"
      if(qmmm_print_debug())
     > write(*,*) "in "//pname

      nqm   = qmmm_get_nqm()
      i_ind = qmmm_get_i_ai()
      i_act = qmmm_get_i_act()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()

c    
c     count number of active atoms
c     ---------------------------
      nact = 0
      do i=1,nqm
       if(log_mb(i_act+i-1)) nact = nact + 1
      end do

      if(.not.ma_push_get(mt_dbl,3*nact,'c_all',h_c_all,i_c_all))
     + call errquit('qmmm: Failed to allocate memory for c_all',
     + 3*nact, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nact,'q_all',h_q_all,i_q_all))
     + call errquit('qmmm: Failed to allocate memory for q_all',nact,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,nact,'m_all',h_m_all,i_m_all))
     + call errquit('qmmm: Failed to allocate memory for q_all',nact,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*nact,'t_all',h_t_all,i_t_all))
     + call errquit('qmmm: Failed to allocate memory for q_all',nact,
     &       MA_ERR)

      icent = 0
      do i=1,nqm
        if(log_mb(i_act+i-1)) then
          icent = icent + 1
          do j=1,16
            byte_mb(i_t_all+16*(icent-1)+j-1)=byte_mb(i_t+16*(i-1)+j-1)
          end do
          dbl_mb(i_m_all+(icent-1))=dbl_mb(i_m+(i-1))
          dbl_mb(i_q_all+(icent-1))=dbl_mb(i_q+(i-1))
          do j=1,3
            dbl_mb(i_c_all+(icent-1)*3+j-1)=dbl_mb(i_c+(i-1)*3+j-1)
          end do
        end if
      end do
 
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set(igeom,nact,byte_mb(i_t_all),
     + dbl_mb(i_c_all),dbl_mb(i_q_all)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,nact,dbl_mb(i_m_all)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

      if(ga_nodeid().eq.0) then
        call util_print_centered(6,
     >      "initial QMMM geometry (active atoms only)",
     >              36, .true.)
        if(.not.geom_print_pdb(igeom,6))
     +    call errquit('qmmm: Failed to print pdb geometry',0)
        call util_print_centered(6,
     >       "initial full QMMM geometry globally indexed",
     >              36, .true.)
        call qmmm_print_pdbi(nqm,int_mb(i_ind),
     >                       dbl_mb(i_c),
     >                       dbl_mb(i_q),
     >                       byte_mb(i_t))
       end if

      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
c     release temporary memory
c     ------------------------
      if(.not.ma_pop_stack(h_t_all))
     & call errquit('qmmm: Failed to deallocate stack t_all',nact,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m_all))
     & call errquit('qmmm: Failed to deallocate stack m_all',nact,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q_all))
     & call errquit('qmmm: Failed to deallocate stack q_all',nact,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c_all))
     & call errquit('qmmm: Failed to deallocate stack c_all',nact,
     &       MA_ERR)

      if(qmmm_print_debug())
     > write(*,*) "out "//pname
 
      end

      subroutine qmmm_load_qm_geom(geom,irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
      character*(*) geom
      integer irtdb
c
      integer igeom
      character*32 pname

      pname = "qmmm_load_qm_geom"
      if(qmmm_print_debug())
     > write(*,*) "in "//pname

      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,geom))
     + call errquit('qmmm: Failed to load geom to rtdb',0, RTDB_ERR)

      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      if(qmmm_print_debug())
     > write(*,*) "out "//pname
 
      end

      subroutine qmmm_create_full_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "bq.fh"
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
      integer i_cbq, h_cbq
      integer i_qbq, h_qbq
      logical nobq
      logical bq_dynamical
      double precision charge
      character*32 pname

      integer hbq

      pname = "qmmm_create_full_geom0"

      if(qmmm_print_debug())
     > write(*,*) " in "//pname

      if (.not. rtdb_get(irtdb,'qmmm:remove_bq',mt_log,1,nobq)) then
        nobq=.false.
      end if

c
c     load bq information (coord and charges)
c     into qmmm common arrays from md module
c     --------------------------------------
      
      bq_dynamical = qmmm_bq_dynamical()
      if(qmmm_print_debug())
     >       write(*,*) pname,"bq_dynamical",bq_dynamical
      if(.not.nobq) then
         if(bq_dynamical) then
           call qmmm_bq_data_load()
           if(qmmm_print_debug())
     >       write(*,*) pname,"updated Bq's"
         else
           if(qmmm_print_debug())
     >       write(*,*) pname,"Bq's are not updated"
         end if
      end if

      nqm   = qmmm_get_nqm()
      nbq   = qmmm_get_nbq()
      ntot  = nbq+nqm

c     nothing to do if no bq's
c     -----------------------
      if(nbq.eq.0) then
      if(qmmm_print_debug())
     > write(*,*) "qmmm: no bq charges were found"
c        return
      end if

      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()

      i_cbq = qmmm_get_i_cbq()
      i_qbq = qmmm_get_i_qbq()
      h_cbq = qmmm_get_h_cbq()
      h_qbq = qmmm_get_h_qbq()

c
c     fix hydrogen link atoms
c     -----------------------
      call qmmm_adjust_link_coord()

c     ---------------------
c     activate bq structure
c     ---------------------
      if(.not.bq_create("qmmm charges",hbq))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR) 
c      if(.not.bq_set(hbq,nbq,dbl_mb(i_qbq),dbl_mb(i_cbq)))
      if(.not.bq_pset(hbq,nbq,h_qbq,h_cbq))
     + call errquit(pname//'Failed bq_pset',0,CALC_ERR)
      if(.not.bq_activate(hbq)) 
     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)

      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set(igeom,nqm,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,nqm,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

      if(ga_nodeid().eq.0) then
        call util_print_centered(6,"QMMM geometry with Bq charges",
     >              36, .true.)
        if(.not.geom_print_pdb(igeom,6))
     +    call errquit('qmmm: Failed to print pdb geometry',0)
      end if

      if(.not.geom_rtdb_store(irtdb,igeom,'full_geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

      if(qmmm_print_debug())
     > write(*,*) pname//" : full geometry is stored in rtdb"


       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      if(.not.rtdb_cput(irtdb,'geometry',1,'full_geometry'))
     + call errquit('qmmm: to store geom name to rtdb',0, RTDB_ERR)

       if(ga_nodeid().eq.0) 
     > call qmmm_print_pdb_bq(nbq,"-bq.pdb",
     >                    dbl_mb(i_cbq),
     >                    dbl_mb(i_qbq))



      if(qmmm_print_debug())
     > write(*,*) pname//" finished"

      end

      subroutine qmmm_create_full_geom0(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
      integer irtdb
c     local variables
      integer i
      integer k
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
      logical bq_dynamical
      double precision charge
      character*32 filename
      character*32 pname

      pname = "qmmm_create_full_geom"

      if(qmmm_print_debug())
     > write(*,*) " in "//pname

      if (.not. rtdb_get(irtdb,'qmmm:remove_bq',mt_log,1,nobq)) then
        nobq=.false.
      end if

c
c     load bq information (coord and charges)
c     into qmmm common arrays from md module
c     --------------------------------------
      
      bq_dynamical = qmmm_bq_dynamical()
      if(qmmm_print_debug())
     >       write(*,*) pname,"bq_dynamical",bq_dynamical
      if(.not.nobq) then
         if(bq_dynamical) then
           call qmmm_bq_data_load()
           if(qmmm_print_debug())
     >       write(*,*) pname,"updated Bq's"
         else
           if(qmmm_print_debug())
     >       write(*,*) pname,"Bq's are not updated"
         end if
      end if

      nqm   = qmmm_get_nqm()
      nbq   = qmmm_get_nbq()
      ntot  = nbq+nqm

c     nothing to do if no bq's
c     -----------------------
      if(nbq.eq.0) then
      if(qmmm_print_debug())
     > write(*,*) "qmmm: no bq charges were found"
c        return
      end if

      i_ind = qmmm_get_i_ai()
      i_num = qmmm_get_i_num()
      i_c   = qmmm_get_i_c()
      i_q   = qmmm_get_i_q()
      i_t   = qmmm_get_i_t()
      i_m   = qmmm_get_i_m()

      i_cbq = qmmm_get_i_cbq()
      i_qbq = qmmm_get_i_qbq()

c
c     fix hydrogen link atoms
c     -----------------------
      call qmmm_adjust_link_coord()
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


       if(ga_nodeid().eq.0) 
     > call qmmm_print_pdb(ntot,"-full.pdb",
     >                    dbl_mb(i_c_all),
     >                    dbl_mb(i_q_all),
     >                    byte_mb(i_t_all))


      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_cart_set(igeom,ntot,byte_mb(i_t_all),
     + dbl_mb(i_c_all),dbl_mb(i_q_all)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,ntot,dbl_mb(i_m_all)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

      if(ga_nodeid().eq.0) then
        call util_print_centered(6,"QMMM geometry with Bq charges",
     >              36, .true.)
        if(.not.geom_print_pdb(igeom,6))
     +    call errquit('qmmm: Failed to print pdb geometry',0)
      end if

      if(.not.geom_rtdb_store(irtdb,igeom,'full_geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

      if(qmmm_print_debug())
     > write(*,*) pname//" : full geometry is stored in rtdb"


       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      if(.not.rtdb_cput(irtdb,'geometry',1,'full_geometry'))
     + call errquit('qmmm: to store geom name to rtdb',0, RTDB_ERR)

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
      

      if(qmmm_print_debug())
     > write(*,*) pname//" finished"

      end

      subroutine qmmm_restore_qm_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "bq.fh"
      integer irtdb
c     local variables
      logical ignore
      integer hbq
      double precision charge
      character*32 pname
c
      pname = "qmmm_restore_qm_geom"


      ignore = rtdb_delete(irtdb,'geometry')
C      if(.not.rtdb_delete(irtdb,'geometry'))
C     + call errquit(pname//' to  delete geom name in rtdb',0, RTDB_ERR)

c      if(.not.bq_get_handle("qmmm charges",hbq))
c     + call errquit(pname//'cannot get qmmm bq handle',0, CALC_ERR)
c      
c      if(.not.bq_destroy(hbq))
c     + call errquit(pname//'cannot destroy qmmm bq handle',0, CALC_ERR)

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
      character*(*) myname
      character*(nw_max_path_len) filename

      integer i
      integer n
      integer nf
      integer ns

      call util_file_prefix("",filename)
      nf=inp_strlen(filename)-1
      ns=inp_strlen(myname)
      filename = filename(1:nf)//myname(1:ns)

      open(unit=46,status="unknown",form="formatted",file=filename)

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

      subroutine qmmm_print_pdb_bq(nt,myname,c,q)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      double precision c(3,nt)
      double precision q(nt)
      character*16 t
      character*(*) myname
      character*(nw_max_path_len) filename

      integer i
      integer n
      integer nf
      integer ns

      call util_file_prefix("",filename)
      nf=inp_strlen(filename)-1
      ns=inp_strlen(myname)
      filename = filename(1:nf)//myname(1:ns)

      open(unit=46,status="unknown",form="formatted",file=filename)

      t = "Bq"
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


      subroutine qmmm_print_pdbi(nt,ai,c,q,t)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
      integer nt
      integer ai(nt)
      double precision c(3,nt)
      double precision q(nt)
      character*16 t(nt)

      integer i
      integer n
      integer nf
      integer ns

      do i=1,nt
       write(6,FMT=9000)
     >           ai(i),t(i),c(1,i)*0.529177249d00,
     >            c(2,i)*0.529177249d00,
     >            c(3,i)*0.529177249d00, q(i)

      end do
9000  FORMAT("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

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

      subroutine qmmm_link_atoms0(nl,n,lind,inum,q,t,c,m)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "qmmm.fh"
      integer nl
      integer n
      integer lind(2*nl)
      integer inum(n)
      double precision q(n)
      double precision c(3,n)
      character*16 t(n)
      double precision m(n)

c     local variables
c     ---------------
      integer tmp(2)
      integer i
      integer i1
      integer il
      integer k
      logical h_link
      integer qmmm_geom_heavy_neighbor
      external qmmm_geom_heavy_neighbor
      double precision r1(3),r2(3)
      double precision qmmm_geom_distance
      external qmmm_geom_distance

      il = 0
      h_link=qmmm_h_link()
      if(h_link) then
      do i=n-nl+1,n
        i1=qmmm_geom_heavy_neighbor(i,n,inum,c)
        il=il+1
        lind(il)=i
        il=il+1
        lind(il)=i1
        do k=1,3
         r1(k)=c(k,i)
         r2(k)=c(k,i1)
        end do
        tmp(1)=i
        tmp(2)=i1
        call qmmm_adjust_link(r1,r2)
        do k=1,3
         c(k,i) =r1(k)
         c(k,i1)=r2(k)
        end do
        t(i)="H_L"
        q(i)=1.0
      end do
      else
      do i=n-nl+1,n
        il=il+1
        lind(il)=i
        il=il+1
        lind(il)=i
        q(i)=9.0
        t(i)="F_L"
        if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do
      end if

      end

      subroutine qmmm_adjust_link_params(nl,n,inum,q,t,m)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "qmmm.fh"
      integer nl
      integer n
      integer inum(n)
      double precision q(n)
      double precision c(3,n)
      character*16 t(n)
      double precision m(n)

c     local variables
c     ---------------
      integer i
      logical h_link

      if(nl.eq.0) return

      h_link=qmmm_h_link()
      if(h_link) then
      do i=n-nl+1,n
        t(i)="H_L"
        q(i)=1.0
        inum(i)=1
        if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do
      else
      do i=n-nl+1,n
        q(i)=9.0
        t(i)="F_L"
        inum(i)=9
        if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do
      end if

      end

      subroutine qmmm_adjust_flinks(n,ai,nps,ips,inum,q,t,m)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
#include "mafdecls.fh"
      integer n
      integer nps
      integer ai(n)
      integer ips(nps)
      integer inum(n)
      double precision q(n)
      character*16 t(n)
      double precision m(n)

c     local variables
c     ---------------
      integer i
      integer i_mqm
      integer iq
      integer h_mqm
      character*32 pname

      pname = "qmmm_adjust_flink"

      if(nps.eq.0) return

      if(.not.ma_push_get(mt_int,nps,'mqm',h_mqm,i_mqm))
     + call errquit( pname//'Failed to allocate memory for mqm',
     + nps, MA_ERR)

      if(.not.qmmm_map(n,ai,
     >                 nps,ips,
     >                 int_mb(i_mqm)))  
     >     call errquit( pname//'Failed to map',
     > 0, 0)

 
      do i=1,nps
        iq = int_mb(i_mqm+i-1)
        q(iq)=9.0
        t(iq)="F_L"
        inum(iq)=9
        if(.not.geom_tag_to_default_mass(t(iq),m(iq)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do

      if(.not.ma_pop_stack(h_mqm))
     + call errquit( pname//'Failed to deallocate memory for mqm',
     + nps, MA_ERR)

      end

      subroutine qmmm_adjust_link_coord()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
#include "mm_utils.fh"
c     local variables
      integer nqm
      integer nlink
      integer i,k
      integer i_c
      integer i_lb
      integer h_cl,i_cl
      integer h_mqm,i_mqm
      integer i_ai
      integer il,iq,ih
      integer ind_lb,ind_q
      double precision l_qh,l_ql
      double precision rh(3)
      double precision rq(3)
      double precision rl(3)
      double precision rql(3)
      double precision a,b
      character*32 pname

      double precision dnrm2
      external dnrm2

      pname = "qmmm_adjust_link_coord"

      if(qmmm_print_debug())
     >   write(*,*) "in "//pname

      nlink = qmmm_get_nlink()

c
c     nothing to do if no links
c     -------------------------
      if(nlink.eq.0) return
c
c     nothing to do if links are not hydrogens
c     ----------------------------------------
      if(.not.qmmm_h_link()) return

      nqm = qmmm_get_nqm()

      if(.not.ma_push_get(mt_dbl,3*nlink,'cl',h_cl,i_cl))
     + call errquit( pname//'Failed to allocate memory for cl',
     + 3*nlink, MA_ERR)

      if(.not.ma_push_get(mt_int,nlink,'mqm',h_mqm,i_mqm))
     + call errquit( pname//'Failed to allocate memory for mqm',
     + nlink, MA_ERR)


c
c     get pointer to index of link pairs
c     ----------------------------------
      i_lb = qmmm_get_i_lb()

c
c     get pointer to global index
c     ---------------------------
      i_ai   = qmmm_get_i_ai()
c
c     get pointer to coord array
c     --------------------------
      i_c   = qmmm_get_i_c()
c
c     get coordinates of mm link atoms
c     -------------------------------
      call mm_get_solute_coord_gen(nlink,
     &                   mm_link,
     &                   int_mb(i_lb+nlink),
     &                   dbl_mb(i_cl))

      if(.not.qmmm_map(nqm,int_mb(i_ai),
     >                 nlink,int_mb(i_lb),
     >                 int_mb(i_mqm)))  
     >     call errquit( pname//'Failed to map',
     > 0, 0)

c
c     loop over link atom pairs
c     -------------------------
      do i=1,nlink
        iq = int_mb(i_mqm+i-1)
        do k=1,3
           rq(k)  = dbl_mb(i_c+(iq-1)*3+k-1)
           rl(k)  = dbl_mb(i_cl+(i-1)*3+k-1)
           rql(k) = rq(k)-rl(k)
        end do
        l_ql=dnrm2(3,rql,1)
c        a = 1.0d0/(l_ql*0.529)
        a = 0.709d0
        ih = nqm-nlink+i
        do k=1,3
           rh(k) = rq(k) - rql(k)*a
           dbl_mb(i_c+(ih-1)*3+k-1) = rh(k)
        end do
      end do

      if(.not.ma_pop_stack(h_mqm))
     + call errquit( pname//'Failed to deallocate memory for cq',
     + 3*nlink, MA_ERR)
 
      if(.not.ma_pop_stack(h_cl))
     + call errquit( pname//'Failed to deallocate memory for cl',
     + 3*nlink, MA_ERR)

      end

      subroutine qmmm_link_index(nl,n,lind,inum,c)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "qmmm.fh"
      integer nl
      integer n
      integer lind(2*nl)
      integer inum(n)
      double precision c(3,n)

c     local variables
c     ---------------
      integer i
      integer i1
      integer il
      logical h_link
      integer qmmm_geom_heavy_neighbor
      external qmmm_geom_heavy_neighbor

      il = 0
      if(qmmm_h_link()) then
      do i=n-nl+1,n
        i1=qmmm_geom_heavy_neighbor(i,n,inum,c)
        il=il+1
        lind(il)=i
        il=il+1
        lind(il)=i1
      end do
      end if

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
      integer i1,j
      double precision r,rmin
      character*32 pname
   
      pname = "qmmm_geom_heavy_neighbor"

      write(*,*) pname, (atn(i1),i1=1,ntot)


      j = -1
      rmin=1000.0d0
      do i1=1,ntot
        if(atn(i1).ne.1) then
           r=qmmm_geom_distance(c(1,i),c(1,i1))
          if(r.lt.rmin) then
           if(r.gt. 1.0e-2) then
             j=i1
             rmin=r
           endif
          endif
        endif
      end do
      if(j.eq.-1) 
     > call errquit(pname//' failed to find  a link',0, 0)

      qmmm_geom_heavy_neighbor=j

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

      subroutine qmmm_push_qm_coord_mm(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      logical link
      integer nqm
      integer i_ai
      integer i_c

      link = .not.qmmm_h_link()

      nqm  = qmmm_get_nqm()
      i_ai = qmmm_get_i_ai()
      i_c = qmmm_get_i_c()

      call mm_set_solute_coord_quant(nqm,
     &                   link,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_c))

      end

      subroutine qmmm_get_qm_coord_mm(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      logical link
      integer nqm
      integer i_ai
      integer i_c

      link = .not.qmmm_h_link()

      nqm  = qmmm_get_nqm()
      i_ai = qmmm_get_i_ai()
      i_c = qmmm_get_i_c()

      call mm_get_solute_coord_quant(nqm,
     &                   link,
     &                   int_mb(i_ai),
     &                   dbl_mb(i_c))
        
c      call qmmm_bq_data_load()

      end

      subroutine qmmm_get_qm_coord_qm(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
      integer irtdb
c     local variables
      integer nqm
      integer ncent
      integer icent
      integer i
      integer k


      integer igeom
      integer i_c
      integer i_act
      integer i_ctmp,h_ctmp
      character*32 pname

      pname = "qmmm_get_qm_coord_qm"

      nqm   = qmmm_get_nqm()
      i_c   = qmmm_get_i_c()
      i_act = qmmm_get_i_act()

c    --------------------------------------
c    get qm coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(.not.ma_push_get(mt_dbl,3*ncent,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ncent, MA_ERR)

      if(.not. geom_cart_coords_get(igeom,dbl_mb(i_ctmp))) 
     >    call errquit("qmmm:geom_cart_coords_get",0,0)

      icent = 0
      do i=1,nqm
        if(log_mb(i_act+i-1)) then
          icent = icent + 1
          do k=1,3
            dbl_mb(i_c+3*(i-1)+k-1) = dbl_mb(i_ctmp+3*(icent-1)+k-1)
          end do
        end if
      end do


      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      end

      subroutine qmmm_push_qm_coord_qm(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
      integer irtdb
c     local variables
      integer nqm
      integer ncent
      integer icent
      integer i
      integer k

      integer igeom
      integer i_c
      integer i_act
      integer i_ctmp,h_ctmp
      character*32 pname

      pname = "qmmm_push_qm_coord_qm"

      nqm   = qmmm_get_nqm()
      i_c   = qmmm_get_i_c()
      i_act = qmmm_get_i_act()

c    --------------------------------------
c    get qm coordinates out of geometry
c    --------------------------------------
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom, ncent) ) 
     >    call errquit("qmmm:geom_ncent",0,0)

      if(.not.ma_push_get(mt_dbl,3*ncent,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ncent, MA_ERR)

      icent = 0
      do i=1,nqm
        if(log_mb(i_act+i-1)) then
          icent = icent + 1
          do k=1,3
            dbl_mb(i_ctmp+3*(icent-1)+k-1) = dbl_mb(i_c+3*(i-1)+k-1)
          end do
        end if
      end do

      if(.not. geom_cart_coords_set(igeom,dbl_mb(i_ctmp))) 
     >    call errquit("qmmm:geom_cart_coords_get",0,0)

      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('Failed to store geometry',0, GEOM_ERR)


       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      end

      subroutine qmmm_print_current_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh" 
#include "geom.fh"
#include "global.fh"
      integer irtdb
      integer igeom

      character*23 pname

      pname = "qmmm_print_current_geom"
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,"geometry"))
     >    call errquit(pname//":load geometry",0,0)

c     print geometry
c     -------------
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,"Current QM/MM geometry", 36, .true.)
      if(.not.geom_print(igeom))
     +    call errquit('qmmm: Failed to print geometry',0)
      end if

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)



      end

