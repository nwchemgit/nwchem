c
c $Id: qmmm_geom.F,v 1.59 2005-05-27 17:59:22 marat Exp $
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
#include "qmmm_geom_data.fh"
#include "inp.fh"
#include "qmmm_params.fh"
      integer irtdb
c     local variables
      integer i
      integer nqm,nlink
      integer i_ind,i_lb
      integer i_tmp,h_tmp
      integer i_act,h_act
      integer ig,ng0
      integer nbq,i_ibq
      character*30 region
      character*32 pname

      pname = "qmmm_geom_init"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      nqm   = qmmm_get_nqm()
      nlink = qmmm_get_nlink()
      i_ind = qmmm_get_i_ai()
      i_lb  = qmmm_get_i_lb()

      nbq = qmmm_get_nbq()
      i_ibq = qmmm_get_i_ibq()

      call mm_get_tot_nact(ng)

      if(.not.ma_alloc_get(mt_int,ng,'qmmm geom index',h_ig,i_ig))
     + call errquit(pname//'Failed to allocate heap',ng,
     &       MA_ERR)

      call mm_get_solute_act_ind(ng,int_mb(i_ig))

      write(*,*) "ng,i_ig",ng,(int_mb(i_ig+i-1),i=1,ng)
 
      call qmmm_create_geom(irtdb,
     >                      ng,
     >                      int_mb(i_ig),
     >                      "geometry",.false.)


      nfg = nqm + nlink

      if(.not.ma_alloc_get(mt_int,nfg,'qmmm fullg ind',h_ifg,i_ifg))
     + call errquit(pname//'Failed to allocate heap',nfg,
     &       MA_ERR)
      call icopy(nqm,int_mb(i_ind),1,int_mb(i_ifg),1)
      call icopy(nlink,int_mb(i_lb+nlink),1,int_mb(i_ifg+nqm),1)

      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_init1(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
#include "qmmm_params.fh"
      integer irtdb
c     local variables
      integer i
      integer nqm,nlink
      integer i_ind,i_lb
      integer i_tmp,h_tmp
      integer i_act,h_act
      integer ig,ng0
      integer nbq,i_ibq
      character*30 region
      character*32 pname

      pname = "qmmm_geom_init"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      nqm   = qmmm_get_nqm()
      nlink = qmmm_get_nlink()
      i_ind = qmmm_get_i_ai()
      i_lb  = qmmm_get_i_lb()

      nbq = qmmm_get_nbq()
      i_ibq = qmmm_get_i_ibq()

c
c     initializing auxiliary geometry index
c     -------------------------------------
      if(.not.ma_push_get(mt_int,max_driver_atoms,
     +               'qmmm geom tmp index',h_tmp,i_tmp))
     + call errquit(pname//'Failed to allocate heap',ng0,
     &       MA_ERR)

      if (.not. rtdb_cget(irtdb, 'qmmm:opt_region', 1, region)) then
         region = "qm"
      end if
     
      if(region.eq."qm") then
        ng0 = nqm
        call icopy(ng0,int_mb(i_ind),1,int_mb(i_tmp),1)
      else if(region.eq."qmlink") then
        ng0 = nqm+nlink
        call icopy(nqm,int_mb(i_ind),1,int_mb(i_tmp),1)
        call icopy(nlink,int_mb(i_lb+nlink),1,int_mb(i_tmp+nqm),1)
      else if(region.eq."bq") then
        ng0 = nbq 
        call icopy(nbq,int_mb(i_ibq),1,int_mb(i_tmp),1)
      else if(region.eq."qmbq") then
        ng0 = nbq+nqm
        call icopy(nqm,int_mb(i_ind),1,int_mb(i_tmp),1)
        call icopy(nbq,int_mb(i_ibq),1,int_mb(i_tmp+nqm),1)
      else if(region.eq."user") then
        ng0 = 1
        int_mb(i_tmp) = 24
      else 
        call errquit("unknown opt_region",0,0)
      end if


      if(.not.ma_push_get(mt_log,ng0,'qmmm_geom temp',h_act,i_act))
     + call errquit(pname//'Failed to allocate heap',ng0,
     &       MA_ERR)
c       
c     get mm constraints if any
c     ------------------------
      call mm_get_solute_cons(ng0,
     >                       int_mb(i_tmp),
     >                       log_mb(i_act))

      ng = 0
      do i=1,ng0
        if( log_mb(i_act+i-1)) ng = ng +1
      end do


      if(.not.ma_alloc_get(mt_int,ng,'qmmm geom index',h_ig,i_ig))
     + call errquit(pname//'Failed to allocate heap',ng,
     &       MA_ERR)

      ig = 0
      do i=1,ng0
        if(log_mb(i_act+i-1)) then
          write(*,*) "found"
          ig = ig + 1
          int_mb(i_ig+ig-1) = int_mb(i_tmp+i-1) 
        end if  
      end do

      if(.not.ma_pop_stack(h_act))
     + call errquit(pname//'Failed to deallocate memory ',ng,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_tmp))
     + call errquit(pname//'Failed to deallocate memory ',ng,
     &       MA_ERR)


      write(*,*) "ng,i_ig",ng,(int_mb(i_ig+i-1),i=1,ng)
 
      call qmmm_create_geom(irtdb,
     >                      ng,
     >                      int_mb(i_ig),
     >                      "geometry",.true.)


      nfg = nqm + nlink

      if(.not.ma_alloc_get(mt_int,nfg,'qmmm fullg ind',h_ifg,i_ifg))
     + call errquit(pname//'Failed to allocate heap',nfg,
     &       MA_ERR)
      call icopy(nqm,int_mb(i_ind),1,int_mb(i_ifg),1)
      call icopy(nlink,int_mb(i_lb+nlink),1,int_mb(i_ifg+nqm),1)



      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_data_release()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_geom_data.fh"
      

      if(nfg.ne.0) then
       if(.not.ma_free_heap(h_ifg))
     &   call errquit('qmmm h_ifg: Failed to deallocate heap',h_ifg,
     &       MA_ERR)
       nfg = 0
      end if


      if(ng.ne.0) then
       if(.not.ma_free_heap(h_ig))
     &   call errquit('qmmm h_ig: Failed to deallocate heap',ng,
     &       MA_ERR)
       ng = 0
      end if

       return
      end

      subroutine qmmm_geom_init0(irtdb)
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
      logical hide_geom
      character*30 operation
      character*32 pname

      pname = "qmmm_geom_init"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      hide_geom  = .true.
      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'
      if (operation.eq.'property') hide_geom = .false.

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


      call mm_get_solute_geom_quant(nqm,
     >                 int_mb(i_ind),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q),
     >                 dbl_mb(i_m))


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

c     ---------------------
c     start geometry hiding 
c     ---------------------
      do i=1,nqm
        log_mb(i_act+i-1)=.true.
      end do

      if(hide_geom) then
        call mm_get_solute_cons_quant(nqm,
     >                         int_mb(i_ind),
     >                         log_mb(i_act))


        if(qmmm_h_link() .and. nlink.gt.0) then
          do i=1,nlink
            i1=nqm-nlink+i
            log_mb(i_act+i1-1)=.false.
          end do
        end if

        if(qmmm_fixlinks() .and. nlink.gt.0) then

          if(.not.ma_push_get(mt_int,nlink,'mqm',h_mqm,i_mqm))
     +       call errquit( pname//'Failed to allocate memory for mqm',
     +       nlink, MA_ERR)
           if(.not.qmmm_map(nqm,int_mb(i_ai),
     >                   nlink,int_mb(i_lb),
     >                   int_mb(i_mqm)))  
     >       call errquit( pname//'Failed to map', 0, 0)


          do i=1,nlink
            i1=int_mb(i_mqm+i-1)
            log_mb(i_act+i1-1)=.false.
          end do

          if(.not.ma_pop_stack(h_mqm))
     +     call errquit( pname//'Failed to deallocate memory for mqm',
     +     nlink, MA_ERR)

        end if
c     ----------------------
c     end of geometry hiding 
c     ----------------------

        if(qmmm_print_debug()) then
          write(*,*) "list of active atoms"
          do i=1,nqm
            write(*,*) i, log_mb(i_act+i-1)
          end do
        end if

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

      subroutine qmmm_create_geom(irtdb,n,ai,aname,link)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
      integer n
      integer ai(n)
      character*(*) aname
      logical link
c
      integer i
      integer igeom
      character*32 pname
      integer i_q,h_q
      integer i_c,h_c
      integer i_t,h_t
      integer i_m,h_m
      integer i_atn,h_atn
      integer i_il
 
      pname = "qmmm_create_geom"
      if(qmmm_print_debug())
     > write(*,*) "in "//pname

      if(.not.ma_push_get(mt_dbl,3*n,'c',h_c,i_c))
     + call errquit('qmmm: Failed to allocate memory for c',
     + 3*n, MA_ERR)
      if(.not.ma_push_get(mt_dbl,n,'q',h_q,i_q))
     + call errquit('qmmm: Failed to allocate memory for q',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,n,'m',h_m,i_m))
     + call errquit('qmmm: Failed to allocate memory for m',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_int,n,'inum',h_atn,i_atn))
     + call errquit('qmmm: Failed to allocate memory for atn',n,
     &       MA_ERR)
      if(.not.ma_push_get(mt_byte,16*n,'t',h_t,i_t))
     + call errquit('qmmm: Failed to allocate memory for t',n,
     &       MA_ERR)

      call mm_get_solute_geom(n,
     >                 ai,
     >                 int_mb(i_atn),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q))

      if(.not.geom_create(igeom,aname))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(link)  
     +   call qmmm_links_adjust(n,ai,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q),int_mb(i_atn))


      do i=1,n
         if(.not.geom_atn_to_default_mass(int_mb(i_atn+i-1),
     >                                    dbl_mb(i_m+i-1)))  
     >    call errquit(pname,0, GEOM_ERR)
      end do

      if(.not.geom_cart_set(igeom,n,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,n,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)

      if(.not.geom_print(igeom))
     + call errquit('qmmm: Failed to print geom',0, RTDB_ERR)

      if(.not.geom_rtdb_store(irtdb,igeom,aname))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
c     release temporary memory
c     ------------------------
      if(.not.ma_pop_stack(h_t))
     & call errquit('qmmm: Failed to deallocate stack t_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_atn))
     & call errquit('qmmm: Failed to deallocate stack t_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_m))
     & call errquit('qmmm: Failed to deallocate stack m_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_q))
     & call errquit('qmmm: Failed to deallocate stack q_all',n,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_c))
     & call errquit('qmmm: Failed to deallocate stack c_all',n,
     &       MA_ERR)

      if(qmmm_print_debug())
     > write(*,*) "out "//pname
 
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
#include "mm_utils.fh"
#include "qmmm_utils.fh"
#include "qmmm_geom_data.fh"
#include "inp.fh"
      integer irtdb
c     local variables
      character*32 pname

      pname = "create full geom"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      write(*,*) "nfg=",nfg
      call qmmm_create_geom(irtdb,
     >                      nfg,
     >                      int_mb(i_ifg),
     >                      "full_geometry",.true.)


      if(.not.rtdb_cget(irtdb,'geometry',1,oldgeom))
     + oldgeom = ' '

      if(.not.rtdb_cput(irtdb,'geometry',1,'full_geometry'))
     + call errquit(pname//' storing geom name to rtdb',0, RTDB_ERR)
     
      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_nbq_act(nbq_act)
      implicit none
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_geom_data.fh"
#include "qmmm.fh"
c
      integer nbq_act
c
      integer i,j
      integer ind
      integer nbq,i_ibq

      nbq = qmmm_get_nbq()
      i_ibq = qmmm_get_i_ibq()
      
      nbq_act = 0
      do i=1,ng
        ind = int_mb(i_ig+i-1)
        do j=1,nbq
          if(ind.eq.int_mb(i_ibq+j-1)) then
             nbq_act = nbq_act + 1
          end if
        end do
      end do

      end

      subroutine qmmm_restore_geom(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_geom_data.fh"
#include "bq.fh"
      integer irtdb
c     local variables
      logical ignore
      integer hbq
      double precision charge
      character*32 pname
c
      pname = "qmmm_restore_qm_geom"

      if(qmmm_print_debug())
     > write(*,*) " in "//pname

      ignore = rtdb_delete(irtdb,'geometry')

      if(oldgeom.ne.' ') then
        if(.not.rtdb_cput(irtdb,'geometry',1,oldgeom))
     +    call errquit(pname//' storing geom name to rtdb',0, RTDB_ERR)
      end if

      if(qmmm_print_debug())
     > write(*,*) " out "//pname

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

      subroutine qmmm_push_coord_mm(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_geom_data.fh"
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

      pname = "qmmm_push_coord_mm"

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

      if(ncent.ne.ng) 
     >    call errquit(pname//"mismatch in geometries",ncent+10000*ng,0)

      call mm_set_solute_coord(ng,
     >                         int_mb(i_ig),
     >                         dbl_mb(i_ctmp))


      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

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

      function qmmm_geom_get_nfg()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_nfg

      qmmm_geom_get_nfg = nfg

      end 

      function qmmm_geom_get_ng()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_ng

      qmmm_geom_get_ng = ng

      end 

      function qmmm_geom_get_i_ig()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_i_ig

      qmmm_geom_get_i_ig = i_ig

      end 

      function qmmm_geom_get_i_ifg()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_i_ifg

      qmmm_geom_get_i_ifg = i_ifg

      end 

