c
c $Id: qmmm_geom.F,v 1.61 2005-06-01 16:27:03 marat Exp $
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

      call mm_get_solute_tot_nact(ngs)
      call mm_get_solvent_tot_nact(ngw)
      ng = ngs+ngw
      write(*,*) "number of active waters",ngw

      if(.not.ma_alloc_get(mt_int,ng,'qmmm geom index',h_ig,i_ig))
     + call errquit(pname//'Failed to allocate heap',ng,
     &       MA_ERR)

      call mm_get_solute_act_ind(ngs,int_mb(i_ig))
      call mm_get_solvent_act_ind(ngw,int_mb(i_ig+ngs))

      call qmmm_geom_create_active(irtdb)
c      call qmmm_geom_create(irtdb,
c     >                      ngs,
c     >                      int_mb(i_ig),
c     >                      ngw,
c     >                      int_mb(i_ig+ngs),
c     >                      "geometry",.false.)

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

      subroutine qmmm_geom_create(irtdb,ns,ais,nw,aiw,aname,link)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
      integer irtdb
      integer ns
      integer ais(ns)
      integer nw
      integer aiw(nw)
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
      integer n
 
      pname = "qmmm_geom_create"
      if(qmmm_print_debug())
     > write(*,*) "in "//pname

      n = ns+nw
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

      if(ns.ne.0) 
     >  call mm_get_solute_geom(ns,
     >                 ais,
     >                 int_mb(i_atn),
     >                 byte_mb(i_t),
     >                 dbl_mb(i_c),
     >                 dbl_mb(i_q))

      if(nw.ne.0) 
     >  call mm_get_solvent_geom(nw,
     >                   aiw,
     >                   int_mb(i_atn+ns),
     >                   byte_mb(i_t+16*ns),
     >                   dbl_mb(i_c+3*ns),
     >                   dbl_mb(i_q+ns))

      if(.not.geom_create(igeom,aname))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

      if(link)  
     +   call qmmm_links_adjust(ns,ais,byte_mb(i_t),
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

      subroutine qmmm_geom_create_full(irtdb)
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

      pname = "qmmm_geom_create_full"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      write(*,*) "nfg=",nfg
      call qmmm_geom_create(irtdb,
     >                      nfg,
     >                      int_mb(i_ifg),
     >                      0,
     >                      0,
     >                      "full_geometry",.true.)


      if(.not.rtdb_cget(irtdb,'geometry',1,oldgeom))
     + oldgeom = ' '

      if(.not.rtdb_cput(irtdb,'geometry',1,'full_geometry'))
     + call errquit(pname//' storing geom name to rtdb',0, RTDB_ERR)
     
      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_create_active(irtdb)
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

      pname = "qmmm_geom_create_active"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      call qmmm_geom_create(irtdb,
     >                      ngs,
     >                      int_mb(i_ig),
     >                      ngw,
     >                      int_mb(i_ig+ngs),
     >                      "geometry",.false.)


      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_push_active(irtdb)
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

      pname = "qmmm_geom_push_active"

      if(qmmm_print_debug())
     >    write(*,*) "in",pname

      call qmmm_geom_push(irtdb,
     >                      ngs,
     >                      int_mb(i_ig),
     >                      ngw,
     >                      int_mb(i_ig+ngs),
     >                      "geometry")


      if(qmmm_print_debug())
     >    write(*,*) "out",pname

      end

      subroutine qmmm_geom_restore(irtdb)
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

      subroutine qmmm_geom_push(irtdb,ns,ais,nw,aiw,aname)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
      integer irtdb
      integer ns
      integer ais(ns)
      integer nw
      integer aiw(nw)
      character*(*) aname    
c     local variables
      integer ntot
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

      pname = "qmmm_geom_push"

      ntot = ns+nw

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

      if(ncent.ne.ntot) 
     >    call errquit(pname//"mismatch in geometries",
     >                 ncent+10000*ntot,0)

      call mm_set_solute_coord(ns,
     >                         ais,
     >                         dbl_mb(i_ctmp))

      call mm_set_solvent_coord(nw,
     >                         aiw,
     >                         dbl_mb(i_ctmp+3*ns))


      if(.not.ma_pop_stack(h_ctmp))
     & call errquit('qmmm: 
     >              Failed to deallocate stack c_tmp',ncent,
     &       MA_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

      end

      subroutine qmmm_geom_print_current(irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh" 
#include "geom.fh"
#include "global.fh"
      integer irtdb
      integer igeom

      character*23 pname

      pname = "qmmm_geom_print_current"
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

      function qmmm_geom_get_ngs()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_ngs

      qmmm_geom_get_ngs = ngs

      end 

      function qmmm_geom_get_ngw()
      implicit none
#include "qmmm_geom_data.fh"

      integer qmmm_geom_get_ngw

      qmmm_geom_get_ngw = ngw

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

