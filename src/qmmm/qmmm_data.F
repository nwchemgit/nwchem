      subroutine qmmm_data_init(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"

      integer irtdb
c
      character*32 tag
      character*32 pname
      character*32 optimize
      character*32 interface
      real*8 eatoms,qatoms
      logical xyz,xyzq,fixlinks,fixqm
      integer linkatm,bq_exclude

      pname = "qmmm_data_init: "
      eatoms = 0.0d0
      fixlinks = .false.
      fixqm = .false.
      linkatm=qmmm_halogen
      bq_exclude=qmmm_hydrogen
      optimize='qm'
      interface='qm'
c
      call util_print_centered(6,"QMMM Parameters",32,.true.)

c     initialize rtdb parameters
      tag = "qmmm:uqmatm"
      if (.not.rtdb_get(irtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) then
        write(*,*) pname//"Warning! qm reference energy is set to", 
     >             eatoms
        if (.not.rtdb_put(irtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) 
     >      call errquit(pname//'failed to store eatoms',0,
     >       RTDB_ERR)
      end if

      tag = "qmmm:fixlinks"
      write(*,*) "doing "//tag
      if (.not.rtdb_get(irtdb,tag,mt_log,1,fixlinks)) then
        if (.not.rtdb_put(irtdb,tag,mt_log,1,fixlinks)) 
     >      call errquit(pname//'failed to store fixlinks',0,
     >       RTDB_ERR)

       end if

      tag = "qmmm:xyz"
      write(*,*) "doing "//tag
      if (.not.rtdb_get(irtdb,tag,mt_log,1,xyz)) then
        if (.not.rtdb_put(irtdb,tag,mt_log,1,xyz)) 
     >      call errquit(pname//'failed to store xyz',0,
     >       RTDB_ERR)

       end if

      tag = "qmmm:xyzq"
      write(*,*) "doing "//tag
      if (.not.rtdb_get(irtdb,tag,mt_log,1,xyzq)) then
        if (.not.rtdb_put(irtdb,tag,mt_log,1,xyzq)) 
     >      call errquit(pname//'failed to store xyzq',0,
     >       RTDB_ERR)

       end if

      tag ="qmmm:fixqm" 
      write(*,*) "doing "//tag
      if (.not.rtdb_get(irtdb,tag,mt_log,1,fixqm)) then
        if (.not.rtdb_put(irtdb,tag,mt_log,1,fixqm)) 
     >      call errquit(pname//'failed to store fixqm',0,
     >       RTDB_ERR)

       end if

      tag ="qmmm:bq_exclude" 
      write(*,*) "doing "//tag
      if (.not.rtdb_get(irtdb,tag,mt_int,1,bq_exclude)) then
        if (.not.rtdb_put(irtdb,tag,mt_int,1,bq_exclude)) 
     >      call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)

       end if

      tag ="qmmm:linkatm" 
      write(*,*) "doing "//tag
      if (.not.rtdb_get(irtdb,tag,mt_int,1,linkatm)) then
        if (.not.rtdb_put(irtdb,tag,mt_int,1,linkatm)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

       end if
     
      tag ="qmmm:interface"
      write(*,*) "doing "//tag
      if (.not.rtdb_cget(irtdb,tag,1,interface)) then
        if (.not.rtdb_cput(irtdb,tag,1,interface)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

       end if

      tag ="qmmm:optimize"
      write(*,*) "doing "//tag
      if (.not.rtdb_cget(irtdb,tag,1,optimize)) then
        if (.not.rtdb_cput(irtdb,tag,1,optimize)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

       end if

20    FORMAT(A32,A16,F12.6)
21    FORMAT(A32,A16,L2)

      call mm_get_tot_nqm(nqm)
      call mm_get_tot_nlink(nlink)
c
c     allocate data arrays
c     --------------------
      call qmmm_data_alloc(nqm)

c 
c     establish persistent mapping between mm and qm atoms
c     first go pure quantum atoms and then links
c     ----------------------------------------------------
      call mm_get_solute_quant_ind(nqm-nlink,int_mb(i_ai))
      if(nlink.gt.0) then
        call mm_get_solute_link_ind(nlink,int_mb(i_ai+nqm-nlink))
      end if

      if(.not.rtdb_get(irtdb,'charge',mt_dbl,1,qcharge)) then
        qcharge = 0.0d0    
      end if

       if (.not.rtdb_get(irtdb,'qmmm:linkatm',mt_int,1,link_atom))
     + call errquit('qmmm_data_init: qmmm:linkatm',link_atom,
     &       RTDB_ERR)


      end

      subroutine qmmm_data_alloc(nat)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
      integer nat

      nqm = nat
c
c     indexing array
c     -------------
      if(.not.ma_alloc_get(MT_INT, nat, 'qmmm index array',
     &      h_ai, i_ai) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call ifill(nqm,0,int_mb(i_ai),1)
c
c     coordinate array
c     -------------
      if(.not.ma_alloc_get(MT_DBL, 3*nat, 'qmmm coord array',
     &      h_c, i_c) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(3*nat,0,dbl_mb(i_c),1)

c
c     charge array
c     -------------
      if(.not.ma_alloc_get(MT_DBL, nat, 'qmmm charge array',
     &      h_q, i_q) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(nat,0,dbl_mb(i_q),1)

c
c     mass array
c     -------------
      if(.not.ma_alloc_get(MT_DBL, nat, 'qmmm mass array',
     &      h_m, i_m) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call dfill(nat,0,dbl_mb(i_m),1)

c
c     tag array
c     -------------
      if(.not.ma_alloc_get(MT_BYTE, 16*nat, 'qmmm tag array',
     &      h_t, i_t) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)

c
c     atomic number array
c     -------------------
      if(.not.ma_alloc_get(MT_INT, nat, 'qmmm inum array',
     &      h_num, i_num) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call ifill(nqm,0,int_mb(i_num),1)
    
       
      end

      subroutine qmmm_data_release()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"

      if(.not.ma_free_heap(h_num))
     & call errquit('qmmm c: Failed to deallocate heap',nqm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_t))
     & call errquit('qmmm c: Failed to deallocate heap',nqm,
     &       MA_ERR)
 
      if(.not.ma_free_heap(h_m))
     & call errquit('qmmm c: Failed to deallocate heap',nqm,
     &       MA_ERR)
 
      if(.not.ma_free_heap(h_q))
     & call errquit('qmmm c: Failed to deallocate heap',nqm,
     &       MA_ERR)
     
      if(.not.ma_free_heap(h_c))
     & call errquit('qmmm c: Failed to deallocate heap',nqm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_ai))
     & call errquit('qmmm ai: Failed to deallocate heap',nqm,
     &       MA_ERR)

       return
      end

      function qmmm_get_nqm()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nqm

      qmmm_get_nqm = nqm

      end 

      function qmmm_get_nlink()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nlink

      qmmm_get_nlink = nlink

      end 

      function qmmm_get_i_ai()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_ai

      qmmm_get_i_ai = i_ai

      end 

      function qmmm_get_i_q()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_q

      qmmm_get_i_q = i_q

      end 

      function qmmm_get_i_c()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_c

      qmmm_get_i_c = i_c

      end 

      function qmmm_get_i_t()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_t

      qmmm_get_i_t = i_t

      end 

      function qmmm_get_i_num()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_num

      qmmm_get_i_num = i_num

      end 

      function qmmm_get_i_m()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_m

      qmmm_get_i_m = i_m

      end 

      function qmmm_get_i_g()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_g

      qmmm_get_i_g = i_g

      end 

      function qmmm_get_i_gmm()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_gmm

      qmmm_get_i_gmm = i_gmm

      end 

      function qmmm_get_qcharge()
      implicit none
#include "qmmm_data.fh"
      double precision qmmm_get_qcharge

      qmmm_get_qcharge = qcharge

      end 

      function qmmm_get_link_atom()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_link_atom

      qmmm_get_link_atom = link_atom

      end 

      function qmmm_h_link()
      implicit none
#include "qmmm_data.fh"
#include "qmmm_params.fh"

      logical qmmm_h_link

      qmmm_h_link=link_atom.eq.qmmm_hydrogen

      end 

