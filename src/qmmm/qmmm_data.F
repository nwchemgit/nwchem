      subroutine qmmm_data_init(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"

      integer irtdb
c
      character*32 tag
      character*32 pname
      character*32 optimization
      character*32 link_ecp
      real*8 eatoms,qatoms
      logical fixqm
      integer linkatm,bq_exclude
      character*30 operation
      integer i_list,h_list
      integer i_nlist,h_nlist
      integer i
      logical qmmm_h_link
      external qmmm_h_link


      pname = "qmmm_data_init: "

      call mm_get_tot_nqm(nqm)
      if(ga_nodeid().eq.0) 
     >  write(*,*) "number of quantum atoms",nqm
c
c     allocate data arrays
c     --------------------
      call qmmm_data_alloc()
c 
c     establish persistent mapping between mm and qm atoms
c     first go pure quantum atoms and then links
c     ----------------------------------------------------
      call mm_get_solute_quant_ind(nqm,int_mb(i_ai))
      if(.not.rtdb_get(irtdb,'charge',mt_dbl,1,qcharge)) then
        qcharge = 0.0d0    
      end if

c       
c     get mm constraints if any
c     ------------------------
      call mm_get_solute_cons(nqm,
     >                       int_mb(i_ai),
     >                       log_mb(i_act))


      end

      subroutine qmmm_get_nqm_act(nqm_act)
      implicit none
#include "mafdecls.fh"
#include "qmmm_data.fh"
c
      integer nqm_act
c
      integer i
      
      nqm_act = 0
      do i=1,nqm
        if(log_mb(i_act+i-1)) nqm_act = nqm_act + 1
      end do

      end

      subroutine qmmm_data_rdinp(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_link_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"

      integer irtdb
c
      character*32 tag
      character*32 pname
      character*32 optimization
c      character*32 interface
      character*32 link_ecp
      double precision  eatoms,qatoms
      logical fixqm
      integer linkatm,bq_exclude
      character*30 operation
      logical ignore
      logical flag

      pname = "qmmm_data_init: "

      master = ga_nodeid().eq.0
      fixlinks = .false.

      eatoms = 0.0d0
      linkatm=qmmm_hydrogen
      bq_exclude=qmmm_no_atoms
      optimization = "qm"
      link_ecp = "auto"
c
c     initialize rtdb parameters
      if(master)
     > call util_print_centered(6,
     >    "QM/MM Interface Parameters",32,.true.)

      tag = "qmmm:uqmatm"
      if (.not.rtdb_get(irtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) then
        if (.not.rtdb_put(irtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) 
     >      call errquit(pname//'failed to store eatoms',0,
     >       RTDB_ERR)
      end if
      if(master) then 
       if(eatoms.eq.0) then
         write(6,19) "reference energy",tag,eatoms,
     >     "  <--- Warning zero value is not advisable !!!"
       else
         write(6,20) "reference energy",tag,eatoms
       end if
      end if

c
c     default value for bqzone from md
      if(.not.rtdb_get(irtdb,'md:rqmmm',mt_dbl,1,bqzone))
     + call errquit(pname//'failed to get bqzone from mm',0)
c
      tag = "qmmm:bqzone"
      if (rtdb_get(irtdb,tag,mt_dbl,1,bqzone)) then
        if (.not.rtdb_put(irtdb,'md:rqmmm',mt_dbl,1,bqzone)) 
     >      call errquit(pname//'failed to store '//tag,0,
     >       RTDB_ERR)
      end if
      if(master) 
     > write(6,20) "bqzone radius",tag,bqzone


      tag = "qmmm:fixlinks"
      if (.not.rtdb_get(irtdb,tag,mt_log,1,fixlinks)) then
        if (.not.rtdb_put(irtdb,tag,mt_log,1,fixlinks)) 
     >      call errquit(pname//'failed to store fixlinks',0,
     >       RTDB_ERR)

       end if
      if(master)
     > write(6,21) "fixed link bonds",tag,fixlinks


      tag ="qmmm:bq_exclude" 
      if (.not.rtdb_get(irtdb,tag,mt_int,1,bq_exclude)) then
        if (.not.rtdb_put(irtdb,tag,mt_int,1,bq_exclude)) 
     >      call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)

       end if
      if(master) then
        if(bq_exclude.eq.qmmm_linkbond_H) then
          write(6,22) "excluded bq's",tag,
     >                 "hydrogens bonded to link atoms"
        else if(bq_exclude.eq.qmmm_linkbond) then
          write(6,22) "excluded bq's",tag,"all atoms bonded to links"
        else if(bq_exclude.eq.qmmm_all_atoms) then
          write(6,22) "excluded bq's",tag,"all"
        else if(bq_exclude.eq.qmmm_no_atoms) then
          write(6,22) "excluded bq's",tag,"none"
        else
          call errquit(pname//'invalid bq_exclude',0,RTDB_ERR)
        end if
      end if

      tag ="qmmm:linkatm" 
      if (.not.rtdb_get(irtdb,tag,mt_int,1,linkatm)) then
        if (.not.rtdb_put(irtdb,tag,mt_int,1,linkatm)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

       end if
       if(master) then
         if(linkatm.eq.qmmm_hydrogen) then
           write(6,22) "link atom type",tag,"hydrogens"
         else if(linkatm.eq.qmmm_halogen) then
           write(6,22) "link atom type",tag,"halogens"
         else
           call errquit(pname//'invalid link atom type',0,RTDB_ERR)
         end if
       end if
     
      tag ="qmmm:optimization"
      if (.not.rtdb_cget(irtdb,tag,1,optimization)) then
        if (.not.rtdb_cput(irtdb,tag,1,optimization)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

       else if (optimization .ne. "all" .and.
     >          optimization .ne. "mm" .and.
     >          optimization .ne. "qm" ) then
            call errquit(pname//'unknown value'//tag,0,
     >       RTDB_ERR)

       end if
      if (master) 
     >  write(6,22) "optimization region",tag,optimization

      tag ="qmmm:link_ecp"
      if (.not.rtdb_cget(irtdb,tag,1,link_ecp)) then
       if (.not.rtdb_cput(irtdb,tag,1,link_ecp)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

       else if (link_ecp .ne. "auto" .and.
     >          link_ecp .ne. "user" ) then
            call errquit(pname//'unknown value'//tag,0,
     >       RTDB_ERR)

       end if
      if (master) 
     >  write(6,22) "ecp on link atoms",tag,link_ecp

      fixqm = .false.
      qmmm_interface='qm'

      tag='task:operation'
      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'
      if (master) 
     >  write(6,22) "operation",tag,operation
c
      if (inp_compare(.false., 'energy', operation)) then
        if (optimization.eq.'all') then
           fixqm = .false.
           qmmm_interface = "mm"
           bq_dynamical=.true.
        else if (optimization.eq.'mm') then
           fixqm = .true.
           qmmm_interface = "mm"
           bq_dynamical=.true.
        else if (optimization.eq.'qm') then
           fixqm = .false.
           qmmm_interface = "qm"
        endif
      else if (inp_compare(.false., 'optimize', operation)) then
        if (optimization.eq.'all') then
           fixqm = .false.
           qmmm_interface = "mm"
           bq_dynamical=.true.
        else if (optimization.eq.'mm') then
           fixqm = .true.
           qmmm_interface = "mm"
           bq_dynamical=.true.
        else if (optimization.eq.'qm') then
           fixqm = .false.
           qmmm_interface = "qm"
           bq_dynamical=.false.
        endif
      else if (inp_compare(.false., 'saddle', operation)) then
        fixqm = .false.
        qmmm_interface = "qm"
      else if (inp_compare(.false.,'dynamics',operation)) then
        qmmm_interface = "mm"
      end if

      tag ="qmmm:fixqm" 
      if (.not.rtdb_put(irtdb,tag,mt_log,1,fixqm)) 
     >      call errquit(pname//'failed to store fixqm',0,
     >       RTDB_ERR)
      if (master) 
     >  write(6,21) "fixed qm region",tag,fixqm


      tag ="qmmm:interface"
      ignore = rtdb_cget(irtdb,tag,1,qmmm_interface)
c      if (.not.rtdb_cget(irtdb,tag,1,interface)) then
c         if (.not.rtdb_cput(irtdb,tag,1,interface)) 
c     >        call errquit(pname//'failed to store'//tag,0,
c     >         RTDB_ERR)
c      end if
      if (master) 
     >  write(6,22) "interface api",tag,qmmm_interface

      tag = "qmmm:bq_dynamical"
      if (.not.rtdb_get(irtdb,tag,mt_log,1,bq_dynamical)) then
        if (.not.rtdb_put(irtdb,tag,mt_log,1,bq_dynamical)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

      end if
      if(master)
     > write(6,21) "dynamical bq update",tag,bq_dynamical

      if (master) 
     >  write(6,23) 
       

19    FORMAT(1X,A,T24,A,T46,F12.6,A)
20    FORMAT(1X,A,T24,A,T46,F12.6)
21    FORMAT(1X,A,T24,A,T46,L3)
22    FORMAT(1X,A,T24,A,T46,A)
23    FORMAT(1X,54("-"),//)

       if (.not.rtdb_get(irtdb,'qmmm:linkatm',mt_int,1,link_atom))
     + call errquit('qmmm_data_init: qmmm:linkatm',link_atom,
     &       RTDB_ERR)

      nqm   = -1 

      end

      subroutine qmmm_data_alloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
c
      integer i
c
c     indexing array
c     -------------
      if(.not.ma_alloc_get(MT_INT, nqm, 'qmmm index array',
     &      h_ai, i_ai) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call ifill(nqm,-1,int_mb(i_ai),1)
c   
c     index list of active atoms  
c     --------------------------
      if(.not.ma_alloc_get(MT_LOG, nqm, 'qmmm active list',
     &      h_act, i_act) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      
       do i=1,nqm
         log_mb(i_act+i-1) = .true.
       end do
c
c
cc     coordinqme array
cc     -------------
c      if(.not.ma_alloc_get(MT_DBL, 3*nqm, 'qmmm coord array',
c     &      h_c, i_c) ) call errquit(
c     &      'qmmm_data_alloc: unable to allocate heap space',
c     &      nqm, MA_ERR)
c      call dfill(3*nqm,0,dbl_mb(i_c),1)
c
cc
cc     charge array
cc     -------------
c      if(.not.ma_alloc_get(MT_DBL, nqm, 'qmmm charge array',
c     &      h_q, i_q) ) call errquit(
c     &      'qmmm_data_alloc: unable to allocate heap space',
c     &      nqm, MA_ERR)
c      call dfill(nqm,0,dbl_mb(i_q),1)
c
cc
cc     mass array
cc     -------------
c      if(.not.ma_alloc_get(MT_DBL, nqm, 'qmmm mass array',
c     &      h_m, i_m) ) call errquit(
c     &      'qmmm_data_alloc: unable to allocate heap space',
c     &      nqm, MA_ERR)
c      call dfill(nqm,0,dbl_mb(i_m),1)
c
cc
cc     tag array
cc     -------------
c      if(.not.ma_alloc_get(MT_BYTE, 16*nqm, 'qmmm tag array',
c     &      h_t, i_t) ) call errquit(
c     &      'qmmm_data_alloc: unable to allocate heap space',
c     &      nqm, MA_ERR)
c
cc
cc     atomic number array
cc     -------------------
c      if(.not.ma_alloc_get(MT_INT, nqm, 'qmmm inum array',
c     &      h_num, i_num) ) call errquit(
c     &      'qmmm_data_alloc: unable to allocate heap space',
c     &      nqm, MA_ERR)
c      call ifill(nqm,0,int_mb(i_num),1)
c
cc     charge array
cc     -------------
c      if(.not.ma_alloc_get(MT_DBL, nqm, 'qmmm esp charge array',
c     &      h_qs, i_qs) ) call errquit(
c     &      'qmmm_data_alloc: unable to allocate heap space',
c     &      nqm, MA_ERR)
c      call dfill(nqm,0,dbl_mb(i_qs),1)

    
      end

      subroutine qmmm_data_release()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"

c      if(.not.ma_free_heap(h_num))
c     & call errquit('qmmm c: Failed to deallocate heap',nqm,
c     &       MA_ERR)
c
c      if(.not.ma_free_heap(h_t))
c     & call errquit('qmmm c: Failed to deallocate heap',nqm,
c     &       MA_ERR)
c 
c      if(.not.ma_free_heap(h_m))
c     & call errquit('qmmm c: Failed to deallocate heap',nqm,
c     &       MA_ERR)
c 
c      if(.not.ma_free_heap(h_q))
c     & call errquit('qmmm c: Failed to deallocate heap',nqm,
c     &       MA_ERR)
c     
c      if(.not.ma_free_heap(h_c))
c     & call errquit('qmmm c: Failed to deallocate heap',nqm,
c     &       MA_ERR)

      if(.not.ma_free_heap(h_act))
     & call errquit('qmmm ai: Failed to deallocate heap',nqm,
     &       MA_ERR)

      if(.not.ma_free_heap(h_ai))
     & call errquit('qmmm ai: Failed to deallocate heap',nqm,
     &       MA_ERR)

c
c      if(.not.ma_free_heap(h_qs))
c     & call errquit('qmmm qs: Failed to deallocate heap',nqm,
c     &       MA_ERR)


       return
      end

      function qmmm_get_nqm()
      implicit none
#include "qmmm_data.fh"

      integer qmmm_get_nqm

      qmmm_get_nqm = nqm

      end 

      function qmmm_get_i_ai()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_ai

      qmmm_get_i_ai = i_ai

      end 

      function qmmm_get_i_act()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_act

      qmmm_get_i_act = i_act

      end 

      function qmmm_get_i_q()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_q

      qmmm_get_i_q = i_q

      end 

      function qmmm_get_i_qs()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_qs

      qmmm_get_i_qs = i_qs

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

      function qmmm_master()
      implicit none
#include "qmmm_data.fh"
#include "qmmm_params.fh"

      logical qmmm_master

      qmmm_master = master

      end 

      function qmmm_bq_dynamical()
      implicit none
#include "qmmm_data.fh"
#include "qmmm_params.fh"

      logical qmmm_bq_dynamical

      qmmm_bq_dynamical = bq_dynamical

      end 

      function qmmm_fixlinks()
      implicit none
#include "qmmm_data.fh"
#include "qmmm_params.fh"

      logical qmmm_fixlinks

      qmmm_fixlinks = fixlinks

      end 

      function qmmm_print_default(name)
      implicit none
#include "qmmm_data.fh"
#include "util.fh"
#include "global.fh"

      logical qmmm_print_default
      character*(*) name
c 
      logical oprint

      oprint = util_print(name, print_default)

      qmmm_print_default = oprint .and. (ga_nodeid().eq.0)

      end 

      function qmmm_print_debug()
      implicit none
#include "qmmm_data.fh"
#include "util.fh"
#include "global.fh"

      logical qmmm_print_debug
c 
      logical oprint

      oprint = util_print("debug_info", print_debug)

      qmmm_print_debug = oprint .and. (ga_nodeid().eq.0)

      end 

      function qmmm_get_bqzone()
      implicit none
#include "qmmm_data.fh"
      double precision qmmm_get_bqzone

      qmmm_get_bqzone = bqzone

      end 

      function qmmm_get_interface()
      implicit none
#include "qmmm_data.fh"
      character*(*) qmmm_get_interface

      qmmm_get_interface = qmmm_interface      

      end 

