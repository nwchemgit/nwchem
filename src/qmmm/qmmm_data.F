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
      character*32 interface
      character*32 link_ecp
      real*8 eatoms,qatoms
      logical fixqm
      integer linkatm,bq_exclude
      character*30 operation

      pname = "qmmm_data_init: "

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


20    FORMAT(1X,A,T24,A,T46,F12.6)
21    FORMAT(1X,A,T24,A,T46,L3)
22    FORMAT(1X,A,T24,A,T46,A)
23    FORMAT(1X,54("-"),//)


      end

      subroutine qmmm_data_rdinp(irtdb)
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
      character*32 interface
      character*32 link_ecp
      double precision  eatoms,qatoms
      logical fixqm
      integer linkatm,bq_exclude
      character*30 operation

      pname = "qmmm_data_init: "

      master = ga_nodeid().eq.0
      fixlinks = .false.

      bqzone = 6.0
      eatoms = 0.0d0
      linkatm=qmmm_halogen
      bq_exclude=qmmm_hydrogen
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

      tag = "qmmm:bqzone"
      if (.not.rtdb_get(irtdb,tag,mt_dbl,1,bqzone)) then
        if (.not.rtdb_put(irtdb,tag,mt_dbl,1,bqzone)) 
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
        if(bq_exclude.eq.qmmm_hydrogen) then
          write(6,22) "excluded bq's",tag,
     >                 "hydrogens bonded to link atoms"
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
      interface='qm'

      tag='task:operation'
      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'
      if (master) 
     >  write(6,22) "operation",tag,operation
c
      if (inp_compare(.false., 'optimize', operation)) then
        if (optimization.eq.'all') then
           fixqm = .false.
           interface = "mm"
        else if (optimization.eq.'mm') then
           fixqm = .true.
           interface = "mm"
        else if (optimization.eq.'qm') then
           fixqm = .false.
           interface = "qm"
        endif
      else if (inp_compare(.false., 'saddle', operation)) then
        fixqm = .false.
        interface = "qm"
      else if (inp_compare(.false.,'dynamics',operation)) then
        interface = "mm"
      end if

      tag ="qmmm:fixqm" 
      if (.not.rtdb_put(irtdb,tag,mt_log,1,fixqm)) 
     >      call errquit(pname//'failed to store fixqm',0,
     >       RTDB_ERR)
      if (master) 
     >  write(6,21) "fixed qm region",tag,fixqm


      tag ="qmmm:interface"
      if (.not.rtdb_cput(irtdb,tag,1,interface)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
      if (master) 
     >  write(6,22) "interface api",tag,interface

      if (master) 
     >  write(6,23) 
       

19    FORMAT(1X,A,T24,A,T46,F12.6,A)
20    FORMAT(1X,A,T24,A,T46,F12.6)
21    FORMAT(1X,A,T24,A,T46,L3)
22    FORMAT(1X,A,T24,A,T46,A)
23    FORMAT(1X,54("-"),//)


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

      function qmmm_halogen_link()
      implicit none
#include "qmmm_data.fh"
#include "qmmm_params.fh"

      logical qmmm_halogen_link

      qmmm_halogen_link=link_atom.eq.qmmm_halogen

      end 

      function qmmm_master()
      implicit none
#include "qmmm_data.fh"
#include "qmmm_params.fh"

      logical qmmm_master

      qmmm_master = master

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

