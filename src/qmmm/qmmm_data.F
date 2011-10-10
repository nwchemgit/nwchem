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
      character*32 pname
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
#include "util.fh"

      integer irtdb
c
      character*84 tag
      character*32 pname
      character*32 optimization
      character*32 link_ecp
      double precision  eatoms
      integer linkatm,bq_exclude
      character*30 operation
      character*30 region(3)
      integer nregion
      logical ignore
      integer i

      pname = "qmmm_data_init: "

      master = ga_nodeid().eq.0
c
      oprint_debug = util_print("debug_info", print_debug)
      oprint_debug = oprint_debug .and. (ga_nodeid().eq.0)
c
      oprint_low = util_print("low", print_low)
c
      eatoms = 0.0d0
      linkatm=qmmm_hydrogen
      bq_exclude=qmmm_no_atoms
      link_ecp = "auto"
c
c     initialize rtdb parameters
      if(master)
     > call util_print_centered(6,
     >    "QM/MM Interface Parameters",32,.true.)
c
c     get the operation
c     -----------------
      tag='task:operation'
      if (.not. rtdb_cget(irtdb, 'task:operation', 1, operation))
     $     operation = 'energy'
      if (master) 
     >  write(6,22) "operation",tag,operation
c
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
c     convert from nm to angstrom
      bqzone = bqzone*cnm2ang
c
      tag = "qmmm:bqzone"
      if (rtdb_get(irtdb,tag,mt_dbl,1,bqzone)) then
        if (.not.rtdb_put(irtdb,'md:rqmmm',mt_dbl,1,bqzone/cnm2ang)) 
     >      call errquit(pname//'failed to store '//tag,0,
     >       RTDB_ERR)
      end if
      if(master) 
     > write(6,20) "bqzone radius",tag,bqzone


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
     
      optimization = "bfgs"
      tag ="qmmm:optimization"
      if (.not.rtdb_cget(irtdb,tag,1,optimization)) then
        if (.not.rtdb_cput(irtdb,tag,1,optimization)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

      end if
      if (master) 
     >  write(6,22) "optimization method",tag,optimization

c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(irtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) then
c        if(operation.ne."energy" .and.
c     >     operation.ne."gradient" .and.
c     >     operation.ne."property" .and.
c     >     operation.ne."fep" .and.
c     >     operation.ne."neb" )
c     >     call errquit(pname//'cannot find '//tag,0,
c     >       RTDB_ERR)
        nregion = 0
      end if      
      if(nregion.ne.0) then
        if (.not.rtdb_cget(irtdb,tag,nregion,region)) 
     >      call errquit(pname//"failed to get"//tag,0,RTDB_ERR)
      end if
      if (.not.rtdb_put(irtdb,tag(1:inp_strlen(tag))//"_n",
     >                  mt_int,1,nregion)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >                   RTDB_ERR)

      if (master) then
        tag = " "
        do i=1,nregion
          write(tag,'("qmmm:region ",I1)') i
          write(6,22) tag,region(i)
        end do
      end if

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

      qmmm_interface='qm'
      tag ="qmmm:interface"
      ignore = rtdb_cget(irtdb,tag,1,qmmm_interface)
      if (.not.rtdb_cget(irtdb,tag,1,qmmm_interface)) then
         qmmm_interface = "qm"
      end if
      if (master) 
     >  write(6,22) "interface api",tag,qmmm_interface

      tag = "qmmm:bq_dynamical"
      if (.not.rtdb_get(irtdb,tag,mt_log,1,bq_dynamical)) 
     >  bq_dynamical = .false.

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
c
c     indexing array
c     -------------
      if(.not.ma_alloc_get(MT_INT, nqm, 'qmmm index array',
     &      h_ai, i_ai) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nqm, MA_ERR)
      call ifill(nqm,-1,int_mb(i_ai),1)
     
      end

      subroutine qmmm_data_release()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"

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

      function qmmm_get_i_ai()
      implicit none
#include "qmmm_data.fh"
      integer qmmm_get_i_ai

      qmmm_get_i_ai = i_ai

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

c      oprint = util_print("debug_info", print_debug)
c      qmmm_print_debug = oprint .and. (ga_nodeid().eq.0)

       qmmm_print_debug = oprint_debug
      end 

      function qmmm_print_low()
      implicit none
#include "qmmm_data.fh"
#include "util.fh"
#include "global.fh"

      logical qmmm_print_low
c 

       qmmm_print_low = oprint_low
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
      character*32 qmmm_get_interface

      qmmm_get_interface = qmmm_interface      

      end 

      subroutine qmmm_set_interface(a_interface)
      implicit none
#include "qmmm_data.fh"
      character*(*) a_interface

      qmmm_interface = a_interface      

      end 

c $Id$
