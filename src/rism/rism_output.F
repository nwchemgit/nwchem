c
c $Id$
c
c
      subroutine rism_print_header()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"

      character*48 string

      if(ga_nodeid().eq.0) then
      
        string = "  NWChem 1D-RISM Module"
        call banner(luout, 
     >      string,
     >      "*", .true. , .false., .true.)
        string = " "
        call banner(luout, 
     >      string,
     >      "*", .false., .false., .true.)

        string = "  Marat Valiev "// 
     >           "and Gennady N. Chuev 2011"
        call banner(luout, 
     >      string,
     >      "*", .false., .true., .true.)
        
        call util_flush(luout)
      end if

      end

       subroutine rism_print_params(rtdb)
       implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "rism_params.fh"
      integer rtdb
      integer nu,nv,nvv,ngr
      integer icl,icr
      integer i_rgrid,i_kgrid
      integer i_tv,i_den,i_isv,i_mv
      integer i_xv,i_yv,i_zv
      integer i_ims
      logical result
      character*32 sname,dname,pname
      double precision t,tau,lambd,tol,solperm
      integer dd
      integer i_sgvv,i_epsvv,i_qvv
      integer i_sigu,i_epsiu,i_qqu,i_wu
      character*72 filename,tag
      character*16 aname
      logical okspace
c       
      pname = "rism_wrapper"
c      
      if(ga_nodeid().ne.0) return
c      
      call util_print_centered(luout,
     A     "1D-RISM Calculation Parameters",32,.true.)
c
c     grid data
c     -------------- 
      sname = "grid"
      dname = "npoints"
      call db_data_get_int(sname,dname,1,ngr,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,19) "Number of grid points",ngr
      tag = "okspace"
      call db_data_get_log(sname,tag,1,okspace,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
c     parameters
c     -------------- 
      sname = "parameters"
      dname = "closure"
      call db_data_get_int(sname,dname,1,icl,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      if(icl.eq.t_hnc) then
        aname = "hnc"
      else if (icl.eq.t_kh) then
        aname = "kh"
      else
        aname = "uknown"
      end if
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,22) "Closure type",aname

      dname = "vdw_rule"
      call db_data_get_int(sname,dname,1,icr,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      if(icr.eq.t_arithmetic) then
        aname = "arithmetic"
      else if (icl.eq.t_geometric) then
        aname = "geometric"
      else
        aname = "uknown?"
      end if
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,22) "VDW combination rule",aname

      dname = "solvent_permittivity"
      call db_data_get_dbl(sname,dname,1,solperm,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,20) "Solvent permittivity",solperm

      dname = "tau"
      call db_data_get_dbl(sname,dname,1,tau,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,20) "Coulomb decay length",solperm


      dname = "tolerance"
      call db_data_get_dbl(sname,dname,1,tol,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,20) "Convergence tolerance",tol

      dname = "mixing"
      call db_data_get_dbl(sname,dname,1,lambd,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,20) "Mixing parameter",lambd

      dname = "temperature"
      call db_data_get_dbl(sname,dname,1,t,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,20) "Temperature",t

      dname = "diis"
      call db_data_get_int(sname,dname,1,dd,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      write(luout,19) "Diis iterations",dd
c      
c     get filename for solvent g(r)
c     -----------------------------
      tag = "rism:solvent:rdf"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) 
     >  call errquit("cannot get "//tag,0,0)
      tag = sname(1:inp_strlen(sname))//":"//
     >      dname(1:inp_strlen(dname))
      if(okspace) then
        aname = "Solvent rdf file in kspace"
      else
        aname = "Solvent rdf file "
      end if
      write(luout,22) aname,filename
c      
c     get solute filename
c     --------------------
      tag = "rism:solute:configuration"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) 
     >  call errquit("cannot get "//tag,0,0)
      write(luout,22) "Solute configuration file",filename
c

c
c      tag ="qmmm:bq_exclude" 
c      if(bq_exclude.eq.qmmm_hydrogen) then
c        write(luout,22) "excluded bq's",tag,
c     V        "hydrogens bonded to link atoms"
c      else if(bq_exclude.eq.qmmm_all_atoms) then
c        write(luout,22) "excluded bq's",tag,"all"
c      else if(bq_exclude.eq.qmmm_no_atoms) then
c        write(luout,22) "excluded bq's",tag,"none"
c      else
c        call errquit(pname//'invalid bq_exclude',0,RTDB_ERR)
c      end if
c      
c      tag ="qmmm:linkatm" 
c      if(linkatm.eq.qmmm_hydrogen) then
c        write(luout,22) "link atom type",tag,"hydrogens"
c      else if(linkatm.eq.qmmm_halogen) then
c        write(luout,22) "link atom type",tag,"halogens"
c      else
c        call errquit(pname//'invalid link atom type',0,RTDB_ERR)
c      end if
c
c      tag ="qmmm:optimization"
c      write(luout,22) "optimization region",tag,optimization
c

19    FORMAT(A,T33,I4)
23    FORMAT(A,T33,F12.6)
20    FORMAT(A,T33,G10.4)
21    FORMAT(A,T33,L3)
22    FORMAT(A,T33,A)


       return
      end 

      subroutine rism_print_solute_configuration(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
      integer rtdb
c
      character*(72) filename
      integer i
      integer fn
      character*72 tag
      character*72 buffer
      character*1 b1
      integer master
      integer nu
      character*32 sname
      logical result
c 
      integer  i_epsiu
      integer  i_qqu
      integer  i_xu
      integer  i_yu
      integer  i_zu
      integer  i_sigu
      integer  i_isu
      integer  i_mu
      integer  i_tu
      integer  i_tt
      integer  i_wu
      integer  i_kgrid,ngrid
      character*80 pname 
c
      pname = "rism_print_solute_configuration"
#ifdef RISM_DEBUG      
      call rism_message("in "//pname)
#endif
c
      sname = "solute"
      tag = "natoms"
      call db_data_get_int(sname,tag,1,nu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"tu",i_tu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"type",i_tt,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"xu",i_xu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"yu",i_yu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"zu",i_zu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"sigma",i_sigu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"epsilon",i_epsiu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
      call db_data_get_index(sname,"charge",i_qqu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               ,0,0)
 
      call rism_print_solute_configuration0(
     >                     nu,
     >                     byte_mb(i_tu),
     >                     byte_mb(i_tt),
     >                     dbl_mb(i_xu),
     >                     dbl_mb(i_yu),
     >                     dbl_mb(i_zu),
     >                     dbl_mb(i_sigu),
     >                     dbl_mb(i_epsiu),
     >                     dbl_mb(i_qqu)
     >                     )
#ifdef RISM_DEBUG      
      call rism_message("out "//pname)
#endif

      return
      end

      subroutine rism_print_solute_configuration0(
     >                         n,tn,tt,x,y,z,sigma,eps,q)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
#include "stdio.fh"
      integer n
      character(db_string_size) tn(n),tt(n)
      double precision x(n)
      double precision y(n)
      double precision z(n)
      double precision eps(n)
      double precision sigma(n)
      double precision q(n)
      character*80 pname 
      integer i
c
      pname = "rism_print_solute_configuration0"
#ifdef RISM_DEBUG      
      call rism_message("in "//pname)
#endif
c
c      
      call util_print_centered(luout,
     A     "Solute configuration",32,.true.)
c
      write(*,9)
      do i=1,n
        write(*,10) tn(i),tt(i),x(i),y(i),z(i),sigma(i),eps(i),q(i)
      end do
9     FORMAT("name",T6,"type",T17,"x",11X,"y",11X,"z",11X,"sigma",
     >        T65,"eps",T77,"charge")
10    FORMAT(2A5,T13,6F12.6) 
#ifdef RISM_DEBUG      
      call rism_message("out "//pname)
#endif

      return
      end
