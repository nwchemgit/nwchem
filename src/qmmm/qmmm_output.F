c
c $Id: qmmm_output.F,v 1.1 2004-10-13 23:38:43 marat Exp $
c
c
      subroutine qmmm_print_header(irtdb)
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
      character*48 string

      
      pname = "qmmm_print_header:"
      string = "     NWChem QM/MM Interface Module"
      call banner(6, 
     >    string,
     >    "*", .true. , .false., .true.)
      string = " "
      call banner(6, 
     >    string,
     >    "*", .false., .false., .true.)

      string = "     Developed Marat Valiev 2004"
      call banner(6, 
     >    string,
     >    "*", .false., .true., .true.)

      end

      subroutine qmmm_print_params(irtdb)
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

      
      pname = "qmmm_print_header:"
      call util_print_centered(6,"QM/MM Interface Parameters",32,.true.)

      tag = "qmmm:uqmatm"
      if (.not.rtdb_get(irtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) then
           call errquit(pname//'failed to get '//tag,0,
     >       RTDB_ERR)
      end if
      write(6,20) "reference energy",tag,eatoms

      tag = "qmmm:fixlinks"
      if (.not.rtdb_get(irtdb,tag,mt_log,1,fixlinks)) then
           call errquit(pname//'failed to get '//tag,0,
     >       RTDB_ERR)
      end if
      write(6,21) "fixed link bonds",tag,fixlinks


      tag ="qmmm:fixqm" 
      if (.not.rtdb_get(irtdb,tag,mt_log,1,fixqm)) then
           call errquit(pname//'failed to store fixqm',0,
     >       RTDB_ERR)
       end if
      write(6,21) "fixed qm region",tag,fixqm

      tag ="qmmm:bq_exclude" 
      if (.not.rtdb_get(irtdb,tag,mt_int,1,bq_exclude)) then
           call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)
      end if
      if(bq_exclude.eq.qmmm_hydrogen) then
        write(6,22) "excluded bq's",tag,"hydrogens bonded to link atoms"
      else if(bq_exclude.eq.qmmm_all_atoms) then
        write(6,22) "excluded bq's",tag,"all"
      else if(bq_exclude.eq.qmmm_no_atoms) then
        write(6,22) "excluded bq's",tag,"none"
      else
        call errquit(pname//'invalid bq_exclude',0,RTDB_ERR)
      end if
      
      tag ="qmmm:linkatm" 
      if (.not.rtdb_get(irtdb,tag,mt_int,1,linkatm)) then
           call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
       end if
      if(linkatm.eq.qmmm_hydrogen) then
        write(6,22) "link atom type",tag,"hydrogens"
      else if(linkatm.eq.qmmm_halogen) then
        write(6,22) "link atom type",tag,"halogens"
      else
        call errquit(pname//'invalid link atom type',0,RTDB_ERR)
      end if

      tag ="qmmm:optimize"
      if (.not.rtdb_cget(irtdb,tag,1,optimize)) then
           call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
      end if
      write(6,22) "optimization region",tag,optimize

20    FORMAT(A,T24,A,T46,F12.6)
21    FORMAT(A,T24,A,T46,L3)
22    FORMAT(A,T24,A,T46,A)


      end

