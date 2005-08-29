c
c $Id: qmmm_output.F,v 1.18 2005-08-29 01:19:57 marat Exp $
c
c
      subroutine qmmm_print_header()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "global.fh"
#include "stdio.fh"

      character*48 string

      if(ga_nodeid().eq.0) then
      
        string = "     NWChem QM/MM Interface Module"
        call banner(luout, 
     >      string,
     >      "*", .true. , .false., .true.)
        string = " "
        call banner(luout, 
     >      string,
     >      "*", .false., .false., .true.)

        string = "     Developed Marat Valiev 2004"
        call banner(luout, 
     >      string,
     >      "*", .false., .true., .true.)
        
        call util_flush(luout)
      end if
      

      end

      subroutine qmmm_print_params(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "qmmm.fh"
#include "stdio.fh"

      integer irtdb
c
      character*32 tag
      character*32 pname
      character*32 optimization
      real*8 eatoms
      integer linkatm,bq_exclude

      
      pname = "qmmm_print_header:"

      tag = "qmmm:uqmatm"
      if (.not.rtdb_get(irtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) then
           call errquit(pname//'failed to get '//tag,0,
     >       RTDB_ERR)
      end if

      tag ="qmmm:bq_exclude" 
      if (.not.rtdb_get(irtdb,tag,mt_int,1,bq_exclude)) then
           call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)
      end if
c
      tag ="qmmm:linkatm" 
      if (.not.rtdb_get(irtdb,tag,mt_int,1,linkatm)) then
           call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
       end if
c
      tag ="qmmm:optimization"
      if (.not.rtdb_cget(irtdb,tag,1,optimization)) then
           call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
      end if

      if(.not.qmmm_master()) return

      call util_print_centered(luout,
     A     "QM/MM Interface Parameters",32,.true.)

      tag = "qmmm:uqmatm"
      write(luout,20) "reference energy",tag,eatoms

      tag ="qmmm:bq_exclude" 
      if(bq_exclude.eq.qmmm_hydrogen) then
        write(luout,22) "excluded bq's",tag,
     V        "hydrogens bonded to link atoms"
      else if(bq_exclude.eq.qmmm_all_atoms) then
        write(luout,22) "excluded bq's",tag,"all"
      else if(bq_exclude.eq.qmmm_no_atoms) then
        write(luout,22) "excluded bq's",tag,"none"
      else
        call errquit(pname//'invalid bq_exclude',0,RTDB_ERR)
      end if
      
      tag ="qmmm:linkatm" 
      if(linkatm.eq.qmmm_hydrogen) then
        write(luout,22) "link atom type",tag,"hydrogens"
      else if(linkatm.eq.qmmm_halogen) then
        write(luout,22) "link atom type",tag,"halogens"
      else
        call errquit(pname//'invalid link atom type',0,RTDB_ERR)
      end if

      tag ="qmmm:optimization"
      write(luout,22) "optimization region",tag,optimization


20    FORMAT(A,T24,A,T46,F12.6)
21    FORMAT(A,T24,A,T46,L3)
22    FORMAT(A,T24,A,T46,A)


      end

      subroutine qmmm_print_energy(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "qmmm.fh"

      integer irtdb
c
      character*32 pname
      double precision energy, mm_energy, qm_energy,eatoms

      
      pname = "qmmm_print_energy:"


      if (.not.rtdb_get(irtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if(.not.qmmm_master()) return

      write(*,19) 
      write(*,21) "QM/MM Energy"
      write(*,19) 
      write(*,21) "quantum energy", qm_energy, qm_energy*cau2kj
      write(*,21) "quantum energy adjusted", qm_energy-eatoms,
     >            (qm_energy-eatoms)*cau2kj
      write(*,21) "classical energy", mm_energy,
     >            mm_energy*cau2kj
      write(*,21) "total qmmm energy", energy, energy*cau2kj
      write(*,19) 

      write(*,*) 
19    FORMAT(72("-"))
20    FORMAT(1X,A,T34,F18.9)
21    FORMAT(1X,A,T34,F18.9," (",E12.6,2X, "kjoule/mol)")


      end

      subroutine qmmm_print_energy3(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "qmmm.fh"

      integer irtdb
c
      character*32 pname
      double precision energy, mm_energy, qm_energy,eatoms

      
      pname = "qmmm_print_energy:"


      if (.not.rtdb_get(irtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if(.not.qmmm_master()) return

      write(*,19) 
      write(*,21) "QM/MM Energy"
      write(*,19) 
      write(*,21) "quantum energy", qm_energy, qm_energy*cau2kj
      write(*,21) "quantum energy adjusted", qm_energy-eatoms,
     >            (qm_energy-eatoms)*cau2kj
      write(*,21) "classical energy", mm_energy,
     >            mm_energy*cau2kj
      write(*,21) "total qmmm energy", energy, energy*cau2kj
      write(*,19) 

      write(*,*) 
19    FORMAT("@",72("-"))
20    FORMAT("@",1X,A,T34,F18.9)
21    FORMAT("@",1X,A,T34,F18.9," (",E12.6,2X, "kjoule/mol)")


      end

      subroutine qmmm_print_energy2(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "qmmm.fh"

      integer irtdb
c
      character*32 pname
      double precision energy, mm_energy, qm_energy,eatoms

      
      pname = "qmmm_print_energy:"


      if (.not.rtdb_get(irtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

c      if (.not. rtdb_get(irtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
c     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if(.not.qmmm_master()) return

      mm_energy = energy-qm_energy+eatoms

      write(*,19) 
      write(*,21) "QM/MM Energy"
      write(*,19) 
      write(*,21) "quantum energy", qm_energy, qm_energy*cau2kj
      write(*,21) "quantum energy adjusted", qm_energy-eatoms,
     >            (qm_energy-eatoms)*cau2kj
      write(*,21) "classical energy", mm_energy,
     >            mm_energy*cau2kj
      write(*,21) "total qmmm energy", energy, energy*cau2kj
      write(*,19) 

      write(*,*) 
19    FORMAT(72("-"))
20    FORMAT(1X,A,T34,F18.9)
21    FORMAT(1X,A,T34,F18.9," (",E12.6,2X, "kjoule/mol)")


      end

      subroutine qmmm_print_energy1(irtdb,tot_energy,qm_energy)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "qmmm_params.fh"
#include "rtdb.fh"
#include "qmmm.fh"

      integer irtdb
      double precision tot_energy
      double precision energy, mm_energy, qm_energy,eatoms
c
      character*32 pname

      
      pname = "qmmm_print_energy1"


      if (.not.rtdb_get(irtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      energy = tot_energy

      qm_energy = qm_energy/cau2kj
      energy = energy/cau2kj
      mm_energy=energy-qm_energy
      qm_energy = qm_energy+eatoms


      if (.not. rtdb_put(irtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if (.not. rtdb_put(irtdb,'qmmm:mm_energy',mt_dbl,1,mm_energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      if (.not. rtdb_put(irtdb,'qmmm:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed get energy', 0, RTDB_ERR)

      call qmmm_print_energy(irtdb)

      end

