      subroutine ccsd_get_param(rtdb, maxit,convi,iprt, theory,
     $   maxdiis)
C     $Id: ccsd_get_param.F,v 1.5 1997-04-14 06:34:55 gg502 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
c     
      integer rtdb, maxit, convi, iprt, maxdiis
      Character*(*) theory
c
      double precision mone
      parameter (mone = -1.0d0)
#include "ccsd_debug.fh"
c     
      if (.not. rtdb_get(rtdb, 'ccsd:maxiter',
     $     MT_INT, 1, maxit)) maxit = 20
      if (.not. rtdb_get(rtdb, 'ccsd:thresh',
     $     MT_INT, 1, convi)) convi = 6
      if (.not. rtdb_get(rtdb, 'ccsd:iprt',
     $     MT_INT, 1, iprt)) iprt = 0
      if (.not. rtdb_get(rtdb, 'ccsd:maxdiis',
     $     MT_INT, 1, maxdiis)) maxdiis = 5
C
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $   call errquit('task:energy: theory not specified',0)
c
C     Set these guys in the /DEBUG/ common block
C
      If ( .NOT. RTDB_Get( RTDB, 'ccsd:doa', MT_Int, LenDbg, DoA) )
     $   Call IFill(LenDbg, 2, DoA, 1)
      If ( .NOT. RTDB_Get( RTDB, 'ccsd:dob', MT_Int, LenDbg, DoB) )
     $   Call IFill(LenDbg, 2, DoB, 1)
      If ( .NOT. RTDB_Get( RTDB, 'ccsd:dog', MT_Int, LenDbg, DoG) )
     $   Call IFill(LenDbg, 1, DoG, 1)
      If ( .NOT. RTDB_Get( RTDB, 'ccsd:doh', MT_Int, LenDbg, Doh) )
     $   Call IFill(LenDbg, 1, Doh, 1)
      If ( .NOT. RTDB_Get( RTDB, 'ccsd:dojk', MT_Int, LenDbg, Dojk) )
     $   Call IFill(LenDbg, 2, Dojk, 1)
      If ( .NOT. RTDB_Get( RTDB, 'ccsd:dos', MT_Int, LenDbg, DoS) )
     $   Call IFill(LenDbg, 1, DoS, 1)
      If ( .NOT. RTDB_Get( RTDB, 'ccsd:dod', MT_Int, LenDbg, DoD) )
     $   Call IFill(LenDbg, 1, DoD, 1)
C
c$$$      do i=1,LenDbg
c$$$       doa(i)=2
c$$$       dob(i)=2
c$$$       dog(i)=1
c$$$       doh(i)=1
c$$$       dojk(i)=2
c$$$       dos(i)=1
c$$$       dod(i)=1
c$$$      enddo
c
      end
