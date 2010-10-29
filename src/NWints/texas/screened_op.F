c
c     Interface and computation routines for screened exchange project
c
c     BGJ - 2/01
c
c     $Id$
c
      function screened_exchange()
      implicit none
      logical screened_exchange, screened_x
      double precision Ks
#include "bgj.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      double precision scrk_Ks
      logical scrk_op_enabled
      common /screened_k/ scrk_Ks, scrk_op_enabled
      logical initd
      save initd
      data initd / .false. /
c
      if (.not. initd) then
         scrk_op_enabled = .false.
         initd = .true.
      endif
c
      screened_x = rtdb_get(bgj_get_rtdb_handle(),
     &     'screened_exchange:Ks',MT_DBL,1,Ks)
      scrk_Ks = Ks
      screened_exchange = screened_x
c
      return
      end
c
      function screened_operator_Ks()
      implicit none
      double precision screened_operator_Ks
#include "bgj.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      double precision scrk_Ks
      logical scrk_op_enabled
      common /screened_k/ scrk_Ks, scrk_op_enabled
c
      screened_operator_Ks = scrk_Ks
c
      return
      end
c
      subroutine enable_screened_operator()
      implicit none
      double precision Ks
#include "bgj.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      double precision scrk_Ks
      logical scrk_op_enabled
      common /screened_k/ scrk_Ks, scrk_op_enabled
c
      if (.not.rtdb_put(bgj_get_rtdb_handle(),
     &     'screeend_operator_enabled',MT_LOG,1,.true.))
     &     call errquit('Error enabling screened operator',0, INT_ERR)
c
      scrk_op_enabled = .true.
      if (.not.rtdb_get(bgj_get_rtdb_handle(),
     &     'screened_exchange:Ks',MT_DBL,1,Ks))
     &     call errquit('Could not get screening constant',0, INT_ERR)
      scrk_Ks = Ks
c
      return
      end
c
      subroutine disable_screened_operator()
      implicit none
#include "bgj.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      double precision scrk_Ks
      logical scrk_op_enabled
      common /screened_k/ scrk_Ks, scrk_op_enabled
c
      if (.not.rtdb_put(bgj_get_rtdb_handle(),
     &     'screened_operator_enabled',MT_LOG,1,.false.))
     &     call errquit('Error disabling screened operator',0, INT_ERR)
c
      scrk_op_enabled = .false.
c
      return
      end
c
      function screened_operator_enabled()
      implicit none
      logical screened_operator_enabled
#include "bgj.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      double precision scrk_Ks
      logical scrk_op_enabled
      common /screened_k/ scrk_Ks, scrk_op_enabled
c
      screened_operator_enabled = scrk_op_enabled
c
      return
      end
c
      subroutine density_average(rho, qwgt, nq, ipol, ifunc)
      implicit none
      integer nq, ipol, ifunc, i
      double precision t(100), rho(nq,*), qwgt(nq), p, rho_avg,
     &     pi, CF, CTF
      save t
#include "screened_op.fh"
      if (ifunc .eq. 0) then
c        initialize
         call dfill(100,0.D0,t,1)
         return
      endif
      if (ifunc .lt. 0) then
c        print
         pi = acos(-1d0)
         CF = (3.d0*pi*pi)**(1d0/3d0)
         CTF = sqrt(4d0*CF/pi)
         p = t(1)
         rho_avg = t(2)/p
#if 0
         write(*,*)'=== Density average begin ==='
         write(*,*)'int rho:         ',t(1)
         write(*,*)'avg rho:         ',rho_avg, rho_avg/p
         write(*,*)'avg kTF:         ',CTF*t(3)/p
         write(*,*)'kTF from avg rho:',CTF*rho_avg**(1d0/6d0)
         write(*,*)'avg kF:          ',CF*t(4)/p
         write(*,*)'kF from avg rho: ',CF*rho_avg**(1d0/3d0)
         write(*,*)'=== Density average end ==='
#endif
         return
      endif
c
c     Update integrated quantities
c
c     1 - rho
c     2 - rho**2 (average of rho)
c     3 - rho**(7/6) (average of rho**(1/6))
c     4 - rho**(4/3) (average of rho**(1/3))     
c
      do i = 1, nq
         p = rho(i,1)
         t(1) = t(1) + qwgt(i) * p
         t(2) = t(2) + qwgt(i) * p*p
         t(3) = t(3) + qwgt(i) * p*p**(1d0/6d0)
         t(4) = t(4) + qwgt(i) * p*p**(1d0/3d0)
      enddo
      return
      end
