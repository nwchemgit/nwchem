c
c     qmd_thermostat: apply thermostat
c
      subroutine qmd_thermostat(n,dt,m,v,KE_in,KE_targ,act_at)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer n                ! number of atoms
      logical act_at(n)        ! active atoms
      double precision dt      ! step size
      double precision m(n)    ! mass
      double precision v(3,n)  ! velocities
      double precision KE_in   ! instantaneous kinetic energy
      double precision KE_targ ! target kinetic energy(only for rescale)
c
      character*30 pname
c
c     preliminaries
      pname='qmd_thermostat: '
c
c     calculate kinetic energy
      call qmd_kinetic(n,m,v,KE_in)
c
c     apply thermostat
      if (do_nosehoover) then
         call qmd_nosehoover(n,dt,m,v,KE_in)
      else if (do_berendsen) then
         call qmd_berendsen(n,dt,m,v,KE_in)
      else if (do_langevin) then
         call qmd_langevin(n,dt,m,v,KE_in,act_at)
      else if (do_rescale) then
         call qmd_rescale(n,m,v,KE_in,KE_targ)
      else if (do_svr) then
         call qmd_svr(n,dt,m,v,KE_in)
      else
         call errquit(pname//'unknown thermostat',0,RTDB_ERR)
      end if
c
      return
      end
c
c     qmd_nosehoover: nose-hoover thermostat
c
      subroutine qmd_nosehoover(n,dt,m,v,KE_in)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer n
      double precision dt
      double precision m(n)    ! mass
      double precision v(3,n)
      double precision KE_in
c
      double precision dt2,dt4,dt8
c
      integer i,j,ios,isy,ims,mdeg
      double precision scalefac,aa,kbt

c
      character*30 pname
c
c     preliminaries   
      pname='qmd_nosehoover: '
      dt2 = 0.5d0*dt
      dt4 = 0.25d0*dt
      dt8 = 0.125d0*dt
      kbt = kb_au*targ_temp

      mdeg = ndeg + 6
      if (do_linear) mdeg = mdeg - 1

      scalefac = 1d0
      m_nh(1) = m_nh(1)*mdeg

      g_nh(1) = (2d0*KE_in - mdeg*kbt)/m_nh(1)

      do ims=1,3
      do isy=1,nnhsy
        v_nh(nnhos) = v_nh(nnhos) + g_nh(nnhos)*dt4*syw_nh(isy)/3
        do ios=1,nnhos-1
          aa = dexp(-dt8*syw_nh(isy)*v_nh(nnhos+1-ios)/3)
          v_nh(nnhos-ios) = v_nh(nnhos-ios)*aa**2 +
     $                      dt4*syw_nh(isy)*g_nh(nnhos-ios)*aa/3
        enddo
        aa = dexp(-dt2*syw_nh(isy)*v_nh(1)/3)
        scalefac = scalefac*aa
        g_nh(1) = (2d0*KE_in*scalefac**2 - mdeg*kbt)/m_nh(1)
        do ios=1,nnhos
          r_nh(ios) = r_nh(ios) + v_nh(ios)*dt2*syw_nh(isy)/3
        enddo
        do ios=1,nnhos-1
          aa = dexp(-dt8*syw_nh(isy)*v_nh(ios+1)/3)
          v_nh(ios) = v_nh(ios)*aa**2 + dt4*syw_nh(isy)*g_nh(ios)*aa/3
          g_nh(ios+1) = (m_nh(ios)*v_nh(ios)**2 - kbt)/m_nh(ios+1)
        enddo
        v_nh(nnhos) = v_nh(nnhos) + g_nh(nnhos)*dt4*syw_nh(isy)/3
      enddo
      enddo

      do i=1,n
        v(1:3,i) = scalefac*v(1:3,i)
      enddo

      m_nh(1) = m_nh(1)/mdeg

c
c     nose-hoover 2 chain formula: Frenkel & Smit

c      g_nh(2) = (m_nh(1)*v_nh(1)*v_nh(1)-kbt)/m_nh(2)
c
c      v_nh(2) = v_nh(2) + g_nh(2)*dt4
c      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
c
c      g_nh(1) = (2.0d0*KE_in - ndeg*kbt)/m_nh(1)
c
c      v_nh(1) = v_nh(1) + g_nh(1)*dt4
c      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
c
c      r_nh(1) = r_nh(1) + v_nh(1)*dt2
c      r_nh(2) = r_nh(2) + v_nh(2)*dt2
c
c      scalefac = dexp(-v_nh(1)*dt2)
c
c     scale kinetic energy
c      KE_in = KE_in*scalefac*scalefac
c
c     nose-hoover 2 chain formula: Frenkel & Smit
c      g_nh(1) = (2.0d0*KE_in - ndeg*kbt)/m_nh(1)
c
c      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
c      v_nh(1) = v_nh(1) + g_nh(1)*dt4
c      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
c
c      g_nh(2) = (m_nh(1)*v_nh(1)*v_nh(1)-kbt)/m_nh(2)
c
c      v_nh(2) = v_nh(2) + g_nh(2)*dt4
c
      return
      end
c
c     qmd_berendsen: berendsen thermostat
c
      subroutine qmd_berendsen(n,dt,m,v,KE_in)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer n
      double precision dt
      double precision m(n)    ! mass
      double precision v(3,n)
      double precision KE_in
c
      integer i,j
      double precision inst_temp
      double precision ratio_temp
      double precision scalefac
c
      character*32 pname
      pname='qmd_berendsen: '
c
c     instantaneous temperature
      inst_temp = 2.0d0*KE_in/kb_au/(1.d0*ndeg)
c
c     velocity scaling factor: Berendsen scaling
      ratio_temp = targ_temp/inst_temp
      scalefac = dsqrt(1.d0 + (dt/(2.d0*tau))*(ratio_temp-1.d0))
c
      do i = 1,n
        v(1,i) = scalefac*v(1,i)
        v(2,i) = scalefac*v(2,i)
        v(3,i) = scalefac*v(3,i)
      end do
c
      return
      end
c
c     qmd_langevin: Langevin thermostat: 
c     Bussi and Parrinello, Phys. Rev. E (2007)
c
      subroutine qmd_langevin(n,dt,m,v,KE_in,act_at)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer n
      logical act_at(n)
      double precision dt
      double precision m(n)    ! mass
      double precision v(3,n)
      double precision KE_in
c
      integer i,j
      double precision inst_temp
      double precision ratio_temp
      double precision c1,c2
c
      external qmd_gasdev
      double precision qmd_gasdev
c
      character*32 pname
      pname='qmd_langevin: '
c
c     calculate c1
      c1 = exp(-friction*dt/2.d0)   
      do i = 1,n
        if (act_at(i)) then
          c2=sqrt((1.d0-c1**2)*(kb_au*targ_temp)/(m(i)*camu2au)) 
          do j = 1,3
           v(j,i) = c1*v(j,i)+c2*qmd_gasdev(idum)
          end do
        end if
      end do
c
      return
      end
c
c
c     qmd_rescale: velocity rescaling
c
      subroutine qmd_rescale(n,m,v,KE_in,KE_targ)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer n
      double precision m(n)
      double precision v(3,n)
      double precision KE_in
      double precision KE_targ
c
      integer i
      double precision ratio_KE
      double precision scalefac
c
      character*32 pname
      pname='qmd_rescale: '
c
c     velocity rescaling
      ratio_KE = KE_targ/KE_in
      scalefac = dsqrt(ratio_KE)
c
      do i = 1,n
        v(1,i) = scalefac*v(1,i)
        v(2,i) = scalefac*v(2,i)
        v(3,i) = scalefac*v(3,i)
      end do
c
      return
      end
c
c
c     qmd_svr: stochastic velocity rescaling
c     Bussi, Donadio, and Parrinello
c     J. Chem. Phys. 126, 014101 (2007)
c     see also
c     Comp. Phys. Comm. 179, 26 (2008)
c
      subroutine qmd_svr(n,dt,m,v,KE_in)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer n
      double precision dt
      double precision m(n)
      double precision v(3,n)
      double precision KE_in
c
      integer i
      double precision scheck
      double precision sigma
      double precision rr, ss, cc
      double precision scalefac
c
      external qmd_gasdev
      double precision qmd_gasdev
      external qmd_sumno
      double precision qmd_sumno
c
      character*32 pname
      pname='qmd_svr: '
c
c     idum, tau, ndeg, targ_temp in common block
c
      if (tau.ne.0.d0) then
        cc = exp(-dt/(2.d0*tau))
      else
        cc = 0.d0
      end if
c
      rr = qmd_gasdev(idum)
      sigma = dble(ndeg)*kb_au*targ_temp/2.d0
      ss = qmd_sumno(ndeg-1,idum)
c     If temperature set to zero, then sigma=0 and scheck would diverge
c      to positive infinity. Instead just set the value to 1 since we
c      only need the sign
      if (sigma.gt.1.d-8) then
        scheck = rr + dsqrt(cc*dble(ndeg)*KE_in/((1.d0-cc)*sigma))
      else
        scheck = 1.d0
      end if
c
      scalefac = cc + (sigma/(dble(ndeg)*KE_in))*(1.d0-cc)*(ss+rr**2) 
     $          + 2.d0*rr*dsqrt(cc*(1.d0-cc)*sigma/(dble(ndeg)*KE_in))
c
      scalefac = dsign(dsqrt(scalefac),scheck)
c
      do i = 1,n
        v(1,i) = scalefac*v(1,i)
        v(2,i) = scalefac*v(2,i)
        v(3,i) = scalefac*v(3,i)
      end do
c
      return
      end
c $Id$
