      Subroutine diis_driver(derr,ispin,ipol,icall,
     &     MXhist,nbf,geom,AO_bas_han,
     &     g_fock,g_dens,g_tmp,g_svecs,g_diis,isvals)

C$Id: diis_driver.F,v 1.5 1996-05-04 21:06:26 og845 Exp $
      implicit none
      
      integer ispin
      integer ipol
      integer MXhist
      integer nbf
      integer g_fock
      integer g_dens
      integer g_tmp
      integer g_svecs
      integer g_diis
      integer isvals
      Integer nhist(2), icall(2)
      integer geom,AO_bas_han
      double precision derr ! Diis error norm [output]
      
      integer lpulay,npulay,ipulay
      Save lpulay,npulay,ipulay,nhist

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"

#if defined(DFT_TIMER)
#include "timer.fh"
#endif
c
      logical LResult
      integer g_fsm12,g_s12p,me,nproc,iw,lenp
      integer lfwght,ifwght,lrow,irow,lcol,icol
      double precision one,zero
      parameter (one = 1.d0,zero = 0.d0)

c     
c     Use the DIIS scheme to accelerate convergence.
c     
c*****************************************************************************
      if( icall(ispin).lt.0 )then
        LResult = MA_Free_Heap(lpulay)
        return
      end if
      me=ga_nodeid()
      nproc=ga_nnodes()
c
#if defined(DFT_TIMER)
      call ga_sync
      call t_stat(10,-1)
#endif
c
c
c     Allocate scratch space.
c
      if( icall(1).eq.0 )then
        npulay  = ipol*(MXhist+1)*(MXhist+2)/2
        LResult = MA_Alloc_Get(MT_Dbl,npulay,'Pulay',lpulay,ipulay)
      end if

      LResult = MA_Push_Get(MT_Dbl,MXhist+1,'Fock wghts',lfwght,ifwght)
c
c     Increment call counter.
c
      icall(ispin) = icall(ispin) + 1
c
#if defined(DFT_TIMER)
      call dft_second(tp1)
#endif
c       
c     Construct a temporary global array.
c       
      call gacrea(geom,AO_bas_han,
     &     nbf,nbf,'ga temp2',g_fsm12,'atom')
c
C     Mat mult Fock * s-1/2
c
      iw=2
      call diis_bld12(dbl_mb(isvals),g_svecs,g_tmp,g_fsm12,nbf,iw)
      call ga_dgemm('N','N',nbf,nbf,nbf,
     *     one,g_fock,g_tmp,zero,g_fsm12)
      call gacrea(geom,AO_bas_han,
     &     nbf,nbf,'ga temp3',g_s12p,'atom')
      iw=3
      call diis_bld12(dbl_mb(isvals),g_svecs,g_tmp,g_s12p,nbf,iw)
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(1)='  s12   '
      call tcllct(t_sub(1,1),tp2,tp1)
      call dft_second(tp1)
#endif
c
C     Mat mult S+1/2*P (transp)
c
      call ga_dgemm('N','T',nbf,nbf,nbf,
     *     one,g_tmp,g_dens,zero,g_s12p)
c
C     Mat mult (S+*P(transp))(transp)*(F*S-)
c
      call ga_dgemm('N','N',nbf,nbf,nbf,
     *     one,g_s12p,g_fsm12,zero,g_tmp)
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(2)='  mxm   '
      call tcllct(t_sub(1,2),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Destroy temp globals.
c
      call gadest(g_s12p)
      call gadest(g_fsm12)

      LResult = MA_Push_Get(MT_Dbl,nbf,'row',lrow,irow)
      LResult = MA_Push_Get(MT_Dbl,nbf,'column',lcol,icol)

      call diis_refvec(nbf,g_tmp,dbl_mb(irow),DBL_MB(icol))

      LResult = MA_Pop_Stack(lcol)
      LResult = MA_Pop_Stack(lrow)

      call ga_sync
c
c
c     Compute Pulays B matrix.
c
      call diis_core(derr,Dbl_MB(ipulay),MXhist,nhist,
     &     icall,nbf,ispin,ipol,
     &     g_diis,g_tmp)
      lenp = (mxhist+1)*(mxhist+2)/2
c
c
c     Solve DIIS equation for Hamiltonian weights.
c
      call diis_bsolve(Dbl_MB(ipulay+lenp*(ispin-1)),Dbl_MB(ifwght),
     &     MXhist,nhist(ispin),icall(ispin))
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(3)='  bsolve   '
      call tcllct(t_sub(1,3),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Compute an extrapolated Hamiltonian using the weights.
c
      call diis_hamwgt(Dbl_MB(ifwght),MXhist,nhist,nbf,icall,
     &     ispin,ipol,
     &     g_fock,g_diis,g_tmp)
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(4)='  hamwght '
      call tcllct(t_sub(1,4),tp2,tp1)
#endif
c
      LResult = MA_Pop_Stack(lfwght)
c
#if defined(DFT_TIMER)
      ltime = 7
      if( me.eq.0 )then
        write(ltime,*)
        write(ltime,9996) '      Analysis of DIIS'
        write(ltime,*)
        write(ltime,9997) ' number of nodes:  ',nproc
        write(ltime,*)
      end if
      subname='DIIS  '
      call t_stat(4,1)
 9996 format(10x,a)
 9997 format(10x,a,i6)
#endif
c
      return
      end
