      Subroutine diis_driver(toll_s, derr, ispin, ipol, icall, MXhist, 
     &                       nbf, geom, AO_bas_han, g_fock, g_dens, 
     &                       g_tmp, g_svecs, isvals, diis, nodiis)
c
C$Id: diis_driver.F,v 1.13 1999-06-16 00:10:27 d3e129 Exp $
c
      implicit none
c      
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
      integer geom, AO_bas_han
      double precision derr ! Diis error norm [output]
      double precision toll_s
      integer lpulay, npulay, ipulay
      Save lpulay, npulay, ipulay, nhist, g_diis
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "msgids.fh"
c
#if defined(DFT_TIMER)
#include "timer.fh"
#endif
c
      integer g_fsm12, g_s12p, me, nproc, iw, lenp
      integer lfwght, ifwght, lrow, irow, lcol, icol
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical diis, nodiis
      double precision avail, dw_agg_avail
      integer n_n2_objs
      double precision one, zero
      parameter (one = 1.d0, zero = 0.d0)
c
c     Use the DIIS scheme to accelerate convergence.
c     
      if (icall(ispin).lt.0 )then
         if (.not.MA_Free_Heap(lpulay))
     &      call errquit('diis_driver: cannot free heap',0)
c         call gadest(g_diis)
         if (.not. ga_destroy(g_diis)) call errquit
     &      ('diis_driver: could not destroy g_diis', 0)
         return
      endif
      me = ga_nodeid()
      nproc = ga_nnodes()
c
#if defined(DFT_TIMER)
      call ga_sync
      call t_stat(10,-1)
#endif
c
c     Allocate scratch space.
c
      if ( icall(1).eq.0 )then
c
c       Create necessary DIIS arrays.
c
        npulay  = ipol*(MXhist+1)*(MXhist+2)/2
        if (.not.MA_Alloc_Get(MT_Dbl,npulay,'Pulay',lpulay,ipulay))
     &     call errquit('diis_driver: cannot allocate Pulay',0)
        call dfill(npulay, 0.0d0, dbl_mb(ipulay), 1)
c
c       Inquire about amount of global memory left.  DIIS uses
c       several n^2 arrays to store fock and error matrices and
c       needs probably at least 3 of each to be effective.
c
        avail = ga_memory_avail()
        call ga_igop(msg_min_ga_avail, avail, 1, 'min')
c
c        - subtract off a couple of n^2 objects for later use,
c
        dw_agg_avail = dble(nproc)*dble(avail)/8.0d0-2.0d0*dble(nbf)**2
        n_n2_objs = nint(dw_agg_avail/dble(nbf)**2)
        if (2*mxhist*ipol.gt.n_n2_objs)then
c
c          Reset mxhist.
c          
           mxhist = n_n2_objs/(2*ipol)
           if (mxhist.lt.3)then
              diis = .false.
              nodiis = .true.
              if (me.eq.0)write(LuOut,2118)
              return
           else
              if (me.eq.0)write(LuOut,2119)mxhist
           endif
         endif
c        call gacrea(geom,AO_bas_han,
c     &    nbf,nbf*(2*mxhist*ipol),
c     &    'diis bufs',g_diis,'row')
*ga:1:0
         if (.not. ga_create(mt_dbl, nbf, nbf*(2*mxhist*ipol), 
     &                       'diis bufs', nbf, 0, g_diis)) 
     &      call errquit('diis_driver: error creating g_diis',0)
         call ga_zero(g_diis)
c        if (me.eq.0)write(LuOut,*)'avail,dw_agg_avail,n_n2_objs,mxhist:',
c     &                            avail,dw_agg_avail,n_n2_objs,mxhist
      endif
 2118 format(/,10x,'Not enough global memory currently available for ',
     &       /,10x,'DIIS.  DIIS option disabled. ')
 2119 format(/,10x,'The amount of global memory currently available',
     &       /,10x,'for DIIS only allows: ',i2,' Fock matrices to be',
     &       /,10x,'stored and used in the extrapolation.')
      if (.not.MA_Push_Get(MT_Dbl,MXhist+1,'Fock wghts',lfwght,ifwght))
     &   call errquit('diis_driver: cannot allocate Fock wghts',0)
      call dfill(MXhist+1, 0.0d0, dbl_mb(ifwght), 1)
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
c      call gacrea(geom,AO_bas_han,
c     &     nbf,nbf,'ga temp2',g_fsm12,'atom')
      g_fsm12 = ga_create_atom_blocked(geom, AO_bas_han, 'ga temp2')
c
c     Mat mult Fock * s-1/2
c
      iw = 2
      call diis_bld12(toll_s, dbl_mb(isvals), g_svecs, g_tmp, g_fsm12,
     &                nbf, iw)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_fock,g_tmp,zero,g_fsm12)
c      call gacrea(geom,AO_bas_han,
c     &     nbf,nbf,'ga temp3',g_s12p,'atom')
      g_s12p = ga_create_atom_blocked(geom, AO_bas_han, 'ga temp3')
      iw = 3
      call diis_bld12(toll_s, dbl_mb(isvals), g_svecs, g_tmp, g_s12p, 
     &                nbf, iw)
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(1) = '  s12   '
      call tcllct(t_sub(1,1),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Mat mult S+1/2*P (transp)
c
      call ga_dgemm('N','T',nbf,nbf,nbf,one,g_tmp,g_dens,zero,g_s12p)
c
C     Mat mult (S+*P(transp))(transp)*(F*S-)
c
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_s12p,g_fsm12,zero,g_tmp)
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(2) = '  mxm   '
      call tcllct(t_sub(1,2),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Destroy temp globals.
c
c      call gadest(g_s12p)
      if (.not. ga_destroy(g_s12p)) call errquit
     &   ('diis_driver: could not destroy g_s12p', 0)
c      call gadest(g_fsm12)
      if (.not. ga_destroy(g_fsm12)) call errquit
     &   ('diis_driver: could not destroy g_fsm12', 0)
c
      if (.not.MA_Push_Get(MT_Dbl,nbf,'row',lrow,irow))
     &   call errquit('diis_driver: cannot allocate row',0)
      if (.not.MA_Push_Get(MT_Dbl,nbf,'column',lcol,icol))
     &   call errquit('diis_driver: cannot allocate',0)
c
      call diis_refvec(nbf,g_tmp,dbl_mb(irow),dbl_mb(icol))
c
      if (.not.ma_pop_stack(lcol))
     &   call errquit('diis_driver: cannot pop stack',0)
      if (.not.ma_pop_stack(lrow))
     &   call errquit('diis_driver: cannot pop stack',0)
c
      call ga_sync
c
c
c     Compute Pulays B matrix.
c
      call diis_core(derr, Dbl_MB(ipulay), MXhist, nhist, icall, nbf,
     &               ispin, ipol, g_diis, g_tmp)
c
      lenp = (mxhist+1)*(mxhist+2)/2
c
c     Solve DIIS equation for Hamiltonian weights.
c
      call diis_bsolve(Dbl_MB(ipulay+lenp*(ispin-1)), Dbl_MB(ifwght),
     &                 MXhist, nhist(ispin), icall(ispin))
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(3) = '  bsolve   '
      call tcllct(t_sub(1,3),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Compute an extrapolated Hamiltonian using the weights.
c
      call diis_hamwgt(Dbl_MB(ifwght), MXhist, nhist, nbf, icall, 
     &                 ispin, ipol, g_fock, g_diis, g_tmp)
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(4) = '  hamwght '
      call tcllct(t_sub(1,4),tp2,tp1)
#endif
c
      if (.not.ma_pop_stack(lfwght))
     &   call errquit('diis_driver: cannot pop stack',0)
c
#if defined(DFT_TIMER)
      ltime = 7
      if ( me.eq.0 )then
         write(ltime,*)
         write(ltime,9996) '      Analysis of DIIS'
         write(ltime,*)
         write(ltime,9997) ' number of nodes:  ',nproc
         write(ltime,*)
      endif
      subname = 'DIIS  '
      call t_stat(4,1)
 9996 format(10x,a)
 9997 format(10x,a,i6)
#endif
c
      return
      end
