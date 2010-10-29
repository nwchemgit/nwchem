      Subroutine diis_driver_so(toll_s, derr, icall, MXhist, 
     &     nbf, g_fock, g_dens, 
     &     g_svecs, svals, diis, nodiis)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
cso      integer ispin
cso      integer ipol
      integer MXhist
      integer nbf, nbf12
      integer g_fock(2)
      integer g_dens(2)
      integer g_tmp,g_scr,g_scr2
      integer g_svecs
      integer g_diis(2)
      double precision svals(*)
      Integer nhist, icall
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
      integer g_fsm12(2), g_s12p(2), me, nproc, iw, lenp
      integer g_pf(2)
      integer lfwght,ifwght
      integer lrowre,irowre,lcolre,icolre,lrowim,irowim,lcolim,icolim
c     
      logical diis, nodiis
      double precision avail, dw_agg_avail
      integer n_n2_objs
      double precision one, mone, zero
      parameter (one = 1.d0, mone = -1.0d0, zero = 0.d0)
c      integer i
cso
cso   initialize some variables 
cso 
      nbf12 = nbf/2
c
c     Use the DIIS scheme to accelerate convergence.
c     
      if (icall.lt.0 )then
         if (.not.MA_Free_Heap(lpulay))
     &      call errquit('diis_driver: cannot free heap',0, MA_ERR)
         if (.not. ga_destroy(g_diis(1))) call errquit
     &      ('diis_driver: could not destroy g_diis(1)', 0, GA_ERR)
         if (.not. ga_destroy(g_diis(2))) call errquit
     &      ('diis_driver: could not destroy g_diis(2)', 0, GA_ERR)
         return
      endif
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     Allocate scratch space.
c
      if ( icall.eq.0 )then
c
c       Create necessary DIIS arrays.
c
        npulay  = (MXhist+1)*(MXhist+2)/2
        if (.not.MA_Alloc_Get(MT_Dbl,npulay,'Pulay',lpulay,ipulay))
     &     call errquit('diis_driver: cannot allocate Pulay',0, MA_ERR)
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
        dw_agg_avail = dble(nproc)*dble(avail)/16.0d0-4.0d0*dble(nbf)**2
        n_n2_objs = nint(dw_agg_avail/(2*dble(nbf)**2))
        if (2*mxhist.gt.n_n2_objs)then
c
c          Reset mxhist.
c          
           mxhist = n_n2_objs/4
           if (mxhist.lt.3)then
              diis = .false.
              nodiis = .true.
              if (me.eq.0)write(LuOut,2118)
              return
           else
              if (me.eq.0)write(LuOut,2119)mxhist
           endif
         endif
         if (.not. ga_create(mt_dbl, nbf, nbf*(2*mxhist), 
     &                       're diis bufs', nbf, 1, g_diis(1))) 
     &      call errquit('diis_driver: error creating g_diis(1)',0,
     &       GA_ERR)
         call ga_zero(g_diis(1))
         if (.not. ga_create(mt_dbl, nbf, nbf*(2*mxhist), 
     &                       'im diis bufs', nbf, 1, g_diis(2))) 
     &      call errquit('diis_driver: error creating g_diis(2)',0,
     &       GA_ERR)
         call ga_zero(g_diis(2))
c        if (me.eq.0)write(LuOut,*)'avail,dw_agg_avail,n_n2_objs,mxhist:',
c     &                            avail,dw_agg_avail,n_n2_objs,mxhist
      endif
 2118 format(/,10x,'Not enough global memory currently available for ',
     &       /,10x,'DIIS.  DIIS option disabled. ')
 2119 format(/,10x,'The amount of global memory currently available',
     &       /,10x,'for DIIS only allows: ',i2,' Fock matrices to be',
     &       /,10x,'stored and used in the extrapolation.')
      if (.not.MA_Push_Get(MT_Dbl,MXhist+1,'Fock wghts',lfwght,ifwght))
     &   call errquit('diis_driver: cannot allocate Fock wghts',0,
     &       MA_ERR)
      call dfill(MXhist+1, 0.0d0, dbl_mb(ifwght), 1)
c
c     Increment call counter.
c
      icall = icall + 1
c
c     Construct temporary global arrays. g_scr is used to construct sm12, 
c     g_fsmre and g_fsmim are for real and imaginary parts of Fock * s-1/2
c
      if(.not.ga_create(mt_dbl,nbf12,nbf12,'scr',0,0,g_scr))
     &     call errquit('diis_driver: error creating g_scr',0, GA_ERR)
      if(.not.ga_create(mt_dbl, nbf12,nbf12,'tmp',0,0,g_tmp))
     &     call errquit('diis_driver: error creating g_tmp',0, GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nbf,'scr2',0,0,g_scr2))
     &     call errquit('diis_driver: error creating g_scr2',0, GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nbf,'fsm12re',0,0,g_fsm12(1)))
     &     call errquit('diis_driver: error creating g_fsm12(1)',0,
     &       GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nbf,'fsm12im',0,0,g_fsm12(2)))
     &     call errquit('diis_driver: error creating g_fsm12(2)',0,
     &       GA_ERR)
c
c     Mat mult Fock * s-1/2
c
      iw = 2
      call diis_bld12_so(toll_s, svals, g_svecs, g_scr, g_tmp,
     &                nbf12, iw)
      call ga_zero(g_scr2)
      call ga_fock_sf(g_scr, g_scr2, nbf12)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_fock(1),g_scr2,
     &     zero,g_fsm12(1))
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_fock(2),g_scr2,
     &     zero,g_fsm12(2))
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'s12pre',0,0,g_s12p(1)))
     &     call errquit('diis_driver: error creating g_s12p(1)',0,
     &       GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nbf,'s12pre',0,0,g_s12p(2)))
     &     call errquit('diis_driver: error creating g_s12p(2)',0,
     &       GA_ERR)
c
c     Mat mult S+1/2*P (transp)
c
      iw = 3
      call diis_bld12_so(toll_s, svals, g_svecs, g_scr, g_tmp, 
     &                nbf12, iw)
      call ga_zero(g_scr2)
      call ga_fock_sf(g_scr, g_scr2, nbf12)

      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_scr2,g_dens(1),
     &     zero,g_s12p(1))
      call ga_dgemm('N','N',nbf,nbf,nbf,mone,g_scr2,g_dens(2),
     &     zero,g_s12p(2))
c
C     Mat mult (S+*P(transp))(transp)*(F*S-)
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'pfre',0,0,g_pf(1)))
     &     call errquit('diis_driver: error creating g_pf(1)',0, GA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nbf,'pfim',0,0,g_pf(2)))
     &     call errquit('diis_driver: error creating g_pf(2)',0, GA_ERR)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_s12p(1),g_fsm12(1),
     &     zero,g_pf(1))
      call ga_dgemm('N','N',nbf,nbf,nbf,mone,g_s12p(2),g_fsm12(2),
     &     one,g_pf(1)) 
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_s12p(1),g_fsm12(2),
     &     zero,g_pf(2))
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_s12p(2),g_fsm12(1),
     &     one,g_pf(2)) 
c
c     Destroy temp globals.
c
      if (.not.ga_destroy(g_s12p(1))) call errquit
     &   ('diis_driver: could not destroy g_s12p(1)', 0, GA_ERR)
      if (.not.ga_destroy(g_s12p(2))) call errquit
     &   ('diis_driver: could not destroy g_s12p(2)', 0, GA_ERR)
      if (.not. ga_destroy(g_fsm12(1))) call errquit
     &   ('diis_driver: could not destroy g_fsm12(1)', 0, GA_ERR)
      if (.not. ga_destroy(g_fsm12(2))) call errquit
     &   ('diis_driver: could not destroy g_fsm12(2)', 0, GA_ERR)
      if (.not. ga_destroy(g_scr)) call errquit
     &   ('diis_driver: could not destroy g_scr', 0, GA_ERR)
      if (.not. ga_destroy(g_scr2)) call errquit
     &   ('diis_driver: could not destroy g_scr2', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp)) call errquit
     &   ('diis_driver: could not destroy g_tmp', 0, GA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl,nbf,'rowim',lrowim,irowim))
     &   call errquit('diis_driver: cannot allocate row',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nbf,'columnim',lcolim,icolim))
     &   call errquit('diis_driver: cannot allocate',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nbf,'rowre',lrowre,irowre))
     &   call errquit('diis_driver: cannot allocate rowre',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nbf,'columnre',lcolre,icolre))
     &   call errquit('diis_driver: cannot allocate',0, MA_ERR)
c
      call diis_refvec_so(nbf,g_pf,dbl_mb(irowre),dbl_mb(icolre),
     &     dbl_mb(irowim),dbl_mb(icolim))
c
      if (.not.ma_pop_stack(lcolre))
     &   call errquit('diis_driver: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lrowre))
     &   call errquit('diis_driver: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcolim))
     &   call errquit('diis_driver: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lrowim))
     &   call errquit('diis_driver: cannot pop stack',0, MA_ERR)
c
      call ga_sync
c
c
c     Compute Pulays B matrix.
c
      call diis_core_so(derr,Dbl_MB(ipulay),MXhist,nhist,icall,nbf,
     &     g_diis, g_pf)
c     
      lenp = (mxhist+1)*(mxhist+2)/2
c
c     Solve DIIS equation for Hamiltonian weights.
c
c      write(*,*)"nhist=", nhist
c      write(*,*)(Dbl_MB(ipulay+i-1),i=1,(nhist+1)*(nhist+2)/2)
      call diis_bsolve_so(Dbl_MB(ipulay), Dbl_MB(ifwght),
     &     MXhist, nhist)
c      write(*,*)(Dbl_MB(ifwght+i-1), i=1,nhist+1)
c
c     Compute an extrapolated Hamiltonian using the weights.
c
c      write(*,*)"before weight"
      call diis_hamwgt_so(Dbl_MB(ifwght), MXhist, nhist, nbf, icall, 
     &     g_fock, g_diis, g_pf)
c      write(*,*)"after diis"
c
      if (.not.ma_pop_stack(lfwght))
     &     call errquit('diis_driver: cannot pop stack',0, MA_ERR)
c
      if(.not.ga_destroy(g_pf(1)))
     &     call errquit('diis_driver:can not destroy pa_pf',0, GA_ERR)
      if(.not.ga_destroy(g_pf(2)))
     &     call errquit('diis_driver:can not destroy pa_pf',0, GA_ERR)
      return
      end
