      Subroutine diis_driver(toll_s, derr, ipol, icall, MXhist, 
     &                       nbf, geom, AO_bas_han, g_focks, g_dens, 
     &                       g_tmp, g_svecs, svals, diis, nodiis,
     .     g_diis,keeps,g_s12m,g_s12p,reset)
c
C$Id: diis_driver.F,v 1.29 2004-01-20 19:08:43 edo Exp $
c
      implicit none
#include "errquit.fh"
c      
      integer ispin
      integer ipol
      integer MXhist
      integer nbf
      integer g_focks(2)
      integer g_dens(2)
      integer g_tmp
      integer g_svecs
      integer g_diis(mxhist)
      double precision svals(*)
      Integer nhist(2), icall(2)
      integer geom, AO_bas_han
      double precision derr(2) ! Diis error norm [output]
      double precision toll_s
      integer lpulay, npulay, ipulay
      logical reset
      Save lpulay, npulay, ipulay, nhist
c
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "dftps.fh"
c
c
      integer g_s12m, g_s12p, me, nproc, g_tmp2,g_rr
      integer lfwght, ifwght
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical diis, nodiis,keeps
      double precision avail, dw_agg_avail
      integer n_n2_objs,i
      double precision one, zero
      parameter (one = 1.d0, zero = 0.d0)
c
c
c     Use the DIIS scheme to accelerate convergence.
c     
      if(reset) then
         if( icall(1).eq.0 ) return
         icall(1) = -1
      endif
      if (icall(1).lt.0 )then
         if (.not.MA_Free_Heap(lpulay))
     &      call errquit('diis_driver: cannot free heap',0, MA_ERR)
         do i=1,2*mxhist*ipol
            if (.not. ga_destroy(g_diis(i))) call errquit
     &           ('diis_driver: could not destroy g_diis', i, GA_ERR)
         enddo
         if(reset) then
            icall(1) = 0
            icall(2) = 0
         endif
         return
      endif
      if (odftps) call pstat_on(ps_diis)
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     Allocate scratch space.
c
      if ( icall(1).eq.0 )then
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
              if (odftps) call pstat_off(ps_diis)
              return
           else
              if (me.eq.0)write(LuOut,2119)mxhist
           endif
         endif
*ga:1:0
         do i=1,mxhist*ipol*2
            g_diis(i) = ga_create_atom_blocked(
     .           geom, AO_bas_han, 'g_diis')
         enddo
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
      do ispin=1,ipol
      call dfill(MXhist+1, 0.0d0, dbl_mb(ifwght), 1)
c
c     Increment call counter.
c
      icall(ispin) = icall(ispin) + 1
c
c     Construct a temporary global array.
c       
      g_tmp2=ga_create_atom_blocked(geom, AO_bas_han, 'gtmp2')
      g_rr = ga_create_atom_blocked(geom, AO_bas_han, 'grr')
      if (odftps) call pstat_on(ps_bld12)
      if(.not.keeps) then
         g_s12m=ga_create_atom_blocked(geom, AO_bas_han, 'gs12m')
c
c     get s-1/2
c
         call diis_bld12(toll_s, svals, g_svecs, g_s12m, 
     .        g_tmp,2)
      endif
c
c     Mat mult Fock * s-1/2
c
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_focks(ispin),
     .     g_s12m,zero,g_rr)
      if (odftps) call pstat_off(ps_dgemm)
      if(.not.keeps) then
         if (.not. ga_destroy(g_s12m)) call errquit
     &        ('diis_driver: could not destroy g_s12m', 0, GA_ERR)
      endif

c
c     get s+1/2
c
      if(.not.keeps) then
         g_s12p = ga_create_atom_blocked(geom, AO_bas_han, 'gs12p')
         call diis_bld12(toll_s, svals, g_svecs,  g_s12p, 
     .     g_tmp2, 3)
      endif
      if (odftps) call pstat_off(ps_bld12)
c
c     Mat mult S+1/2*P (transp)
c
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_s12p,g_dens(ispin),
     ,     zero,g_tmp2)
      if (odftps) call pstat_off(ps_dgemm)
      if(.not.keeps) then
         if (.not. ga_destroy(g_s12p)) call errquit
     &        ('diis_driver: could not destroy g_s12p', 0, GA_ERR)
      endif
c
C     Mat mult (S+*P(transp))*(F*S-)
c
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_tmp2,g_rr,zero,g_tmp)
      if (odftps) call pstat_off(ps_dgemm)
      if (.not. ga_destroy(g_rr)) call errquit
     &     ('diis_driver: could not destroy g_rr', 0, GA_ERR)
      
c
      call ga_transpose(g_tmp,g_tmp2)
      call ga_add(-1d0,g_tmp,1d0,g_tmp2,g_tmp)

      if (.not. ga_destroy(g_tmp2)) call errquit
     &     ('diis_driver: could not destroy g_tmp2', 0, GA_ERR)
c
c
c     Compute Pulays B matrix.
c
      call diis_core(derr(ispin), Dbl_MB(ipulay), MXhist, nhist,
     .     icall(ispin), 
     &               ispin, ipol, g_diis, g_tmp)
c
      enddo
c
c     Solve DIIS equation for Hamiltonian weights.
c
      call diis_bsolve(Dbl_MB(ipulay), Dbl_MB(ifwght),
     &                 MXhist, nhist(1))
c
c     Compute an extrapolated Hamiltonian using the weights.
c
      call diis_hamwgt(Dbl_MB(ifwght), MXhist, nhist, icall, 
     &                 ipol,  g_focks, g_diis)
c
      if (.not.ma_pop_stack(lfwght))
     &   call errquit('diis_driver: cannot pop stack',0, MA_ERR)
c
      if (odftps) call pstat_off(ps_diis)
      return
      end
