      double precision function dft_dencvg(g_dens,g_denso,nbf)
C$Id: dft_dencvg.F,v 1.4 1996-07-03 21:47:26 d3h449 Exp $
      implicit none

      integer g_dens ! GA handle for Densit Matrix [input/output]
      integer g_denso ! GA handle for old DM [input]
      integer nbf

      logical LResult
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
c
      integer idone,me,nproc,i,j
      integer ldm,idm,ldmo,idmo
      double precision rms
c
      idone=0
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
c--   >  Evaluate change in energy.
c     
c--   >  Check for density matrix convergence.
c     

      rms = 0.0d0
      if(.not.MA_Push_Get(MT_Dbl,nbf,'dmat',ldm,idm))
     &  call errquit('dft_dencvg: cannot allocate dmat',0)
      if(.not.MA_Push_Get(MT_Dbl,nbf,'dmatold',ldmo,idmo))
     &  call errquit('dft_dencvg: cannot allocate dmatold',0)
      do i=1+me,nbf,nproc
        call get_col(g_dens,nbf,i,Dbl_MB(idm))
        call get_col(g_denso,nbf,i,Dbl_MB(idmo))
        do j=0,nbf-1
          rms=rms+(DBL_MB(idm+j)-DBL_MB(idmo+j))**2
        enddo
      enddo
      if(.not.ma_pop_stack(ldmo))
     &  call errquit('dft_dencvg: cannot pop stack',0)
      if(.not.ma_pop_stack(ldm))
     &  call errquit('dft_dencvg: cannot pop stack',0)
      call ga_sync
      call ga_dgop(Msg_rms_dens,rms,1,'+')
      dft_dencvg = rms/dfloat(nbf*nbf)

      return
      end
