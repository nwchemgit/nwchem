      double precision function dft_dencvg(g_dens,g_denso,nbf)
C$Id: dft_dencvg.F,v 1.7 2003-10-17 22:55:44 carlfahl Exp $
      implicit none
#include "errquit.fh"
c
      integer g_dens ! GA handle for Densit Matrix [input/output]
      integer g_denso ! GA handle for old DM [input]
      integer nbf
c
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
     &  call errquit('dft_dencvg: cannot allocate dmat',0, MA_ERR)
      if(.not.MA_Push_Get(MT_Dbl,nbf,'dmatold',ldmo,idmo))
     &  call errquit('dft_dencvg: cannot allocate dmatold',0, MA_ERR)
      do i=1+me,nbf,nproc
        call get_col(g_dens,nbf,i,Dbl_MB(idm))
        call get_col(g_denso,nbf,i,Dbl_MB(idmo))
        do j=0,nbf-1
          rms=rms+(DBL_MB(idm+j)-DBL_MB(idmo+j))**2
        enddo
      enddo
      if(.not.ma_pop_stack(ldmo))
     &  call errquit('dft_dencvg: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(ldm))
     &  call errquit('dft_dencvg: cannot pop stack',0, MA_ERR)
      call ga_sync
      call ga_dgop(Msg_rms_dens,rms,1,'+')
      dft_dencvg = rms/dble(nbf*nbf)

      return
      end
