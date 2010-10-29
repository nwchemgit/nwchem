      double precision function dft_dencvg(g_dens,g_denso,nbf)
C$Id$
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
      integer me,nproc,i,j
      integer ldm,idm,ldmo,idmo
      integer ilo,jlo,ihi,jhi,numi,numj
      double precision rms
c
      call ga_sync
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
c--   >  Evaluate change in energy.
c     
c--   >  Check for density matrix convergence.
c     

      rms = 0.0d0
#if 1
      call ga_distribution(g_dens,
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. ilo.le.ihi) then
         numi =  ihi-ilo+1
         numj =  jhi-jlo+1
         if (numi.gt.0 .and. numj.gt.0) then
            if (.not. MA_Push_Get(MT_Dbl,numi*numj,'dmat',ldm, idm))
     &           call errquit('dft_dencvg: failed to alloc tmpm',0,0)
            if (.not. MA_Push_Get(MT_Dbl, numi*numj,'dmato',ldmo,idmo))
     &           call errquit('dft_dencvg: failed to alloc tmpm',0,0)
            call ga_get(g_dens,ilo,ihi,jlo,jhi,
     .           dbl_mb(idm),numi)
            call ga_get(g_denso,ilo,ihi,jlo,jhi,
     .           dbl_mb(idmo),numi)
            do j=0,numi*numj-1
               rms=rms+(DBL_MB(idm+j)-DBL_MB(idmo+j))**2
            enddo
            if(.not.ma_chop_stack(ldm))
     &           call errquit('dft_dencvg: cannot pop stack',0, MA_ERR)
            endif
         endif
#else
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
#endif
      call ga_dgop(Msg_rms_dens,rms,1,'+')
      dft_dencvg = rms/dble(nbf*nbf)

      return
      end
