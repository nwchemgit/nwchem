      double precision function dft_dencvg(iga_dens,igdenso,nbf)
C$Id: dft_dencvg.F,v 1.2 1995-10-17 08:54:52 d3g681 Exp $
      implicit double precision(a-h,o-z)
      logical LResult
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
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
      LResult = MA_Push_Get(MT_Dbl,nbf,'dmat',ldm,idm)
      LResult = MA_Push_Get(MT_Dbl,nbf,'dmatold',ldmo,idmo)
      do i=1+me,nbf,nproc
        call get_col(iga_dens,nbf,i,Dbl_MB(idm))
        call get_col(igdenso,nbf,i,Dbl_MB(idmo))
        do j=0,nbf-1
          rms=rms+(DBL_MB(idm+j)-DBL_MB(idmo+j))**2
        enddo
      enddo
      LResult = MA_Pop_Stack(ldmo)
      LResult = MA_Pop_Stack(ldm)
      call ga_sync
      call ga_dgop(Msg_rms_dens,rms,1,'+')
      dft_dencvg = rms/dfloat(nbf*nbf)

      return
      end
