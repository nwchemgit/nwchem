      subroutine bse_wmn(pars,w,nmo,nri,ipol)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bse.fh"

      type(bse_params_t) :: pars
      integer nmo, nri, ipol

      character(*), parameter :: pname = 'bse_wmn: '

      double precision w(*)

      integer isp,ilo,ihi,jlo,jhi,keri,ld

      do isp=1,ipol
        ! copy and transform occ-vir block
        if(.not.ga_duplicate(pars%g_eriov(isp),pars%g_wov(isp),'Wov'))
     &    call errquit(pname//'could not duplicate ERI OV',0,GA_ERR)
        call ga_copy(pars%g_eriov(isp),pars%g_wov(isp))
        if(pars%mynpoles(isp).eq.0) cycle
        call ga_access(pars%g_wov(isp),1,nri,pars%ovlo(isp),
     &                 pars%ovhi(isp),keri,ld)
        call dtfsm('n','l','l','n','n',nri,pars%mynpoles(isp),1d0,
     &              w,dbl_mb(keri),nri)
        call ga_release_update(pars%g_wov(isp),1,nri,pars%ovlo(isp),
     &                         pars%ovhi(isp))    

        ! transform occ-occ and vir-vir blocks in place
        call ga_distribution(pars%g_erioo(isp),pars%me,ilo,ihi,jlo,jhi)
        if (jlo.gt.jhi) goto 101
        call ga_access(pars%g_erioo(isp),ilo,ihi,jlo,jhi,keri,ld)
        call dtfsm('n','l','l','n','n',nri,jhi-jlo+1,1.0d0,w,
     &              dbl_mb(keri),nri)
        call ga_release_update(pars%g_erioo(isp),1,nri,jlo,jhi)

  101   continue
        call ga_distribution(pars%g_erivv(isp),pars%me,ilo,ihi,jlo,jhi)
        if (jlo.gt.jhi) cycle
        call ga_access(pars%g_erivv(isp),ilo,ihi,jlo,jhi,keri,ld)
        call dtfsm('n','l','l','n','n',nri,jhi-jlo+1,1.0d0,w,
     &              dbl_mb(keri),nri)
        call ga_release_update(pars%g_erivv(isp),1,nri,jlo,jhi)
      enddo   

      end
