      subroutine gw_analytic_rpa(pars,eia,omega,maxpoles,maxnpoles,nri)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "errquit.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: maxpoles, maxnpoles, nri

      double precision, intent(in)  :: eia(maxpoles,*)
      double precision, intent(out) :: omega(maxnpoles,*)

      integer g_rpa, g_xpy
      integer npoles,mynpoles
      integer kRPA,kOV1,kOV2,kXPY,ld

      character(len=17),parameter :: pname = 'gw_analytic_rpa: '

      logical,external :: ga_create, ga_destroy


      call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $               kOV1,ld)
      if (pars%ipol.gt.1) then
        call ga_access(pars%g_eriov(2),1,nri,pars%ovlo(2),pars%ovhi(2),
     $                 kOV2,ld)
      else
        kOV2 = 1
      endif

      if (pars%ipol.eq.1) then

        npoles = pars%npoles(1)
        mynpoles = pars%mynpoles(1)

        if (.not.ga_create(mt_dbl,npoles,npoles,'rpa',0,npoles,g_rpa))
     $    call errquit(pname//'can''t create array',0,GA_ERR)
        if (.not.ga_create(mt_dbl,npoles,npoles,'x+y',0,npoles,g_xpy))
     $    call errquit(pname//'can''t create array',0,GA_ERR) 

        call ga_access(g_rpa,pars%ovlo(1),pars%ovhi(1),1,npoles,
     $                 kRPA,ld)   
        call ga_access(g_xpy,pars%ovlo(1),pars%ovhi(1),1,npoles,
     $                 kXPY,ld)   

        call gw_analytic_ovscale(1,dbl_mb(kOV1),eia(1,1),nri,mynpoles)
        call ga_dgemm('t','n',pars%npoles(1),pars%npoles(1),nri,4d0,
     $                 pars%g_eriov(1),pars%g_eriov(1),0d0,g_rpa)
        call gw_analytic_addrpa(dbl_mb(kRPA),eia(1,1),mynpoles,
     $                          pars%ovlo(1))

        call ga_diag_std(g_rpa,g_xpy,omega)
        omega(:,1) = dsqrt(omega(:,1))


        call gw_analytic_scalxpy(dbl_mb(kXPY),omega,mynpoles,npoles)
        call ga_dgemm('n','n',nri,npoles,npoles,1d0,pars%g_eriov(1),
     $                 g_xpy,0d0,pars%g_erim(1))
        call gw_analytic_ovscale(2,dbl_mb(kOV1),eia(1,1),nri,mynpoles)

        if (.not.(ga_destroy(g_rpa) .and. ga_destroy(g_xpy)))
     $    call errquit(pname//'can''t destroy arrays',0,GA_ERR)


      else

      endif

      end subroutine
