      subroutine gw_analytic_graph(pars,Sigmax,Sigma,Eold,Enew,Omega,
     $                             efermi,nri,maxpoles)

      implicit none
#include "errquit.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri, maxpoles
      double precision,intent(in) :: efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Enew(pars%nmo,pars%ipol)
      double precision,intent(in) :: Omega(maxpoles,pars%ipol)


      end
