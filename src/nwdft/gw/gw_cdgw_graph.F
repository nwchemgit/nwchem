      subroutine gw_cdgw_graph(pars,Sigmax,Sigma,Eold,Enew,Eia,iW,
     $                         iWidx,efermi,nri,mynmo,maxqp,maxpoles)

      implicit none
#include "errquit.fh"
#include "gw.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: nri, mynmo, maxqp, maxpoles
      integer,intent(in) :: iWidx(0:mynmo,maxqp,pars%ipol)
      double precision,intent(in) :: efermi(pars%ipol)
      double precision,intent(in) :: Sigmax(pars%nmo,pars%ipol)
      double precision,intent(in) :: Sigma(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eold(pars%nmo,pars%ipol)
      double precision,intent(in) :: Enew(pars%nmo,pars%ipol)
      double precision,intent(in) :: Eia(maxpoles,pars%ipol)
      double precision,intent(in) :: iW(pars%nmo,pars%ngl,maxqp,
     $                                  pars%ipol)


      end
