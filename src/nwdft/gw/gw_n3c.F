      double precision function gw_n3c(pars)
      implicit none

#include "errquit.fh"
#include "cdft.fh"
#include "bas.fh"
#include "basP.fh"
#include "global.fh"
#include "schwarz.fh"
#include "stdio.fh"
#include "util.fh"
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
      character(len=8), parameter :: pname = 'gw_n3c: '
c
      integer ao_han
      integer nsh_ao
      integer ishc,nshbfc,ilastc,ifirstc
      integer ishd,nshbfd,ilastd,ifirstd
      double precision ischw, eri_est, n3c, n3csym
c
      ao_han = pars%ao_bas_han
c
      do ishc=1,pars%nshells_ao
        if (.not.bas_cn2bfr(ao_han,ishc,ifirstc,ilastc))
     $    call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
        nshbfc = ilastc - ifirstc + 1
c
        do ishd=1,ishc
c
          eri_est = schwarz_shell(ishc,ishd)
          if (eri_est.lt.pars%tol2e*0.001d0) cycle
c
          if (.not.bas_cn2bfr(ao_han,ishd,ifirstd,ilastd))
     $      call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
          nshbfd = ilastd - ifirstd + 1
c
          if (ishc.eq.ishd) then
            ischw = ischw + dble(nshbfd*(nshbfd+1)/2)
          else
            ischw = ischw + dble(nshbfd*nshbfc)
          endif
        enddo
      enddo
c
      ischw = ischw*pars%nbf_ri
c
      n3c = dble(pars%nbf_ao**2 * pars%nbf_ri)
      n3csym = (n3c+dble(pars%nbf_ao*pars%nbf_ri))/2.0
c
      write(luout,*)
      write(luout,9000) n3c, n3csym, ischw
c
 9000 format(10x,'  Maximum number of 3-center AO ERIs:',f17.0,
     $     /,10x,'Reduced by permutational symmetry to:',f17.0,
     $     /,10x,'     Reduced by Schwarz screening to:',f17.0)
c
      return
      end
