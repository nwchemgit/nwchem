      Subroutine dft_getvc( rtdb,geom,AO_bas_han,CD_bas_han,
     &                     CD_coef,g_vc,Ecoul,nbf_cd,
     &                     iVcoul_opt,g_2ceri,g_dens,IOLGC,itol2e)

C$Id: dft_getvc.F,v 1.5 1995-12-27 23:33:18 og845 Exp $
      implicit none
      

c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han
      integer nbf_cd
      integer g_vc
      integer g_2ceri
      integer g_dens(2)
      integer itol2e
      double precision Ecoul
      integer iVcoul_opt
c
c     Charge Density Expansion Coefficients
c
      double precision CD_coef(nbf_cd)
#include "bas.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"

      logical IOLGC
      double precision tol2e,Ecoul2

c     
c     This driver routine solves for the Coulomb potential (Vcoul) by one of the
c     following methods:
c     
c     1) evaluating the set of 4-ctr ERIs formed from the AO basis set,
c     2) fitting the charge density via Dunlap method and, then, evaluating 
c     the set of 3-ctr ERIs formed from the AO and CD basis sets,
c     3) solving Poisson's equation via (Delley,Becke,Feibelman,?)'s method. 
c     
c******************************************************************************
      if( iVcoul_opt.eq.0 )then

c       Compute the Coulomb potential via evaluation of the AO set of 4-ctr ERIs.
        tol2e=10.d0**(-itol2e)*100.
        call ga_zero(g_vc)
        call rhf_fock_2e(geom, AO_bas_han, g_dens(1), g_vc, tol2e,
     &       .true., .false., .false.)
        call ga_sync
c
c       Compute the Coulomb energy.
c
        Ecoul = 0.5d0*ga_ddot(g_dens(1),g_vc)

      else if( iVcoul_opt.eq.1 )then
c       Compute the Coulomb potential via a LSQ fitting procedure.
        call dft_fitvc( rtdb, geom, AO_bas_han, CD_bas_han,
     &       CD_coef,g_vc,Ecoul2,nbf_cd,g_2ceri,IOLGC,itol2e)
c
c       Compute the Coulomb energy.
c
        Ecoul = Ecoul - Ecoul2 

      else if( iVcoul_opt.eq.2 )then

        write(LuOut,*) 'GETVCOUL:  Poisson Eq. solver (iVcoul_opt.eq.2)'
        write(LuOut,*) '           has yet to be implemented.'
        call errquit('Exiting in getvcoul.',1)

      else

        write(LuOut,*) 'GETVCOUL:  Unknown option indicated by the'
        write(LuOut,*) '           variable iVcoul_opt = ',iVcoul_opt
        call errquit('Exiting in getvcoul.',2)

      end if

      end
