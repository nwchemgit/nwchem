      Subroutine dft_getvc( rtdb,geom,AO_bas_han,CD_bas_han,
     &                     CD_coef,iga_vc,Ecoul,nbf_cd,
     &                     iVcoul_opt,iga_2ceri,iga_dens,IOLGC,itol2e)

C$Id: dft_getvc.F,v 1.2 1995-04-26 18:09:51 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      logical IOLGC

c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han
c**** 
#include <bas.fh>
#include <rtdb.fh>
#include <stdio.fh>
#include <mafdecls.h>
#include <global.fh>
#include <tcgmsg.fh>

c
c     Charge Density Expansion Coefficients
c
      Dimension CD_coef(nbf_cd),iga_dens(2)

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
        tol2e=10.d0**(-itol2e)
        call ga_zero(iga_vc)
        call rhf_fock_2e(geom, AO_bas_han, iga_dens(1), iga_vc, tol2e,
     &       .true., .false., .false.)
        call ga_sync
c
c       Compute the Coulomb energy.
c
        Ecoul = 0.5d0*ga_ddot(iga_dens(1),iga_vc)

      else if( iVcoul_opt.eq.1 )then
c       Compute the Coulomb potential via a LSQ fitting procedure.
        call dft_fitvc( rtdb, geom, AO_bas_han, CD_bas_han,
     &       CD_coef,iga_vc,Ecoul2,nbf_cd,iga_2ceri,IOLGC,itol2e)
c
c       Compute the Coulomb energy.
c
c        Ecoul1 = 0.5D0*ga_ddot(iga_dens,iga_vc)
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
