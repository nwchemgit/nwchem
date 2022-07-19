      Subroutine dft_getvc(CD_coef, i3c_ERI, Ecoul,  g_vc,
     &                     iVcoul_opt, n_batch,
     &                     n_bmax, iwhat, n3c_dbl,
     &                     IOLGC, 
     &                     derfit, nmat)
c
C$Id$
      implicit none
#include "errquit.fh"
c
#include "bas.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "util.fh"
c
      integer g_vc(*)
      integer iVcoul_opt
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl
      integer nmat           ! number of fitted J matrices to make [in]
      logical IOLGC, derfit
c
      double precision CD_coef(nbf_cd,nmat)
      double precision i3c_ERI(n3c_dbl)
      double precision Ecoul
      double precision Ecoul2
c     
c     This driver routine solves for the Coulomb potential (Vcoul) by one of the
c     following methods:
c     
c     1) evaluating the set of 4-ctr ERIs formed from the AO basis set,
c     2) fitting the charge density via Dunlap method and, then, evaluating 
c     the set of 3-ctr ERIs formed from the AO and CD basis sets,
c     3) solving Poisson's equation via (Delley,Becke,Feibelman,?)'s method. 
c     
      call ga_sync
c
      if( iVcoul_opt.eq.0 )then
        call ga_zero(g_vc(1))
        if (ipol.gt.1) call ga_zero(g_vc(2))
c
      else if( iVcoul_opt.eq.1 )then
c       Compute the Coulomb potential via a LSQ fitting procedure.

        call dft_fitvc(CD_coef, i3c_ERI, Ecoul2,  g_vc,
     &                 n_batch, n_bmax, iwhat, n3c_dbl,
     &                 IOLGC, derfit, nmat,
     &                 .true.)
c        call ga_sync
c
c       Compute the Coulomb energy.
c
        Ecoul = Ecoul - Ecoul2 
c
      else if( iVcoul_opt.eq.2 )then
c
        write(LuOut,*) 'GETVCOUL:  Poisson Eq. solver (iVcoul_opt.eq.2)'
        write(LuOut,*) '           has yet to be implemented.'
        call errquit('Exiting in getvcoul.',1, CALC_ERR)
c
      else
c
        write(LuOut,*) 'GETVCOUL:  Unknown option indicated by the'
        write(LuOut,*) '           variable iVcoul_opt = ',iVcoul_opt
        call errquit('Exiting in getvcoul.',2, CALC_ERR)
c
      end if
c
      return
      end

