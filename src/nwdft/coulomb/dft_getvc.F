      Subroutine dft_getvc(CD_coef, i3c_ERI, Ecoul, rtdb, g_vc,
     &                     iVcoul_opt, g_dens, n_batch,
     &                     n_bmax, iwhat, n3c_dbl, iwhat_max, 
     &                     n_semi_bufs, fd, IOLGC, rdens_atom,
     &                     cntoce, cntobfr, cetobfr, natoms,
     &                     derfit, nmat)
c
C$Id: dft_getvc.F,v 1.19 1999-05-10 22:43:45 edo Exp $
      implicit none
c
#include "bas.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "util.fh"
c
      integer natoms
      integer cntoce(nshells_ao), cntobfr(2,nshells_ao), 
     &        cetobfr(2,natoms)
      double precision rdens_atom(natoms,natoms,ipol)
      integer rtdb
      integer g_vc
      integer iVcoul_opt
      integer g_dens(2)
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl
      integer iwhat_max
      integer n_semi_bufs
      integer fd
      integer nmat           ! number of fitted J matrices to make [in]
      logical IOLGC, derfit
c
      double precision CD_coef(nbf_cd,nmat)
      double precision i3c_ERI(n3c_dbl)
      double precision Ecoul
      double precision Ecoul2
      double precision tol2e
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
c$$$c
c$$$c       Compute the Coulomb potential via evaluation of the AO set of 4-ctr ERIs.
c$$$c
c$$$        tol2e=10.d0**(-itol2e)
         call ga_zero(g_vc)
c$$$        call rhf_fock_2e(geom, AO_bas_han, g_dens(1), g_vc, tol2e,
c$$$     &       .true., .false., oskel)
c$$$        call ga_sync
c$$$c
c$$$c       Symmetrize Vcoul?
c$$$c
c$$$c        if (oskel)
c$$$c     &     call sym_symmetrize(geom, AO_bas_han, .false., g_vc)
c$$$c
c$$$c       Compute the Coulomb energy.
c$$$c
c$$$        Ecoul = 0.5d0*ga_ddot(g_dens(1),g_vc)
            Ecoul = 0.0d0
c
      else if( iVcoul_opt.eq.1 )then
c       Compute the Coulomb potential via a LSQ fitting procedure.

c        call dft_tstamp(' Before call to dft_fitvc')
c        if (.not.ma_verify_allocator_stuff()) stop ' ma ck 7'
        call dft_fitvc(CD_coef, i3c_ERI, Ecoul2, rtdb, g_vc,
     &                 n_batch, n_bmax, iwhat, n3c_dbl,
     &                 iwhat_max, n_semi_bufs, fd, IOLGC, rdens_atom,
     &                 cntoce, cntobfr, cetobfr, natoms, derfit, nmat,
     &                 .true.)
c        call dft_tstamp(' After  call to dft_fitvc')
c        if (.not.ma_verify_allocator_stuff()) stop ' ma ck 7'
c
        call ga_sync
c
c       Symmetrize Vcoul?
c
c        if (oskel)
c     &     call sym_symmetrize(geom, AO_bas_han, .false., g_vc)
c
c       Compute the Coulomb energy.
c
        Ecoul = Ecoul - Ecoul2 
c
      else if( iVcoul_opt.eq.2 )then
c
        write(LuOut,*) 'GETVCOUL:  Poisson Eq. solver (iVcoul_opt.eq.2)'
        write(LuOut,*) '           has yet to be implemented.'
        call errquit('Exiting in getvcoul.',1)
c
      else
c
        write(LuOut,*) 'GETVCOUL:  Unknown option indicated by the'
        write(LuOut,*) '           variable iVcoul_opt = ',iVcoul_opt
        call errquit('Exiting in getvcoul.',2)
c
      end if
c
      return
      end

