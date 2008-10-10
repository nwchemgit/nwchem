c#define NBACC 1
c
      subroutine shell_fock_buildxc(geom,basis, nder, ndens,
     $     jfac, kfac, tol2e, g_dens, g_fock, osym)
c
c     == adapted from shell_fock_build ==
c     == just builds the exchange-correlation part of the fock matrix ==

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"

      integer basis, nder, ndens, nfock
      integer geom
      integer g_dens, g_fock
      double precision jfac(*), kfac(*), tol2e

      integer nsh, nbf, max_at_bf, max_sh_bf, nshblocks
      integer ityp, ndim, dims(3)
c
      logical osym
      logical xc_gotxc
      external xc_gotxc
c
c     == Get the dimensions of g_fock and eventually g_dens. ==
      call nga_inquire(g_fock, ityp, ndim, dims)
      nfock = dims(1)
c
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('new_fock: could not get nbf',0, BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('new_fock: could not get nsh',0, BASIS_ERR)
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('new_fock: could not get max_at_bf',0,
     &       BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('new_fock: could not get max_sh_bf',0,
     &       BASIS_ERR)
c
c     == Get the exchange-correlation contribution ==
      if(xc_gotxc().and.nder.eq.0)then
         call xc_newfock(geom,basis,jfac,kfac,nbf,
     &        nder, ndens, nfock, g_fock,g_dens)
      endif
c
      if (osym) then
         call ga_symmetrize(g_fock)
      else
         call ga_antisymmetrize(g_fock)
      endif
c
      return
      end
