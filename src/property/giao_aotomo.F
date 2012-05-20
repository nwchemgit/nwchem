c
c     Transform the matrix/ces in g_mat from AO to MO basis
c
      subroutine giao_aotomo(gmat,gvecs,nocc,nvirt,ndens,nblk,nbf)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer ndens, nblk, nbf
      integer gmat, gvecs(ndens), nocc(ndens), nvirt(ndens)
c
      integer gtmp
      integer ind, nmo, in, ndns
      integer alo(3),ahi(3),blo(3),bhi(3),clo(3),chi(3),dlo(2),dhi(2)
c
      if (.not.ga_create(MT_DBL,nbf,nbf,'TMP',-1,-1,gtmp))
     &    call errquit('giao_aotomo: ga_create failed gtmp',0,GA_ERR)
c
      alo(1) = 1
      alo(2) = 1
      ahi(1) = nbf
      dlo(1) = 1
      dlo(2) = 1
      dhi(2) = nbf
      blo(1) = 1
      blo(2) = 1
      bhi(1) = nbf
      bhi(2) = nbf
      clo(1) = 1
      clo(2) = 1
c
      ind = 0
      do ndns = 1, ndens
        nmo = nocc(ndns) + nvirt(ndns)
        ahi(2) = nmo  
        dhi(1) = nmo
        chi(1) = nmo
        chi(2) = nmo
        do in = 1, nblk 
          ind = ind + 1
          blo(3) = ind
          bhi(3) = ind
          clo(3) = ind
          chi(3) = ind
          call ga_zero(gtmp)
          call nga_matmul_patch ('t','n', 1.D0, 0.D0,
     $                           gvecs(ndns), dlo, dhi,
     $                           gmat       , blo, bhi,
     $                           gtmp       , dlo, dhi)
          call nga_matmul_patch ('n','n', 1.D0, 0.D0,
     $                           gtmp       , dlo, dhi,
     $                           gvecs(ndns), alo, ahi,
     $                           gmat       , clo, chi)
c
        enddo
      enddo
c
      if (.not.ga_destroy(gtmp)) call 
     &    errquit('giao_aotomo: ga_destroy failed gtmp',0,GA_ERR)
c
      return
      end
c $Id$
