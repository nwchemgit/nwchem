c $Id: hnd_dd1emem.F,v 1.5 1999-12-29 23:30:24 windus Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 1e- derivative memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_dd1emem(bases,nbas)
      implicit none
#include "apiP.fh"
#include "nwc_const.fh"
#include "bas.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
c
      logical api_is_atomic_basis
      external api_is_atomic_basis
c
      integer nbas          ! [input] number of basis sets to be used
      integer bases(nbas)   ! [input] basis set handles
c
      integer scrmem, maxbuf, ibas, nbfmax, nbf
      integer basis, geom, nat
c
c     set memory requirements for 1 e- 2nd derivative integrals
c     The maxbuf size is the size of a shell of integrals
c     The scrmem is the size of the scratch area needed.
c
      scrmem = 0
      maxbuf = 0
      nbfmax = 0
c
c get the number of atoms, note that this does not do any checking
c to see if the different basis sets are on different geometries
c
      basis = bases(1) + BASIS_HANDLE_OFFSET
      geom = ibs_geom(basis)
      nat = ncenter(geom)
c
c figure out the maximum number of basis functions for any shell
c
      do 00100 ibas = 1,nbas
        if (api_is_atomic_basis(bases(ibas))) then
          call int_nbf_max(bases(ibas),nbf)
          nbfmax = max(nbfmax,nbf)
        endif
00100 continue
c
c This is what I would like it to be, but for testing....
c The extra 3 holds the T terms and the nonHelFey terms
c     maxbuf = nbfmax*nbfmax*3*3*3*nat
      maxbuf = nbfmax*nbfmax*3*3*(3*nat+3)
      isz_1e = max(isz_1e,maxbuf)
c
c scrmem should eventually either go away or be used for other
c arrays that are now statically allocated
c
      scrmem = nbfmax*nbfmax*3*3*4
      mem_1e = max(mem_1e,scrmem) 
c
      return
      end
