c $Id: hnd_dd1emem.F,v 1.2 1999-05-06 17:25:49 windus Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 1e- derivative memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_dd1emem(bases,nbas)
      implicit none
#include "apiP.fh"
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
      integer scrmem, maxbuf, ibas, jbas, nbfmax, nbf
      integer basis, geom, nat
c
c     set memory requirements for 1 e- 2nd derivative integrals
c     The maxbuf size is the size of a shell of integrals
c     The scrmem is the size of the scratch area needed.  This will
c     need to be fixed when I get to passing in temporary space instead
c     allocating everything by dimension (statically in executable)
c
c     scrmem = 3*3*4*255
c     maxbuf = scrmem*2
      scrmem = 0
      maxbuf = 0
      nbfmax = 0
c
c get the number of atoms, note that this does not do any checking
c to see if the different basis sets are on different geometries
c
      write(6,*) 'before the setup stuff'
      call flush(6)
      basis = bases(1) + BASIS_HANDLE_OFFSET
      write(6,*) 'after the basis definition'
      call flush(6)
      geom = ibs_geom(basis)
      write(6,*) 'after the geom definition'
      call flush(6)
      nat = ncenter(geom)
      write(6,*) 'after the nat definition'
      call flush(6)
c
c figure out the maximum number of basis functions for any shell
c
      write(6,*) 'number of basis sets ',nbas
      call flush(6)
      do 00100 ibas = 1,nbas
        if (api_is_atomic_basis(bases(ibas))) then
          write(6,*) 'before the call to int_nbf_max'
          call flush(6)
          call int_nbf_max(bases(ibas),nbf)
          write(6,*) 'after the call to int_nbf_max'
          call flush(6)
          nbfmax = max(nbfmax,nbf)
          write(6,*) 'nbfmax ',nbfmax
          call flush(6)
        endif
00100 continue
c
c The 4 should be changed to a 3 when I have everything worked out
c
      maxbuf = nbfmax*nbfmax*3*3*4*nat
      isz_1e = max(isz_1e,maxbuf)
c
c scrmem is for the didjvij temporary space and for holding the T second
c ders during the calculation of H1
c
      scrmem = nbfmax*nbfmax*3*3*nat + nbfmax*nbfmax*3*3*3
      mem_1e = max(mem_1e,scrmem) 
c
      return
      end
