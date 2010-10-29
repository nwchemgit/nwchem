c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 1e- derivative memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_dd1emem(bases,nbas,maxang)
      implicit none
#include "apiP.fh"
#include "nwc_const.fh"
#include "bas.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "hnd_pointers.fh"
c
      logical api_is_atomic_basis
      external api_is_atomic_basis
c
      integer nbas          ! [input] number of basis sets to be used
      integer bases(nbas)   ! [input] basis set handles
c
      integer scrmem, maxbuf, ibas, nbfmax, nbf, ll
      integer basis, geom, nat, maxang
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
c     Largest scratch space needed comes from hnd_tvd_ij.F:
c     3(L+3)(L+5)+(3(L+3)(L+3)+12(L+1)^2 for T 1st der
c     3(L+3)^2(L+2)+6(L+1)^2(L+2)        for V 1st der
c     18(L+1)^2                          for T 2nd der
c     9(L+1)^2(L+2)                      for V 2nd der
c
c     We are also setting up the pointers into the scratch block
c     We still need 1st der pointers
c
      tvp (1) = 1                                 ! T (L+3)(L+5)
      tvp (2) = (maxang+3)*(maxang+5)+1           ! T (L+3)(L+5)
      tvp (3) = 2*tvp (2)+1                       ! T (L+3)(L+5)
      tvp (4) = 3*tvp (2)+1                       ! T (L+3)(L+3)
      tvp (5) = tvp (4)+(maxang+3)**2             ! T (L+3)(L+3)
      tvp (6) = tvp (5)+(maxang+3)**2             ! T (L+3)(L+3)
      tvp (7) = tvp (6)+(maxang+3)**2             ! T (L+1)(L+3)
      do ll = 8, 19
         tvp (ll) = tvp (ll-1)+(maxang+1)**2      ! T 11*(L+1)(L+3) + V (L+2)(L+3)^2
      enddo
      tvp (20) = tvp(19)+(maxang+3)**2*(maxang+2) ! V (L+2)(L+3)^2
      tvp (21) = tvp(20)+(maxang+3)**2*(maxang+2) ! V (L+2)(L+3)^2
      tvp (22) = tvp(21)+(maxang+3)**2*(maxang+2) ! V (L+2)(L+1)^2
      do ll = 23, 28                              ! V 5*(L+2)(L+1)^2 + T 2nd (L+1)^2
         tvp(ll) = tvp(ll-1)+(maxang+1)**2*(maxang+2)
      enddo
c
c     And now we define the 2nd der pointers
c
      do ll = 29, 46                             
         tvp(ll) = tvp(ll-1)+(maxang+1)**2        ! V 2nd 17*(L+1)^2 + (L+3)^2(L+2)
      enddo
      do ll = 47, 54                              ! V 2nd 8*(L+1)^2(L+2)
         tvp(ll) = tvp(ll-1)+(maxang+1)**2*(maxang+2)
      enddo
      scrmem = tvp(54)+(maxang+3)**2*(maxang+2)
      mem_1e = max(mem_1e,scrmem)
c
c     Pointers defined for hnd_tvd_ij.F can and are reused by hnd_dipder.F, 
c     hnd_hlf_ij.F, and hnd_sd_ij.F
c
c     hnd_hlfd_ij.F needs its own pointers as its blocks don't match those
c     here
c
      hdp(1) = 1
      do ll = 2, 7                                ! 5*(L+2)^2(L+2) + (L+2)^2(L+2)
         hdp(ll) = hdp(ll-1)+(maxang+2)**2*(maxang+2)
      enddo
      do ll = 8, 18                               ! 11*(L+1)^2(L+2)
         hdp(ll) = hdp(ll-1)+(maxang+1)**2*(maxang+2)
      enddo
c
      return
      end
