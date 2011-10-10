       subroutine dens_prod(nbf,ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &                     ish,jsh,ksh,lsh, norm,GG,ng,blen,
     &                     pa_ij,pa_kl,pa_ji,pa_lk,pa_il,pa_kj,
     &                     pa_jl,pa_ki,pa_ik,pa_lj,pa_jk,pa_li, 
     &                     pb_ij,pb_kl,pb_ji,pb_lk,pb_il,pb_kj,
     &                     pb_jl,pb_ki,pb_ik,pb_lj,pb_jk,pb_li) 
       implicit none
c
c  this routine multiplies the densities and forms the array GG
c  (a product of densities) and calculates the maximum value, which
c  is used as a screening parameter.
c
c      ilo,ihi,jlo,jhi,klo,khi,llo,lhi - bf index ranges [input]
c      ish,jsh,ksh,lsh - shell indixes                   [input]
c      pa,pb - alpha and beta densities                  [input]
c      norm                                              [output]
c      GG - densit products                              [output]
c      ng - number of elements in GG                     [output]
c 
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
c
       integer nbf,ng
       integer i,j,k,l, ii,jj,kk,ll
       integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi 
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl, blen
       double precision Gijkl,Gjikl,Gijlk,Gjilk
       double precision Gilkj,Giklj,Gjlki,Gjkli
       double precision GGijkl,GGilkj
       double precision dij,dkl,dijkl
       double precision delta
       double precision GG(10000), dabsmax, norm
       double precision pa_ij(blen,blen),pa_kl(blen,blen)
       double precision pa_ji(blen,blen),pa_lk(blen,blen)
       double precision pa_il(blen,blen),pa_kj(blen,blen)
       double precision pa_jl(blen,blen),pa_ki(blen,blen)
       double precision pa_ik(blen,blen),pa_lj(blen,blen)
       double precision pa_jk(blen,blen),pa_li(blen,blen)
       double precision pb_ij(blen,blen),pb_kl(blen,blen) 
       double precision pb_ji(blen,blen),pb_lk(blen,blen)
       double precision pb_il(blen,blen),pb_kj(blen,blen)
       double precision pb_jl(blen,blen),pb_ki(blen,blen)
       double precision pb_ik(blen,blen),pb_lj(blen,blen)
       double precision pb_jk(blen,blen),pb_li(blen,blen)
       external dabsmax
c
c
c  initialize ng
c  =============
c
       ng = 0
c
c  loop over bf indices
c  ====================
c
      do i = ilo, ihi
c
          maxj = jhi
          if(jsh.eq.ish) maxj = i
          do j = jlo, maxj
c
            maxk = khi
            if(ksh.eq.ish.AND.jsh.eq.lsh) maxk = i
            do k = klo, maxk
c
              maxl = lhi
              if(lsh.eq.ksh) maxl = k
              if(ish.eq.ksh.AND.jsh.eq.lsh.AND.k.eq.i) maxl=j
              do l = llo, maxl
c
              ii = i - ilo+1
              jj = j - jlo+1
              kk = k - klo+1
              ll = l - llo+1
c
                Gijkl =   pa_ij(ii,jj)*pa_kl(kk,ll)
     &                  + pa_ij(ii,jj)*pb_kl(kk,ll)
     &                  + pb_ij(ii,jj)*pa_kl(kk,ll)
     &                  + pb_ij(ii,jj)*pb_kl(kk,ll)
c
                Gjikl =   pa_ji(jj,ii)*pa_kl(kk,ll)
     &                  + pa_ji(jj,ii)*pb_kl(kk,ll)
     &                  + pb_ji(jj,ii)*pa_kl(kk,ll)
     &                  + pb_ji(jj,ii)*pb_kl(kk,ll)
c
                Gijlk =   pa_ij(ii,jj)*pa_lk(ll,kk)
     &                  + pa_ij(ii,jj)*pb_lk(ll,kk)
     &                  + pb_ij(ii,jj)*pa_lk(ll,kk)
     &                  + pb_ij(ii,jj)*pb_lk(ll,kk)
c
                Gjilk =   pa_ji(jj,ii)*pa_lk(ll,kk)
     &                  + pa_ji(jj,ii)*pb_lk(ll,kk)
     &                  + pb_ji(jj,ii)*pa_lk(ll,kk)
     &                  + pb_ji(jj,ii)*pb_lk(ll,kk)
c
                Gilkj = - pa_il(ii,ll)*pa_kj(kk,jj)
     &                  - pb_il(ii,ll)*pb_kj(kk,jj)
c
                Gjlki = - pa_jl(jj,ll)*pa_ki(kk,ii)
     &                  - pb_jl(jj,ll)*pb_ki(kk,ii)
c
                Giklj = - pa_ik(ii,kk)*pa_lj(ll,jj)
     &                  - pb_ik(ii,kk)*pb_lj(ll,jj)
c
                Gjkli = - pa_jk(jj,kk)*pa_li(ll,ii)
     &                  - pb_jk(jj,kk)*pb_li(ll,ii)
c
                dij = 0.0d0 
                if (i.eq.j) dij = 1.0d0
c
                dkl = 0.0d0 
                if (k.eq.l) dkl = 1.0d0
c
                dijkl = 0.0d0
                if (i.eq.k.AND.j.eq.l) dijkl = 1.0d0
c
                delta = (1.0d0+dij)*(1.0d0+dkl)
c
                GGijkl = Gijkl + Gjikl + Gijlk + Gjilk
                GGilkj = Gjlki + Gilkj + Gjkli + Giklj
c
                ng = ng + 1
                GG(ng)  = 2.0d0*((GGijkl + GGilkj)/delta)/(1.0d0+dijkl)
c
               enddo
            enddo
         enddo
      enddo
c
      norm = dabsmax(ng,GG)
c
      return
      end


c $Id$
