       subroutine dens_prod(nbf,ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &                   ish,jsh,ksh,lsh,pa,pb,norm,GG,ng)
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
       integer i,j,k,l
       integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi 
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl
       double precision pa(nbf,nbf)
       double precision pb(nbf,nbf)
       double precision Gijkl,Gjikl,Gijlk,Gjilk
       double precision Gilkj,Giklj,Gjlki,Gjkli
       double precision GGijkl,GGilkj
       double precision dij,dkl,dijkl
       double precision delta
       double precision GG(10000), dabsmax, norm
       external dabsmax
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
                Gijkl =   pa(i,j)*pa(k,l) 
     &                  + pa(i,j)*pb(k,l)
     &                  + pb(i,j)*pa(k,l)
     &                  + pb(i,j)*pb(k,l)
c
                Gjikl =   pa(j,i)*pa(k,l) 
     &                  + pa(j,i)*pb(k,l)
     &                  + pb(j,i)*pa(k,l)
     &                  + pb(j,i)*pb(k,l)
c
                Gijlk =   pa(i,j)*pa(l,k) 
     &                  + pa(i,j)*pb(l,k)
     &                  + pb(i,j)*pa(l,k)
     &                  + pb(i,j)*pb(l,k)
c
                Gjilk =   pa(j,i)*pa(l,k) 
     &                  + pa(j,i)*pb(l,k)
     &                  + pb(j,i)*pa(l,k)
     &                  + pb(j,i)*pb(l,k)
c
c
                Gilkj = - pa(i,l)*pa(k,j) 
     &                  - pb(i,l)*pb(k,j)
c
                Gjlki = - pa(j,l)*pa(k,i) 
     &                  - pb(j,l)*pb(k,i)
c
                Giklj = - pa(i,k)*pa(l,j) 
     &                  - pb(i,k)*pb(l,j)
c
                Gjkli = - pa(j,k)*pa(l,i) 
     &                  - pb(j,k)*pb(l,i)
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
                GGijkl = GGijkl/delta
                GGilkj = GGilkj/delta
c
                ng = ng + 1
                GG(ng)  =          2.0d0*GGijkl/(1.0d0+dijkl)
                GG(ng) = GG(ng) +  2.0d0*GGilkj/(1.0d0+dijkl)
c
               enddo
            enddo
         enddo
      enddo
c
      norm = dabsmax(ng,GG)
c
 1111 format (a31,1x,4F16.12)
       return
       end


