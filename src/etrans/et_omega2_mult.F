       subroutine omega2_mult(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     &                   ish,jsh,ksh,lsh,G,ng,gJ,omega2)
       implicit none
c
c  this routine multiplies the integrals and density product, G,
c  and returns the contribution to the 2e part of Vab.
c
c  variables are defined in et_dens_mult.F
c
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
c
       integer ng
       integer i,j,k,l
       integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi 
       integer ish, jsh, ksh, lsh
       integer maxj, maxk, maxl
       double precision gJ(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
       double precision G(ng)
       double precision omega2
c
c  initialize omega2 and ng
c  =========================
c
       omega2 = 0.0d0
       ng = 0
c
c  loop over bf indices
c  ===================
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
             ng = ng + 1
c
               omega2 = omega2 + gJ(l,k,j,i)*G(ng)
c
               enddo
            enddo
         enddo
      enddo
c
 1111 format (2(F20.16),4x,2(F20.16))
       return
       end


c $Id$
