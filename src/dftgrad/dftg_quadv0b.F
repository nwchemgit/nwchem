      Subroutine dftg_quadv0_b(ictr, ldew, AO_bas_han, ipol, GRAD, 
     &                         nbf_ao, mbf_ao, nat, nq, chi, delchi, 
     &                         heschi, Amat, Cmat, ibf, iniz, ifin, 
     &                         PPP, F, zzz, hh, g_dens, force, 
     &                         oactive)
*
* $Id: dftg_quadv0b.F,v 1.6 2000-07-18 01:13:20 edo Exp $
*
      Implicit none
c
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
c     
      logical ldew   ! true if weight derivatives are included [input]
      integer ictr
      integer ipol, nbf_ao, mbf_ao, nat, nq
      Integer AO_bas_han
      logical GRAD, oactive(nat)
      double precision zzz(3,mbf_ao)
      double precision PPP(nbf_ao), F(mbf_ao)
      integer ibf(mbf_ao), iniz(nat), ifin(nat)
      integer g_dens(2)
c     
c     scratch array for hessian product
c     
      double precision hh(3,nq)
c     
c     functions and grad
c     
      double precision chi(nq,mbf_ao), delchi(nq,3,mbf_ao),
     &                 heschi(nq,6,mbf_ao)
c     
c     Vxc and its GCA part
c     
      double precision Amat(nq,ipol), Cmat(nq,3,ipol)
c     
c     force vector
c     
      double precision force(3,nat)
c
c     local
c
      logical ldo
      integer iat, n, ispin, mu, mmu, nu, ilast, ifirst
      double precision C1, C2, C3, z1, z2, z3, DC, fatto
      double precision out(3)
      double precision two
      parameter(two=2.d0)
c
      do iat = 1, nat
         ldo = .true.
         if (ldew)ldo = ldo.and.(iat.ne.ictr)
         if (ldo)then
            if (iniz(iat).ne.0)then
               if (.not.bas_ce2bfr(AO_bas_han,iat,ifirst,ilast))
     &            call errquit
     &            (' dftg_quadv0_b: number of basis with atom ?',1)
               do mu = iniz(iat),ifin(iat)
                  mmu = ibf(mu)
                  do ispin = 1, ipol
                     call ga_get(g_dens(ispin),1,nbf_ao,mmu,mmu,PPP,1)
                     call dft_gather(mbf_ao,F,PPP,ibf,0)
c               
c                    LDA part
c               
                     do n = 1, nq
                        fatto = Amat(n,ispin)
                        hh(1,n) = fatto*delchi(n,1,mu)
                        hh(2,n) = fatto*delchi(n,2,mu)
                        hh(3,n) = fatto*delchi(n,3,mu)
                     enddo
                     if (GRAD)then
                        do n = 1, nq
                           C1 = Cmat(n,1,ispin)
                           C2 = Cmat(n,2,ispin)
                           C3 = Cmat(n,3,ispin)
                           hh(1,n) = hh(1,n) + C1*heschi(n,1,mu) +
     &                                         C2*heschi(n,2,mu) +
     &                                         C3*heschi(n,3,mu)
                           hh(2,n) = hh(2,n) + C1*heschi(n,2,mu) +
     &                                         C2*heschi(n,4,mu) +
     &                                         C3*heschi(n,5,mu)
                           hh(3,n) = hh(3,n) + C1*heschi(n,3,mu) +
     &                                         C2*heschi(n,5,mu) +
     &                                         C3*heschi(n,6,mu)
                        enddo
                     endif
                     call dgemm('n','n',3,mbf_ao,nq,1.d0,hh,3,chi,
     &                          nq,0.d0,zzz,3)
                     if (GRAD)then
                        do nu = 1, mbf_ao
                           z1 = 0.d0
                           z2 = 0.d0
                           z3 = 0.d0
                           do n = 1, nq
c                  
c                             GCA part
c                  
                              DC = Cmat(n,1,ispin)*delchi(n,1,nu) +
     &                             Cmat(n,2,ispin)*delchi(n,2,nu) +
     &                             Cmat(n,3,ispin)*delchi(n,3,nu)
                              z1 = z1 + delchi(n,1,mu)*DC
                              z2 = z2 + delchi(n,2,mu)*DC
                              z3 = z3 + delchi(n,3,mu)*DC
                           enddo
                           zzz(1,nu) = zzz(1,nu) + z1
                           zzz(2,nu) = zzz(2,nu) + z2
                           zzz(3,nu) = zzz(3,nu) + z3
                        enddo
                     endif
                     call dgemv('n',3,mbf_ao,-two,zzz,3,F,1,0d0,
     &                             out,1)
                     if (oactive(iat))then
                        force(1,iat)=force(1,iat)+out(1)
                        force(2,iat)=force(2,iat)+out(2)
                        force(3,iat)=force(3,iat)+out(3)
                     endif
                     if (ldew.and.oactive(ictr)) then
                        force(1,ictr)=force(1,ictr)-out(1)
                        force(2,ictr)=force(2,ictr)-out(2)
                        force(3,ictr)=force(3,ictr)-out(3)
                     endif
                  enddo
               enddo
            endif
         endif
      enddo
c
      return 
      end 
c     
