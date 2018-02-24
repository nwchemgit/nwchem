      subroutine dftg_quadv0_b(ictr, ldew, AO_bas_han, ipol, GRAD, kske,
     &                         kslap, nbf_ao, mbf_ao, nat, nq, chi, 
     &                         delchi, heschi, d3chi, Amat, Cmat, Mmat, 
     &                         Lmat, ccrap, 
     I     ibf, iniz, ifin, cetobfr,
     P     PPP, F, zzz, hh, g_dens, force, 
     &                         oactive)
*
* $Id$
*
      Implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "dftpara.fh"
c     
      logical ldew   ! true if weight derivatives are included [input]
      integer ictr
      integer ipol, nbf_ao, mbf_ao, nat, nq
      Integer AO_bas_han
      logical GRAD, oactive(nat), kske, kslap
      double precision zzz(*)
      double precision PPP(*), F(*)
      integer ibf(mbf_ao), iniz(nat), ifin(nat)
      integer cetobfr(2,*)
      integer g_dens(2)
c     
c     scratch array for hessian product
c     
      double precision hh(nq,3,*)
c     
c     functions and grad
c     
      double precision chi(nq,mbf_ao), delchi(nq,3,mbf_ao),
     &                 heschi(nq,6,mbf_ao), d3chi(nq,10,mbf_ao)
c     
c     Vxc and its GCA part
c     
      double precision Amat(nq,ipol), Cmat(nq,3,ipol)
      double precision Mmat(nq,ipol), Lmat(nq,ipol)
c     
c     force vector
c     
      double precision force(3,nat)
c
c     local
c
      logical ldo
      integer iat, n, ispin, mu,  nu, mu1, zzzptr
      integer ifinia,inizia,nnia,nbfia,ifirst,ilast
      integer nbhand
      double precision C1, C2, C3, z1, z2, z3, fatto
      double precision out(3)
      double precision ccrap(nq,mbf_ao,*)
      double precision two
      parameter(two=2.d0)
c
      if (GRAD)then
         do ispin=1,ipol
            do nu = 1,mbf_ao
               do n = 1, nq
                  ccrap(n,nu,ispin) = Cmat(n,1,ispin)*delchi(n,1,nu) +
     &                 Cmat(n,2,ispin)*delchi(n,2,nu) +
     &                 Cmat(n,3,ispin)*delchi(n,3,nu)
               enddo
            enddo
         enddo
      endif
      do iat = 1, nat
         ldo = .true.
         if (ldew)ldo = ldo.and.(iat.ne.ictr)
         if (ldo)then
            if (iniz(iat).ne.0)then

         inizia = iniz(iat)
         ifinia = ifin(iat)
         ifirst = cetobfr(1,iat)
         ilast  = cetobfr(2,iat)
         nbfia = ilast-ifirst+1
         nnia = ifinia-inizia+1
c
         out(1)=0d0
         out(2)=0d0
         out(3)=0d0
         do ispin = 1, ipol
            if(dftnbget) then
               call ga_nbget(g_dens(ispin), 
     %              ifirst, ilast, 1, nbf_ao, PPP,nbfia, nbhand)
            else
               if(truerepdm) then
                  call dftg_dmget(dbl_mb(k_repdm(ispin)),nbf_ld,
     %                 ifirst, ilast, PPP,nbfia)
               else
                  call ga_get(g_dens(ispin), 
     %                 ifirst, ilast, 1, nbf_ao, PPP,nbfia)
               endif
            endif
               do mu = inizia,ifinia
                  mu1=mu-inizia+1
c               
c                    LDA part
c               
                     do n = 1, nq
                        fatto = Amat(n,ispin)
                        hh(n,1,mu1) =  fatto*delchi(n,1,mu)
                        hh(n,2,mu1) =  fatto*delchi(n,2,mu)
                        hh(n,3,mu1) =  fatto*delchi(n,3,mu)
                     enddo
                     if (GRAD)then
                        do n = 1, nq
                           C1 = Cmat(n,1,ispin)
                           C2 = Cmat(n,2,ispin)
                           C3 = Cmat(n,3,ispin)
                           hh(n,1,mu1)=hh(n,1,mu1)+C1*heschi(n,1,mu) +
     &                                         C2*heschi(n,2,mu) +
     &                                         C3*heschi(n,3,mu)
                           hh(n,2,mu1)=hh(n,2,mu1)+C1*heschi(n,2,mu) +
     &                                         C2*heschi(n,4,mu) +
     &                                         C3*heschi(n,5,mu)
                           hh(n,3,mu1)=hh(n,3,mu1)+C1*heschi(n,3,mu) +
     &                                         C2*heschi(n,5,mu) +
     &                                         C3*heschi(n,6,mu)
                        enddo
                     endif
                     if (kslap) then
                       do n=1,nq
                         hh(n,1,mu1) = hh(n,1,mu1) + (d3chi(n,1,mu) +
     &                   d3chi(n,4,mu) + d3chi(n,6,mu))*Lmat(n,ispin)
                         hh(n,2,mu1) = hh(n,2,mu1) + (d3chi(n,2,mu) +
     &                   d3chi(n,7,mu) + d3chi(n,9,mu))*Lmat(n,ispin)
                         hh(n,3,mu1) = hh(n,3,mu1) + (d3chi(n,3,mu) +
     &                   d3chi(n,8,mu) + d3chi(n,10,mu))*Lmat(n,ispin)
                       end do
                     end if
                     enddo ! mu
cold zzz(3,mbf)   =hh(3,q)   * chi(q,mbf)
c    zzz(3,nnia,mbf)=hh(q,3,nnia)* chi(q,mbf)
                     call dgemm('t','n',3*nnia,mbf_ao,nq,1.d0,hh,nq,
     C                    chi, nq, 0d0,zzz,3*nnia)
                     if (GRAD)then
c     zzz(3,nnia,mbf)=dchi(q,3,nnia)* cc(q,mbf)
                        call dgemm('t','n',3*nnia,mbf_ao,nq,1.d0,
     D                       delchi(1,1,inizia),nq,
     C                       ccrap(1,1,ispin), nq, 1.d0,zzz,3*nnia)
                     endif
                      if (kske)then
                         do mu = inizia,ifinia
                            mu1=mu-inizia+1
                            do nu = 1, mbf_ao
                               z1 = 0.d0
                               z2 = 0.d0
                               z3 = 0.d0
                               do n = 1, nq
                                  z1 = z1 + 
     A                                 (delchi(n,1,nu)*heschi(n,1,mu) +
     &                                 delchi(n,2,nu)*heschi(n,2,mu) +
     &                                 delchi(n,3,nu)*heschi(n,3,mu))
     &                                 *Mmat(n,ispin)
                                  z2 = z2 + 
     A                                 (delchi(n,1,nu)*heschi(n,2,mu) +
     &                                 delchi(n,2,nu)*heschi(n,4,mu) +
     &                                 delchi(n,3,nu)*heschi(n,5,mu))
     &                                 *Mmat(n,ispin)
                                  z3 = z3 + 
     A                                 (delchi(n,1,nu)*heschi(n,3,mu) +
     &                                 delchi(n,2,nu)*heschi(n,5,mu) +
     &                                 delchi(n,3,nu)*heschi(n,6,mu))
     &                                 *Mmat(n,ispin)
                               enddo
c     zzz(3,nnia,mbf)
c     zzz(1,mu,nu)
                               zzzptr=1+(mu1-1)*3+(nu-1)*3*nnia
                               zzz(zzzptr)   = zzz(zzzptr)   + z1
                               zzz(zzzptr+1) = zzz(zzzptr+1) + z2
                               zzz(zzzptr+2) = zzz(zzzptr+2) + z3
                            enddo
                         enddo
                      endif
                      if (kslap)then
                         do mu = inizia,ifinia
                            mu1=mu-inizia+1
                            do nu = 1, mbf_ao
                               z1 = 0.d0
                               z2 = 0.d0
                               z3 = 0.d0
                               do n = 1, nq
                                  z1 = z1 + (2d0*
     A                                 (delchi(n,1,nu)*heschi(n,1,mu) +
     &                                 delchi(n,2,nu)*heschi(n,2,mu) +
     &                                 delchi(n,3,nu)*heschi(n,3,mu)) +
     &                                 (heschi(n,1,nu)+heschi(n,4,nu) +
     &                                  heschi(n,6,nu))*delchi(n,1,mu))
     &                                 *Lmat(n,ispin)
                                  z2 = z2 + (2d0*
     A                                 (delchi(n,1,nu)*heschi(n,2,mu) +
     &                                 delchi(n,2,nu)*heschi(n,4,mu) +
     &                                 delchi(n,3,nu)*heschi(n,5,mu)) +
     &                                 (heschi(n,1,nu)+heschi(n,4,nu) +
     &                                  heschi(n,6,nu))*delchi(n,2,mu))
     &                                 *Lmat(n,ispin)
                                  z3 = z3 + (2d0*
     A                                 (delchi(n,1,nu)*heschi(n,3,mu) +
     &                                 delchi(n,2,nu)*heschi(n,5,mu) +
     &                                 delchi(n,3,nu)*heschi(n,6,mu)) +
     &                                 (heschi(n,1,nu)+heschi(n,4,nu) +
     &                                  heschi(n,6,nu))*delchi(n,3,mu))
     &                                 *Lmat(n,ispin)
                               enddo
c     zzz(3,nnia,mbf)
c     zzz(1,mu,nu)
                               zzzptr=1+(mu1-1)*3+(nu-1)*3*nnia
                               zzz(zzzptr)   = zzz(zzzptr)   + z1
                               zzz(zzzptr+1) = zzz(zzzptr+1) + z2
                               zzz(zzzptr+2) = zzz(zzzptr+2) + z3
                            enddo
                         enddo
                      endif
                      if (dftnbget) then
                         call ga_nbwait(nbhand)
                      endif
                      call gat_mat(F, PPP,  nbfia, nbf_ao,  nnia,mbf_ao,
     1                     ifirst, 1, ibf(inizia), ibf)
cold      zzz(3,mbf)*F(mbf)
c     zzz(3,nnia,mbf)*F(nnia,mbf)
                     call dgemv('n',3,nnia*mbf_ao,-two,zzz,3,
     F                    F,1,1d0,out,1)
                  enddo         !ispin
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
            endif
         endif
      enddo
c
      return 
      end 
c     
      subroutine dftg_dmget(repdm, nbf_ld,
     %                 ilo, ihi,  Pmat,nbfia)
      implicit none
#include "mafdecls.fh"
      integer  ilo, ihi,  nbfia,nbf_ld
      double precision pmat(nbfia,*)
      double precision repdm(*)
      integer ij,nnn
      integer i,j
c
c     3blocks
c
c     full block ilo,ihi;1,ilo-1
      if(ilo.gt.1) then
         do j=1,ilo-1
            ij=((j-1)*(2*(nbf_ld+1)-j)+1)/2+ilo-j+1
            nnn=nbfia
            call dcopy(nnn,repdm(ij),1,pmat(1,j),1)
         enddo
      endif
c     diag block ilo,ihi;ilo,ihi
      nnn=nbfia
      do j=ilo,ihi
         ij=((j-1)*(2*(nbf_ld+1)-j)+1)/2+1
         call dcopy(nnn,repdm(ij),1,pmat(j-ilo+1,j),1)
         nnn=nnn-1
      enddo
c     copy offdiag terms (aka transp)
      if(ihi.gt.ilo) then
         do j=ilo,ihi
            do i=j+1,ihi
               pmat(j-ilo+1,i)=pmat(i-ilo+1,j)
            enddo
         enddo
      endif
c     full block ilo,ihi;1,ilo-1 to be transposed
c     transpose ihi+1,nbf;ilo,ihi into ilo,ihi;ihi+1;nbf
      if(ihi.lt.nbf_ld) then
         nnn=nbf_ld-ihi
         do j=ilo,ihi
            ij=((j-1)*(2*(nbf_ld+1)-j)+1)/2+(ihi+1)-j+1
            call dcopy(nnn,repdm(ij),1,pmat(j-ilo+1,ihi+1),nbfia)
         enddo
      endif
      return
      end
