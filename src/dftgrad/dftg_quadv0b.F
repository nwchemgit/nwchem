      Subroutine dftg_quadv0_b(
     &     AO_bas_han,ipol,GRAD,
     &     nbf_ao,mbf_ao,nat,nq,
     &     chi,delchi,heschi,
     &     Amat,Cmat,
     &     ibf,iniz,ifin,
     &     PPP,F,zzz,hh,
     &     iga_dens,force,
     &     oactive, nactive )
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
c     
c     handles
c     
      Integer AO_bas_han
      Integer nactive
      logical GRAD, oactive(nat)
      dimension zzz(3,mbf_ao)
C     
C     DM's
C     
      dimension PPP(nbf_ao),F(mbf_ao)
      dimension ibf(mbf_ao),iniz(nat),ifin(nat)
      dimension iga_dens(2)
c     
c     scratch array for hessian product
c     
      dimension hh(3,nq)
C     
C     functions and grad
C     
      
      dimension chi(nq,mbf_ao),delchi(3,nq,mbf_ao),
     &     heschi(6,nq,mbf_ao)
C     
C     Vxc and its GCA part
C     
      dimension Amat(ipol,nq),Cmat(ipol,3,nq)

c     
c     force vector
c     
      Real*8 force(3,nat)
c--   >
      twoi=2.d0
C     
C     
C     
      do iat=1,nat
        if( oactive(iat) .and. iniz(iat).ne.0) then
c****     if( iniz(iat).ne.0) then
          if(.not.bas_ce2bfr(AO_bas_han,iat,ifirst,ilast))
     &         call errquit(
     &         ' dftg_quadv0_b: number of basis with atom ?',1)
          do mu=iniz(iat),ifin(iat)
            mmu=ibf(mu)
            do ispin=1,ipol
              call ga_get(iga_dens(ispin),1,nbf_ao,mmu,mmu,PPP,1)
              call dft_gather(mbf_ao,
     &             F,PPP,ibf,0)

C             
C             LDA part
C             
              do n=1,nq
                fatto=Amat(ispin,n)
                hh(1,n)=fatto*delchi(1,n,mu)
                hh(2,n)=fatto*delchi(2,n,mu)
                hh(3,n)=fatto*delchi(3,n,mu)
              enddo

              if (GRAD) then

                do n=1,nq
                  C1=Cmat(ispin,1,n)
                  C2=Cmat(ispin,2,n)
                  C3=Cmat(ispin,3,n)
                  
                  hh(1,n)=hh(1,n)+C1*heschi(1,n,mu)+
     *                 C2*heschi(2,n,mu)+C3*heschi(3,n,mu)
                  hh(2,n)=hh(2,n)+ C1*heschi(2,n,mu)+
     *                 C2*heschi(4,n,mu)+C3*heschi(5,n,mu)
                  hh(3,n)=hh(3,n)+C1*heschi(3,n,mu)+
     *                 C2*heschi(5,n,mu)+C3*heschi(6,n,mu)
                enddo
              endif

              call dgemm('n','n',3,mbf_ao,nq,
     &             1.d0,hh,3,chi,nq,0.d0,
     &             zzz,3)
              if (GRAD) then

                do nu=1,mbf_ao
                  z1=0.d0
                  z2=0.d0
                  z3=0.d0
                  do n=1,nq
C                   
C                   GCA part
C                   
                    DC=  Cmat(ispin,1,n)*delchi(1,n,nu)+
     *                   Cmat(ispin,2,n)*delchi(2,n,nu)+
     *                   Cmat(ispin,3,n)*delchi(3,n,nu)

                    z1=z1+delchi(1,n,mu)*DC
                    z2=z2+delchi(2,n,mu)*DC
                    z3=z3+delchi(3,n,mu)*DC

                  enddo
                  zzz(1,nu)=zzz(1,nu)+z1
                  zzz(2,nu)=zzz(2,nu)+z2
                  zzz(3,nu)=zzz(3,nu)+z3
                enddo
              endif

              call dgemv('n',3,mbf_ao,-twoi,zzz,3,F,1,1.d0,
     &             force(1,iat),1)
            enddo
          enddo
        endif

      enddo

      return 
      end 
c     
