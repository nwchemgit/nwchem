c
c     Construct blocks of the 2-electron density matrix given
c     blocks of the various 1-electron density matrices.
c     Given are the ranges for i, j, k, and l basis function labels.
c     HvD 7/2006, NG 11/2012
c
      subroutine tddft_grad_make_2pdm(tda,otriplet,ipol,jfac,kfac,
     +           dija,dijb,dila,dilb,dika,dikb,
     +           dkla,dklb,djla,djlb,dkja,dkjb,
     +           pija,pijb,pila,pilb,pika,pikb,
     +           pkla,pklb,pjla,pjlb,pkja,pkjb,
     +           xpyija,xpyijb,xpykla,xpyklb,xpyila,xpyilb,
     +           xpykja,xpykjb,xpyika,xpyikb,xpyjla,xpyjlb,
     +           xmyila,xmyilb,xmykja,xmykjb,
     +           xmyika,xmyikb,xmyjla,xmyjlb,
     +           blen,
     +           ilo,jlo,klo,llo,
     +           ihi,jhi,khi,lhi,pdm2d,psum)
c
      implicit none
c
#include "errquit.fh"
c
c     Input:
c
      logical tda      ! True if Tamm-Dancoff approximation
      logical otriplet ! True if triplet excited state
      integer ipol     ! =1 (restricted), =2 (unrestricted)
      integer blen     ! the dimension of the density matrix blocks
      integer ilo, ihi ! limits on the i-label
      integer jlo, jhi ! limits on the j-label
      integer klo, khi ! limits on the k-label
      integer llo, lhi ! limits on the l-label
c
      double precision jfac ! the fraction of Coulomb integrals
      double precision kfac ! the fraction of Hartree-Fock exchange
c
      double precision dija(blen,blen),dijb(blen,blen)
      double precision dila(blen,blen),dilb(blen,blen)
      double precision dika(blen,blen),dikb(blen,blen)
      double precision dkla(blen,blen),dklb(blen,blen)
      double precision djla(blen,blen),djlb(blen,blen)
      double precision dkja(blen,blen),dkjb(blen,blen)
c
      double precision pija(blen,blen),pijb(blen,blen)
      double precision pila(blen,blen),pilb(blen,blen)
      double precision pika(blen,blen),pikb(blen,blen)
      double precision pkla(blen,blen),pklb(blen,blen)
      double precision pjla(blen,blen),pjlb(blen,blen)
      double precision pkja(blen,blen),pkjb(blen,blen)
c
      double precision xpyija(blen,blen), xpyijb(blen,blen)
      double precision xpykla(blen,blen), xpyklb(blen,blen)
      double precision xpyila(blen,blen), xpyilb(blen,blen)
      double precision xpykja(blen,blen), xpykjb(blen,blen)
      double precision xpyika(blen,blen), xpyikb(blen,blen)
      double precision xpyjla(blen,blen), xpyjlb(blen,blen)
c
      double precision xmyila(blen,blen), xmyilb(blen,blen)
      double precision xmykja(blen,blen), xmykjb(blen,blen)
      double precision xmyika(blen,blen), xmyikb(blen,blen)
      double precision xmyjla(blen,blen), xmyjlb(blen,blen)
c
c     Output:
c
      double precision pdm2d(blen,blen,blen,blen) ! 2-el density matrix
      double precision psum ! norm of this block of 2-el density matrix
c
c     Local:
c
      logical ohfexch ! do we have Hartree-Fock exchange?
c
      integer i,j,k,l ! counters over orbital labels
c
      double precision denc ! coulomb accumulator
      double precision denx ! exchange accumulator
      double precision half, one
      parameter (half = 0.5d0, one = 1.0d0)
c
c     Code:
c
      ohfexch = (kfac.gt.0.0d0)
cDEBUG
c     write(*,*)'*** tddft_grad_make_2pdm: ohfexch,kfac=',ohfexch,kfac
c     write(*,*)'*** Pij(',ilo,':',ihi,',',jlo,':',jhi,')'
c     call print_matrix_patch(pija,ilo,ihi,jlo,jhi,blen)
c     write(*,*)'*** Dkl(',klo,':',khi,',',llo,':',lhi,')'
c     call print_matrix_patch(dkla,klo,khi,llo,lhi,blen)
cDEBUG
c
      if (ipol.eq.1) then
cDEBUG
c       do i = ilo, ihi
c         do l = llo, lhi
c           write(*,"('*** pila =',2i4,f22.8)")i,l,pila(i,l)
c           write(*,"('*** dila =',2i4,f22.8)")i,l,dila(i,l)
c         enddo
c       enddo
c       do k = klo, khi
c         do j = jlo, jhi
c           write(*,"('*** pkja =',2i4,f22.8)")k,j,pkja(k,j)
c           write(*,"('*** dkja =',2i4,f22.8)")k,j,dkja(k,j)
c         enddo
c       enddo
c       do i = ilo, ihi
c         do k = klo, khi
c           write(*,"('*** pika =',2i4,f22.8)")i,k,pika(i,k)
c           write(*,"('*** dika =',2i4,f22.8)")i,k,dika(i,k)
c         enddo
c       enddo
c       do j = jlo, jhi
c         do l = llo, lhi
c           write(*,"('*** pjla =',2i4,f22.8)")j,l,pjla(j,l)
c           write(*,"('*** djla =',2i4,f22.8)")j,l,djla(j,l)
c         enddo
c       enddo
cDEBUG
        do i = ilo, ihi
          do j = jlo, jhi
            do k = klo, khi
              do l = llo, lhi
                denc = pija(i,j)*dkla(k,l)
     +               + pkla(k,l)*dija(i,j)
                if (.not.otriplet) then
c                  denc = denc + 2*xpyija(i,j)*xpykla(k,l)
                  denc = denc + 2.0d0*xpyija(i,j)*xpykla(k,l)
                endif
c Daniel (4-9-13): Need the Coulomb integral factor here.
                denc = denc*jfac
                denx = 0.0d0
                if (ohfexch) then
                  denx = - half*half*kfac*(
     +                   pila(i,l)*dkja(k,j)+dila(i,l)*pkja(k,j)
     +                 + pika(i,k)*djla(j,l)+dika(i,k)*pjla(j,l))
cnext line are correct 
                  denx = denx - half*kfac*(
     +                 + xpyila(i,l)*xpykja(k,j)
     +                 + xpyika(i,k)*xpyjla(j,l)
     +                   )
cDEBUG
c     write(*,*)'*** tda=',tda
cDEBUG
c Daniel (1-7-13): We can avoid this part if pure DFT is used.
c Otherwise, we need the (X-Y) contribution, even if CIS/TDA is
c being performed.
                  denx = denx - half*kfac*(
     +                 - xmyila(i,l)*xmykja(k,j)
     +                 + xmyika(i,k)*xmyjla(j,l)
     +                   )
                endif
cDEBUG
c               den2 = pija(i,j)*dkla(k,l)
c    +               + xpyija(i,j)*xpykla(k,l)
c               den2 = 2*den2
c               if (ohfexch) then
c                 den2 = den2 - half*kfac*(
c    +                   pila(i,l)*dkja(k,j)
c    +                 + pika(i,k)*djla(j,l)
c    +                 + xpyila(i,l)*xpykja(k,j)
c    +                 + xpyika(i,k)*xpyjla(j,l)
c    +                   )
c                 if (.not.tda) then
c                   den2 = den2 - half*kfac*(
c    +                   - xmyila(i,l)*xmykja(k,j)
c    +                   + xmyika(i,k)*xmyjla(j,l))
c                 endif
c               endif
cDEBUG
                pdm2d(l,k,j,i) = denc+denx
cDEBUG
c               pdm2d(l,k,j,i) = denx
cDEBUG
              enddo
            enddo
          enddo
        enddo
      else if (ipol.eq.2) then
cDEBUG
c       do i = ilo, ihi
c         do l = llo, lhi
c           write(*,"('*** pila =',2i4,f22.8)")i,l,pila(i,l)*2
c           write(*,"('*** pilb =',2i4,f22.8)")i,l,pilb(i,l)*2
c           write(*,"('*** dila =',2i4,f22.8)")i,l,dila(i,l)*2
c           write(*,"('*** dilb =',2i4,f22.8)")i,l,dilb(i,l)*2
c         enddo
c       enddo
c       do k = klo, khi
c         do j = jlo, jhi
c           write(*,"('*** pkja =',2i4,f22.8)")k,j,pkja(k,j)*2
c           write(*,"('*** pkjb =',2i4,f22.8)")k,j,pkjb(k,j)*2
c           write(*,"('*** dkja =',2i4,f22.8)")k,j,dkja(k,j)*2
c           write(*,"('*** dkjb =',2i4,f22.8)")k,j,dkjb(k,j)*2
c         enddo
c       enddo
c       do i = ilo, ihi
c         do k = klo, khi
c           write(*,"('*** pika =',2i4,f22.8)")i,k,pika(i,k)*2
c           write(*,"('*** pikb =',2i4,f22.8)")i,k,pikb(i,k)*2
c           write(*,"('*** dika =',2i4,f22.8)")i,k,dika(i,k)*2
c           write(*,"('*** dikb =',2i4,f22.8)")i,k,dikb(i,k)*2
c         enddo
c       enddo
c       do j = jlo, jhi
c         do l = llo, lhi
c           write(*,"('*** pjla =',2i4,f22.8)")j,l,pjla(j,l)*2
c           write(*,"('*** pjlb =',2i4,f22.8)")j,l,pjlb(j,l)*2
c           write(*,"('*** djla =',2i4,f22.8)")j,l,djla(j,l)*2
c           write(*,"('*** djlb =',2i4,f22.8)")j,l,djlb(j,l)*2
c         enddo
c       enddo
cDEBUG
        do i = ilo, ihi
          do j = jlo, jhi
            do k = klo, khi
              do l = llo, lhi
                denc = pija(i,j)*dkla(k,l)+pkla(k,l)*dija(i,j)
     +               + pijb(i,j)*dklb(k,l)+pklb(k,l)*dijb(i,j)
     +               + pija(i,j)*dklb(k,l)+pkla(k,l)*dijb(i,j)
     +               + pijb(i,j)*dkla(k,l)+pklb(k,l)*dija(i,j)
     +               + xpyija(i,j)*xpykla(k,l)
     +               + xpyijb(i,j)*xpyklb(k,l)
     +               + xpyija(i,j)*xpyklb(k,l)
     +               + xpyijb(i,j)*xpykla(k,l)
c Daniel (4-9-13): Need the Coulomb integral factor here.
                denc = denc*jfac
                denx = 0.0d0
                if (ohfexch) then
                  denx = denx - half*one*kfac*(
     +                   pila(i,l)*dkja(k,j)+dila(i,l)*pkja(k,j)
     +                 + pilb(i,l)*dkjb(k,j)+dilb(i,l)*pkjb(k,j)
     +                 + pika(i,k)*djla(j,l)+dika(i,k)*pjla(j,l)
     +                 + pikb(i,k)*djlb(j,l)+dikb(i,k)*pjlb(j,l)
     +                   )
                  denx = denx - half*kfac*(
     +                 + xpyila(i,l)*xpykja(k,j)
     +                 + xpyilb(i,l)*xpykjb(k,j)
     +                 + xpyika(i,k)*xpyjla(j,l)
     +                 + xpyikb(i,k)*xpyjlb(j,l))
c Daniel (1-7-13): We can avoid this part if pure DFT is used.
c Otherwise, we need the (X-Y) contribution, even if CIS/TDA is
c being performed.
                  denx = denx - half*kfac*(
     +                 - xmyila(i,l)*xmykja(k,j)
     +                 - xmyilb(i,l)*xmykjb(k,j)
     +                 + xmyika(i,k)*xmyjla(j,l)
     +                 + xmyikb(i,k)*xmyjlb(j,l))
                endif
                pdm2d(l,k,j,i) = denc+denx
cDEBUG
c               pdm2d(l,k,j,i) = denx
cDEBUG
              enddo
            enddo
          enddo
        enddo
      else
        call errquit("tddft_grad_make_2pdm: illegal value for ipol",0,
     +               UERR)
      endif
c
      psum = 0.0d0
      do i = ilo, ihi
        do j = jlo, jhi
          do k = klo, khi
            do l = llo, lhi
              psum = psum + pdm2d(l,k,j,i)**2
            enddo
          enddo
        enddo
      enddo
      psum = dsqrt(psum)
cDEBUG
c     write(*,"('*** psum =',f22.8)") psum
cDEBUG
c
      end
c $Id$
