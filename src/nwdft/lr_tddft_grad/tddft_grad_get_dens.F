c     Get the blocks of the density matrices needed to construct
c     the 2-electron density matrix for the gradient evaluation.
c     Given are the ranges for i, j, k, and l basis function labels.
c     HvD 2006, NG 2012
c
      subroutine tddft_grad_get_dens(tda,ipol,kfac,
     +           dija,dijb,dila,dilb,dika,dikb,
     +           dkla,dklb,djla,djlb,dkja,dkjb,
     +           pija,pijb,pila,pilb,pika,pikb,
     +           pkla,pklb,pjla,pjlb,pkja,pkjb,
     +           xpyija,xpyijb,xpykla,xpyklb,xpyila,xpyilb,
     +           xpykja,xpykjb,xpyika,xpyikb,xpyjla,xpyjlb,
     +           xmyila,xmyilb,xmykja,xmykjb,
     +           xmyika,xmyikb,xmyjla,xmyjlb,
     +           blen,
     +           ilo,ihi,jlo,jhi,klo,khi,llo,lhi,g_d)
c
c
      implicit none
c
c     Input:
c
      logical tda      ! True if Tamm-Dancoff approximation
      integer ipol     ! =1 (restricted), =2 (unrestricted)
      integer blen     ! the dimension of the density matrix blocks
      integer ilo, ihi ! limits on the i-label
      integer jlo, jhi ! limits on the j-label
      integer klo, khi ! limits on the k-label
      integer llo, lhi ! limits on the l-label
      integer g_d(8)   ! the TDDFT density matrices
c
      double precision kfac ! the fraction of Hartree-Fock exchange
c
c     Output:
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
c     Local:
c
      logical ohfexch ! do we have Hartree-Fock exchange?
c
c     Code:
c
      ohfexch = (kfac.gt.0.0d0)
cDEBUG
c     write(*,*)'*** tddft_grad_get_dens: enter'
cDEBUG
c
      call ga_get(g_d(1),ilo,ihi,jlo,jhi,dija,blen)
      call ga_get(g_d(1),klo,khi,llo,lhi,dkla,blen)
      call ga_get(g_d(3),ilo,ihi,jlo,jhi,pija,blen)
      call ga_get(g_d(3),klo,khi,llo,lhi,pkla,blen)
      call ga_get(g_d(5),ilo,ihi,jlo,jhi,xpyija,blen)
      call ga_get(g_d(5),klo,khi,llo,lhi,xpykla,blen)
cDEBUG
c     call ga_get(g_d(7),ilo,ihi,jlo,jhi,xmyija,blen)
c     call ga_get(g_d(7),klo,khi,llo,lhi,xmykla,blen)
cDEBUG
c
      if (ohfexch) then
        call ga_get(g_d(1),ilo,ihi,llo,lhi,dila,blen)
        call ga_get(g_d(1),ilo,ihi,klo,khi,dika,blen)
        call ga_get(g_d(1),klo,khi,jlo,jhi,dkja,blen)
        call ga_get(g_d(1),jlo,jhi,llo,lhi,djla,blen)
        call ga_get(g_d(3),ilo,ihi,llo,lhi,pila,blen)
        call ga_get(g_d(3),ilo,ihi,klo,khi,pika,blen)
        call ga_get(g_d(3),klo,khi,jlo,jhi,pkja,blen)
        call ga_get(g_d(3),jlo,jhi,llo,lhi,pjla,blen)
        call ga_get(g_d(5),ilo,ihi,llo,lhi,xpyila,blen)
        call ga_get(g_d(5),ilo,ihi,klo,khi,xpyika,blen)
        call ga_get(g_d(5),klo,khi,jlo,jhi,xpykja,blen)
        call ga_get(g_d(5),jlo,jhi,llo,lhi,xpyjla,blen)
c Daniel (1-7-13): We can avoid this part if pure DFT is used.
c Otherwise, we need the (X-Y) contribution, even if CIS/TDA is
c being performed.
c        if (.not.tda) then
        call ga_get(g_d(7),ilo,ihi,llo,lhi,xmyila,blen)
        call ga_get(g_d(7),ilo,ihi,klo,khi,xmyika,blen)
        call ga_get(g_d(7),klo,khi,jlo,jhi,xmykja,blen)
        call ga_get(g_d(7),jlo,jhi,llo,lhi,xmyjla,blen)
c        endif
      endif
c
      if (ipol.eq.2) then
        call ga_get(g_d(2),ilo,ihi,jlo,jhi,dijb,blen)
        call ga_get(g_d(2),klo,khi,llo,lhi,dklb,blen)
        call ga_get(g_d(4),ilo,ihi,jlo,jhi,pijb,blen)
        call ga_get(g_d(4),klo,khi,llo,lhi,pklb,blen)
        call ga_get(g_d(6),ilo,ihi,jlo,jhi,xpyijb,blen)
        call ga_get(g_d(6),klo,khi,llo,lhi,xpyklb,blen)
cDEBUG
c       call ga_get(g_d(8),ilo,ihi,jlo,jhi,xmyijb,blen)
c       call ga_get(g_d(8),klo,khi,llo,lhi,xmyklb,blen)
cDEBUG
c
        if (ohfexch) then
          call ga_get(g_d(2),ilo,ihi,llo,lhi,dilb,blen)
          call ga_get(g_d(2),ilo,ihi,klo,khi,dikb,blen)
          call ga_get(g_d(2),klo,khi,jlo,jhi,dkjb,blen)
          call ga_get(g_d(2),jlo,jhi,llo,lhi,djlb,blen)
          call ga_get(g_d(4),ilo,ihi,llo,lhi,pilb,blen)
          call ga_get(g_d(4),ilo,ihi,klo,khi,pikb,blen)
          call ga_get(g_d(4),klo,khi,jlo,jhi,pkjb,blen)
          call ga_get(g_d(4),jlo,jhi,llo,lhi,pjlb,blen)
          call ga_get(g_d(6),ilo,ihi,llo,lhi,xpyilb,blen)
          call ga_get(g_d(6),ilo,ihi,klo,khi,xpyikb,blen)
          call ga_get(g_d(6),klo,khi,jlo,jhi,xpykjb,blen)
          call ga_get(g_d(6),jlo,jhi,llo,lhi,xpyjlb,blen)
c Daniel (1-7-13): We can avoid this part if pure DFT is used.
c Otherwise, we need the (X-Y) contribution, even if CIS/TDA is
c being performed.
c          if (.not.tda) then
          call ga_get(g_d(8),ilo,ihi,llo,lhi,xmyilb,blen)
          call ga_get(g_d(8),ilo,ihi,klo,khi,xmyikb,blen)
          call ga_get(g_d(8),klo,khi,jlo,jhi,xmykjb,blen)
          call ga_get(g_d(8),jlo,jhi,llo,lhi,xmyjlb,blen)
c          endif
        endif
      endif
cDEBUG
c     write(*,*)'*** tddft_grad_get_dens: leave'
cDEBUG
c
      end
c $Id$
