      subroutine ccsd_sortmo(ncor,nocc,nvir,nact,nbf,scr,
     $                       lnobj,g_objv,g_coul,g_exch,g_ncoul,g_nexch,
     $                       iprt)
C     $Id: ccsd_sortmo.F,v 2.4 1997-05-15 22:15:41 rg240 Exp $
      implicit none
      integer ncor,nocc,nvir,nact,nbf,g_objv,g_coul,g_exch,
     $        g_ncoul,g_nexch,iprt,lnobj
      double precision scr(2*lnobj)
c
      integer iam,ga_nodeid
      integer ihi,ilo,jhi,jlo,icnt,ij,kl,i,j,k,a,b,av,ad1,ad2,ad3,ad4
      double precision xx
c
      if (iprt.gt.5)print *,' sort_ccsd '
c      call ga_summarize(.true.)
c      call ga_print(g_coul)
c      call ga_print(g_exch)
c
      iam = ga_nodeid()
c
      call ga_distribution(g_objv,iam,jlo,jhi,ilo,ihi)
      av=0
      do a=nocc+1,nocc+nvir
       av=av+1
       if (av.ge.ilo.and.av.le.ihi)then
        icnt=0
c
c dint
        do b=nocc+1,nocc+nvir
         do i=1,nocc
          do j=1,nocc
           icnt=icnt+1
           if (i.ge.j)then
            ij = i*(i-1)/2+j
            kl = (a-1)*nact+b
           else
            ij = j*(j-1)/2+i
            kl = (b-1)*nact+a
           endif
           call ga_get(g_exch,kl,kl,ij,ij,xx,1)
           scr(icnt)=xx
crk       write(6,*)'icnt,dint',icnt,scr(icnt)
          enddo
         enddo
        enddo
c eint
        do i=1,nocc
         do j=1,nocc
          do k=1,nocc
           icnt=icnt+1
           if (i.ge.j)then
            ij = i*(i-1)/2+j
            kl = (k-1)*nact+a
           else
            ij = j*(j-1)/2+i
            kl = (a-1)*nact+k
           endif
           call ga_get(g_exch,kl,kl,ij,ij,xx,1)
           scr(icnt)=xx
crk       write(6,*)'icnt,eint',icnt,scr(icnt)
          enddo
         enddo
        enddo
        if (iprt.gt.50)then
         print *,'mo object ',av,icnt
         write(*,435)(scr(i),i=1,icnt)
 435     format(4e15.5)
        endif
        call ga_put(g_objv,1,icnt,av,av,scr,icnt)
c
       endif
      enddo
c
      call ga_distribution(g_coul,iam,jlo,jhi,ilo,ihi)
      ij=0
      do i=1,nocc
        do j=1,i
          ij=ij+1
          if (ij.ge.ilo.and.ij.le.ihi)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
            call ga_get(g_coul,1,nact*nact,ij,ij,scr,nact*nact)
            call dcopy((nact*nact), scr, 1, scr(lnobj+1), 1)
c$$$            do a=1,nact
c$$$            do b=1,nact
c$$$            ad3=(b-1)*nact+a
c$$$            ad4=(b-1)*nact+a
c$$$            scr(lnobj+ad4)=scr(ad3)
c$$$            enddo
c$$$            enddo
            call ga_put(g_ncoul,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scr(lnobj+1),nact)
            call mxtrn(scr(lnobj+1),scr,nact,nact)
            call ga_put(g_ncoul,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scr,nact)
c
            call ga_get(g_exch,1,nbf*nbf,ij,ij,scr,nbf*nbf)
c$$$            do a=1,nact
c$$$            do b=1,nact
c$$$            ad3=(b-1)*nact+a
c$$$            ad4=(b-1)*nact+a
c$$$            scr(lnobj+ad4)=scr(ad3)
c$$$            enddo
c$$$            enddo
            call dcopy((nact*nact), scr, 1, scr(lnobj+1), 1)
            call ga_put(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scr(lnobj+1),nact)
            call mxtrn(scr(lnobj+1),scr,nact,nact)
            call ga_put(g_nexch,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scr,nact)
          endif
        enddo
      enddo
      return
      end

