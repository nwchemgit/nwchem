      subroutine ccsd_sortmo(nocc,nvir,nbf,aint,ndima,scr,
     $                       lnobj,g_objv,g_coul,g_exch,iprt,ncor)
      implicit none
      integer nocc,nvir,nbf,g_objv,g_coul,g_exch,iprt,ndima,lnobj,ncor
      double precision aint(ndima),scr(lnobj)
c
      integer iam,ga_nodeid
      integer ihi,ilo,jhi,jlo,icnt,ij,kl,i,j,k,l,a,b,av
      double precision xx,gmoint2
      double precision zip
      data zip/0.0d00/
c
      if (iprt.gt.5)print *,' sort_ccsd '
c      call ga_summarize(.true.)
c      call ga_print(g_coul)
c      call ga_print(g_exch)
c
      write(6,*)'ndima,lnobj',ndima,lnobj
      write(6,*)'(ncor+nocc)*nbf',(ncor+nocc)*nbf
      iam = ga_nodeid()
      call ga_distribution(g_coul,iam,jlo,jhi,ilo,ihi)
crk   ndima=nocc*nocc*nocc*(nocc+1)/2
      call dfill(ndima,zip,aint,1)
      do ij=1,nocc*(nocc+1)/2
       if (ij.ge.ilo.and.ij.le.ihi)then
        call ga_get(g_exch,1,(ncor+nocc)*nbf,ij,ij,scr,(ncor+nocc)*nbf)
        icnt=(ij-1)*nocc*nocc
        do k=ncor+1,ncor+nocc
         do l=ncor+1,ncor+nocc
          icnt=icnt+1
          aint(icnt)=scr((k-1)*nbf+l)
crk       write(6,*)'icnt,aint',icnt,aint(icnt)
         enddo
        enddo
       endif
      enddo
      call ga_dgop(804,aint,ndima, '+')
c
      call ga_distribution(g_objv,iam,jlo,jhi,ilo,ihi)
      av=0
      do a=ncor+nocc+1,ncor+nocc+nvir
       av=av+1
       if (av.ge.ilo.and.av.le.ihi)then
        icnt=0
c cint
        do b=ncor+nocc+1,ncor+nocc+nvir
         do i=1,nocc
          do j=1,nocc
           icnt=icnt+1
           if (i.ge.j)then
            ij = i*(i-1)/2+j
            kl = (a-1)*nbf+b
           else
            ij = j*(j-1)/2+i
            kl = (b-1)*nbf+a
           endif
           call ga_get(g_coul,kl,kl,ij,ij,xx,1)
           scr(icnt)=xx
crk       write(6,*)'icnt,cint',icnt,scr(icnt)
          enddo
         enddo
        enddo
c dint
        do b=ncor+nocc+1,ncor+nocc+nvir
         do i=1,nocc
          do j=1,nocc
           icnt=icnt+1
           if (i.ge.j)then
            ij = i*(i-1)/2+j
            kl = (a-1)*nbf+b
           else
            ij = j*(j-1)/2+i
            kl = (b-1)*nbf+a
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
          do k=ncor+1,ncor+nocc
           icnt=icnt+1
           if (i.ge.j)then
            ij = i*(i-1)/2+j
            kl = (k-1)*nbf+a
           else
            ij = j*(j-1)/2+i
            kl = (a-1)*nbf+k
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
      return
      end

