      subroutine ccsd_sortmo(nocc,nvir,nbf,aint,ndima,scr,
     $                       lnobj,g_objv,g_coul,g_exch,iprt)
      implicit none
      integer nocc,nvir,nbf,g_objv,g_coul,g_exch,iprt,ndima,lnobj
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
      iam = ga_nodeid()
      call ga_distribution(g_coul,iam,jlo,jhi,ilo,ihi)
      ndima=nocc*nocc*nocc*(nocc+1)/2
      call dfill(ndima,zip,aint,1)
      do ij=1,nocc*(nocc+1)/2
       if (ij.ge.ilo.and.ij.le.ihi)then
        call ga_get(g_exch,1,nocc*nbf,ij,ij,scr,nocc*nbf)
        icnt=(ij-1)*nocc*nocc
        do k=1,nocc
         do l=1,nocc
          icnt=icnt+1
          aint(icnt)=scr((k-1)*nbf+l)
         enddo
        enddo
       endif
      enddo
      call ga_dgop(804,aint,ndima, '+')
c
      call ga_distribution(g_objv,iam,jlo,jhi,ilo,ihi)
      av=0
      do a=nocc+1,nocc+nvir
       av=av+1
       if (av.ge.ilo.and.av.le.ihi)then
        icnt=0
c cint
        do b=nocc+1,nocc+nvir
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
c           print *,'a,b,i,j ',a,b,i,j,xx,gmoint2(a,b,i,j)
c           scr(icnt)=gmoint(a,b,i,j)
          enddo
         enddo
        enddo
c dint
        do b=nocc+1,nocc+nvir
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
c           print *,'a,i,b,j ',a,i,b,j,xx,gmoint2(a,i,b,j)
c           scr(icnt)=gmoint(a,i,b,j)
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
            kl = (k-1)*nbf+a
           else
            ij = j*(j-1)/2+i
            kl = (a-1)*nbf+k
           endif
           call ga_get(g_exch,kl,kl,ij,ij,xx,1)
           scr(icnt)=xx
c           print *,'a,j,i,k ',a,j,i,k,xx,gmoint2(a,j,i,k)
c           scr(icnt)=gmoint(a,j,i,k)
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

