      subroutine ccsd_sortmo(ncor,nocc,nvir,nact,nbf,scra,scrb,
     $                       g_coul,g_exch,g_ncoul,g_nexch,
     $                       iprt)
C     $Id: ccsd_sortmo.F,v 2.5 1997-05-29 16:42:54 gg502 Exp $
      implicit none
      integer ncor,nocc,nvir,nact,nbf,g_coul,g_exch,
     $        g_ncoul,g_nexch,iprt
      double precision scra(nbf*nbf),scrb(nact*nact)
c
      integer iam,ga_nodeid
      integer ihi,ilo,jhi,jlo,ij,kl,i,j,k,l,a,b,ad1,ad2,ad3,ad4
c
      if (iprt.gt.5)print *,' sort_ccsd '
c      call ga_summarize(.true.)
c      call ga_print(g_coul)
c      call ga_print(g_exch)
c
      iam = ga_nodeid()
c
      call ga_distribution(g_coul,iam,jlo,jhi,ilo,ihi)
      ij=0
      do i=1,nocc
        do j=1,i
          ij=ij+1
          if (ij.ge.ilo.and.ij.le.ihi)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
            call ga_get(g_coul,1,nbf*nbf,ij,ij,scra,nbf*nbf)
            do a=1,nact
            do b=1,nact
            ad3=(ncor+b-1)*nbf+ncor+a
            ad4=(b-1)*nact+a
            scrb(ad4)=scra(ad3)
            enddo
            enddo
            call ga_put(g_ncoul,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scrb,nact)
            call mxtrn(scrb,scra,nact,nact)
            call ga_put(g_ncoul,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scra,nact)
c
            call ga_get(g_exch,1,nbf*nbf,ij,ij,scra,nbf*nbf)
            do a=1,nact
            do b=1,nact
            ad3=(ncor+b-1)*nbf+ncor+a
            ad4=(b-1)*nact+a
            scrb(ad4)=scra(ad3)
            enddo
            enddo
            call ga_put(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scrb,nact)
            call mxtrn(scrb,scra,nact,nact)
            call ga_put(g_nexch,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scra,nact)
          endif
        enddo
      enddo
      return
      end
