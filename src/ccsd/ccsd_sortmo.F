      subroutine ccsd_sortmo(ncor,nocc,nvir,nact,nbf,scra,scrb,
     $                       g_coul,g_exch,g_ncoul,g_nexch,
     $                       iprt)
C     $Id: ccsd_sortmo.F,v 2.6 1997-09-19 17:52:17 rg240 Exp $
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

c$$$            call ga_get(g_coul,1,nbf*nbf,ij,ij,scra,nbf*nbf)
c$$$            do a=1,nact
c$$$            do b=1,nact
c$$$            ad3=(ncor+b-1)*nbf+ncor+a
c$$$            ad4=(b-1)*nact+a
c$$$            scrb(ad4)=scra(ad3)
c$$$            enddo
c$$$            enddo
c$$$            call ga_put(g_ncoul,ad2+1,ad2+nact,ad1+1,ad1+nact,
c$$$     &                  scrb,nact)
c$$$            call mxtrn(scrb,scra,nact,nact)
c$$$            call ga_put(g_ncoul,ad1+1,ad1+nact,ad2+1,ad2+nact,
c$$$     &                  scra,nact)

            call ga_get(g_coul,1,nact*nact,ij,ij,scra,1)
            call ga_put(g_ncoul,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scra,nact)
            call mxtrn(scra,scrb,nact,nact)
            call ga_put(g_ncoul,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scrb,nact)

c
c$$$            call ga_get(g_exch,1,nbf*nbf,ij,ij,scra,nbf*nbf)
c$$$            do a=1,nact
c$$$            do b=1,nact
c$$$            ad3=(ncor+b-1)*nbf+ncor+a
c$$$            ad4=(b-1)*nact+a
c$$$            scrb(ad4)=scra(ad3)
c$$$            enddo
c$$$            enddo
c$$$            call ga_put(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
c$$$     &                  scrb,nact)
c$$$            call mxtrn(scrb,scra,nact,nact)
c$$$            call ga_put(g_nexch,ad1+1,ad1+nact,ad2+1,ad2+nact,
c$$$     &                  scra,nact)

            call ga_get(g_exch,1,nact*nact,ij,ij,scra,1)
            call ga_put(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scra,nact)
            call mxtrn(scra,scrb,nact,nact)
            call ga_put(g_nexch,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scrb,nact)
            
          endif
        enddo
      enddo
      return
      end
