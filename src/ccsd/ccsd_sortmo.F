      subroutine ccsd_sortmo(ncor,nocc,nvir,nact,nbf,scra,scrb,
     $                       g_coul,g_exch,g_ncoul,g_nexch,
     $                       iprt)
C     $Id$
      implicit none
      integer ncor,nocc,nvir,nact,nbf,g_coul,g_exch,
     $        g_ncoul,g_nexch,iprt
      double precision scra(nact*nact),scrb(nact*nact)
c
      integer iam,ga_nodeid
      integer ihi,ilo,jhi,jlo,ij,i,j,ad1,ad2
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
            call ga_get(g_coul,1,nact*nact,ij,ij,scra,nact*nact)
            call ga_put(g_ncoul,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scra,nact)
            call mxtrn(scra,scrb,nact,nact)
            call ga_put(g_ncoul,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scrb,nact)
c
            call ga_get(g_exch,1,nact*nact,ij,ij,scra,nact*nact)
            call ga_put(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scra,nact)
            call mxtrn(scra,scrb,nact,nact)
            call ga_put(g_nexch,ad1+1,ad1+nact,ad2+1,ad2+nact,
     &                  scrb,nact)
          endif
        enddo
      enddo
c      call ga_print(g_ncoul)
c      call ga_print(g_nexch)
      return
      end
