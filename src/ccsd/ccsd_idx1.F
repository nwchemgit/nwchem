      subroutine ccsd_idx1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,idx,idim,nsh,shinf)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "schwarz.fh"
#include "eaf.fh"
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,idim,nsh
      integer idx(idim),shinf(nsh,3)
      double precision eri1(*),
     &                 eri2(*),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
c
      integer i,j,k,l,ii,ncnt,ia,ka
      integer nfi,nfk,ilo,ihi,klo,khi,ish,ksh
      integer ict,off1,off2
      integer mlkj,mkl,mkj
     
      double precision int1, int2, scale
c
      if (occsdps) then
         call pstat_on(ps_idx1)
      else
         call qenter('idx1',0)
      endif
c
      ict=1
      do ncnt=3,idim
        ksh=(idx(ncnt)-1)/nsh+1
        ish=idx(ncnt)-(ksh-1)*nsh
        nfk=shinf(ksh,1)
        klo=shinf(ksh,2)
        khi=shinf(ksh,3)
        nfi=shinf(ish,1)
        ilo=shinf(ish,2)
        ihi=shinf(ish,3)
        if (ish.eq.ksh) then
          call ccsd_idx1_wrk2(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                        nfl,nocc,nbf,ilo,ihi,klo,khi)
        else
          call ccsd_idx1_wrk1(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                        nfl,nocc,nbf,ilo,ihi,klo,khi)
        endif
        ict=ict+nfk*nfl*nfi*nfj
      enddo
      if (occsdps) then
         call pstat_off(ps_idx1)
      else
         call qexit('idx1',0)
      endif
      return
      end
#if defined(USE_OPENMP)
      subroutine ccsd_idx1_wrk1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,ilo,ihi,klo,khi)
      implicit none
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,ilo,ihi,klo,khi
      integer i,j,k,l,ii
      double precision eri1(nfl,klo:khi,nfj,ilo:ihi),
     &                 eri2(nfj,klo:khi,nfl,ilo:ihi),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
      double precision int1,int2
      call pstat_on(ps_idx11)
!$omp parallel do schedule(static)
!$omp& default(none)
!$omp& private(int1,int2)
!$omp& firstprivate(nfj,nfl,nocc,nbf,ilo,ihi,klo,khi)
!$omp& shared(eri1,eri2,snsi,sisn,t1ao)
      do ii=1,nocc
        do l=1,nfl
          do j=1,nfj
            do i=ilo,ihi
              do k=klo,khi
                int1=eri1(l,k,j,i)
                int2=eri2(j,k,l,i)
                snsi(i,j,l,ii)=snsi(i,j,l,ii)+t1ao(k,ii)*int1
                snsi(k,j,l,ii)=snsi(k,j,l,ii)+t1ao(i,ii)*int2
                sisn(k,l,j,ii)=sisn(k,l,j,ii)+t1ao(i,ii)*int1
                sisn(i,l,j,ii)=sisn(i,l,j,ii)+t1ao(k,ii)*int2
              enddo
            enddo
          enddo
        enddo
      enddo
!$omp end parallel do
      call pstat_off(ps_idx11)
      return
      end
      subroutine ccsd_idx1_wrk2(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,ilo,ihi,klo,khi)
      implicit none
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,ilo,ihi,klo,khi
      integer i,j,k,l,ii
      double precision eri1(nfl,klo:khi,nfj,ilo:ihi),
     &                 eri2(nfj,klo:khi,nfl,ilo:ihi),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
      double precision int1,int2
      call pstat_on(ps_idx12)
!$omp parallel do schedule(static)
!$omp& default(none)
!$omp& private(int1,int2)
!$omp& firstprivate(nfj,nfl,nocc,nbf,ilo,ihi,klo,khi)
!$omp& shared(eri1,eri2,snsi,sisn,t1ao)
      do ii=1,nocc
        do l=1,nfl
          do j=1,nfj
            do i=ilo,ihi
              do k=klo,khi
                int1=eri1(l,k,j,i)*0.5d0
                int2=eri2(j,k,l,i)*0.5d0
                snsi(i,j,l,ii)=snsi(i,j,l,ii)+t1ao(k,ii)*int1
                snsi(k,j,l,ii)=snsi(k,j,l,ii)+t1ao(i,ii)*int2
                sisn(k,l,j,ii)=sisn(k,l,j,ii)+t1ao(i,ii)*int1
                sisn(i,l,j,ii)=sisn(i,l,j,ii)+t1ao(k,ii)*int2
              enddo
            enddo
          enddo
        enddo
      enddo
!$omp end parallel do
      call pstat_off(ps_idx12)
      return
      end
#else // USE_OPENMP
      subroutine ccsd_idx1_wrk1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,ilo,ihi,klo,khi)
      implicit none
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,ilo,ihi,klo,khi
      integer i,j,k,l,ii
      double precision eri1(nfl,klo:khi,nfj,ilo:ihi),
     &                 eri2(nfj,klo:khi,nfl,ilo:ihi),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
      double precision int1,int2
      call pstat_on(ps_idx11)
      do l=1,nfl
        do j=1,nfj
          do i=ilo,ihi
            do k=klo,khi
              int1=eri1(l,k,j,i)
              int2=eri2(j,k,l,i)
              do ii=1,nocc
                snsi(i,j,l,ii)=snsi(i,j,l,ii)+t1ao(k,ii)*int1
                snsi(k,j,l,ii)=snsi(k,j,l,ii)+t1ao(i,ii)*int2
                sisn(k,l,j,ii)=sisn(k,l,j,ii)+t1ao(i,ii)*int1
                sisn(i,l,j,ii)=sisn(i,l,j,ii)+t1ao(k,ii)*int2
              enddo
            enddo
          enddo
        enddo
      enddo
      call pstat_off(ps_idx11)
      return
      end
      subroutine ccsd_idx1_wrk2(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,ilo,ihi,klo,khi)
      implicit none
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,ilo,ihi,klo,khi
      integer i,j,k,l,ii
      double precision eri1(nfl,klo:khi,nfj,ilo:ihi),
     &                 eri2(nfj,klo:khi,nfl,ilo:ihi),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
      double precision int1,int2
      call pstat_on(ps_idx12)
      do l=1,nfl
        do j=1,nfj
          do i=ilo,ihi
            do k=klo,khi
              int1=eri1(l,k,j,i)*0.5d0
              int2=eri2(j,k,l,i)*0.5d0
              do ii=1,nocc
                snsi(i,j,l,ii)=snsi(i,j,l,ii)+t1ao(k,ii)*int1
                snsi(k,j,l,ii)=snsi(k,j,l,ii)+t1ao(i,ii)*int2
                sisn(k,l,j,ii)=sisn(k,l,j,ii)+t1ao(i,ii)*int1
                sisn(i,l,j,ii)=sisn(i,l,j,ii)+t1ao(k,ii)*int2
              enddo
            enddo
          enddo
        enddo
      enddo
      call pstat_off(ps_idx12)
      return
      end
#endif // USE_OPENMP
