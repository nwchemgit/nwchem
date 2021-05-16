      subroutine ccsd_idx1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,idx,idim,nsh,shinf,use_ccsd_omp)
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
      integer, intent(in) :: nfj,nfl,nocc,nbf,idim,nsh
      double precision, intent(inout) :: snsi(nbf,nfj,nfl,nocc)
      double precision, intent(inout) :: sisn(nbf,nfl,nfj,nocc)
      double precision, intent(in) :: eri1(*)
      double precision, intent(in) :: eri2(*)
      double precision, intent(in) :: t1ao(nbf,nocc)
      integer, intent(in) :: idx(idim),shinf(nsh,3)
      logical, optional, intent(in) :: use_ccsd_omp
      !
      integer             :: ncnt,ict
      integer             :: ilo,ihi,ish,nfi
      integer             :: klo,khi,ksh,nfk
      double precision    :: factor
      if (.not.present(use_ccsd_omp)) then
          call errquit('ccsd_idx1: use_ccsd_omp not present!',0,0)
      endif
      if (occsdps) then
         call pstat_on(ps_idx1)
      else
         call qenter('idx1',0)
      endif
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
          factor = 0.5d0
        else
          factor = 1.0d0
        endif
        call ccsd_idx1_wrk(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                     nfl,nocc,nbf,ilo,ihi,klo,khi,factor,
     &                     use_ccsd_omp)
        ict=ict+nfk*nfl*nfi*nfj
      enddo
      if (occsdps) then
         call pstat_off(ps_idx1)
      else
         call qexit('idx1',0)
      endif
      return
      end
      subroutine ccsd_idx1_wrk(snsi,sisn,eri1,eri2,t1ao,
     &                         nfj,nfl,nocc,nbf,ilo,ihi,klo,khi,
     &                         factor,use_ccsd_omp)
      implicit none
      integer, intent(in)             :: nfj,nfl,nocc,nbf
      integer, intent(in)             :: ilo,ihi,klo,khi
      double precision, intent(inout) :: snsi(nbf,nfj,nfl,nocc)
      double precision, intent(inout) :: sisn(nbf,nfl,nfj,nocc)
      double precision, intent(in)    :: eri1(nfl,klo:khi,nfj,ilo:ihi)
      double precision, intent(in)    :: eri2(nfj,klo:khi,nfl,ilo:ihi)
      double precision, intent(in)    :: t1ao(nbf,nocc)
      double precision, intent(in)    :: factor
      logical, intent(in)             :: use_ccsd_omp
      !
      double precision                :: int1,int2
      integer                         :: i,j,k,l,ii
      if (use_ccsd_omp) then
!$omp parallel do schedule(static)
!$omp& default(none)
!$omp& private(int1,int2)
!$omp& shared(nfj,nfl,nocc,nbf,ilo,ihi,klo,khi)
!$omp& shared(eri1,eri2,snsi,sisn,t1ao,factor)
        do ii=1,nocc
          do l=1,nfl
            do j=1,nfj
              do i=ilo,ihi
                do k=klo,khi
                  int1=eri1(l,k,j,i)*factor
                  int2=eri2(j,k,l,i)*factor
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
      else
        do l=1,nfl
          do j=1,nfj
            do i=ilo,ihi
              do k=klo,khi
                int1=eri1(l,k,j,i)*factor
                int2=eri2(j,k,l,i)*factor
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
      endif
      return
      end
