      subroutine diis_gets12pm(geom,bas,toll_s,svals, iga_svecs, 
     .     iga_s12m, iga_s12p)
      implicit none
#include "errquit.fh"
#include "global.fh"
      integer geom               ! [in] handle
      integer bas               ! [in] handle
      double precision toll_s   ! [in]
      double precision svals(*) ! [in]
      integer iga_svecs         ! [in]
      integer iga_s12m,iga_s12p ! [out]
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer iga_stmp
c
      iga_s12m = ga_create_atom_blocked(geom, bas, 'gs12m')
      iga_s12p = ga_create_atom_blocked(geom, bas, 'gs12p')
      iga_stmp = ga_create_atom_blocked(geom, bas, 'gtmp')
      call diis_bld12(toll_s, svals, iga_svecs, iga_s12m, 
     &                      iga_stmp, 2)
      call diis_bld12(toll_s, svals, iga_svecs, iga_s12p, 
     &                      iga_stmp, 3)
      if (.not. ga_destroy(iga_stmp)) call errquit
     &     ('diis_gets: could not destroy g_tmp', 0, GA_ERR)
c
      return
      end
      subroutine diis_bld12(toll_s, svals, iga_svecs, iga_sout, 
     &                      iga_stmp,  iwhich)
c***********************************************************************
c     sbuild:  Build S, S**(-1/2) or S**(+1/2) from S eigenvalues and
c     S eigenvectors.
c     iwhich = 1; build S
c     iwhich = 2; build S**(-1/2)
c     iwhich = 3; build S**(+1/2)
c     iwhich = 4; build S**(-1)
c***********************************************************************
C$Id: diis_bld12.F,v 1.19 2003-10-17 22:55:44 carlfahl Exp $
c
      implicit none
#include "errquit.fh"
c
      double precision svals(*) !  S evals [input]
      integer iga_svecs ! GA handle for S evecs [input]
      integer iga_sout  ! GA handle for S^(n) [ouput]
      integer iga_stmp  ! GA handle for scratch [input]
      integer iwhich    ! 2->S(-1/2) 3->S(1/2) 
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "dftps.fh"
c     
      logical LResult
      integer me, nproc, i
      integer ltmpm, itmpm
      integer ilo,ihi,jlo,jhi,numi,numj,l_s,k_s,icount
      double precision toll_s, aaz
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      call ga_sync
c
c
c     
c        Build f(S)
c     
      call ga_copy(iga_svecs,iga_stmp)
      call ga_distribution(iga_stmp,
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. ilo.le.ihi) then
         numi =  ihi-ilo+1
         numj =  jhi-jlo+1
         if (.not. MA_Push_Get(MT_Dbl, numi*numj, 'matS', l_s, k_s))
     &        call errquit('diis_bld12: failed to alloc tmpm',0)
         if (numi.gt.0 .and. numj.gt.0) then
            call ga_get(iga_stmp,ilo,ihi,jlo,jhi,
     .           dbl_mb(k_s),numi)
            icount=k_s
            do i=jlo,jhi
               aaz = svals(i)
               if (aaz.ge.toll_s)  then
                  if(iwhich.eq.2) then
                     aaz = 1.d0/sqrt(aaz)
                  elseif(iwhich.eq.3) then
                     aaz = dsqrt(aaz)
                  elseif(iwhich.eq.4) then
                     aaz=1d0/aaz
                  endif
                  call dscal(numi, aaz, DBL_MB(icount), 1)
               endif
               icount=icount+numi
            enddo
            call ga_put(iga_stmp,ilo,ihi,jlo,jhi,
     .           dbl_mb(k_s),numi)
         endif
         if(.not.MA_Pop_Stack(l_s)) call 
     .        errquit('diisbld12:postack failed',0)
      endif
#ifdef OLD
      if (.not. MA_Push_Get(MT_Dbl, nbf_ao, 'tmpm', ltmpm, itmpm))
     &   call errquit('diis_bld12: failed to alloc tmpm',0)
      do i = me+1, nbf_ao, nproc
         aaz = svals(i)
         if (aaz.ge.toll_s)  then
            call get_col(iga_svecs, nbf_ao, i, DBL_MB(itmpm))
            if(iwhich.eq.2) then
               aaz = 1.d0/sqrt(aaz)
            elseif(iwhich.eq.3) then
               aaz = dsqrt(aaz)
            elseif(iwhich.eq.4) then
               aaz=1d0/aaz
            endif
            call dscal(nbf_ao, aaz, DBL_MB(itmpm), 1)
            call put_col(iga_stmp, nbf_ao, i, DBL_MB(itmpm))
         else
c     write(LuOut,*)' Small eig of S in build of S^-1/2',
c     &           i, aaz, toll_s
            call get_col(iga_svecs, nbf_ao, i, DBL_MB(itmpm))
            call dfill(nbf_ao,0d0,dbl_mb(itmpm),1)
            call put_col(iga_stmp, nbf_ao, i, DBL_MB(itmpm))
         endif
      enddo
      LResult = MA_Pop_Stack(ltmpm)
#endif
c
c
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N', 'T', nbf_ao, nbf_ao, nbf_ao, 1.d0, 
     &              iga_stmp, iga_svecs, 0.d0, iga_sout)
      if (odftps) call pstat_off(ps_dgemm)
c
      return
      end
