      subroutine rohf_canon(oaufbau, oprint)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "cscfps.fh"
#include "cscf.fh"
#include "crohf.fh"
      logical oaufbau
      logical oprint
c
      integer g_u, g_fock, g_tmp
      double precision one, zero
      data one, zero/1.d0, 0.d0/
c
c     This routine assumes that rohf_energy/rohf_fock have been called
c     so that the contents of /crohf/ are current.
c
c     Diagonalize the ROHF 'Fock' matrix
c     
c     If (oaufbau) 
c        diagonalize the whole thing and allow mixing of closed-open-virt
c     else
c        diagonalize separately the cloed-closed, open-open, and 
c        virt-virt parts
c
c     Transform Fock matrices and MO coefficients into the new canonical basis
c
      if (oscfps) call pstat_on(ps_diag)
c
      if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_canon: u',
     $     32, 32, g_u)) call errquit('rohf_canon: ga failed for u', 0,
     &       GA_ERR)
      call ga_zero(g_u)         ! ESSENTIAL FOR SYMM BLOCKED SUBSPACE CANON
      if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_canon: fock',
     $     32, 32, g_fock)) call errquit
     $     ('rohf_canon: ga failed for fock', 0, GA_ERR)
c
      call rohf_get_fock(g_fock)
c
      if (oaufbau) then
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
         call ga_pdsyev(g_fock, g_u, dbl_mb(k_eval), 0)
#else
         call ga_diag_std(g_fock, g_u, dbl_mb(k_eval))
#endif
#else 
         call ga_diag_std_seq(g_fock, g_u, dbl_mb(k_eval))
#endif
      else
c
c     closed-closed piece
c
         if (nclosed .gt. 0) call rohf_canon_subspace
     $        (g_fock, g_u, dbl_mb(k_eval), 1, nclosed,int_mb(k_irs))
c
c     open-open piece
c
         if (nopen .gt. 0) call rohf_canon_subspace
     $        (g_fock, g_u, dbl_mb(k_eval+nclosed),
     $        nclosed+1, nclosed+nopen, int_mb(k_irs))
c     
c     virt-virt piece
c     
         if (nmo-nclosed-nopen .gt. 0)
     $        call rohf_canon_subspace(g_fock, g_u,
     $        dbl_mb(k_eval+nclosed+nopen),
     $        nclosed+nopen+1, nmo, int_mb(k_irs))
      endif
c
      if (oscfps) call pstat_off(ps_diag)
c
      call movecs_fix_phase(g_u)
c
c     Apply rotation to orbitals and fock matrix
c
*ga:1:0
      if (.not. ga_duplicate(g_movecs, g_tmp, 'rohf_canon: tmp'))
     $     call errquit
     $     ('rohf_canon: ga_dup for tmp', 0, GA_ERR)
c
      call ga_copy(g_movecs, g_tmp)
      call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_tmp, g_u,
     $     zero, g_movecs)
c
      if (nbf .ne. nmo) then
         if (.not. ga_destroy(g_tmp)) 
     $        call errquit('rohf_canon: ga_destroy?',0, GA_ERR)
         if (.not. ga_create(MT_DBL, nmo, nmo, 'rohf_canon: tmp',
     $        32, 32, g_tmp)) call errquit
     $        ('rohf_canon: ga failed for tmp', 0, GA_ERR)
      endif
c
      call ga_dgemm('n', 'n', nmo, nmo, nmo, one, crohf_g_fcv, g_u,
     $     zero, g_tmp)
      call ga_dgemm('t', 'n', nmo, nmo, nmo, one, g_u, g_tmp, 
     $     zero, crohf_g_fcv)
c
      if (nopen .gt. 0) then
         call ga_dgemm('n', 'n', nmo, nmo, nmo, one, crohf_g_fcp, g_u,
     $        zero, g_tmp)
         call ga_dgemm('t', 'n', nmo, nmo, nmo, one, g_u, g_tmp, 
     $        zero, crohf_g_fcp)
         call ga_dgemm('n', 'n', nmo, nmo, nmo, one, crohf_g_fpv, g_u,
     $        zero, g_tmp)
         call ga_dgemm('t', 'n', nmo, nmo, nmo, one, g_u, g_tmp, 
     $        zero, crohf_g_fpv)
      endif
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,*)
         write(6,*)
         call util_print_centered(6, 'Eigenvalues', 20, .true.)
         call output(dbl_mb(k_eval), 1, min(nclosed+nopen+5,nmo), 
     $        1, 1, nmo, 1, 1)
         call util_flush(6)
      endif
c
      if (.not. ga_destroy(g_u))
     $     call errquit('rohf_canon: destroy', 0, GA_ERR)
      if (.not. ga_destroy(g_fock))
     $     call errquit('rohf_canon: destroy', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp))
     $     call errquit('rohf_canon: destroy', 0, GA_ERR)
c
      end
      subroutine rohf_canon_subspace(g_fock, g_u, evals, ilo, ihi,
     $     irs)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_fock, g_u, ilo, ihi, irs(*)
      double precision evals(*)
c
      integer k_map,l_map
      integer k_work,l_work,k_tmp,l_tmp
      integer sizew
c
c     map dim ihi-ilo+1
c
!     sizew=ihi-ilo+1
      sizew=ihi
      if (.not.ma_push_get(MT_int,sizew,'map',l_map,k_map))
     &   call errquit('rohf_can: cannot allocate map',0, MA_ERR)
      if (.not.ma_push_get(MT_Dbl,sizew,'work',l_work,k_work))
     &   call errquit('rohf_can: cannot allocate work',0, MA_ERR)
      if (.not.ma_push_get(MT_Dbl,sizew,'tmp',l_tmp,k_tmp))
     &   call errquit('rohf_can: cannot allocate tmp',0, MA_ERR)
      call rohf_canon_sub0(g_fock, g_u, evals, ilo, ihi,
     $     irs, int_mb(k_map),
     ,     dbl_mb(k_work), dbl_mb(k_tmp))
      if (.not.ma_chop_stack(l_map))
     &   call errquit('rohf_can: cannot pop stack',1, MA_ERR)
      return
      end
c     
      subroutine rohf_canon_sub0(g_fock, g_u, evals, ilo, ihi,
     $     irs, map, work, tmp)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_fock, g_u, ilo, ihi, irs(*)
      double precision evals(*)
c
      integer map(*)        ! Should be dynamically allocated
      double precision work(*) 
      double precision tmp(*) 
c
      integer n, i, j, nirs, irrep
      integer g_v
      integer g_tmp
c     
      nirs = 0
      do i = ilo, ihi
         nirs = max(irs(i),nirs)
      enddo
c     
      do irrep = 1, nirs
         n = 0
         do i = ilo, ihi
            if (irs(i) .eq. irrep) then
               n = n + 1
               map(n) = i
            endif
         enddo
         if (n .gt. 0) then
c
*ga:1:0
            if (.not. ga_create(MT_DBL, n, n, 'rohf_canon: tmp',
     $           0, 0, g_tmp)) call errquit
     $           ('rohf_canon: ga failed for tmp', 0, GA_ERR)
*ga:1:0
            if (.not. ga_create(MT_DBL, n, n, 'rohf_canon: v',
     $           n, 0, g_v)) call errquit
     $           ('rohf_canon: ga failed for v', 0, GA_ERR)
c
            call ga_sync
            do i = 1+ga_nodeid(),n,ga_nnodes()
               call ga_get(g_fock,ilo,ihi,map(i),map(i),work(ilo),
     I              ihi-ilo+11)
               do j = 1, n
                  tmp(j) = work(map(j))
               enddo
               call ga_put(g_tmp,1, n, i, i, tmp, n)
            enddo
            call ga_sync
c
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
      call ga_pdsyev(g_tmp, g_v, tmp, 0) 
#else
      call ga_diag_std(g_tmp, g_v, tmp)
#endif
#else 
      call ga_diag_std_seq(g_tmp, g_v, tmp)
#endif
            call ga_sync
            do i = 1, n
               evals(map(i)-ilo+1) = tmp(i)
            enddo
            do i = 1+ga_nodeid(),n,ga_nnodes()
               call ga_get(g_v,1, n, i, i, tmp, n)
               call dfill((ihi-ilo+1), 0.0d0, work(ilo), 1)
               do j = 1, n
                  work(map(j)) = tmp(j)
               enddo
               call ga_acc(g_u,ilo,ihi,map(i),map(i),work(ilo),
     I              ihi-ilo+1,1.0d0)
            enddo
            call ga_sync
            if (.not. (ga_destroy(g_tmp) .and. ga_destroy(g_v)))
     $           call errquit('rohf_canon: ga_destroy ?', 0, GA_ERR)
         endif
      enddo
c
      if (nirs .eq. 1) return
c
c     Now we have diagonalized each symmetry block we need to reorder
c     so that have Auf Bau ordering within the entire block
c
      do i = ilo, ihi
	 map(i) = i
	 do j = ilo, i-1
	    if (evals(map(i)-ilo+1).lt.evals(map(j)-ilo+1)) then
	       n = map(i)
	       map(i) = map(j)
	       map(j) = n
            endif
         enddo
      enddo
      do i = ilo, ihi
	 if (map(i) .ne. i) goto 100
      enddo
      return
c
100   continue
*      do i = ilo, ihi
*         if (map(i) .ne. i) write(6,*) ' FLIP ', i, map(i)
*      enddo
      n = ihi - ilo + 1
*ga:1:0
      if (.not. ga_create(MT_DBL, n, n, 'rohf_canon: tmp',
     $        n, 0, g_tmp)) call errquit
     $        ('rohf_canon: ga failed for tmp', 0, GA_ERR)
      do i = ilo+ga_nodeid(),ihi,ga_nnodes()
         call ga_get(g_u,ilo,ihi, map(i), map(i), tmp(ilo), ihi-ilo+1)
         call ga_put(g_tmp,1,n,i-ilo+1,i-ilo+1,tmp(ilo),n)
      enddo
      call ga_sync
      do i = ilo+ga_nodeid(),ihi,ga_nnodes()
         call ga_get(g_tmp,1,n,i-ilo+1,i-ilo+1,work(ilo),n)
         call ga_put(g_u,ilo,ihi, i, i, work(ilo), ihi-ilo+1)
      enddo
*      write(6,*) ' input '
*      call output(evals,ilo,ihi,1,1,ihi,1,1)
      do i = ilo, ihi
	 tmp(i) = evals(map(i)-ilo+1)
	 work(i)= irs(map(i))
      enddo
      do i = ilo, ihi
	 evals(i-ilo+1) = tmp(i)
	 irs(i)   = work(i)
      enddo
*      write(6,*) ' output '
*      call output(evals,ilo,ihi,1,1,ihi,1,1)
c
      if (.not. ga_destroy(g_tmp)) call errquit('rohf_canon: ga?',0,
     &       GA_ERR)
c
      end
