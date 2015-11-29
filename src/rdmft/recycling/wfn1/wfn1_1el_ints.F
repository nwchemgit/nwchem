C> \ingroup wfn1
C> @{
C>
C> \brief Compute the 1-electron Hamiltonian
C>
C> This routine computes the 1-electron Hamiltonian, i.e. it computes
C> the kinetic energy plus the nuclear attraction terms and stores
C> it in the corresponding matrix.
C>
      subroutine wfn1_1e_tv(ibasis, nbf, h1)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
c     
c     arguments
c     
      integer ibasis !< [Input] The basis set handle
      integer nbf    !< [Input] The number of basis functions
c
      double precision h1(nbf,nbf) !< [Output] The 1-electron Hamiltonian
c     
c     local variables
c     
      integer nshell
      integer ishell, jshell, mem1, max1e
      integer ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
      integer jj
      logical odoit
      double precision q2
      integer i, noffset,g_loc
c
      logical ocache_save
c
      logical odbug
      logical osome
c
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibasis, nshell)) call errquit
     $     ('wfn1_1e_tv: bas_numcont failed for ibasis', ibasis,
     &       BASIS_ERR)
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c     
c     k_* are the offsets corrsponding to the l_* handles
c     
cso
      call int_mem_1e(max1e, mem1)
cso
      if (.not. MA_push_get(MT_DBL,max1e,'wfn1_1e_tv:buf',l_buf,k_buf))
     $     call errquit('wfn1_1e_tv: ma failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'wfn1_1e_tv:scr',l_scr,k_scr))
     $     call errquit('wfn1_1e_tv: ma failed', mem1, MA_ERR)
c     
c     Loop thru shells
c     
      q2 = 1.0d0
      do jshell = 1, nshell
         if (.not. bas_cn2bfr(ibasis, jshell, jlo, jhi))
     $      call errquit('wfn1_1e_tv: bas_cn2bfr', ibasis, BASIS_ERR)
         do ishell = 1, nshell
c
            if (.not. bas_cn2bfr(ibasis, ishell, ilo, ihi)) 
     $         call errquit('wfn1_1e_tv: bas_cn2bfr', ibasis, BASIS_ERR)
            idim = ihi - ilo + 1
c     
c           Generate the integrals
c     
            call dfill(idim*(jhi-jlo+1),0.0d0,dbl_mb(k_buf),1)
            call int_1eke (ibasis, jshell, ibasis, ishell, 
     $           mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            do jj = jlo, jhi
              call dcopy(idim,dbl_mb(k_buf+idim*(jj-jlo)),1,
     $                   h1(ilo,jj),1)
            enddo
            call dfill(idim*(jhi-jlo+1),0.0d0,dbl_mb(k_buf),1)
            call int_1epe (ibasis, jshell, ibasis, ishell, 
     $           mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            do jj = jlo, jhi
              call daxpy(idim,1.0d0,dbl_mb(k_buf+idim*(jj-jlo)),1,
     $                   h1(ilo,jj),1)
            enddo
         end do
      end do
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('wfn1_1e_tv: pop failed', 0, MA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('wfn1_1e_tv: pop failed', 0, MA_ERR)
c     
      end
C>
C> \brief Compute the overlap matrix
C>
C> This routine computes the overlap matrix.
C>
      subroutine wfn1_1e_s(ibasis, nbf, ov)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
c     
c     arguments
c     
      integer ibasis !< [Input] The basis set handle
      integer nbf    !< [Input] The number of basis functions
c
      double precision ov(nbf,nbf) !< [Output] The overlap matrix
c     
c     local variables
c     
      integer nshell
      integer ishell, jshell, mem1, max1e
      integer ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
      logical odoit
      integer jj
      double precision q2
      integer i, noffset,g_loc
c
      logical ocache_save
c
      logical odbug
      logical osome
c
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibasis, nshell)) call errquit
     $     ('wfn1_1e_s: bas_numcont failed for ibasis', ibasis,
     &       BASIS_ERR)
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c     
c     k_* are the offsets corrsponding to the l_* handles
c     
cso
      call int_mem_1e(max1e, mem1)
cso
      if (.not. MA_push_get(MT_DBL,max1e,'wfn1_1e_s:buf',l_buf,k_buf))
     $     call errquit('wfn1_1e_s: ma failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'wfn1_1e_s:scr',l_scr,k_scr))
     $     call errquit('wfn1_1e_s: ma failed', mem1, MA_ERR)
c     
c     Loop thru shells
c     
      q2 = 1.0d0
      do jshell = 1, nshell
         if (.not. bas_cn2bfr(ibasis, jshell, jlo, jhi))
     $      call errquit('wfn1_1e_s: bas_cn2bfr', ibasis, BASIS_ERR)
         do ishell = 1, nshell
c
            if (.not. bas_cn2bfr(ibasis, ishell, ilo, ihi)) 
     $         call errquit('wfn1_1e_s: bas_cn2bfr', ibasis, BASIS_ERR)
            idim = ihi - ilo + 1
c     
c           Generate the integrals
c     
            call int_1eov (ibasis, jshell, ibasis, ishell, 
     $           mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            do jj = jlo, jhi
              call dcopy(idim,dbl_mb(k_buf+idim*(jj-jlo)),1,
     $                   ov(ilo,jj),1)
            enddo
c
         end do
      end do
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('wfn1_1e_s: pop failed', 0, MA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('wfn1_1e_s: pop failed', 0, MA_ERR)
c     
      end
C>
C> @}
