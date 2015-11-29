C> \ingroup wfn1
C> @{
C>
C> \brief Compute the 2-electron integrals
C>
C> This routine computes the 2-electron integrals, i.e. it computes
C> the 2-electron integrals and stores them in the corresponding
C> matrix/tensor.
C>
      subroutine wfn1_2e_rep(ibasis, nbf, eri, erix)
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
      double precision eri(nbf,nbf,nbf,nbf) !< [Output] The 2-electron
                                            !< integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Output] The 2-electron
                                             !< integrals sorted for
                                             !< Hartree-Fock exchange
c     
c     local variables
c     
      integer nshell
      integer ishell, jshell, kshell, lshell, mem1, max1e
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
      integer type
      integer max2e, mem2
      logical odoit
      double precision q2
      integer i, noffset,g_loc, loff
      integer jj, kk, ll
c
      logical ocache_save
c
      logical odbug
      logical osome
c
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibasis, nshell)) call errquit
     $     ('wfn1_2e_rep: bas_numcont failed for ibasis', ibasis,
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
      call int_mem_2e4c(max2e, mem2)
cso
      if (.not. MA_push_get(MT_DBL,max2e,'wfn1_2e_rep:buf',l_buf,k_buf))
     $     call errquit('wfn1_2e_rep: ma failed', max2e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem2,'wfn1_2e_rep:scr',l_scr,k_scr))
     $     call errquit('wfn1_2e_rep: ma failed', mem2, MA_ERR)
c     
c     Loop thru shells
c     
      q2 = 1.0d0
      do lshell = 1, nshell
         if (.not. bas_cn2bfr(ibasis, lshell, llo, lhi)) 
     $      call errquit('wfn1_2e_rep: bas_cn2bfr', ibasis, BASIS_ERR)
         do kshell = 1, nshell
            if (.not. bas_cn2bfr(ibasis, kshell, klo, khi)) 
     $         call errquit('wfn1_2e_rep: bas_cn2bfr', ibasis,
     $                      BASIS_ERR)
            do jshell = 1, nshell
               if (.not. bas_cn2bfr(ibasis, jshell, jlo, jhi)) 
     $            call errquit('wfn1_2e_rep: bas_cn2bfr', ibasis,
     $                         BASIS_ERR)
               do ishell = 1, nshell
                  if (.not. bas_cn2bfr(ibasis, ishell, ilo, ihi)) 
     $               call errquit('wfn1_2e_rep: bas_cn2bfr', ibasis,
     $                            BASIS_ERR)
c
                  idim = ihi - ilo + 1
c     
c                 Generate the integrals
c     
                  call int_2e4c(ibasis, lshell, kshell,
     $                          ibasis, jshell, ishell,
     $                          mem2, dbl_mb(k_scr),
     $                          max2e, dbl_mb(k_buf))
c
c                 Put the integrals in the right place
c
                  loff = 0
                  do ll = llo, lhi
                    do kk = klo, khi
                      do jj = jlo, jhi
                        call dcopy(idim,dbl_mb(k_buf+loff),1,
     $                             eri(ilo,jj,kk,ll),1)
                        call dcopy(idim,dbl_mb(k_buf+loff),1,
     $                             erix(ilo,ll,kk,jj),1)
                        loff = loff + idim
                      enddo
                    enddo
                  enddo
               end do
            end do
         end do
      end do
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_1e_ga: pop failed', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_1e_ga: pop failed', 0, GA_ERR)
c     
      end
C> @}
