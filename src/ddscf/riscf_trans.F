      subroutine riscf_trans_int (g_three, ao_basis, nsh, nbf, nsqhalf,
     $     bfsquare, npair, ri_basis, nsh_F, Nff, max_bf_ri, g_Vsqrt, 
     $     tol2e, g_indx, V_row, T_tmp, buf, lbuf, scr, lscr, odisk)

C$Id: riscf_trans.F,v 1.3 1996-02-07 00:41:32 vg038 Exp $
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "tcgmsg.fh"
#ifdef RISCF_USE_DA
#include "da.fh"
#endif

      integer g_three, ao_basis, ri_basis, g_Vsqrt, g_indx, s, sfirst, 
     $     slast, bfsquare, nsh, nbf, nsh_F, Nff, max_bf_ri, lbuf, lscr,
     $     npair

      double precision tol2e, V_row, T_tmp, buf, scr

      logical odisk

      dimension V_row (Nff, max_bf_ri)
      dimension T_tmp (Nff, bfsquare)
      dimension buf (lbuf)
      dimension scr (lscr)

      integer icount, i, j, ifirst, ilast, jfirst, jlast, nfi, nfj, nfs, 
     $     ij, nij, next, indx, nproc, ic, is, nsqhalf, ndone, g_tmp

      logical status, IeqJ

      double precision ddot
      external ddot

      nproc = ga_nnodes()

      if ( odisk ) then
C       first collect them in a temporary global array
        status = ga_create (MT_DBL, nproc*max_bf_ri*max_bf_ri, Nff, 
     $       '3c tmp', Nff, 1, g_tmp)
        if (.not. status ) call errquit ('couldn''t create 3c tmp ga',
     $       110)
      endif

      icount = -1
      ndone = 0
      next = nxtval(nproc)
      do i=1, nsh
        do j=1, i
          call ga_get (g_indx, i, i, j, j, indx, 1)
          if ( indx .ne. 0) then
            icount = icount + 1
            if ( next .eq. icount ) then
              
              IeqJ = (i .eq. j)
              status = bas_cn2bfr(ao_basis, i, ifirst, ilast)
              nfi = ilast - ifirst + 1
              if (IeqJ) then
                nfj = nfi
                nij = ((nfi+1) * nfi) / 2
              else
                status = bas_cn2bfr(ao_basis, j, jfirst, jlast)
                nfj = jlast - jfirst + 1
                nij = nfi * nfj
              endif

              do s=1, nsh_F
                status = bas_cn2bfr(ri_basis, s, sfirst, slast)
                nfs = slast - sfirst + 1
                call int_2e3c (ri_basis, s, ao_basis, i, j, lscr, scr, 
     $                         lbuf, buf)
                call put_int (buf, sfirst, nfs, nfi, nfj, 
     $               T_tmp, Nff, IeqJ)
              enddo

C             transform three center integrals
              ic = 0
              do s=1, nsh_F
                status = bas_cn2bfr(ri_basis, s, sfirst, slast)
                call ga_get (g_Vsqrt, 1, Nff, sfirst, slast, V_row, Nff)
                nfs = slast - sfirst + 1
                do is=1, nfs
                  do ij = 1, nij
                    ic = ic + 1
                    scr(ic) = ddot (Nff, V_row(1,is), 1, T_tmp(1,ij),1)
                  enddo
                enddo
              enddo ! s
              if ( odisk ) then
#ifdef RISCF_USE_DA
                call riscf_put_da ( g_three, scr, indx, nij, Nff, npair, 
     $               ndone, g_tmp )
#endif
              else
                call ga_put (g_three, indx, indx+nij-1,  1, Nff, scr, 
     $               nij)
              endif
              next = nxtval(nproc)
            endif ! next
          endif ! sparsity
        enddo   ! j
      enddo     ! i

#ifdef RISCF_USE_DA
      if ( odisk ) then
C       riscf_put_da is collective
C       nodes that already exited the loop have to take part
        if ( ndone .ne. npair ) then
          call riscf_put_da ( g_three, scr, indx, 0, Nff, npair, 
     $         ndone, g_tmp )
        endif
        status = ga_destroy ( g_tmp )
      endif
#endif

      next = nxtval(-nproc)

CDEBUG
C     call riscf_da_print ( g_three )
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine put_int (buf, sfirst, nfs, nfi, nfj, T_tmp, 
     $                    Nff, IeqJ)
      implicit none
      integer sfirst, nfs, nfi, nfj, Nff
      double precision buf, T_tmp
      logical IeqJ
      dimension T_tmp (Nff, *)
      dimension buf (nfj, nfi, nfs)

      integer ic, i, j, s

      ic = 0
      if (IeqJ) then
        do i=1, nfi
          do j=1, i
            ic = ic + 1
            do s=1, nfs
              T_tmp(sfirst-1+s, ic) = buf(j,i,s)
            enddo
          enddo
        enddo
      else
        do i=1, nfi
          do j=1, nfj
            ic = ic + 1
            do s=1, nfs
              T_tmp(sfirst-1+s, ic) = buf(j,i,s)
            enddo
          enddo
        enddo
      endif
      return
      end
