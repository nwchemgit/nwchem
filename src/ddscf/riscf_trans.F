      subroutine riscf_trans_int (g_three, ao_basis, nsh, nbf, nsqhalf,
     $     bfsquare, ri_basis, nsh_F, Nff, max_bf_ri, g_Vsqrt, 
     $     tol2e, g_indx, V_row, T_tmp, buf, lbuf, scr, lscr)

C$Id: riscf_trans.F,v 1.2 1995-12-29 02:04:15 vg038 Exp $
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "tcgmsg.fh"

      integer g_three, ao_basis, ri_basis, g_Vsqrt, g_indx, s, sfirst, 
     $     slast, bfsquare, nsh, nbf, nsh_F, Nff, max_bf_ri, lbuf, lscr

      double precision tol2e, V_row, T_tmp, buf, scr

      dimension V_row (Nff, max_bf_ri)
      dimension T_tmp (Nff, bfsquare)
      dimension buf (lbuf)
      dimension scr (lscr)

      integer icount, i, j, ifirst, ilast, jfirst, jlast, nfi, nfj, nfs, 
     $     ij, nij, next, indx, nproc, ic, is, nsqhalf

      logical status, IeqJ

      double precision ddot
      external ddot

      nproc = ga_nnodes()

      icount = -1
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
              do s=1, Nsh_f
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
              call ga_put (g_three, indx, indx+nij-1,  1, Nff, scr, nij)
              
              next = nxtval(nproc)
            endif ! next
          endif ! sparsity
        enddo   ! j
      enddo     ! i

      next = nxtval(-nproc)

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
