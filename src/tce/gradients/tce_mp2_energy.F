      subroutine tce_mp2_energy(d_f1,d_v2,d_u1,d_u2,d_r1,d_r2,d_e,
     1           k_f1_offset,k_v2_offset,k_u1_offset,k_u2_offset,
     1           k_r1_offset,k_r2_offset,k_e_offset,
     1           size_e,size_r1,size_r2,
     1           size_u1,size_u2,ref,corr)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
c
      integer d_f1,d_v2,d_u1,d_u2,d_r1,d_r2,d_e 
      integer k_f1_offset,k_v2_offset,k_u1_offset,k_u2_offset
      integer k_r1_offset,k_r2_offset,k_e_offset
      integer size_e,size_r1,size_r2,size_u1,size_u2
      double precision ref,corr
c
c      integer iter
      double precision cpu, wall
      double precision r1,r2
      double precision residual,mbpt2
      logical nodezero
      integer dummy
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
      call tce_diis_init()
      do iter=1,maxiter
         cpu=-util_cpusec()
         wall=-util_wallsec()
         if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "MBPT(2)"
         call tce_filename('e',filename)
         call createfile(filename,d_e,size_e)
         call mbpt2_e(d_e,d_u2,d_v2,
     1                k_e_offset,k_u2_offset,k_v2_offset)
         call reconcilefile(d_e,size_e)
         call tce_filename('r1',filename)
         call createfile(filename,d_r1,size_r1)
         call mbpt1_t1(d_f1,d_r1,d_u1,d_u2,k_f1_offset,k_r1_offset,
     1                 k_u1_offset,k_u2_offset)
         call reconcilefile(d_r1,size_r1)
         call tce_filename('r2',filename)
         call createfile(filename,d_r2,size_r2)
         call mbpt1_t2(d_f1,d_r2,d_u1,d_u2,d_v2,k_f1_offset,
     1                 k_r2_offset,k_u1_offset,k_u2_offset,
     2                 k_v2_offset)
         call reconcilefile(d_r2,size_r2)
         call tce_residual_t1(d_r1,k_r1_offset,r1)
         call tce_residual_t2(d_r2,k_r2_offset,r2)
         residual = max(r1,r2)
         call get_block(d_e,mbpt2,1,0)
         corr = mbpt2
         cpu=cpu+util_cpusec()
         wall=wall+util_wallsec()
         if (nodezero) write(LuOut,9100) iter,residual,mbpt2,cpu,wall
         if (residual .lt. thresh) then
           if (nodezero) then
             write(LuOut,9060)
             write(LuOut,9070) "MBPT(2)",mbpt2
             write(LuOut,9080) "MBPT(2)",ref + corr
             call util_flush(LuOut)
           endif
           call deletefile(d_r2)
           call deletefile(d_r1)
           call deletefile(d_e)
           call tce_diis_tidy()
           goto 120
         endif
         call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                  d_r1,d_u1,k_u1_offset,size_u1,
     2                  d_r2,d_u2,k_u2_offset,size_u2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
         call deletefile(d_r2)
         call deletefile(d_r1)
         call deletefile(d_e)
         if (nodezero) call util_flush(LuOut)
      enddo
      call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 120  continue
 9050 format(/,1x,A,' iterations',/,
     1  1x,'--------------------------------------------------------',/
     2  1x,'Iter          Residuum       Correlation     Cpu    Wall',/
     3  1x,'--------------------------------------------------------')
 9060 format(
     1  1x,'--------------------------------------------------------',/
     2  1x,'Iterations converged')
 9070 format(1x,A,' correlation energy / hartree = ',f25.15)
 9080 format(1x,A,' total energy / hartree       = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
      end

c $Id$
