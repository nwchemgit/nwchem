      subroutine tce_cisd_energy(d_f1,d_v2,d_t1,d_t2,d_r1,d_r2,d_e,
     1           k_f1_offset,k_v2_offset,k_t1_offset,k_t2_offset,
     1           k_r1_offset,k_r2_offset,k_e_offset,
     1           size_e,size_r1,size_r2,
     1           size_t1,size_t2,ref,corr)
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
      integer d_f1,d_v2,d_t1,d_t2,d_r1,d_r2,d_e 
      integer k_f1_offset,k_v2_offset,k_t1_offset,k_t2_offset
      integer k_r1_offset,k_r2_offset,k_e_offset
      integer size_e,size_r1,size_r2,size_t1,size_t2
      double precision ref,corr
c
c      integer iter
      double precision cpu, wall
      double precision r1,r2
      double precision residual
      logical nodezero
      integer dummy
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
      call tce_diis_init()
      do iter=1,maxiter
         cpu=-util_cpusec()
         wall=-util_wallsec()
         if (nodezero.and.(iter.eq.1)) write(LuOut,9050) "CISD"
         call tce_filename('e',filename)
         call createfile(filename,d_e,size_e)
         call cisd_e(d_f1,d_e,d_t1,d_t2,d_v2,
     1               k_f1_offset,k_e_offset,
     2               k_t1_offset,k_t2_offset,k_v2_offset)
         call reconcilefile(d_e,size_e)
         call tce_filename('r1',filename)
         call createfile(filename,d_r1,size_r1)
         call cisd_c1(d_e,d_f1,d_r1,d_t1,d_t2,d_v2,
     1                k_e_offset,k_f1_offset,k_r1_offset,
     2                k_t1_offset,k_t2_offset,k_v2_offset)
         call reconcilefile(d_r1,size_r1)
         call tce_filename('r2',filename)
         call createfile(filename,d_r2,size_r2)
         call cisd_c2(d_e,d_f1,d_r2,d_t1,d_t2,d_v2,
     1                k_e_offset,k_f1_offset,k_r2_offset,
     2                k_t1_offset,k_t2_offset,k_v2_offset)
         call reconcilefile(d_r2,size_r2)
         call tce_residual_t1(d_r1,k_r1_offset,r1)
         call tce_residual_t2(d_r2,k_r2_offset,r2)
         residual = max(r1,r2)
         call get_block(d_e,corr,1,0)
         cpu=cpu+util_cpusec()
         wall=wall+util_wallsec()
         if (nodezero) write(LuOut,9100) iter,residual,corr,cpu,wall
         if (residual .lt. thresh) then
            if (nodezero) then
               write(LuOut,9060)
               write(LuOut,9070) "CISD",corr
               write(LuOut,9080) "CISD",ref + corr
            endif
            call deletefile(d_r2)
            call deletefile(d_r1)
            call deletefile(d_e)
            call tce_diis_tidy()
            goto 100
         endif
         call tce_diis(.false.,iter,.true.,.true.,.false.,.false.,
     1                 d_r1,d_t1,k_t1_offset,size_t1,
     2                 d_r2,d_t2,k_t2_offset,size_t2,
     3                 dummy,dummy,dummy,dummy,
     4                 dummy,dummy,dummy,dummy)
         call deletefile(d_r2)
         call deletefile(d_r1)
         call deletefile(d_e)
         if (nodezero) call util_flush(LuOut)
      enddo
      call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 100  continue
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
