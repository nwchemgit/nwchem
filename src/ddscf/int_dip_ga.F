      subroutine int_dip_ga(ibas, jbas, g_x, g_y, g_z)
*
* $Id: int_dip_ga.F,v 1.3 2003-10-17 22:54:33 carlfahl Exp $
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
      integer ibas, jbas
      integer g_x, g_y, g_z
      integer nshell_i, nshell_j
      integer ishell, jshell, iproc, nproc, mscratch, max1e
      integer ijshell, ilo, ihi, jlo, jhi, idim
      integer l_buf, l_scr
      integer k_buf, k_scr
c
      double precision center(3)
c
      call dfill(3, 0.0d0, center, 1)
c
      call ga_sync()
      call ga_zero(g_x)
      call ga_zero(g_y)
      call ga_zero(g_z)
      if(oscfps) call pstat_on(ps_int_dip)
c
c     grab basis set info type stuff
c
      if(.not.bas_numcont(ibas,nshell_i))call errquit
     $  ('int_dip_ga: bas_numcont failed for ibas',ibas, BASIS_ERR)
      if(.not.bas_numcont(jbas,nshell_j))call errquit
     $  ('int_dip_ga: bas_numcont failed for jbas',jbas, BASIS_ERR)
c
c allocate temp. arrays
c
      call int_mem_dipole(max1e,mscratch,ibas,jbas,1)
c      
      if(.not.MA_push_get(MT_DBL,max1e,'int_dip_ga:buf',l_buf,k_buf))
     $  call errquit('int_dip_ga: ma failure',max1e, MA_ERR)
      if(.not.MA_push_get(MT_DBL,mscratch,'int_dip_ga:scr',
     $   l_scr,k_scr))
     $  call errquit('int_dip_ga: ma failure',mscratch, MA_ERR)
c
c     loop over shells
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      do jshell = 1,nshell_j
         do ishell = 1,nshell_i
            if(mod(ijshell,nproc).eq.iproc)then
               if(.not.bas_cn2bfr(ibas,ishell,ilo,ihi))
     $              call errquit('int_dip_ga:bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
               if(.not.bas_cn2bfr(jbas,jshell,jlo,jhi))
     $              call errquit('int_dip_ga:bas_cn2bfr ?', jbas,
     &       BASIS_ERR)
               idim=ihi-ilo+1
c
c     get the integrals we want
c
               call int_mpole(ibas, ishell, jbas, jshell, 
     $              1,
     $              center,
     $              mscratch, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
c
c     bung integrals into global array
c
               call int_dip_put_in_g(g_x, g_y, g_z, ilo,ihi,jlo,jhi,
     $              dbl_mb(k_buf+(ihi-ilo+1)*(jhi-jlo+1)))
            endif
            ijshell=ijshell+1
         enddo
      enddo
c
      if(.not.MA_pop_stack(l_scr))
     $   call errquit('int_dip_ga:pop failed',0, MA_ERR)
      if(.not.MA_pop_stack(l_buf)) 
     $   call errquit('int_dip_ga:pop failed',0, MA_ERR)
c
      call ga_sync()
c
c      call ga_print(g_x)
c      call ga_print(g_y)
c      call ga_print(g_z)
c
      if(oscfps) call pstat_off(ps_int_dip)
c
      end          
      subroutine int_dip_put_in_g(g_x, g_y, g_z,
     $     ilo, ihi, jlo, jhi, buf)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
      integer g_x, g_y, g_z, ilo, ihi, jlo, jhi
      double precision buf(jlo:jhi,3,ilo:ihi)
      integer i, j
c
      do i = ilo,ihi
         do j = jlo, jhi
            call ga_put(g_x,i,i,j,j,buf(j,1,i),1)
            call ga_put(g_y,i,i,j,j,buf(j,2,i),1)
            call ga_put(g_z,i,i,j,j,buf(j,3,i),1)
         end do
      end do
c
      end
