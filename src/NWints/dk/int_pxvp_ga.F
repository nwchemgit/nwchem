C> \ingroup nwint
C> @{
C>
C> \brief Gather pxVp integrals and put them in a global array
C>
C> Gather pxVp integrals from int_1epxvpe on a shell by shell basis
C> and put them per component in a GA.
C>
      subroutine int_pxvp_ga(ibas, g_x, g_y, g_z, oskel)
*
* $Id$
*
c
c     Gather pxVp integrals from int_1epxvpe on a shell by shell basis
c     and put them per component in a GA.
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
#include "errquit.fh"
c
c     arguments
c
      integer ibas !< [Input] The basis set handle
      integer g_x  !< [Output] The \f$x\f$ components
      integer g_y  !< [Output] The \f$y\f$ components
      integer g_z  !< [Output] The \f$z\f$ components
      logical oskel !< [Input] Skeletonize the integrals
c
c     local variables
c
      integer nshell_i, idim
      integer ishell, jshell, iproc, nproc, mem1, max1e
      integer ijshell, ilo, ihi, jlo, jhi
      integer l_buf, l_scr
      integer k_buf, k_scr
      logical odoit
      double precision q2
c
      call ga_sync()
      call ga_zero(g_x)
      call ga_zero(g_y)
      call ga_zero(g_z)
c     if(oscfps) call pstat_on(ps_int_pxvp)
c
c     grab basis set info type stuff
c
      if(.not.bas_numcont(ibas,nshell_i))call errquit
     $  ('int_pxvp_ga: bas_numcont failed for ibas',ibas, BASIS_ERR)
c
c      Allocate temporary arrays
c
      call int_mem_1e(max1e, mem1)
c     call emem_1e_pvp(ibas, ibas, max1e, mem1, 3)
      if (.not. MA_push_get(MT_DBL,max1e*3,'int_pxvp_ga:buf',l_buf,
     $    k_buf)) call errquit('int_pxvp_ga: ma failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL, mem1,'int_pxvp_ga:scr',l_scr,k_scr))
     $     call errquit('int_pxvp_ga: ma failed', mem1, MA_ERR)
c      
c     loop over shells
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      q2 = 1.0d0
      do jshell = 1,nshell_i
         do ishell = 1,nshell_i
            if(mod(ijshell, nproc) .eq. iproc) then
               odoit = .true.
               if (oskel)
     $              odoit = sym_shell_pair(ibas, ishell, jshell, q2)
               if (odoit) then
                  if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi))
     $                 call errquit('int_pxvp_ga: bas_cn2bfr ?', ibas,
     &       GA_ERR)
                  if (.not. bas_cn2bfr(ibas, jshell, jlo, jhi))
     $                 call errquit('int_pxvp_ga: bas_cn2bfr ?', ibas,
     &       BASIS_ERR)
                  idim = ihi - ilo + 1
c
c     Generate the integrals
c
                  call int_1epxvpe (ibas, jshell, ibas, ishell,
     $                    mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
c
c     Add integrals into the global arrays
c
                  call ga_acc(g_z, ilo, ihi, jlo, jhi, 
     $                        dbl_mb(k_buf), idim, q2)
                  call ga_acc(g_y, ilo, ihi, jlo, jhi, 
     $                        dbl_mb(k_buf+max1e), idim, q2)
                  call ga_acc(g_x, ilo, ihi, jlo, jhi, 
     $                        dbl_mb(k_buf+(max1e*2)), idim, q2)
               endif
            endif
            ijshell = ijshell + 1
         enddo
      enddo
      call ga_sync()
c
      if(.not.MA_pop_stack(l_scr))
     $   call errquit('int_pxvp_ga:pop failed',0, GA_ERR)
      if(.not.MA_pop_stack(l_buf)) 
     $   call errquit('int_pxvp_ga:pop failed',0, GA_ERR)
c
      call ga_sync()
c
c     if(oscfps) call pstat_off(ps_int_pxvp)
c
      end          
C>
C> @}
