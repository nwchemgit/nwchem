C
C     zfock_os_dftxc.F
C
C     Calculates DFT exchange-correlation part of a open shell complex
C     Fock matrix.  Note that input is the real part of the density
C     matrix, and output is real part of the Fock matrix, and note that
C     unlike other zfock routines this one accumulates the result in Exc
C     and g_dftxc.  Also, note alpha beta parts.
C
      subroutine zfock_os_dftxc (params, g_densre, Exc, g_dftvxc)
      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre(2)       !alpha, beta

      
C     == Outputs ==
      double precision, intent(out) :: Exc(2)       !exch, corr (if separate) *not* alpha,beta
      integer, intent(in)           :: g_dftvxc(2)  !alpha, beta

      
C     == Parameters ==
      character*15, parameter     :: pname = "zfock_os_dftxc:"

      
C     == Variables ==
C      integer g_densre_clones(2)
      integer lTmat, iTmat, nTrows, nTcols
      double precision elapsed


C      call rt_tddft_os_confirm (params,'junk')


      if (params%prof) call prof_start (elapsed)


C     
C     DFT depends XXXX
C

C      g_densre_clones(1) = g_densre
C      g_densre_clones(2) = g_densre


      
C     == Set up stuff if we are fitting XC ==
      if(xcfit) then

         call errquit (pname//"xcfit might not be working", 0, 0)
         
         nTrows = nbf_xc
         nTcols = ipol
         
         if (.not.ma_push_get(MT_Dbl,nTrows*nTcols,
     $        'Tmat',lTmat, iTmat))
     $        call errquit(pname//'cannot allocate Tmat',0,MA_ERR)
         
         call dfill (nTrows*nTcols, 0d0, dbl_mb(iTmat), 1)
      endif

      
C
C     Need to zero out GA as the DFT routine accumulates result rather
C     than replaces input.
C
      call ga_zero (g_dftvxc(1))
      call ga_zero (g_dftvxc(2))
      
C     (technically g_dftvxc should have 3,4 elements but zero in this
C     case, so we just send the length two array to grid_quadv0)

      
C     == Actual call to DFT solver ==
      call grid_quadv0 (params%rtdb, g_densre, g_dftvxc,
     $     params%nExc, params%rho_n, Exc, dbl_mb(itmat))


C     == If fitting XC ==
      if(xcfit) then

c     symmetrize the "T" vector
         if (oskel)
     $        call sym_vec_symmetrize(geom,xc_bas_han,Dbl_MB(iTmat))
         
         call xc_fitv(params%rtdb, Dbl_MB(iTmat), nTrows, nTcols,
     $        params%g_vxc, params%g_xcinv, params%iolgc)
         
         if (.not.ma_pop_stack(lTmat))
     $        call errquit(pname//"cannot pop stack",0, MA_ERR)
      endif

      if (params%prof)  call prof_end (elapsed, "Fock OS dft xc")
      
      end subroutine

c $Id$
