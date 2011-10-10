      subroutine dft_pstrace(g_dens, basis, nbf, oskel)
      
c     calculate tr[P S] = number of electrons, from an input
c     AO density matrix

c     Note: integrals are assumed to be initialized when this routine
c     is executed

      implicit none

#include "errquit.fh"
#include "bas.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"

c     subroutine arguments:
      integer g_dens ! density matrix GA handle
      integer basis  ! basis set handle
      integer nbf    ! number of AO basis fcts.
      logical oskel  ! whether to use skeleton matrices

c     local variables and arrays
      integer g_smat, g_temp
      double precision sum
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      double precision ga_trace_diag
      external ga_trace_diag

      integer me

c     ========================================================================

      me = ga_nodeid()

      if (me.eq.0) write (luout,*) 'dft_pstrace: checking P-matrix'

c     call int_init(rtdb,1,ao_bas_han)
c     call schwarz_init(geom,ao_bas_han)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 1
      if (.not.nga_create(MT_DBL,2,ahi,'smat',alo,g_smat)) call 
     &   errquit('dft_pstrace: nga_create failed g_smat',0,GA_ERR)
      call ga_zero(g_smat)
      call int_1e_ga(basis, basis ,g_smat,'overlap',oskel)

c     debug: multiply density matrix with overlap matrix
c     and extract diagonal elements => N(electrons)
        
      if (.not.nga_create(MT_DBL,2,ahi,'temp',alo,g_temp)) call 
     &   errquit('dft_pstrace: nga_create failed g_temp',0,GA_ERR)
      
      alo(1) = 1                ! range of A, first dimension 
      ahi(1) = nbf
      alo(2) = 1                ! range of A, second dimension
      ahi(2) = nbf
      blo(1) = 1                ! range of B, first dimension, etc.
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      clo(1) = 1
      chi(1) = nbf
      clo(2) = 1
      chi(2) = nbf

      call ga_zero(g_temp)
      call nga_matmul_patch('n','n',1d0,0d0,
     &   g_dens,alo,ahi,
     &   g_smat,blo,bhi,
     &   g_temp,clo,chi)
      sum = ga_trace_diag(g_temp)

      if (me.eq.0) write (luout,'(1x,a,1x,e15.7)')
     &       'tr(P*S): ',sum

      if (.not.ga_destroy(g_temp)) call 
     &   errquit('dft_scf: ga_destroy failed g_temp',0,GA_ERR)
      if (.not.ga_destroy(g_smat)) call 
     &   errquit('dft_scf: ga_destroy failed g_smat',0,GA_ERR)
c     call schwarz_tidy()
c     call int_terminate()
      
      return
      
      end
c $Id$
