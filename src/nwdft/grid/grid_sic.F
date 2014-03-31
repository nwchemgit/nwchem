c     -----------------------------------------------------------------
      logical function grid_sicinit(nqpts,
     ,    l_vect1,i_vect1, ldelrho_sig,idelrho_sig)
c
C$Id$
c
      implicit none
c
#include "mafdecls.fh"
#include "grid_sic.fh"
      integer nqpts ! [in]
c
      integer len_delrho_sig
c
      grid_sicinit= MA_Push_Get(mt_dbl,nqpts,'vec aux',
     &                       l_vect1,i_vect1)
      call dfill(nqpts, 0.0d0, dbl_mb(i_vect1), 1)
      len_delrho_sig=3*nqpts
      grid_sicinit=grid_sicinit.and.
     .     MA_Push_Get(MT_Dbl,len_delrho_sig,'delrho_sig',
     &     ldelrho_sig,idelrho_sig)
      call dfill(len_delrho_sig, 0.0d0, dbl_mb(idelrho_sig), 1)
      return
      end
c
c     -----------------------------------------------------------------
      logical function grid_sicend(l_vect1,ldelrho_sig)
c
      implicit none
c
#include "mafdecls.fh"
#include "grid_sic.fh"
c
      grid_sicend= ma_pop_stack(ldelrho_sig)
      grid_sicend=grid_sicend.and.ma_pop_stack(l_vect1)
      return
      end
c
c     -----------------------------------------------------------------
      subroutine grid_sicpot(nqpts,qxyz,qwght,rho,delrho,amat,cmat,
     ,     delrho_sig,vect1)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
!#include "grid_sic.fh"
c
      integer nqpts
      double precision rho(*),delrho(*),qxyz(*),qwght(*)![in]
      double precision amat(*),cmat(*) ! [in/out]
      double precision delrho_sig(*),vect1(*)
c
      double precision tol,occup_i
      integer lAmat_coul,iAmat_coul
c
      if (.not.MA_Push_Get(mt_dbl,ipol*nqpts,'Amat_coul',
     &     lAmat_coul,iAmat_coul))
     &     call errquit('xc_quadv0_a: cannot allocate Amat_coul',0,
     &       MA_ERR)
      call dfill(nqpts, 0.d0, Dbl_MB(iAmat_coul), 1)
      tol = 1d-8
      occup_i = dble(sic_orb_occ)
      call potential_list(ao_bas_han, g_sic_dens(1), nqpts,
     &     qxyz, Dbl_MB(iAmat_coul), tol)
      call total_Amat(nqpts, ipol, Amat, Dbl_MB(iAmat_coul),
     &     qwght)
      if (.not.ma_pop_stack(lAmat_coul))
     &     call errquit('xc_quadv0_a: cannot pop stack',0, MA_ERR)
      call part_rho_rel(tol_rho, nqpts, Amat, Cmat, rho,
     &     vect1, ipol, occup_i,
     &     delrho, delrho_sig, 
     &     sic_orb_spin,aux_ik)
      return
      end
c
c     -----------------------------------------------------------------
      subroutine grid_sicrho(nqpts,mbf_ao,cetobfr,wmax,
     ,     rho,delrho,grad,ibf_ao, iniz, ifin,
     ,     chi_ao, delchi_ao,   rchi_atom, rdelchi_atom,
     ,     rdens_atom,
     ,     Fmat, Pmat, Bmat,
     ,     vect1,delrho_sig,
     &     ttau, kske)
c
      implicit none
c
#include "cdft.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer nqpts,mbf_ao ! [in]
      logical grad ! [in]
      double precision wmax ! [in]
      double precision cetobfr(*)
      double precision Fmat(*), Pmat(*), Bmat(*)
      double precision rho(*),delrho(*) ! [in]

      double precision ttau(*)
      double precision lapdum
      logical kske

      double precision chi_ao(*), delchi_ao(*),
     ,     rchi_atom(*), rdelchi_atom(*),rdens_atom(*)
      integer iniz(*),ifin(*),ibf_ao(*)
      double precision vect1(*),delrho_sig(*)
c
      integer idmat,ldmat
c
      call read_vect(rho, vect1, nqpts, 2, 
     &     sic_orb_spin)
      call read_delrho(delrho, delrho_sig, nqpts, ipol,
     &     sic_orb_spin)
!      if (do_scf) then
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &     'Dmat', lDmat,iDmat))
     &     call errquit('grid_sicr: cannot allocate Dmat',0,
     &       MA_ERR)
         call xc_rhogen(0,
     T     tol_rho, AO_bas_han, g_sic_dens, nbf_ao_mxnbf_ce, 
     &     ncenters, 1, 1, 1,
     I     ipol, nqpts, nbf_ao, mbf_ao, GRAD, 
     P     (ipol*(ipol+1))/2,
     F     Fmat, Pmat, Bmat,dbl_mb(idmat),
     &        chi_ao, delchi_ao,0d0,
     &        ibf_ao, iniz, ifin, 
     R     rho, delrho, lapdum, rchi_atom, 
     R     rdelchi_atom,  rdens_atom, cetobfr,wmax,
     &         ttau,kske, .false. )
         if (.not.ma_pop_stack(lDmat))
     &        call errquit('grid_iscrho: cannot pop stack', 1, MA_ERR)
!      endif
      return
      end
