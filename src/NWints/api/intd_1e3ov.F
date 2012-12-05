c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute the 3 center overlap integral derivatives
C>
C> Compute the 3 center overlap integral derivatives as defined by
C> \f{eqnarray*}{
C>   \frac{\partial}{\partial X_x}({\mu}{\nu}{\lambda}) &=& \int_{-\infty}^{\infty}
C>   \frac{\partial g_{\mu}(X_\mu,r_{1})g_{\nu}(X_\nu,r_{1})g_{\lambda}(X_\lambda,r_{1})}{\partial X_x}dr_{1}
C> \f}
C> Where \f$X_x\f$ refers to a nuclear coordinate that can be any of 
C> \f$X_\mu\f$, \f$X_\nu\f$ or \f$X_\lambda\f$.
C>
C> The integral derivatives are returned in `dOV3`. The order in which the
C> integrals are stored is equivalent to the array being dimensioned as:
C> `dOV3(nint,ncoord,natom)` where `ncoord` is equal to 3 for the Cartesian
C> coordinates, and `natom` is equal to 3 for the atoms on which the Gaussians
C> are centered. The lexical indeces of the atoms are returned in `idatom`.
C> This relationship can be stated as: the integrals `dOV3(*,*,i)` belong to
C> atom `idatom(i)`.
C>
C> Finally, if any of the entries in `idatom` equals 0 then the data in `dOV3`
C> is invalid.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{intd\_1e3ov}
c:tex-This routine computes the 3 center overlap integral derivatives:
c:tex-\[
c:tex-\frac{\partial}{\partial q}({\mu}{\nu}{\lamda}) = \int_{-\infty}^{\infty}
c:tex-g_{\mu}(X,r_{1})g_{\nu}(X-R,r_{1})g_{\lamda}(X-R,r_{1})dr_{1}
c:tex-\]
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine intd_1e3ov (i_basis, ish, j_basis, jsh, k_basis, ksh,
     &       lscr, scr, ldov3, dOV3, idatom)
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
c:: external subroutines used
c..  errquit
c::functions
      logical int_chk_sh, int_chk_init
      external int_chk_sh, int_chk_init
c::passed
c:tex-\begin{verbatim}
      integer i_basis            !< [Input] basis set handle for ish
      integer ish                !< [Input] i shell/contraction
      integer j_basis            !< [Input] basis set handle for jsh
      integer jsh                !< [Input] j shell/contraction
      integer k_basis            !< [Input] basis set handle for ksh
      integer ksh                !< [Input] k shell/contraction
      integer lscr               !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer ldov3              !< [Input] length of 3c overlap integrals array
      double precision dOV3(*)   !< [Output] 3c overlap integrals
      integer idatom(3) !< [Output] array identifying centers for derivatives
c                       !  e.g., the first nint*3  derivatives go to center idatom(1)
c                       !        the second nint*3 derivatives go to center idatom(2)
c:tex-\end{verbatim}
c
c Order is...   nint*3*3 (3=> xyz, 3=atoms)
c
c  /                   |   
c | nint,d <ij>        |
c |      --------------|
c  \     d[idatom(1),x]|
c                          |
c       nint,d <ij>        |
c            --------------|
c            d[idatom(1),y]|
c                              |
c           nint,d <ij>        |
c                --------------|
c                d[idatom(1),z]|
c                                  |
c               nint,d <ij>        |
c                    --------------|
c                    d[idatom(2),x]|
c                                      |
c                   nint,d <ij>        |
c                        --------------|
c                        d[idatom(2),y]|
c                                          |
c                       nint,d <ij>        |
c                            --------------|
c                            d[idatom(2),z]|
c                                              |
c                           nint,d <ij>        |
c                                --------------|
c                                d[idatom(3),x]|
c                                                  |
c                               nint,d <ij>        |
c                                    --------------|
c                                    d[idatom(3),y]|
c                                                      \
c                                   nint,d <ij>         |
c                                        -------------- |
c                                        d[idatom(3),z]/
c If idatom(?) = 0 then do not use this 3*nint block of integral derivatives
c::local
      logical any_spherical
      logical shells_ok
      integer nint_intrnl
      integer ucont
      integer ibas, jbas, kbas
      integer i_geom, j_geom, k_geom
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
      integer Lk, k_prim, k_gen, k_iexp, k_icfp, k_cent
c
      integer i_pov3, j_pov3, k_pov3
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('intd_1e3ov'))
     &       call errquit('intd_1e3ov: int_init was not called' ,0,
     &            INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      shells_ok = shells_ok .and. int_chk_sh(k_basis,ksh)
      if (.not.shells_ok)
     &       call errquit('intd_1e3ov: invalid contraction/shell',0,
     &            BASIS_ERR)
c
      call int_nogencont_check(i_basis,'intd_1e3ov:i_basis')
      call int_nogencont_check(j_basis,'intd_1e3ov:j_basis')
      call int_nogencont_check(k_basis,'intd_1e3ov:k_basis')
      call int_nospshell_check(i_basis,'intd_1e3ov:i_basis')
      call int_nospshell_check(j_basis,'intd_1e3ov:j_basis')
      call int_nospshell_check(k_basis,'intd_1e3ov:k_basis')
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      kbas = k_basis + BASIS_HANDLE_OFFSET
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
c
      ucont   = (sf_ibs_cn2ucn(ksh,kbas))
      Lk      = infbs_cont(CONT_TYPE ,ucont,kbas)
      k_prim  = infbs_cont(CONT_NPRIM,ucont,kbas)
      k_gen   = infbs_cont(CONT_NGEN ,ucont,kbas)
      k_iexp  = infbs_cont(CONT_IEXP ,ucont,kbas)
      k_icfp  = infbs_cont(CONT_ICFP ,ucont,kbas)
      k_cent  = (sf_ibs_cn2ce(ksh,kbas))
      k_geom  = ibs_geom(kbas)
c
      nint_intrnl = int_nbf_x(Li)*int_nbf_x(Lj)*int_nbf_x(Lk)
      if ((nint_intrnl*9).gt.ldov3) then
        write(luout,*)' buffer size too small '
        write(luout,*)' buffer size : ',ldov3
        write(luout,*)' needed      : ',(nint_intrnl*9)
        stop ' intd_1e3ov: fatal error '
      endif
c
      if ((i_cent.eq.j_cent).and.(j_cent.eq.k_cent)) then
        idatom(1) = 0
        idatom(2) = 0
        idatom(3) = 0
        call dcopy((9*nint_intrnl),0.0d00,0,dOV3,1)
        return
      endif
c
      if ((i_geom.ne.j_geom.or.j_geom.ne.k_geom).and.WarnP.eq.0) then
        write(luout,*)
     &      'intd_1e3ov: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        write(luout,*)'k_basis geometry handle:',k_geom
        WarnP = 1
      endif
c
      call hfd3ois(
     &    coords(1,i_cent,i_geom),
     &    dbl_mb(mb_exndcf(i_iexp,ibas)),
     &    dbl_mb(mb_exndcf(i_icfp,ibas)),
     &    i_prim, 1, Li,
     &    coords(1,j_cent,j_geom),
     &    dbl_mb(mb_exndcf(j_iexp,jbas)),
     &    dbl_mb(mb_exndcf(j_icfp,jbas)),
     &    j_prim, 1, Lj,
     &    coords(1,k_cent,k_geom),
     &    dbl_mb(mb_exndcf(k_iexp,kbas)),
     &    dbl_mb(mb_exndcf(k_icfp,kbas)),
     &    k_prim, 1, Lk,
     &    dOV3,nint_intrnl,
c.........DryRun
     &    .false.,scr,lscr)
c
      i_pov3 = 1                    ! pointer to i block of derivs
      j_pov3 = 3*nint_intrnl+i_pov3 ! pointer to j block of derivs
      k_pov3 = 3*nint_intrnl+j_pov3 ! pointer to k block of derivs
c
      any_spherical = bas_spherical(ibas).or.
     &    bas_spherical(jbas).or.bas_spherical(kbas)
      if ((i_cent.ne.j_cent).and.(i_cent.ne.k_cent).and.
     &                                       (j_cent.ne.k_cent)) then
        idatom(1) = i_cent
        idatom(2) = j_cent
        idatom(3) = k_cent
        if (any_spherical) then
          if (Li.eq.-1) i_gen = 1
          if (Lj.eq.-1) j_gen = 1
          if (Lk.eq.-1) k_gen = 1
          call spcart_3cBtran(dOV3,scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        int_nbf_x(Lk),int_nbf_s(Lk),Lk,k_gen,bas_spherical(kbas),
     &        9,.false.)
        endif
      else if (i_cent.eq.j_cent) then
        call daxpy(3*nint_intrnl,1.0d00,dOV3(j_pov3),1,dOV3(i_pov3),1)
        call dcopy(3*nint_intrnl,dOV3(k_pov3),1,dOV3(j_pov3),1)
        call dcopy(3*nint_intrnl,0.0d00,0,dOV3(k_pov3),1)
        idatom(1) = i_cent
        idatom(2) = k_cent
        idatom(3) = 0
        if (any_spherical) then
          if (Li.eq.-1) i_gen = 1
          if (Lj.eq.-1) j_gen = 1
          if (Lk.eq.-1) k_gen = 1
          call spcart_3cBtran(dOV3,scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        int_nbf_x(Lk),int_nbf_s(Lk),Lk,k_gen,bas_spherical(kbas),
     &        6,.false.)
        endif
      else if (i_cent.eq.k_cent) then
        call daxpy(3*nint_intrnl,1.0d00,dOV3(k_pov3),1,dOV3(i_pov3),1)
        call dcopy(3*nint_intrnl,0.0d00,0,dOV3(k_pov3),1)
        idatom(1) = i_cent
        idatom(2) = j_cent
        idatom(3) = 0
        if (any_spherical) then
          if (Li.eq.-1) i_gen = 1
          if (Lj.eq.-1) j_gen = 1
          if (Lk.eq.-1) k_gen = 1
          call spcart_3cBtran(dOV3,scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        int_nbf_x(Lk),int_nbf_s(Lk),Lk,k_gen,bas_spherical(kbas),
     &        6,.false.)
        endif
      else if (j_cent.eq.k_cent) then
        call daxpy(3*nint_intrnl,1.0d00,dOV3(k_pov3),1,dOV3(j_pov3),1)
        call dcopy(3*nint_intrnl,0.0d00,0,dOV3(k_pov3),1)
        idatom(1) = i_cent
        idatom(2) = j_cent
        idatom(3) = 0
        if (any_spherical) then
          if (Li.eq.-1) i_gen = 1
          if (Lj.eq.-1) j_gen = 1
          if (Lk.eq.-1) k_gen = 1
          call spcart_3cBtran(dOV3,scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        int_nbf_x(Lk),int_nbf_s(Lk),Lk,k_gen,bas_spherical(kbas),
     &        6,.false.)
        endif
      else
        write(luout,*)'ijk->centers',i_cent,j_cent,k_cent
        write(luout,*)' fatal error '
        stop 'intd_1e3ov: how did I get here'
      endif
      end
C> @}
