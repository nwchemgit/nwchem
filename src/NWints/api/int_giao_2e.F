c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Calculate the 4-center 2-electron integrals between GIAO basis functions
C>
C> The GIAO basis functions are important in the treatment of magnetic
C> properties. These properties are not automatically gauge invariant but
C> by multiplying the usual Atomic Orbitals (AO) with a magnetic field dependent
C> factor Gauge Invariant Atomic Orbitals (GIAO) can be defined. Using these
C> basis functions the properties can be cast in an gauge invariant form (see
C> [1,2]).
C> 
C> In the presence of a magnetic field a vector potential is generated given by
C> \f{eqnarray*}{
C>   \vec{A}_A = \vec{B} \otimes \vec{R}_A
C> \f}
C> where \f$ \vec{B} \f$ is the magnetic field, \f$ \vec{R}_A \f$ is the
C> nuclear position relative to the origin of the systems,
C> \f$ \vec{A} \f$ the generated vector potential.
C> A GIAO is now given by
C> \f{eqnarray*}{
C>    g^{GIAO}_\mu(X_\mu,r) = g_\mu(X_\mu,r)e^{-i \vec{A}_A\cdot\vec{r}}
C> \f}
C> The integrals computed here can now be represented as
C> \f{eqnarray*}{
C>    ({\mu}{\rho}|{\nu}{\lambda}) = \int_{-\infty}^{\infty}
C>    g^{GIAO}_{\mu}(X_{\mu},r_{1})g^{GIAO}_{\rho}(X_{\rho},r_{1})\frac{1}{r_{12}}
C>    g^{GIAO}_{\nu}(X_{\nu},r_{2})g^{GIAO}_{\lambda}(X_{\lambda},r_{2})dr_{1}dr_{2}
C> \f}
C>
C> [1] M. Dupuis,
C>     "New integral transforms for molecular properties and application
C>     to a massively parallel GIAO-SCF implementation",
C>     Comp. Phys. Comm. <b>134</b>, 150-166 (2001), DOI:
C>     <a href="https://doi.org/10.1016/S0010-4655(00)00195-8">
C>     10.1016/S0010-4655(00)00195-8</a>.
C>
C> [2] K. Wolinski, J.F. Hinton, P. Pulay,
C>     "Efficient implementation of the gauge-independent atomic orbital
C>      method for NMR chemical shift calculations",
C>     J. Am. Chem. Soc. <b>112</b>, 8251-8260 (1990), DOI:
C>     <a href="https://doi.org/10.1021/ja00179a005">10.1021/ja00179a005</a>.
C>
      subroutine int_giao_2e(brain, ish, jsh, ketin, ksh, lsh,
     &       lscr, scr, leri, eri)
      implicit none
c
c basic api routine to generate a block of two electron integrals
c eri = <bra_g(ish).bra_g(jsh) | ket_g(ksh).ket_g(lsh)>
c
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "apiP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "rel_nwc.fh"
#include "util.fh"
#include "global.fh"
      common/testdata/timing(20),irepeat
      double precision timing
      integer irepeat
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical cando_sp
      logical cando_txs
      external cando_nw
      external cando_sp
      external cando_txs
c:: passed
      integer brain !< [Input] bra basis set handle
      integer ish   !< [Input] shell/contraction index
      integer jsh   !< [Input] shell/contraction index
      integer ketin !< [Input] ket basis set handle
      integer ksh   !< [Input] shell/contraction index
      integer lsh   !< [Input] shell/contraction index
      integer lscr  !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] array
      integer leri  !< [In|Output] length of integral array
      double precision eri(*) !< [Output] 2e4c integrals
c:: local
      integer bra, ket 
      integer ab_geom, cd_geom, ucont, ityp
      integer La, a_prim, a_gen, a_iexp, a_icfp, a_cent, a_icfps
      integer Lb, b_prim, b_gen, b_iexp, b_icfp, b_cent, b_icfps
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent, c_icfps
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent, d_icfps
c
      double precision roff(3)
      double precision q4
      integer nint
      logical dum_log, do_bra, do_ket
      logical status_sp, status_nw, status_txs, status_gen
      integer texas_ang_limit
      logical ieqj, keql
      integer sbas, abas, bras, kets
c
      logical any_spherical
      integer a_nbf,   b_nbf,   c_nbf,   d_nbf
      integer a_nbf_s, b_nbf_s, c_nbf_s, d_nbf_s
      integer ab_gen, ab_cmp, cd_gen, cd_cmp,i,j
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c     timing(10)=timing(10)-util_wallsec()
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      bras = bra
      kets = ket
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
      a_cent  = (sf_ibs_cn2ce(ish,bra))
      b_cent  = (sf_ibs_cn2ce(jsh,bra))
      c_cent  = (sf_ibs_cn2ce(ksh,ket))
      d_cent  = (sf_ibs_cn2ce(lsh,ket))
c
c
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)
c
      status_sp = cando_sp(brain,ish,jsh).and.cando_sp(ketin,ksh,lsh)
      if (.not.status_sp) then
c
        ieqj = ish.eq.jsh
        keql = ksh.eq.lsh
c
        ucont   = sf_ibs_cn2ucn(ish,bra)
        La      = infbs_cont(CONT_TYPE ,ucont,bra)
        a_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        a_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        a_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
        a_icfps = infbs_cont(CONT_ICFP ,ucont,bras)
c
        ucont   = sf_ibs_cn2ucn(jsh,bra)
        Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
        b_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        b_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        b_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
        b_icfps = infbs_cont(CONT_ICFP ,ucont,bras)
c      
        ucont   = sf_ibs_cn2ucn(ksh,ket)
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        c_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c
        ucont   = sf_ibs_cn2ucn(lsh,ket)
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
        d_icfps = infbs_cont(CONT_ICFP ,ucont,kets)
c
        a_nbf = int_nbf_x(La)
        b_nbf = int_nbf_x(Lb)
        c_nbf = int_nbf_x(Lc)
        d_nbf = int_nbf_x(Ld)
        leri=a_nbf*b_nbf*c_nbf*d_nbf*a_gen*b_gen*c_gen*d_gen
        call hnd_giahnd(
     &      coords(1,a_cent,ab_geom), dbl_mb(mb_exndcf(a_iexp,bra)),
     &      dbl_mb(mb_exndcf(a_icfp,bra)), a_prim, a_gen, La,
     &      coords(1,b_cent,ab_geom), dbl_mb(mb_exndcf(b_iexp,bra)),
     &      dbl_mb(mb_exndcf(b_icfp,bra)), b_prim, b_gen, Lb,
     &      coords(1,c_cent,cd_geom), dbl_mb(mb_exndcf(c_iexp,ket)),
     &      dbl_mb(mb_exndcf(c_icfp,ket)), c_prim, c_gen, Lc,
     &      coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &      dbl_mb(mb_exndcf(d_icfp,ket)), d_prim,d_gen,Ld,
     &      ieqj, keql, eri, leri, scr, lscr)
c
c eri has cartesian block of integrals (llo:lhi,klo:khi,jlo:jhi,ilo:ihi,6)
c
        if (any_spherical) then
          a_nbf_s = int_nbf_s(La)
          b_nbf_s = int_nbf_s(Lb)
          c_nbf_s = int_nbf_s(Lc)
          d_nbf_s = int_nbf_s(Ld)
          cd_gen = c_gen*d_gen
          ab_gen = a_gen*b_gen
          ab_cmp = a_nbf*b_nbf
          cd_cmp = c_nbf*d_nbf
          do_bra=bas_spherical(bra)
          do_ket=bas_spherical(ket)
          call giao_to_sph(eri,leri,scr,La,Lb,Lc,Ld,a_nbf,b_nbf,c_nbf,
     &                     d_nbf,a_nbf_s,b_nbf_s,c_nbf_s,d_nbf_s,
     &                     a_gen,b_gen,c_gen,d_gen,ab_gen,ab_cmp,
     &                     cd_gen,cd_cmp,do_bra,do_ket)
        endif
      else
        write(luout,*)'int_giao_2e: cannot do sp integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        write(luout,*)' ksh   :',ksh
        write(luout,*)' lsh   :',lsh
        call errquit('int_giao_2e: fatal error ',0, INT_ERR)
      endif
      end
c
      subroutine giao_to_sph(eri,leri,scr,La,Lb,Lc,Ld,a_nbf,b_nbf,c_nbf,
     &                     d_nbf,a_nbf_s,b_nbf_s,c_nbf_s,d_nbf_s,
     &                     a_gen,b_gen,c_gen,d_gen,ab_gen,ab_cmp,
     &                     cd_gen,cd_cmp,do_bra,do_ket)
      implicit none
c
      integer leri
      double precision eri(leri,6), scr(*)
      integer La,Lb,Lc,Ld,a_nbf,b_nbf,c_nbf,d_nbf
      integer a_nbf_s,b_nbf_s,c_nbf_s,d_nbf_s,a_gen,b_gen,c_gen,d_gen
      integer ab_gen,ab_cmp,cd_gen,cd_cmp
      logical do_bra,do_ket
c
      integer ityp
c
      ab_cmp = a_nbf_s*b_nbf_s
      do ityp = 1, 6
         call spcart_bra2etran(eri(1,ityp),scr,
     &        b_nbf,a_nbf,b_nbf_s,a_nbf_s,
     &        Lb, La, b_gen, a_gen,
     &        cd_gen*cd_cmp,.false.)
         call spcart_ket2etran(eri(1,ityp),scr,
     &        d_nbf,c_nbf,d_nbf_s,c_nbf_s,
     &        Ld, Lc, d_gen, c_gen,
     &        ab_gen*ab_cmp,.false.)
      enddo
      a_nbf = a_nbf_s
      b_nbf = b_nbf_s
      c_nbf = c_nbf_s
      d_nbf = d_nbf_s
      end
c
C> @}
