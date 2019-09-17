C> \ingroup nwint
C> @{
C>
C> \brief Compute 2-electron 2-center Gaussian periodic electron
C> repulsion integral derivatives
C>
C> The code in question always uses density fitting approaches for
C> the 2-center 2-electron integrals are needed. See [1] for details.
C>
C> [1] JE Jaffe, AC Hess,
C>     <i>"Gaussian basis density functional theory for systems
C>     periodic in two or three dimensions: Energy and forces"</i>,
C>    J.Chem.Phys. <b>105</b>, 10983-10998 (1996), DOI:
C>    <a href="https://doi.org/10.1063/1.472866">
C>    10.1063/1.472866</a>
C>
      subroutine intpd_2e2c(brain, ish, ketin, jsh, Rj,
     &       lscr, scr, leri, eri, idatom2)
*
* $Id$
*
      implicit none
c
c basic api routine to generate a block of 2 center two electron integrals
c eri = <bra_g(ish)|ket_g(jsh)>
c
#include "apiP.fh"
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical int_chk_sh
      logical int_chk_init
      external cando_nw
      external int_chk_sh
      external int_chk_init
      integer int_nint_cart
      external int_nint_cart
c:: passed
      integer brain  !< [Input] basis set handle for bra function
      integer ketin  !< [Input] basis set handle for ket function
      integer ish    !< [Input] lexical contraction/shell index
      integer jsh    !< [Input] lexical contraction/shell index
      integer lscr   !< [Input] length of scratch array
      integer leri   !< [Input] length of ERI array
      integer idatom2(*) !< [Output]  array identifying centers for derivatives externalized
      double precision scr(lscr) !< [Scratch] scratch array 
      double precision eri(leri) !< [Output] ERI array
      double precision Rj(3)  !< [Input] translation vector for jsh center (fractional coordinates)
c:: local
      integer idatom(4)
      logical shells_ok
      integer bra, ket 
      integer p_geom, q_geom, ucont, mynint
      integer Lp, p_prim, p_gen, p_iexp, p_icfp, p_cent
      integer Lq, q_prim, q_gen, q_iexp, q_icfp, q_cent
      double precision xyz_new_j(3)  ! new coordinates for jsh function center
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization
c
      if (.not.int_chk_init('intpd_2e2c'))
     &       call errquit('intpd_2e2c: int_init was not called' ,0,
     &           INT_ERR)
c
c     check input shell ranges
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
c
      if (.not.shells_ok)
     &       call errquit('intpd_2e2c: invalid contraction/shell',0,
     &            BASIS_ERR)
c
c  check if gencon/sp shell
c
      call int_nogencont_check(brain,'intpd_2e2c:bra')
      call int_nogencont_check(ketin,'intpd_2e2c:ket')
      call int_nospshell_check(brain,'intpd_2e2c:bra')
      call int_nospshell_check(ketin,'intpd_2e2c:ket')
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      p_geom = ibs_geom(bra)
      q_geom = ibs_geom(ket)
c
      if (p_geom.ne.q_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'intpd_2e2c: WARNING: possible geometry inconsistency'
        write(luout,*)'bra geometry handle:',p_geom
        write(luout,*)'ket geometry handle:',q_geom
        WarnP = 1
      endif
      mynint = int_nint_cart(brain,ish,ketin,jsh,0,0,0,0)
      if (mynint*3*4.gt.leri) then
        write(luout,*) 'nint*3*4 = ',mynint*3*4
        write(luout,*) 'leri     = ',leri
        call errquit('intpd_2e2c: nint>leri error',911, INT_ERR)
      endif
c
      p_cent  = (sf_ibs_cn2ce(ish,bra))
      q_cent  = (sf_ibs_cn2ce(jsh,ket))
      if (p_cent.eq.q_cent) then
        call dcopy((mynint*3*4),0.0d00,0,eri,1)
        call ifill(2,-1,idatom2,1)
        return
      endif
c
c.. translate jsh center coordinates based on Rj
      call intp_txyz(q_cent,q_geom,Rj,xyz_new_j)
      if(cando_nw(brain,ish,0).and.cando_nw(ketin,jsh,0)) then
c
        ucont   = (sf_ibs_cn2ucn(ish,bra))
        Lp      = infbs_cont(CONT_TYPE ,ucont,bra)
        p_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        p_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        p_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        p_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = (sf_ibs_cn2ucn(jsh,ket))
        Lq      = infbs_cont(CONT_TYPE ,ucont,ket)
        q_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        q_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        q_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        q_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c
        call hf2d(
     &      coords(1,p_cent,p_geom), dbl_mb(mb_exndcf(p_iexp,bra)),
     &      dbl_mb(mb_exndcf(p_icfp,bra)), p_prim, p_gen, Lp, p_cent,
     &      coords(1,p_cent,p_geom), DCexp,
     &      DCcoeff           ,      1,     1, 0, p_cent,
     &      xyz_new_j, dbl_mb(mb_exndcf(q_iexp,ket)),
     &      dbl_mb(mb_exndcf(q_icfp,ket)), q_prim, q_gen, Lq, q_cent,
     &      xyz_new_j, DCexp,
     &      DCcoeff           ,      1,     1, 0, q_cent,
c...................... canAB    canCD    canPQ
     &      eri, mynint, .false., .false., .false.,
c........... dryrun
     &      .false., scr, lscr)
c
        call intd_sum2c(eri,mynint,idatom,p_cent,q_cent)
        idatom2(1) = idatom(1)
        idatom2(2) = idatom(3)
        if (bas_spherical(bra).or.bas_spherical(ket)) then
          if (Lp.eq.-1) p_gen = 1
          if (Lq.eq.-1) q_gen = 1
          call spcart_2cBtran(eri,scr,lscr,
     &        int_nbf_x(Lp),int_nbf_s(Lp),Lp,p_gen,bas_spherical(bra),
     &        int_nbf_x(Lq),int_nbf_s(Lq),Lq,q_gen,bas_spherical(ket),
     &        (3*2),.false.)
        endif
      else
        write(luout,*)'intpd_2e2c: could not do nw integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        call errquit('intpd_2e2c: fatal error ',0, INT_ERR)
      endif
      end
C> @}
