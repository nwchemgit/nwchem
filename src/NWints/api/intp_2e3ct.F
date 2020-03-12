C> \ingroup nwint
C> @{
C>
C> \brief Compute transposed 2-electron 3-center Gaussian periodic
C> electron repulsion integrals
C>
C> The code in question always uses density fitting approaches for
C> the 3-center 2-electron integrals are needed. See [1] for details.
C>
C> [1] JE Jaffe, AC Hess,
C>     <i>"Gaussian basis density functional theory for systems 
C>     periodic in two or three dimensions: Energy and forces"</i>,
C>    J.Chem.Phys. <b>105</b>, 10983-10998 (1996), DOI:
C>    <a href="https://doi.org/10.1063/1.472866">
C>    10.1063/1.472866</a>
C>
      subroutine intp_2e3ct(brain, ish, jsh, ketin, ksh, Rj, Rk, 
     &       lscr, scr, leri, eri)
*
* $Id$
*
      implicit none
c
c basic api routine to generate a block of 3 center two 
c  electron integrals
c
c eri = <bra_g(ish)bra_g(jsh)|ket_g(ksh)>
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
c:: passed
      integer brain  !< [Input] basis set handle for bra function
      integer ketin  !< [Input] basis set handle for ket function
      integer ish    !< [Input] lexical contraction/shell index
      integer jsh    !< [Input] lexical contraction/shell index
      integer ksh    !< [Input] lexical contraction/shell index
      integer lscr   !< [Input] length of scratch array
      integer leri   !< [Input] length of ERI array
      double precision scr(lscr) !< [Scratch] scratch array 
      double precision eri(leri) !< [Output] ERI array
c translation vectors are in fractional coordinates !
      double precision Rk(3)  !< [Input] translation vector for ksh center 
      double precision Rj(3)  !< [Input] translation vector for jsh center 
c:: local
      logical shells_ok
      integer bra, ket 
      integer q_geom, ab_geom, ucont
      integer Lq, q_prim, q_gen, q_iexp, q_icfp, q_cent
      integer La, a_prim, a_gen, a_iexp, a_icfp, a_cent
      integer Lb, b_prim, b_gen, b_iexp, b_icfp, b_cent
      integer nintx
      double precision xyz_new_k(3)  ! new coordinates for ksh function
      double precision xyz_new_j(3)  ! new coordinates for jsh function
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
      if (.not.int_chk_init('intp_2e3ct'))
     &       call errquit('intp_2e3ct: int_init was not called' ,0,
     &            INT_ERR)
c
c     check input shell ranges
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(brain,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      if (.not.shells_ok)
     &       call errquit('intp_2e3ct: invalid contraction/shell',0,
     &           BASIS_ERR)
c
c  check if spherical/gencon/sp shell
c
      call int_nogencont_check(brain,'intp_2e3ct:bra')
      call int_nogencont_check(ketin,'intp_2e3ct:ket')
      call int_nospshell_check(brain,'intp_2e3ct:bra')
      call int_nospshell_check(ketin,'intp_2e3ct:ket')
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      ab_geom = ibs_geom(bra)
      q_geom  = ibs_geom(ket)
c
      if (q_geom.ne.ab_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'intp_2e3ct: WARNING: possible geometry inconsistency'
        write(luout,*)'bra geometry handle:',ab_geom
        write(luout,*)'ket geometry handle:',q_geom
        WarnP = 1
      endif
c
      a_cent  = (sf_ibs_cn2ce(ish,bra))
      b_cent  = (sf_ibs_cn2ce(jsh,bra))
      q_cent  = (sf_ibs_cn2ce(ksh,ket))
c.. translate ksh center coordinates based on Rk
      call intp_txyz(q_cent,q_geom,Rk,xyz_new_k)
c.. translate jsh center coordinates based on Rj
      call intp_txyz(b_cent,ab_geom,Rj,xyz_new_j)
c
      if(cando_nw(brain,ish,jsh).and.cando_nw(ketin,ksh,0)) then
c
        ucont   = (sf_ibs_cn2ucn(ish,bra))
        La      = infbs_cont(CONT_TYPE ,ucont,bra)
        a_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        a_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        a_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = (sf_ibs_cn2ucn(jsh,bra))
        Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
        b_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        b_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        b_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = (sf_ibs_cn2ucn(ksh,ket))
        Lq      = infbs_cont(CONT_TYPE ,ucont,ket)
        q_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        q_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        q_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        q_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c      
        call hf2(
     &         coords(1,a_cent,ab_geom), dbl_mb(mb_exndcf(a_iexp,bra)),
     &         dbl_mb(mb_exndcf(a_icfp,bra)), a_prim,a_gen,La,
     &         xyz_new_j, dbl_mb(mb_exndcf(b_iexp,bra)),
     &         dbl_mb(mb_exndcf(b_icfp,bra)), b_prim, b_gen, Lb,
     &         xyz_new_k, dbl_mb(mb_exndcf(q_iexp,ket)),
     &         dbl_mb(mb_exndcf(q_icfp,ket)), q_prim, q_gen, Lq,
     &         xyz_new_k, DCexp,
     &         DCcoeff           ,      1,     1, 0,
c......................... canAB    canCD    canPQ
     &         eri, leri, .false., .false., .false.,
c............ dryrun
     &         .false., scr, lscr)
c
        if (bas_spherical(ket).or.bas_spherical(bra)) then
          nintx = int_nbf_x(Lq)*int_nbf_x(La)*int_nbf_x(Lb)
          if (nintx.lt.lscr)
     &        call errquit
     &        ('intp_2e3ct: not enough scratch for spherical transform',
     &        911, MEM_ERR)
          call spcart_3ctran(eri,scr,lscr,
     &        int_nbf_x(La),int_nbf_s(La),La,a_gen,bas_spherical(bra),
     &        int_nbf_x(Lb),int_nbf_s(Lb),Lb,b_gen,bas_spherical(bra),
     &        int_nbf_x(Lq),int_nbf_s(Lq),Lq,q_gen,bas_spherical(ket),
     &        .false.)          
        endif
      else
        write(luout,*)'intp_2e3ct: could not do nw integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        write(luout,*)' ksh   :',ksh
        call errquit('intp_2e3ct: fatal error ',0, INT_ERR)
      endif
      end
C> @}
