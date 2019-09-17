C> \ingroup nwint
C> @{
C>
C> \brief Compute 2-electron 3-center Gaussian periodic electron
C> repulsion integrals
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
      subroutine intp_2e3c(brain, ish, ketin, jsh, ksh, Ri, Rj, 
     &       lscr, scr, leri, eri)
c $Id$
      implicit none
c
c basic api routine to generate a block of 3 center two 
c  electron integrals
c
c eri = <bra_g(ish)|ket_g(jsh).ket_g(ksh)>
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
      double precision Ri(3)  !< [Input] translation vector for ish center 
      double precision Rj(3)  !< [Input] translation vector for jsh center 
c:: local
      logical shells_ok
      integer bra, ket 
      integer p_geom, cd_geom, ucont
      integer Lp, p_prim, p_gen, p_iexp, p_icfp, p_cent
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent
      integer nintx
      double precision xyz_new_i(3)  ! new coordinates for ish function
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
      if (.not.int_chk_init('intp_2e3c'))
     &       call errquit('intp_2e3c: int_init was not called' ,0,
     &          INT_ERR)
c
c     check input shell ranges
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      if (.not.shells_ok)
     &       call errquit('intp_2e3c: invalid contraction/shell',0,
     &          BASIS_ERR)
c
c  check if spherical/gencon/sp shell
c
      call int_nogencont_check(brain,'intp_2e3c:bra')
      call int_nogencont_check(ketin,'intp_2e3c:ket')
      call int_nospshell_check(brain,'intp_2e3c:bra')
      call int_nospshell_check(ketin,'intp_2e3c:ket')
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      p_geom  = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
c
      if (p_geom.ne.cd_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'intp_2e3c: WARNING: possible geometry inconsistency'
        write(luout,*)'bra geometry handle:',p_geom
        write(luout,*)'ket geometry handle:',cd_geom
        WarnP = 1
      endif
c
      p_cent  = (sf_ibs_cn2ce(ish,bra))
      c_cent  = (sf_ibs_cn2ce(jsh,ket))
      d_cent  = (sf_ibs_cn2ce(ksh,ket))
c.. translate ish center coordinates based on Ri
      call intp_txyz(p_cent,p_geom,Ri,xyz_new_i)
c.. translate jsh center coordinates based on Rj
      call intp_txyz(c_cent,cd_geom,Rj,xyz_new_j)
c
      if(cando_nw(brain,ish,0).and.cando_nw(ketin,jsh,ksh)) then
c
        ucont   = (sf_ibs_cn2ucn(ish,bra))
        Lp      = infbs_cont(CONT_TYPE ,ucont,bra)
        p_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        p_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        p_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        p_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = (sf_ibs_cn2ucn(jsh,ket))
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c
        ucont   = (sf_ibs_cn2ucn(ksh,ket))
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c      
        call hf2(
     &         xyz_new_i, dbl_mb(mb_exndcf(p_iexp,bra)),
     &         dbl_mb(mb_exndcf(p_icfp,bra)), p_prim, p_gen, Lp,
     &         xyz_new_i, DCexp,
     &         DCcoeff           ,      1,     1, 0,
     &         xyz_new_j, dbl_mb(mb_exndcf(c_iexp,ket)),
     &         dbl_mb(mb_exndcf(c_icfp,ket)), c_prim, c_gen, Lc,
     &         coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &         dbl_mb(mb_exndcf(d_icfp,ket)), d_prim,d_gen,Ld,
c......................... canAB    canCD    canPQ
     &         eri, leri, .false., .false., .false.,
c............ dryrun
     &         .false., scr, lscr)
c
        if (bas_spherical(ket).or.bas_spherical(bra)) then
          nintx = int_nbf_x(Lp)*int_nbf_x(Lc)*int_nbf_x(Ld)
          if (nintx.lt.lscr)
     &        call errquit
     &        ('intp_2e3c: not enough scratch for spherical transform',
     &        911, INT_ERR)
          call spcart_3ctran(eri,scr,lscr,
     &        int_nbf_x(Lp),int_nbf_s(Lp),Lp,p_gen,bas_spherical(bra),
     &        int_nbf_x(Lc),int_nbf_s(Lc),Lc,c_gen,bas_spherical(ket),
     &        int_nbf_x(Ld),int_nbf_s(Ld),Ld,d_gen,bas_spherical(ket),
     &        .false.)          
        endif
      else
        write(luout,*)'intp_2e3c: could not do nw integrals'
        write(luout,*)' brain :',brain
        write(luout,*)' ketin :',ketin
        write(luout,*)' ish   :',ish
        write(luout,*)' jsh   :',jsh
        write(luout,*)' ksh   :',ksh
        call errquit('intp_2e3c: fatal error ',0, INT_ERR)
      endif
      end
C> @}
