      subroutine intd_2e2c(brain, ish, ketin, jsh, 
     &       lscr, scr, leri, eri, idatom2)
c $Id: intd_2e2c.F,v 1.10 2004-12-03 22:34:31 bert Exp $
      implicit none
c
c basic api routine to generate 4 center two electron integral derivatives 
c 
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "apiP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
c
c::external subroutines used
c dcopy
c daxpy
c errquit
c::functions
      integer int_nint_cart
      external int_nint_cart
c::passed
      integer brain ! [input] basis set handle for bra basis
      integer ish   ! [input] lexical contraction index 
      integer jsh   ! [input] lexical contraction index 
      integer ketin ! [input] basis set handle for ket basis
      integer lscr  ! [input] length of scratch array
      integer leri  ! [input] length of eri array
      double precision scr(lscr) ! [scratch] scratch array for integral code.
      double precision eri(leri) ! [output]  array for two electron integral derivatives.
      integer idatom2(2)         ! [output]  array identifying centers for derivatives externalized
c NOTE: length of idatom is always 4 because there can be at most 4 centers involved
      integer idatom(4)          ! local array identifying centers for derivatives
c
c                       ! e.g., the first  nint*3 derivatives go to center idatom2(1)
c                       !       the second nint*3 derivatives go to center idatom2(2)
c
c Order is...   nint*3*2 (3=> xyz, 2=atoms)
c
c  /                   |   
c | nint,  d <ij|kl>   |
c |      --------------|
c  \     d[idatom(1),x]|
c                          |
c       nint,  d <ij|kl>   |
c            --------------|
c            d[idatom(1),y]|
c                              |
c           nint,  d <ij|kl>   |
c                --------------|
c                d[idatom(1),z]|
c                                  |
c               nint,  d <ij|kl>   |
c                    --------------|
c                    d[idatom(2),x]|
c                                      |
c                   nint,  d <ij|kl>   |
c                        --------------|
c                        d[idatom(2),y]|
c                                           \
c                       nint,  d <ij|kl>     |
c                            --------------  |
c                            d[idatom(2),z] / 
c                     
c::local
      logical any_spherical
      integer nint, ucont
      integer bra, ket, ab_geom, cd_geom
      integer inp, igen, iexp, icf, itype, iatom
      integer jnp, jgen, jexp, jcf, jtype, jatom
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      nint = int_nint_cart(brain,ish,ketin,jsh,brain,0,ketin,0)
      if (nint*3*4.gt.leri) then
        write(luout,*) 'nint*3*4 = ',nint*3*4
        write(luout,*) 'leri     = ',leri
        call errquit('intd_2e2c: nint>leri error',911, INT_ERR)
      endif
c
c  check if gencont/sp shell
c
      call int_nogencont_check(brain,'intd_2e2c:bra')
      call int_nogencont_check(ketin,'intd_2e2c:ket')
      call int_nospshell_check(brain,'intd_2e2c:bra')
      call int_nospshell_check(ketin,'intd_2e2c:ket')
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
      if (ab_geom.ne.cd_geom) then
        write(luout,*)
     &      'intd_2e2c.F: two different geometries for',
     &      ' derivatives?'
        call errquit('intd_2e2c: geom error ',911, GEOM_ERR)
      endif
      
c
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)
c
      ucont = (sf_ibs_cn2ucn(ish,bra))
      inp   = infbs_cont(CONT_NPRIM,ucont,bra)
      igen  = infbs_cont(CONT_NGEN,ucont,bra)
      iexp  = infbs_cont(CONT_IEXP,ucont,bra)
      icf   = infbs_cont(CONT_ICFP,ucont,bra)
      itype = infbs_cont(CONT_TYPE,ucont,bra) 
      iatom = (sf_ibs_cn2ce(ish,bra))
c
      ucont = (sf_ibs_cn2ucn(jsh,ket))
      jnp   = infbs_cont(CONT_NPRIM,ucont,ket)
      jgen  = infbs_cont(CONT_NGEN,ucont,ket)
      jexp  = infbs_cont(CONT_IEXP,ucont,ket)
      jcf   = infbs_cont(CONT_ICFP,ucont,ket)
      jtype = infbs_cont(CONT_TYPE,ucont,ket) 
      jatom = (sf_ibs_cn2ce(jsh,ket))
c
c... new logic
      if (iatom.eq.jatom) then
        call dcopy((nint*3*4),0.0d00,0,eri,1)
        call ifill(2,-1,idatom2,1)
        return
      endif
c
c
      call hf2d(
     &       coords(1,iatom,ab_geom),dbl_mb(mb_exndcf(iexp,bra)),
     &       dbl_mb(mb_exndcf(icf,bra)),inp,igen,itype,iatom,
c
     &       coords(1,iatom,ab_geom),DCexp,
     &       DCcoeff,1,1,0,iatom,
c
     &       coords(1,jatom,cd_geom),dbl_mb(mb_exndcf(jexp,ket)),
     &       dbl_mb(mb_exndcf(jcf,ket)),jnp,jgen,jtype,jatom,
c
     &       coords(1,jatom,cd_geom),DCexp,
     &       DCcoeff,1,1,0,jatom,
c
     &       eri,nint,.false.,.false.,.false.,.false.,
     &       scr,lscr)
c
      call intd_sum2c(eri,nint,idatom,iatom,jatom)
      idatom2(1) = idatom(1)
      idatom2(2) = idatom(3)
c
      if (any_spherical) then
        if (itype.eq.-1) igen = 1
        if (jtype.eq.-1) jgen = 1
        call spcart_2cBtran(eri,scr,lscr,
     &      int_nbf_x(itype),int_nbf_s(itype),itype,igen,
     &      bas_spherical(bra),
     &      int_nbf_x(jtype),int_nbf_s(jtype),jtype,jgen,
     &      bas_spherical(ket),
     &      6,.false.)
      endif
c
      end
      subroutine intd_sum2c(eri,nint,idatom,iatom,jatom)
      implicit none
c
      integer nint
      double precision eri(nint,3,4)
      integer idatom(4)
      integer iatom,jatom
c
      integer iduse
c
      call intd_logic_atom(idatom,iatom,iatom,jatom,jatom)
c
      do 00100 iduse = 2,4
        if (idatom(iduse).gt.0) then
          continue
        else
          call daxpy(nint*3,1.0d00,eri(1,1,iduse),1,
     &           eri(1,1,abs(idatom(iduse))),1)
        endif
00100 continue
c
      call dcopy(nint*3,eri(1,1,3),1,eri(1,1,2),1)
c
      end
