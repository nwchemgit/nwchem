      subroutine intd_2e3c(brain, ish, ketin, jsh, ksh, 
     &       lscr, scr, leri, eri, idatom)
c $Id: intd_2e3c.F,v 1.2 1996-01-22 21:46:45 d3e129 Exp $
      implicit none
c
c basic api routine to generate 4 center two electron integral derivatives 
c 
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "apiP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
*:need to fix:#if defined(INT_TIMING)
*:need to fix:#include "int_tim.fh"
*:need to fix:#include "pstat.fh"
*:need to fix:#endif
c
c::external subroutines used
c errquit
c::functions
      integer int_nint
      external int_nint
c::passed
      integer brain ! [input] basis set handle for bra basis
      integer ish   ! [input] lexical contraction index 
      integer ketin ! [input] basis set handle for ket basis
      integer jsh   ! [input] lexical contraction index 
      integer ksh   ! [input] lexical contraction index 
      integer lscr  ! [input] length of scratch array
      integer leri  ! [input] length of eri array
      double precision scr(lscr) ! [scratch] scratch array for integral code.
      double precision eri(leri) ! [output]  array for two electron integral derivatives.
c NOTE: length of idatom is always 4 because there can be at most 4 centers involved
      integer idatom(4)          ! [output]  array identifying centers for derivatives
c
c                       ! e.g., the first  nint*3 derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
c                       !       the third  nint*3 derivatives go to center idatom(3)
c
c Order is...   nint*3*3 (3=> xyz, 3=atoms)
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
c                                           |
c                       nint,  d <ij|kl>    |
c                            -------------- |
c                            d[idatom(2),z] |
c                                              |   
c                           nint,  d <ij|kl>   |
c                                --------------|
c                                d[idatom(3),x]|
c                                                  |
c                               nint,  d <ij|kl>   |
c                                    --------------|
c                                    d[idatom(3),y]|
c                                                       \
c                                   nint,  d <ij|kl>     |
c                                        --------------  |
c                                        d[idatom(3),z] / 
c                     
c::local
      integer nint, ucont
      integer bra, ket, ab_geom, cd_geom
      integer inp, igen, iexp, icf, itype, iatom
      integer jnp, jgen, jexp, jcf, jtype, jatom
      integer knp, kgen, kexp, kcf, ktype, katom
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      nint = int_nint(brain,ish,ketin,jsh,ketin,ksh,ketin,0)
      if (nint*3*4.gt.leri) then
        write(6,*) 'nint*3*4 = ',nint*3*4
        write(6,*) 'leri     = ',leri
        call errquit('intd_2e3c: nint>leri error',911)
      endif
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
      if (ab_geom.ne.cd_geom) then
        write(6,*)'intd_2e3c.F: two different geometries for',
     &         ' derivatives?'
        call errquit('intd_2e3c: geom error ',911)
      endif
      
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
      ucont = (sf_ibs_cn2ucn(ksh,ket))
      knp   = infbs_cont(CONT_NPRIM,ucont,ket)
      kgen  = infbs_cont(CONT_NGEN,ucont,ket)
      kexp  = infbs_cont(CONT_IEXP,ucont,ket)
      kcf   = infbs_cont(CONT_ICFP,ucont,ket)
      ktype = infbs_cont(CONT_TYPE,ucont,ket) 
      katom = (sf_ibs_cn2ce(ksh,ket))
c
*:need to fix:#if defined(INT_TIMING)
*:need to fix:      if(o_int_time.and.(int_time_level.ge.2))
*:need to fix:     &   call pstat_on(intt_arr(intt_hf2d))
*:need to fix:#endif
c... new logic
      if (iatom.eq.jatom.and.jatom.eq.katom)then
        call dfill((nint*3*4),0.0d00,eri,1)
        call ifill(4,-1,idatom,1)
        return
      endif
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
     &       coords(1,katom,cd_geom),dbl_mb(mb_exndcf(kexp,ket)),
     &       dbl_mb(mb_exndcf(kcf,ket)),knp,kgen,ktype,katom,
c
     &       eri,nint,.false.,.false.,.false.,.false.,
     &       scr,lscr)
c
      call intd_sum(eri,nint,idatom,iatom,iatom,jatom,katom)
c
*:need to fix:#if defined(INT_TIMING)
*:need to fix:      if(o_int_time.and.(int_time_level.ge.2))
*:need to fix:     &   call pstat_off(intt_arr(intt_hf2d))
*:need to fix:#endif
      end
