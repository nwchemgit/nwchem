      subroutine intd_2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       lscr, scr, leri, eri, idatom)
c $Id: intd_2e4c.F,v 1.8 1996-07-11 22:21:18 d3e129 Exp $
      implicit none
c
c basic api routine to generate 4 center two electron integral derivatives 
c 
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
c::external subroutines used
c errquit
c::functions
      integer int_nint
      external int_nint
c::passed
      integer brain ! [input] basis set handle for bra basis
      integer ish   ! [input] lexical contraction index 
      integer jsh   ! [input] lexical contraction index 
      integer ketin ! [input] basis set handle for ket basis
      integer ksh   ! [input] lexical contraction index 
      integer lsh   ! [input] lexical contraction index 
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
c                       !       the fourth nint*3 derivatives go to center idatom(4)
c
c Order is...   nint*3*4 (3=> xyz, 4=atoms)
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
c                                                      |
c                                   nint,  d <ij|kl>   |
c                                        --------------|
c                                        d[idatom(3),z]|
c                                                          |
c                                       nint,  d <ij|kl>   |
c                                            --------------|
c                                            d[idatom(4),x]|
c                                                              |
c                                           nint,  d <ij|kl>   |
c                                                --------------|
c                                                d[idatom(4),y]|
c                                                                   \
c                                               nint,  d <ij|kl>     |
c                                                    --------------  |
c                                                    d[idatom(4),z] / 
c                     
c::local
      integer nint, ucont
      integer bra, ket, ab_geom, cd_geom
      integer inp, igen, iexp, icf, itype, iatom
      integer jnp, jgen, jexp, jcf, jtype, jatom
      integer knp, kgen, kexp, kcf, ktype, katom
      integer lnp, lgen, lexp, lcf, ltype, latom
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      nint = int_nint(brain,ish,brain,jsh,ketin,ksh,ketin,lsh)
      if (nint*3*4.gt.leri) then
        write(6,*) 'nint*3*4 = ',nint*3*4
        write(6,*) 'leri     = ',leri
        call errquit('intd_2e4c: nint>leri error',911)
      endif
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
      if (ab_geom.ne.cd_geom) then
        write(6,*)'intd_2e4c.F: two different geometries for',
     &         ' derivatives?'
        call errquit('intd_2e4c: geom error ',911)
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
      ucont = (sf_ibs_cn2ucn(jsh,bra))
      jnp   = infbs_cont(CONT_NPRIM,ucont,bra)
      jgen  = infbs_cont(CONT_NGEN,ucont,bra)
      jexp  = infbs_cont(CONT_IEXP,ucont,bra)
      jcf   = infbs_cont(CONT_ICFP,ucont,bra)
      jtype = infbs_cont(CONT_TYPE,ucont,bra) 
      jatom = (sf_ibs_cn2ce(jsh,bra))
c
      ucont = (sf_ibs_cn2ucn(ksh,ket))
      knp   = infbs_cont(CONT_NPRIM,ucont,ket)
      kgen  = infbs_cont(CONT_NGEN,ucont,ket)
      kexp  = infbs_cont(CONT_IEXP,ucont,ket)
      kcf   = infbs_cont(CONT_ICFP,ucont,ket)
      ktype = infbs_cont(CONT_TYPE,ucont,ket) 
      katom = (sf_ibs_cn2ce(ksh,ket))
c
      ucont = (sf_ibs_cn2ucn(lsh,ket))
      lnp   = infbs_cont(CONT_NPRIM,ucont,ket)
      lgen  = infbs_cont(CONT_NGEN,ucont,ket)
      lexp  = infbs_cont(CONT_IEXP,ucont,ket)
      lcf   = infbs_cont(CONT_ICFP,ucont,ket)
      ltype = infbs_cont(CONT_TYPE,ucont,ket) 
      latom = (sf_ibs_cn2ce(lsh,ket))
c
c... new logic
      if (iatom.eq.jatom.and.jatom.eq.katom.and.katom.eq.latom) then
        call dfill((nint*3*4),0.0d00,eri,1)
        call ifill(4,-1,idatom,1)
        return
      endif
c
      call hf2d(
     &       coords(1,iatom,ab_geom),dbl_mb(mb_exndcf(iexp,bra)),
     &       dbl_mb(mb_exndcf(icf,bra)),inp,igen,itype,iatom,
c
     &       coords(1,jatom,ab_geom),dbl_mb(mb_exndcf(jexp,bra)),
     &       dbl_mb(mb_exndcf(jcf,bra)),jnp,jgen,jtype,jatom,
c
     &       coords(1,katom,cd_geom),dbl_mb(mb_exndcf(kexp,ket)),
     &       dbl_mb(mb_exndcf(kcf,ket)),knp,kgen,ktype,katom,
c
     &       coords(1,latom,cd_geom),dbl_mb(mb_exndcf(lexp,ket)),
     &       dbl_mb(mb_exndcf(lcf,ket)),lnp,lgen,ltype,latom,
c
     &       eri,nint,.false.,.false.,.false.,.false.,
     &       scr,lscr)
c
      call intd_sum(eri,nint,idatom,iatom,jatom,katom,latom)
c
      end
      subroutine intd_sum(eri,nint,idatom,iatom,jatom,katom,latom)
      implicit none
c
      integer nint
      double precision eri(nint,3,4)
      integer idatom(4)
      integer iatom,jatom,katom,latom
c
      integer iduse
c
      call intd_logic_atom(idatom,iatom,jatom,katom,latom)
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
      end
      subroutine intd_logic_atom(idat,iat,jat,kat,lat)
      implicit none
      integer iat,jat,kat,lat
      integer idat(4)
c
      idat(1) = iat
      idat(2) = jat
      idat(3) = kat
      idat(4) = lat
      if (iat.eq.jat) idat(2) = -1
      if (iat.eq.kat) idat(3) = -1
      if (iat.eq.lat) idat(4) = -1
      if (jat.eq.kat) then
        if(idat(2).gt.0) then
          idat(3) = -2
        else
          idat(3) = idat(2)
        endif
      endif
      if (jat.eq.lat) then
        if(idat(2).gt.0) then
          idat(4) = -2
        else
          idat(4) = idat(2)
        endif
      endif
      if (kat.eq.lat) then
        if(idat(3).gt.0) then
          idat(4) = -3
        else
          idat(4) = idat(3)
        endif
      endif
      end
