      subroutine intd_1eov(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lOva,Ova,idatom)
c $Id: intd_1eov.F,v 1.3 1996-01-22 21:46:41 d3e129 Exp $
      implicit none
c
c returns 1e overlap derivative integrals for shells/contractions ish,jsh
c
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c::external subroutines used
c... errquit
c::functions
      integer int_nint
      external int_nint
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] "i" contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] "j" contraction index
      integer lscr      ! [input] length of scratch space "scr"
      integer lOva      ! [input] number of overlap integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      integer idatom(*) ! [output] array identifying centers for derivatives
c                       ! e.g., the first nint*3  derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
c
c Order is...   nint*3*2 (3=> xyz, 2=atoms)
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
c                                          \
c                       nint,d <ij>         |
c                            -------------- |
c                            d[idatom(2),z]/
c
      double precision scr(lscr) ! [input] scratch array
      double precision Ova(*)    ! [output] derivative integrals
c::local
      integer nint, ucont
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      nint = int_nint(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*2.gt.lOva) then
        write(6,*) 'nint*3*2 = ',nint*3*2
        write(6,*) 'lOva     = ',lOva
        call errquit('intd_1eov: nint>lOva error',911)
      endif
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      iatom = (sf_ibs_cn2ce(ish,ibas))
      jatom = (sf_ibs_cn2ce(jsh,jbas))
      if (iatom.eq.jatom) then
        idatom(1) = 0
        idatom(2) = 0
        return
      endif
c
      idatom(1) = iatom
      idatom(2) = jatom
c
      ucont = (sf_ibs_cn2ucn(ish,ibas))
      inp   = infbs_cont(CONT_NPRIM,ucont,ibas)
      igen  = infbs_cont(CONT_NGEN,ucont,ibas)
      iexp  = infbs_cont(CONT_IEXP,ucont,ibas)
      icf   = infbs_cont(CONT_ICFP,ucont,ibas)
      itype = infbs_cont(CONT_TYPE,ucont,ibas) 
      igeom = ibs_geom(ibas)
c
      ucont = (sf_ibs_cn2ucn(jsh,jbas))
      jnp   = infbs_cont(CONT_NPRIM,ucont,jbas)
      jgen  = infbs_cont(CONT_NGEN,ucont,jbas)
      jexp  = infbs_cont(CONT_IEXP,ucont,jbas)
      jcf   = infbs_cont(CONT_ICFP,ucont,jbas)
      jtype = infbs_cont(CONT_TYPE,ucont,jbas) 
      jgeom = ibs_geom(jbas)
c
      if (igeom.ne.jgeom) then
        write(6,*)'intd_1eov.F: two different geometries for',
     &         ' derivatives?'
        call errquit('intd_1eov: geom error ',911)
      endif
c
      call hf1d(
     &       coords(1,iatom,igeom),
     &       dbl_mb(mb_exndcf(iexp,ibas)),
     &       dbl_mb(mb_exndcf(icf,ibas)),
     &       inp,igen,itype,iatom,
c
     &       coords(1,jatom,jgeom),
     &       dbl_mb(mb_exndcf(jexp,jbas)),
     &       dbl_mb(mb_exndcf(jcf,jbas)),
     &       jnp,jgen,jtype,jatom,
c
     &       coords(1,1,igeom),charge(1,igeom),ncenter(igeom),
     &       Ova,scr,scr,nint,
c............overlap, k-e,     pot-e,    canab,  dryrun
     &       .true.,  .false., .false., .false., .false.,
     &       scr,lscr)
c
      end
