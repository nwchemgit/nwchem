      subroutine intd_1eke(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lTa,Ta,idatom)
c $Id: intd_1eke.F,v 1.1 1999-04-10 00:14:50 windus Exp $
c
c This routine has NOT been extensively tested. - TLW
c
      implicit none
#include "stdio.fh"
c::functions
      integer int_nint_cart
      external int_nint_cart
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] "i" contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] "j" contraction index
      integer lscr      ! [input] length of scratch space "scr"
      integer lTa       ! [input] number of kinetic integral derivatives in shells ish and jsh
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
      double precision scr(lscr)
      double precision Ta(lTa)
*
      integer nint
*
      nint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*2.gt.lTa) then
        write(luout,*) 'nint*3*2 = ',nint*3*2
        write(luout,*) 'lTa     = ',lTa
        call errquit('intd_1eke: nint>lTa error',911)
      endif
      call intd_1ekeP(i_basis,ish,j_basis,jsh,lscr,scr,nint,
     &    Ta,idatom)
      end
      
      subroutine intd_1ekeP(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,Ta,idatom)
      implicit none
c
c returns 1e kinetic derivative integrals for shells/contractions ish,jsh
c
#include "stdio.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw
      external cando_hnd_1e
      external cando_nw
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] "i" contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] "j" contraction index
      integer lscr      ! [input] length of scratch space "scr"
      integer nint      ! [input] number of kinetic integral derivatives in shells ish and jsh
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
      double precision scr(lscr)     ! [input] scratch array
      double precision Ta(nint,3,2) ! [output] derivative integrals
c::local
      integer ucont
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer nint_x, nint_s
      integer zatom, zyx
c
c  Temporary variable that needs to be taken out after testing!
c
      integer itemp
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c
c  check if gencont, or spshell
c
      call int_nogencont_check(i_basis,'intd_1eke:i_basis')
      call int_nogencont_check(j_basis,'intd_1eke:j_basis')
      call int_nospshell_check(i_basis,'intd_1eke:i_basis')
      call int_nospshell_check(j_basis,'intd_1eke:j_basis')
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
        write(luout,*)'intd_1eke: two different geometries for',
     &         ' derivatives?'
        call errquit('intd_1eke: geom error ',911)
      endif
c
      if (cando_nw(i_basis,ish,0).and.
     &        cando_nw(j_basis,jsh,0)) then
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
     &       scr,Ta,scr,nint,
c............overlap, k-e,     pot-e,    canab,  dryrun
     &       .false.,  .true., .false., .false., .false.,
     &       scr,lscr)
      write(luout,*) 'K 1e- derivative integrals from nwchem'
      do itemp = 1,nint
        write(luout,*) itemp,' 1,1,',Ta(itemp,1,1)
        write(luout,*) itemp,' 2,1,',Ta(itemp,2,1)
        write(luout,*) itemp,' 3,1,',Ta(itemp,3,1)
        write(luout,*) itemp,' 1,2,',Ta(itemp,1,2)
        write(luout,*) itemp,' 2,2,',Ta(itemp,2,2)
        write(luout,*) itemp,' 3,2,',Ta(itemp,3,2)
      enddo
      elseif (cando_hnd_1e(i_basis,ish,0).and.
     &  cando_hnd_1e(j_basis,jsh,0)) then
        call hnd_stvintd(
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
     &       scr,Ta,scr,nint,
c............overlap, k-e,     pot-e,
     &       .false.,  .true., .false., 
     &       scr,lscr)
      write(luout,*) 'K 1e- derivative integrals from hondo'
      do itemp = 1,nint
        write(luout,*) itemp,' 1,1,',Ta(itemp,1,1)
        write(luout,*) itemp,' 2,1,',Ta(itemp,2,1)
        write(luout,*) itemp,' 3,1,',Ta(itemp,3,1)
        write(luout,*) itemp,' 1,2,',Ta(itemp,1,2)
        write(luout,*) itemp,' 2,2,',Ta(itemp,2,2)
        write(luout,*) itemp,' 3,2,',Ta(itemp,3,2)
      enddo
      endif
c
*     Ta now has the cartesian integral block  (jlo:jhi,ilo:ihi,3,2)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4.
c
      if (itype.eq.-1) igen = 1
      if (jtype.eq.-1) jgen = 1
c
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*... transform both i and j integrals
        i_nbf_x = int_nbf_x(Itype)
        i_nbf_s = int_nbf_s(Itype)
        j_nbf_x = int_nbf_x(Jtype)
        j_nbf_s = int_nbf_s(Jtype)
c        
        do zatom = 1,2
          do zyx = 1,3
            call spcart_tran1e(Ta(1,zyx,zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,Itype,igen,
     &          .false.)
          enddo
        enddo
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Itype)
        i_nbf_s = int_nbf_s(Itype)
        j_nbf_x = int_nbf_x(Jtype)
        j_nbf_s = j_nbf_x
        do zatom = 1,2
          do zyx = 1,3
            call spcart_tran1e(Ta(1,Zyx,Zatom),scr,
     &          j_nbf_x,i_nbf_x,0,jgen,
     &          j_nbf_s,i_nbf_s,Itype,igen,
     &          .false.)
          enddo
        enddo
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Itype)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Jtype)
        j_nbf_s = int_nbf_s(Jtype)
        do zatom = 1,2
          do zyx = 1,3
            call spcart_tran1e(Ta(1,Zyx,Zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,0,igen,
     &          .false.)
          enddo
        enddo
      else
        call errquit(
     &        'intd_1eke: should never reach transform blocked else',
     &        911)
      endif
c
c now shuffle transformed buffers to contiguous space
c
      nint_x = i_nbf_x*j_nbf_x
      nint_s = i_nbf_s*j_nbf_s
      if (nint_s.gt.nint_x) then
        call errquit
     &      ('intd_1eke: nint_s >.nint_x diff=',(nint_s-nint_x))
      elseif (nint_s.eq.nint_x) then
        return
      else
        call int_c2s_mv
     &      (Ta,nint_x,nint_s,(3*2),scr,lscr,'intd_1eke')
      endif
c
      end
