C> Old intd_1eov is now intd_1eovQ (below).  This intd_1eov is a wrapper
C> for backward compatibility.  intd_1eovQ has a new input argument,
C> hlfdrv, required for computing the half-derivative overlap matrix used
C> by COLUMBUS to compute non-adiabatic derivative couplings.  Routines
C> which need to calculate the HDOM should call intd_1eovQ directly,
C> specifying hlfdrv=.true.
C>
C> Lachlan Belcher 14 Feb 2023

      subroutine intd_1eov(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lOva,Ova,idatom)
c $Id$
      implicit none
      integer i_basis   !< [Input] ish basis set handle
      integer ish       !< [Input] "i" contraction index
      integer j_basis   !< [Input] jsh basis set handle
      integer jsh       !< [Input] "j" contraction index
      integer lscr      !< [Input] length of scratch space "scr"
      integer lOva      !< [Input] number of overlap integral derivatives in shells ish and jsh
      double precision scr(lscr) !< [Scratch] Work space
      double precision Ova(lOva) !< [Output] The overlap integral derivatives
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      integer idatom(*) !< [Output] array identifying centers for derivatives
c                       ! e.g., the first nint*3  derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
      logical hlfdrv    ! [input] if true, calculate same atom as well
                        ! as different atom integrals. This is necessary
                        ! when calculating the half-derivative overlap
                        ! matrix for NADCTs
c
      hlfdrv=.false.
      CALL intd_1eovQ(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lOva,Ova,idatom,hlfdrv)

      end
C--------------------------------------------------------------------------------------------

C> \ingroup nwint
C> @{
C>
C> \brief Compute the overlap integral derivatives
C>
C> Compute the overlap integral derivatives defined as
C> \f{eqnarray*}{
C>   \frac{\partial (\mu|\nu)}{\partial X_x} 
C>     &=& \int \frac{\partial [g_\mu(X_\mu,r_1)g_\nu(X_\nu,r_1)]}{\partial X_x}dr_1
C> \f}
C> On output the overlap integral derivatives are stored as if the array
C> is dimensioned as `Ova(nint,ncoord,natom)`, where `nint` is the number 
C> of integrals, `ncoord` are the Cartesian coordinates x, y and z, and
C> is 2 for the number of expansion centers on which the basis functions are
C> positioned.
C>
      subroutine intd_1eovQ(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lOva,Ova,idatom,hlfdrv)
c $Id$
      implicit none
#include "stdio.fh"
#include "errquit.fh"
c::functions
      integer int_nint_cart
      external int_nint_cart
c::passed
      integer i_basis   !< [Input] ish basis set handle
      integer ish       !< [Input] "i" contraction index
      integer j_basis   !< [Input] jsh basis set handle
      integer jsh       !< [Input] "j" contraction index
      integer lscr      !< [Input] length of scratch space "scr"
      integer lOva      !< [Input] number of overlap integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      integer idatom(*) !< [Output] array identifying centers for derivatives
c                       ! e.g., the first nint*3  derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
      logical hlfdrv    ! [input] if true, calculate same atom as well
                        ! as different atom integrals. This is necessary
                        ! when calculating the half-derivative overlap
                        ! matrix for NADCTs
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
      double precision scr(lscr) !< [Scratch] Work space
      double precision Ova(lOva) !< [Output] The overlap integral derivatives
*
      integer nint
*
      nint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*2.gt.lOva) then
        write(luout,*) 'nint*3*2 = ',nint*3*2
        write(luout,*) 'lOva     = ',lOva
        call errquit('intd_1eov: nint>lOva error',911,
     &          INT_ERR)
      endif
      call intd_1eovP(i_basis,ish,j_basis,jsh,lscr,scr,nint,
     &    OVa,idatom,hlfdrv)
      end
      
      subroutine intd_1eovP(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,Ova,idatom,hlfdrv)
      implicit none
c
c returns 1e overlap derivative integrals for shells/contractions ish,jsh
c
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geom.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw
      logical cando_sim
      external cando_hnd_1e
      external cando_nw
      external cando_sim
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] "i" contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] "j" contraction index
      integer lscr      ! [input] length of scratch space "scr"
      integer nint      ! [input] number of overlap integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      integer idatom(*) ! [output] array identifying centers for derivatives
c                       ! e.g., the first nint*3  derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
      logical hlfdrv    ! [input] if true, calculate same atom as well
                        ! as different atom integrals. This is necessary
                        ! when calculating the half-derivative overlap
                        ! matrix for NADCTs
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
      double precision Ova(nint,3,2) ! [output] derivative integrals
c::local
      integer ucont,uconts
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom,isbas,icfS
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom,jsbas,jcfS
      logical any_spherical
      logical orel, oirel, ojrel, oNR
      logical ohnd_ok, onw_ok, osim_ok
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer nint_x, nint_s
      integer zatom, zyx
      integer lbas, sbas, abas
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c
c  check if gencont, or spshell
c
      call int_nogencont_check(i_basis,'intd_1eovP:i_basis')
      call int_nogencont_check(j_basis,'intd_1eovP:j_basis')
      call int_nospshell_check(i_basis,'intd_1eovP:i_basis')
      call int_nospshell_check(j_basis,'intd_1eovP:j_basis')
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
c
      iatom = (sf_ibs_cn2ce(ish,ibas))
      jatom = (sf_ibs_cn2ce(jsh,jbas))
      if ((iatom.eq.jatom).and.(.not.hlfdrv))then
        idatom(1) = 0
        idatom(2) = 0
        return
      endif
c
      idatom(1) = iatom
      idatom(2) = jatom
c
      if (dyall_mod_dir) then
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        orel = ibas .eq. abas .and. jbas .eq. abas
      end if
c
c   i shell
c
      ucont = (sf_ibs_cn2ucn(ish,ibas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'intd_1eov: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          ibas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        isbas = ibas
      end if
c
      inp   = infbs_cont(CONT_NPRIM,ucont,ibas)
      igen  = infbs_cont(CONT_NGEN,ucont,ibas)
      iexp  = infbs_cont(CONT_IEXP,ucont,ibas)
      icf   = infbs_cont(CONT_ICFP,ucont,ibas)
      itype = infbs_cont(CONT_TYPE,ucont,ibas) 
      igeom = ibs_geom(ibas)
      icfS  = infbs_cont(CONT_ICFP,uconts,isbas)
c
c   j shell
c
      ucont = (sf_ibs_cn2ucn(jsh,jbas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,jbas) .ne. 0)) then
        ojrel = .true.
        jsbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'intd_1eov: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          jbas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        jsbas = jbas
      end if
c
      jnp   = infbs_cont(CONT_NPRIM,ucont,jbas)
      jgen  = infbs_cont(CONT_NGEN,ucont,jbas)
      jexp  = infbs_cont(CONT_IEXP,ucont,jbas)
      jcf   = infbs_cont(CONT_ICFP,ucont,jbas)
      jtype = infbs_cont(CONT_TYPE,ucont,jbas) 
      jgeom = ibs_geom(jbas)
      jcfS  = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (igeom.ne.jgeom) then
        write(luout,*)'intd_1eovP: two different geometries for',
     &         ' derivatives?'
        call errquit('intd_1eovP: geom error ',911, INT_ERR)
      endif
c
      ohnd_ok = cando_hnd_1e(i_basis,ish,0)
     &    .and. cando_hnd_1e(j_basis,jsh,0)
     &    .and. (.not.geom_any_finuc (igeom))
     &    .and. (.not.geom_any_finuc (jgeom))
      onw_ok = cando_nw(i_basis,ish,0) .and. cando_nw(j_basis,jsh,0)
      osim_ok = cando_sim(i_basis,ish,0) .and. cando_sim(j_basis,jsh,0)
c
      if (orel) then
        call rel_oneld (
     &      coords(1,iatom,igeom),
     &      dbl_mb(mb_exndcf(iexp,ibas)),
     &      dbl_mb(mb_exndcf(icf,ibas)),
     &      dbl_mb(mb_exndcf(icfS,isbas)),inp,igen,itype,iatom,
     &      coords(1,jatom,jgeom),
     &      dbl_mb(mb_exndcf(jexp,jbas)),
     &      dbl_mb(mb_exndcf(jcf,jbas)),
     &      dbl_mb(mb_exndcf(jcfS,jsbas)),jnp,jgen,jtype,jatom,
     &      coords(1,1,igeom),charge(1,igeom),
     &      geom_invnucexp(1,igeom),ncenter(igeom),
c........................       doS   doT    doV    canAB     
     &      Ova,scr,scr,nint,.true.,.false.,.false.,.false.,onw_ok,
c...........       nonrel dryrun
     &      ohnd_ok,oNR,.false.,scr,lscr,rel_dbg,rel_typ)
      else if (onw_ok.or.osim_ok) then
        call hf1dQ(
     &      coords(1,iatom,igeom),
     &      dbl_mb(mb_exndcf(iexp,ibas)),
     &      dbl_mb(mb_exndcf(icf,ibas)),
     &      inp,igen,itype,iatom,
c
     &      coords(1,jatom,jgeom),
     &      dbl_mb(mb_exndcf(jexp,jbas)),
     &      dbl_mb(mb_exndcf(jcf,jbas)),
     &      jnp,jgen,jtype,jatom,
c
     &      coords(1,1,igeom),charge(1,igeom),
     &      geom_invnucexp(1,igeom),ncenter(igeom),
     &      Ova,scr,scr,nint,
c...........overlap, k-e,     pot-e,    canab,  dryrun
     &      .true.,  .false., .false., .false., .false.,
     &      scr,lscr,hlfdrv)
      elseif (ohnd_ok) then
        call hnd_stvintd(
     &      coords(1,iatom,igeom),
     &      dbl_mb(mb_exndcf(iexp,ibas)),
     &      dbl_mb(mb_exndcf(icf,ibas)),
     &      inp,igen,itype,iatom,
c
     &      coords(1,jatom,jgeom),
     &      dbl_mb(mb_exndcf(jexp,jbas)),
     &      dbl_mb(mb_exndcf(jcf,jbas)),
     &      jnp,jgen,jtype,jatom,
c
     &      coords(1,1,igeom),charge(1,igeom),ncenter(igeom),
     &      Ova,scr,scr,nint,
c............overlap, k-e,     pot-e,
     &      .true.,  .false., .false., 
     &      scr,lscr)
      else
        call errquit('intd_1eov: could not do hnd or nw integrals',
     &                0, INT_ERR)
      endif
c
*     Ova now has the cartesian integral block  (jlo:jhi,ilo:ihi,3,2)
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
            call spcart_tran1e(Ova(1,zyx,zatom),scr,
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
            call spcart_tran1e(Ova(1,Zyx,Zatom),scr,
     &          j_nbf_x,i_nbf_x,0,jgen,
     &          j_nbf_s,i_nbf_s,Itype,igen,
     &          .true.)
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
            call spcart_tran1e(Ova(1,Zyx,Zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,0,igen,
     &          .false.)
          enddo
        enddo
      else
        call errquit(
     &        'intd_1eovP: should never reach transform blocked else',
     &        911, INT_ERR)
      endif
c
c now shuffle transformed buffers to contiguous space
c
      nint_x = i_nbf_x*j_nbf_x
      nint_s = i_nbf_s*j_nbf_s
      if (nint_s.gt.nint_x) then
        call errquit
     &      ('intd_1eov: nint_s >.nint_x diff=',(nint_s-nint_x),
     &         INT_ERR)
      elseif (nint_s.eq.nint_x) then
        return
      else
        call int_c2s_mv
     &      (Ova,nint_x,nint_s,(3*2),scr,lscr,'intd_1eov')
      endif
c
      end
C> @}
