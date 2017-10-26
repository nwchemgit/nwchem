C> \ingroup nwint
C> @{
C>
C> \brief Compute the 1-electron derivative integrals
C>
C> This routine compute the 1-electron derivative integrals. The routine
C> computes the nuclear attraction integral derivatives
C> \f{eqnarray*}{
C>   \sum_A\frac{\partial (\mu|Z_AR_A^{-1}|\nu)}{\partial R_x} &=&
C>      \sum_A\int \frac{\partial [g_\mu(X_\mu,r_1)Z_AR^{-1}_{A1}g_\nu(X_\nu,r_1)]}
C>                      {\partial X_x}dr_1
C> \f}
C> It also adds the kinetic energy integral derivatives
C> \f{eqnarray*}{
C>   -\frac{1}{2}\frac{\partial (\mu|\nabla_1^2|\nu)}{\partial R_x} &=&
C>   -\frac{1}{2}\int \frac{\partial [g_\mu(X_\mu,r_1)\nabla_1^2g_\nu(X_\nu,r_1)]}
C>                         {\partial X_x}dr_1
C> \f}
C> The results are stored in `H1a` as if that variable was declared as
C> `H1a(nint,ncoord,natom)` where `nint` is the number of integrals for a given
C> pair of shells, `ncoord` is the number of Cartesian coordinates, and `natom`
C> is the number of atoms.
C>
      subroutine intd_1eh1(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lH1a,H1a)
C $Id$
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
c layer routine to compute the derivative 1 electron hamiltonian integrals
c for shells/contractions ish,jsh
c
c Order is...   nint*3*nat (3=> xyz, nat=number of atoms)
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
c                                           |
c                       nint,d <ij>         |
c                            -------------- |
c                            d[idatom(2),z] |
c
c                                  . . . 
c                                                            |
c                                         nint,d <ij>        |
c                                              --------------|
c                                            d[idatom(nat),x]|
c                                                                |
c                                             nint,d <ij>        |
c                                                  --------------|
c                                                d[idatom(nat),y]|
c                                                                    \
c                                                 nint,d <ij>         |
c                                                      -------------- |
c                                                    d[idatom(nat),z]/
c      
c::functions
      integer int_nint_cart
      external int_nint_cart
c::passed
      integer i_basis   !< [Input] ish basis set handle
      integer ish       !< [Input] ``i'' contraction index
      integer j_basis   !< [Input] jsh basis set handle
      integer jsh       !< [Input] ``j'' contraction index
      integer lscr      !< [Input] length of scratch space
      integer lH1a      !< [Input] number of h1 integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      double precision scr(lscr) !< [Input] scratch array
      double precision H1a(*)    !< [Output] derivative integrals
c
c::local
      integer nint, offset, scrsize, nat
c
      nat = ncenter(ibs_geom((i_basis + Basis_Handle_Offset)))
c
      nint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*nat.gt.lH1a) then
        write(luout,*) 'nint*3*nat = ',nint*3*nat
        write(luout,*) 'lH1a       = ',lH1a
        call errquit('intd_1eh1: nint>lH1a error',911, INT_ERR)
      endif
c
      offset = nint*3*2       ! scratch for Ta array in intd_1eh1P routine
      scrsize = lscr - offset ! new scratch array size
      offset = offset + 1     ! increment for passing to intd_1eh1P
c
      call intd_1eh1P(i_basis,ish,j_basis,jsh,
cc AJL/Begin/SPIN ECPs
     &       scrsize,scr(offset),nint,H1a,scr,1)
c                                             ^ ECP Integral (ALPHA)
c
      end
      subroutine intd_1eh1P(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,H1a,Ta,ecp_channel)
cc AJL/End
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "apiP.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geom.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "rel_nwc.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw
      logical cando_sim
      external cando_hnd_1e
      external cando_sim
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] ``i'' contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] ``j'' contraction index
      integer lscr      ! [input] length of scratch space 
      integer nint      ! [input] number of integrals in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      double precision scr(lscr) ! [input] scratch array
      double precision H1a(nint,3,*)    ! [output] derivative integrals (nint,3,n_atoms)
      double precision Ta(nint,3,2)     ! [scratch] space for kinetic integrals
cc AJL/Begin/SPIN ECPs
      integer ecp_channel !< [Input] ECP Channel (1 = Alpha, 2 = Beta)
cc AJL/End
c::local
      logical doT
      integer ucont,uconts
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom,isbas,icfS
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom,jsbas,jcfS
      integer nat
      integer nintV
      integer offset
c
      logical any_spherical
      logical orel, oirel, ojrel, oNR
      logical ohnd_ok, onw_ok,osim_ok
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer nint_x, nint_s
      integer zatom, zyx
      integer lbas, sbas, abas
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c  check if gencon/sp shells
c
      call int_nogencont_check(i_basis,'intd_1eh1P:i_basis')
      call int_nogencont_check(j_basis,'intd_1eh1P:j_basis')
      call int_nospshell_check(i_basis,'intd_1eh1P:i_basis')
      call int_nospshell_check(j_basis,'intd_1eh1P:j_basis')
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
c
      if (dyall_mod_dir) then
c
c     get basis set handles; relativistic integral option only valid
c     if both ibas and jbas are the ao basis.
c
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
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'intd_1eh1: no relativistic pointer',911, INT_ERR)
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
      iatom = (sf_ibs_cn2ce(ish,ibas))
      icfS  = infbs_cont(CONT_ICFP ,uconts,isbas)
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
     &      'intd_1eh1: no relativistic pointer',911, INT_ERR)
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
      jatom = (sf_ibs_cn2ce(jsh,jbas))
      jcfS  = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (igeom.ne.jgeom) then
        write(luout,*)'intd_1eh1P.F: two different geometries for',
     &         ' derivatives?'
        call errquit('intd_1eh1P: geom error ',911, INT_ERR)
      endif
c
      if (iatom.eq.jatom) then
        doT = .false.
      else
        doT = .true.
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
c........................     doS   doT  doV    canAB  
     &      scr,Ta,H1a,nint,.false.,doT,.true.,.false.,onw_ok,
c...........       nonrel dryrun
     &      ohnd_ok,oNR,.false.,scr,lscr,rel_dbg,rel_typ)
c     send it here for SImint bit
      else if (onw_ok.or.osim_ok) then
        call hf1d(
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
     &      scr,Ta,H1a,nint,
c..............overlap, k-e,  pot-e,  canab,   dryrun
     &      .false., doT, .true., .false., .false.,
     &      scr,lscr)
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
     &      scr,Ta,H1a,nint,
c............overlap, k-e,     pot-e,
     &      .false.,  doT, .true.,
     &      scr,lscr)
      else
        call errquit('intd_1eh1: could not do hnd or nw integrals',
     &                0, INT_ERR)
      endif
c
c if needed add in Ta derivative integrals
c
      if (doT) then
        call daxpy(nint*3,1.0d00,Ta(1,1,1),1,H1a(1,1,iatom),1)
        call daxpy(nint*3,1.0d00,Ta(1,1,2),1,H1a(1,1,jatom),1)
      endif
c
c check for ecp 
c
*
* this should move to hf1dsp when sp is enabled.
*
      nat = ncenter(igeom)  ! needed for both ecp and spherical
      if (any_ecp) then
        nintV = int_nbf_x(itype)*int_nbf_x(jtype)
        offset = nintV*3*nat + 1
*       write(luout,*)' lscr to ecp_hf1:',(lscr-nintV)
*       call dcopy(nintV*3*nat,0.0d00,0,scr,1) ! buffer zeroed in ecp_gradient
c
cc AJL/Begin/SPIN ECPs
        if (ecp_channel.eq.1) then
          call intd_ecp_hf1(
     &        coords(1,iatom,igeom),
     &        dbl_mb(mb_exndcf(iexp,ibas)),
     &        dbl_mb(mb_exndcf(icf,ibas)),
     &        inp,igen,itype,iatom,
c
     &        coords(1,jatom,jgeom),
     &        dbl_mb(mb_exndcf(jexp,jbas)),
     &        dbl_mb(mb_exndcf(jcf,jbas)),
     &        jnp,jgen,jtype,jatom,
c
     &        scr,nintV,nat,
     &        scr(offset),(lscr-offset-1),
     &        .false.)
c
        else
          call intd_ecp_hf1_beta(
     &        coords(1,iatom,igeom),
     &        dbl_mb(mb_exndcf(iexp,ibas)),
     &        dbl_mb(mb_exndcf(icf,ibas)),
     &        inp,igen,itype,iatom,
c
     &        coords(1,jatom,jgeom),
     &        dbl_mb(mb_exndcf(jexp,jbas)),
     &        dbl_mb(mb_exndcf(jcf,jbas)),
     &        jnp,jgen,jtype,jatom,
c
     &        scr,nintV,nat,
     &        scr(offset),(lscr-offset-1),
     &        .false.) 
c
        endif
c*... sum ecp into derivative H1 block
        call daxpy(nintV*3*nat,1.0d00,scr,1,H1a,1)
cc AJL/End
c 
      endif
c
*     H1a now has the cartesian integral block  (jlo:jhi,ilo:ihi)
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
        do zatom = 1,nat
          do zyx = 1,3
            call spcart_tran1e(H1a(1,zyx,zatom),scr,
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
        do zatom = 1,nat
          do zyx = 1,3
            call spcart_tran1e(H1a(1,zyx,zatom),scr,
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
        do zatom = 1,nat
          do zyx = 1,3
            call spcart_tran1e(H1a(1,zyx,zatom),scr,
     &        j_nbf_x,i_nbf_x,Jtype,jgen,
     &        j_nbf_s,i_nbf_s,0,igen,
     &        .false.)
          enddo
        enddo
      else
        call errquit(
     &        'int_1eh1: should never reach transform blocked else',
     &        911, INT_ERR)
      endif
c
c now shuffle transformed buffers to contiguous space
c
      nint_x = i_nbf_x*j_nbf_x
      nint_s = i_nbf_s*j_nbf_s
      if (nint_s.gt.nint_x) then
        call errquit
     &      ('intd_1eh1: nint_s >.nint_x diff=',(nint_s-nint_x),
     &      INT_ERR)
      elseif (nint_s.eq.nint_x) then
        return
      else
        call int_c2s_mv
     &      (H1a,nint_x,nint_s,(3*nat),scr,lscr,'intd_1eh1')
      endif
c

      end
cc AJL/Begin/SPIN ECPs
cc This is a wrapper, identical to intd_1epot except we have 
cc separated out for alpha and beta channels
      subroutine intd_1eh1_beta(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lH1a,H1a)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c::functions
      integer int_nint_cart
      external int_nint_cart
c::passed
      integer i_basis   !< [Input] ish basis set handle
      integer ish       !< [Input] ``i'' contraction index
      integer j_basis   !< [Input] jsh basis set handle
      integer jsh       !< [Input] ``j'' contraction index
      integer lscr      !< [Input] length of scratch space
      integer lH1a      !< [Input] number of h1 integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      double precision scr(lscr) !< [Input] scratch array
      double precision H1a(*)    !< [Output] derivative integrals
c
c::local
      integer nint, offset, scrsize, nat
c
      nat = ncenter(ibs_geom((i_basis + Basis_Handle_Offset)))
c
      nint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*nat.gt.lH1a) then
        write(luout,*) 'nint*3*nat = ',nint*3*nat
        write(luout,*) 'lH1a       = ',lH1a
        call errquit('intd_1eh1: nint>lH1a error',911, INT_ERR)
      endif
c
      offset = nint*3*2       ! scratch for Ta array in intd_1eh1P routine
      scrsize = lscr - offset ! new scratch array size
      offset = offset + 1     ! increment for passing to intd_1eh1P
c
      call intd_1eh1P(i_basis,ish,j_basis,jsh,
     &       scrsize,scr(offset),nint,H1a,scr,2)
c                                             ^ ECP Integral (BETA)
      end 
cc AJL/End
C> @}
