C> \ingroup nwint
C> @{
C>
C> \brief Computes the electron - COSMO charge interactions
C>
C> Compute derivatives of the integrals
C> \f{eqnarray*}{
C>    (\mu\rho|q) = \int_{-\infty}^{\infty}
C>    \chi_\mu(R_A;r)\chi_\rho(R_B;r)|r-R_C|^{-1}q_C dr
C> \f}
C> In this case we need all the derivatives with respect to centers
C> \f$A\f$, \f$B\f$, and \f$C\f$. The reason is that each COSMO charge
C> is associated with a particular atom to which it is fixed. So there
C> are a number of scenarios:
C>
C> - The COSMO charge is associated with C which is neither A or B.
C>   The COSMO charge contributes a force to all 3 centers
C>
C> - The COSMO charge is associated with A.
C>   The COSMO charge imposes a force on the basis functions centered
C>   on A and the basis functions generate a force on the COSMO charge,
C>   but these now cancel. Only the forces on the functions centered
C>   at A and the COSMO charge associated with A due to the basis 
C>   functions centered at B survive.
C>
C> - The COSMO charge is associated with B.
C>   The is the same as above but now the role of centers A and B is 
C>   interchanged.
C>
C> - The COSMO charge is associated with A and A equals B.
C>   In this case all forces must cancel.
C>
C> In practice the derivative integrals can be calculated in a number
C> of ways. Because we need the derivatives with respect to all three
C> centers a suitable approach is to calculate
C> \f{eqnarray*}{
C>    \frac{d(\mu\rho|q)}{dA} && \\\\
C>    \frac{d(\mu\rho|q)}{dB} && \\\\
C>    \frac{d(\mu\rho|q)}{dC} &=& -\frac{d(\mu\rho|q)}{dA}
C>                                -\frac{d(\mu\rho|q)}{dB}
C> \f}
C> The advantage is that the top two lines are easy to evaluate. The
C> derivative of a Gaussian function with respect to the coordinates
C> of its center is again a Gaussian function (just with a higher 
C> angular momentum). So we can use the standard integral codes to
C> evaluate these expressions. The derivative of center C can then be
C> constructed using translational invariance, as specified in the
C> last equation. 
c
      subroutine intd_1epot_cosmo(i_basis,ish,j_basis,jsh,lscr,scr,
     &           lH1a,H1a,efc_c,efc_q,nefc)
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
c                       !< NOTE: nint*3 integral derivatives returned per unique center
      integer nefc      !< [Input] the number of COSMO charges
      double precision efc_c(3,nefc) !< [Input] COSMO charge coordinates
      double precision efc_q(nefc) !< [Input] COSMO charges
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
      call intd_1epotP_cosmo(i_basis,ish,j_basis,jsh,
     &       scrsize,scr(offset),nint,H1a,efc_c,efc_q,nefc)
c
      end
      subroutine intd_1epotP_cosmo(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,H1a,efc_c,efc_q,nefc)
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
#include "prop.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw
      external cando_hnd_1e
      external cando_nw
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] ``i'' contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] ``j'' contraction index
      integer lscr      ! [input] length of scratch space 
      integer nint      ! [input] number of integrals in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      integer nefc      ! [input] the number of COSMO charges
      double precision efc_c(3) ! [input] the COSMO charges coordinates
      double precision efc_q ! [input] the COSMO charges
      double precision H1a(nint,3,3)    ! [output] derivative integrals (nint,3,n_atoms)
      double precision scr(lscr) ! [scratch] 
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
      logical ohnd_ok, onw_ok
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer nint_x, nint_s
      integer zatom, zyx
      integer lbas, sbas, abas
c 
      integer ib,bq_ncent,nb
      integer i_qbq,i_cbq
      integer i_qbq0,i_cbq0,i_xbq0
      integer h_qbq0,h_cbq0,h_xbq0
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
      ohnd_ok = .false.
      onw_ok = cando_nw(i_basis,ish,0) .and. cando_nw(j_basis,jsh,0)
#ifdef USE_SIMINT
c hack to get cosmo gradient running when Simints ints are used      
      if(.not.orel) onw_ok=.true.
#endif      
c
c     get external charges here (MV)
c     ------------------------------
c     if(.not.geom_extbq_on())
c    >   call errquit('int_1eefc:no active bqs',0,0) 
c     nefc = ncenter(cosmo_geom_efc)
c     i_cbq = geom_extbq_coord()
c     i_qbq = geom_extbq_charge()
c     nb = nefc+ncenter(igeom)
c     nb = nefc+2
      nb = 3

c     if(.not.ma_push_get(mt_dbl,3*nb,'cbq',h_cbq0,i_cbq0))
c    +     call errquit( 'intd1_epot_cosmo',0,0)
c     if(.not.ma_push_get(mt_dbl,nb,'qbq',h_qbq0,i_qbq0))
c    +     call errquit( 'intd1_epot_cosmo',0,0)
c     if(.not.ma_push_get(mt_dbl,nb,'xbq',h_xbq0,i_xbq0))
c    +     call errquit( 'intd1_epot_cosmo',0,0)

c     call dfill(nb,0.d0,dbl_mb(i_xbq0),1)
c
c     Create joined geometry data for the QM centers + COSMO centers
c
c     call dcopy(3,coords(1,iatom,igeom),1,dbl_mb(i_cbq0),1)
c     call dcopy(3,coords(1,jatom,igeom),1,dbl_mb(i_cbq0+3),1)
c     call dcopy(3*nefc,efc_c,1,dbl_mb(i_cbq0+6),1)
c
c     call dcopy(1,efc_q,1,dbl_mb(i_qbq0),1)
c     call dcopy(1,efc_q,1,dbl_mb(i_qbq0+1),1)
c     call dcopy(nefc,efc_q,1,dbl_mb(i_qbq0+2),1)
c
c     call dcopy(1,geom_invnucexp(iatom,igeom),1,dbl_mb(i_xbq0),1)
c     call dcopy(1,geom_invnucexp(jatom,igeom),1,dbl_mb(i_xbq0+1),1)
c
c     call dcopy(3*ncenter(igeom),coords(1,1,igeom),1,
c    >           dbl_mb(i_cbq0),1)
c     call dcopy(3*nefc,coords(1,1,cosmo_geom_efc),1,
c    >           dbl_mb(i_cbq0+3*ncenter(igeom)),1)
c
c     call dcopy(ncenter(igeom),charge(1,igeom),1,
c    >           dbl_mb(i_qbq0),1)
c     call dcopy(nefc,charge(1,cosmo_geom_efc),1,
c    >           dbl_mb(i_qbq0+ncenter(igeom)),1)
c
c     call dcopy(ncenter(igeom),geom_invnucexp(1,igeom),1,
c    >           dbl_mb(i_xbq0),1)
c
c     Junk
c
c     call dcopy(3,coords(1,iatom,igeom),1,dbl_mb(i_cbq0),1)
c     call dcopy(3,coords(1,jatom,igeom),1,dbl_mb(i_cbq0+3),1)
c     call dcopy(3*nefc,coords(1,1,cosmo_geom_efc),1,
c    >           dbl_mb(i_cbq0+6),1)
c
c     call dcopy(1,1.0d0,1,dbl_mb(i_qbq0),1)
c     call dcopy(1,1.0d0,1,dbl_mb(i_qbq0+1),1)
c     call dcopy(nefc,charge(1,cosmo_geom_efc),1,dbl_mb(i_qbq0+2),1)
c
c     call dcopy(1,geom_invnucexp(iatom,igeom),1,dbl_mb(i_xbq0),1)
c     call dcopy(1,geom_invnucexp(jatom,igeom),1,dbl_mb(i_xbq0+1),1)
c        
      if (orel) then
        call rel_oneld_cosmo (
     &      coords(1,iatom,igeom),
     &      dbl_mb(mb_exndcf(iexp,ibas)),
     &      dbl_mb(mb_exndcf(icf,ibas)),
     &      dbl_mb(mb_exndcf(icfS,isbas)),inp,igen,itype,iatom,
     &      coords(1,jatom,jgeom),
     &      dbl_mb(mb_exndcf(jexp,jbas)),
     &      dbl_mb(mb_exndcf(jcf,jbas)),
     &      dbl_mb(mb_exndcf(jcfS,jsbas)),jnp,jgen,jtype,jatom,
c     &      coords(1,1,cosmo_geom_efc),charge(1,cosmo_geom_efc),
c     &      geom_invnucexp(1,cosmo_geom_efc),ncenter(cosmo_geom_efc),
     &      efc_c,efc_q,
     &      H1a,nint,.false.,onw_ok,
c    &      dbl_mb(i_cbq0),dbl_mb(i_qbq0),
c    &      dbl_mb(i_xbq0),nb,
c........................     doS   doT  doV    canAB  
c    &      scr,Ta,H1a,nint,.false.,doT,.true.,.false.,onw_ok,
c...........       nonrel dryrun
     &      ohnd_ok,oNR,.false.,scr,lscr,rel_dbg,rel_typ)
      else if (onw_ok) then
        call hf1d_cosmo(
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
c    &      iat_efc,
c    &      coords(1,1,cosmo_geom_efc),charge(1,cosmo_geom_efc),
c    &      geom_invnucexp(1,cosmo_geom_efc),ncenter(cosmo_geom_efc),
     &      efc_c,efc_q,
     &      H1a,nint,
c........... canab,   dryrun
     &      .false., .false.,
     &      scr,lscr)
c     elseif (ohnd_ok) then
c       call hnd_stvintd_cosmo(
c    &      coords(1,iatom,igeom),
c    &      dbl_mb(mb_exndcf(iexp,ibas)),
c    &      dbl_mb(mb_exndcf(icf,ibas)),
c    &      inp,igen,itype,iatom,
c
c    &      coords(1,jatom,jgeom),
c    &      dbl_mb(mb_exndcf(jexp,jbas)),
c    &      dbl_mb(mb_exndcf(jcf,jbas)),
c    &      jnp,jgen,jtype,jatom,
c
c    &      efc_c,efc_q,nefc,
c    &      H1a,nint,
c    &      scr,lscr)
      else
        call errquit('intd_1epotP_cosmo: could not do nw integrals',
     &                0, INT_ERR)
      endif
c
c     if(.not.ma_pop_stack(h_xbq0))
c    +     call errquit( 'intd1_epot',0,0)
c     if(.not.ma_pop_stack(h_qbq0))
c    +     call errquit( 'intd1_epot',0,0)
c     if(.not.ma_pop_stack(h_cbq0))
c    +     call errquit( 'intd1_epot',0,0)

c
c check for ecp 
c
*
* this should move to hf1dsp when sp is enabled.
*
c     nat = ncenter(igeom)  ! needed for both ecp and spherical
c     if (any_ecp) then
c       nintV = int_nbf_x(itype)*int_nbf_x(jtype)
c       offset = nintV*3*nat + 1
c       call errquit("intd_1epot_cosmo: not done ecps yet",0,CAPMIS_ERR)
c       call intd_ecp_hf1(
c    &      coords(1,iatom,igeom),
c    &      dbl_mb(mb_exndcf(iexp,ibas)),
c    &      dbl_mb(mb_exndcf(icf,ibas)),
c    &      inp,igen,itype,iatom,
c
c    &      coords(1,jatom,jgeom),
c    &      dbl_mb(mb_exndcf(jexp,jbas)),
c    &      dbl_mb(mb_exndcf(jcf,jbas)),
c    &      jnp,jgen,jtype,jatom,
c
c    &      scr,nintV,nat,
c    &      scr(offset),(lscr-offset-1),
c    &      .false.)
c 
*... sum ecp into derivative H1 block
c       call daxpy(nintV*3*nat,1.0d00,scr,1,H1a,1)
c     endif
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
        do zatom = 1,3
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
        do zatom = 1,3
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
        do zatom = 1,3
          do zyx = 1,3
            call spcart_tran1e(H1a(1,zyx,zatom),scr,
     &        j_nbf_x,i_nbf_x,Jtype,jgen,
     &        j_nbf_s,i_nbf_s,0,igen,
     &        .false.)
          enddo
        enddo
      else
        call errquit(
     &        'int_1epot_cosmo: should never reach transform blocked '//
     &        'else', 911, INT_ERR)
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
     &      (H1a,nint_x,nint_s,9,scr,lscr,'intd_1eh1')
      endif
c

      end
C> @}
