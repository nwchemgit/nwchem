C> \ingroup nwint
C> @{
C>
C> \brief Compute the integral 2nd derivatives of the 1-electron Hamiltonian
C>
C> Compute the integral 2nd derivatives of the 1-electron Hamiltonian
C> consisting of the kinetic energy integrals and the nuclear attraction
C> integrals. The kinetic energy integral 2nd derivatives are given by
C> \f{eqnarray*}{
C>   -\frac{1}{2}\frac{\partial^2 (\mu|\nabla_1^2|\nu)}{\partial R_x\partial R_y} &=&
C>   -\frac{1}{2}\int \frac{\partial^2 [g_\mu(X_\mu,r_1)\nabla_1^2g_\nu(X_\nu,r_1)]}
C>                         {\partial X_x\partial X_y}dr_1
C> \f}
C> The nuclear attraction integral 2nd derivatives are given by
C> \f{eqnarray*}{
C>   \sum_A\frac{\partial^2 (\mu|Z_AR_A^{-1}|\nu)}{\partial R_x\partial R_y} &=&
C>      \sum_A\int \frac{\partial^2 [g_\mu(X_\mu,r_1)Z_AR^{-1}_{A1}g_\nu(X_\nu,r_1)]}
C>                      {\partial X_x\partial R_y}dr_1
C> \f}
C> The integral 2nd derivatives are returned in `H1a` in an order that is
C> equivalent to the declaration `H1a(nint,ncoordu,ncoordv,natom,ncoorda)`,
C> where `nint` refers to the number of integrals associated with the shell
C> pair, `ncoordu` refers the number of Cartesian coordinates of the atom 
C> associated shell `ish` likewise `ncoordv` refers to atomic coordinates of the
C> `jsh` atom, `natom` is the number of atoms and `ncoorda` is the number
C> of coordinates of each nucleus.
C>
      subroutine intdd_1eh1(i_basis,ish,j_basis,jsh,lscr,scr,
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
      integer ish       !< [Input] `i` contraction index
      integer j_basis   !< [Input] jsh basis set handle
      integer jsh       !< [Input] `j` contraction index
      integer lscr      !< [Input] length of scratch space
      integer lH1a      !< [Input] number of h1 integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      double precision scr(lscr) !< [Input] scratch array
      double precision H1a(*)    !< [Output] derivative integrals
c
c::local
      integer nint, offset, nat
c
      nat = ncenter(ibs_geom((i_basis + Basis_Handle_Offset)))
c
      nint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*3*(nat*3+3).gt.lH1a) then
        write(luout,*) 'nint*3*3*(nat*3+3) = ',nint*3*3*(nat*3+3)
        write(luout,*) 'lH1a       = ',lH1a
        call errquit('intdd_1eh1: nint>lH1a error',911, INT_ERR)
      endif
c
      offset = nint*3*3*nat*3 + 1  ! T is held seperately in H1a
c
      call intdd_1eh1P(i_basis,ish,j_basis,jsh,
     &       lscr,scr,nint,H1a,H1a(offset),nat)
c
      end
      subroutine intdd_1eh1P(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,H1a,Ta,nat)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "apiP.fh"
#include "nwc_const.fh"
#include "int_nbf.fh"
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
      logical cando_hnd_1edd
      logical cando_nw
      external cando_hnd_1edd
      external cando_nw
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] ``i'' contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] ``j'' contraction index
      integer lscr      ! [input] length of scratch space 
      integer nat       ! [input] number of atoms
      integer nint      ! [input] number of integrals in shells ish and jsh
c                       ! NOTE: nint*3*3 integral derivatives returned per unique center
      double precision scr(lscr) ! [input] scratch array
      double precision H1a(nint,3,3,*)    ! [output] derivative integrals (nint,3,3,n_atoms,3)
      double precision Ta(nint,3,3,3)     ! [scratch] space for kinetic integrals
c::local
      logical doT
      integer ucont
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer nint_x, nint_s
      integer zatom, zyx1, zyx2
c
      logical any_spherical
c
c  Temporary variable that needs to be taken out after testing!
c
c     integer itemp,jtemp,ktemp,ltemp
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
c
      ucont = (sf_ibs_cn2ucn(ish,ibas))
      inp   = infbs_cont(CONT_NPRIM,ucont,ibas)
      igen  = infbs_cont(CONT_NGEN,ucont,ibas)
      iexp  = infbs_cont(CONT_IEXP,ucont,ibas)
      icf   = infbs_cont(CONT_ICFP,ucont,ibas)
      itype = infbs_cont(CONT_TYPE,ucont,ibas) 
      igeom = ibs_geom(ibas)
      iatom = (sf_ibs_cn2ce(ish,ibas))
c
      ucont = (sf_ibs_cn2ucn(jsh,jbas))
      jnp   = infbs_cont(CONT_NPRIM,ucont,jbas)
      jgen  = infbs_cont(CONT_NGEN,ucont,jbas)
      jexp  = infbs_cont(CONT_IEXP,ucont,jbas)
      jcf   = infbs_cont(CONT_ICFP,ucont,jbas)
      jtype = infbs_cont(CONT_TYPE,ucont,jbas) 
      jgeom = ibs_geom(jbas)
      jatom = (sf_ibs_cn2ce(jsh,jbas))
c
      if (igeom.ne.jgeom) then
        write(luout,*)'intdd_1eh1P.F: two different geometries for',
     &         ' derivatives?'
        call errquit('intdd_1eh1P: geom error ',911, GEOM_ERR)
      endif
c
      if (iatom.eq.jatom) then
        doT = .false.
      else
        doT = .true.
      endif
c       
      if (cando_hnd_1edd(i_basis,ish,0).and.
     &    cando_hnd_1edd(j_basis,jsh,0)) then
        call hnd_stvintdd(
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
     &       scr,Ta,H1a,nint,
c............overlap, k-e,     pot-e,
     &       .false.,  doT, .true.,
     &       scr,lscr)
      else
        call errquit('intdd_1eh1: could not do hnd integrals',
     &                0, INT_ERR)
      endif
c
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
        do zatom = 1,nat*3
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(H1a(1,zyx1,zyx2,zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,Itype,igen,
     &          .false.)
          enddo
         enddo
        enddo
        do zatom = 1,3
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(Ta(1,zyx1,zyx2,zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,Itype,igen,
     &          .false.)
          enddo
         enddo
        enddo
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Itype)
        i_nbf_s = int_nbf_s(Itype)
        j_nbf_x = int_nbf_x(Jtype)
        j_nbf_s = j_nbf_x
        do zatom = 1,nat*3
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(H1a(1,zyx1,zyx2,zatom),scr,
     &          j_nbf_x,i_nbf_x,0,jgen,
     &          j_nbf_s,i_nbf_s,Itype,igen,
     &          .false.)
          enddo
         enddo
        enddo
        do zatom = 1,3
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(Ta(1,zyx1,zyx2,zatom),scr,
     &          j_nbf_x,i_nbf_x,0,jgen,
     &          j_nbf_s,i_nbf_s,Itype,igen,
     &          .false.)
          enddo
         enddo
        enddo
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Itype)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Jtype)
        j_nbf_s = int_nbf_s(Jtype)
        do zatom = 1,nat*3
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(H1a(1,zyx1,zyx2,zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,0,igen,
     &          .false.)
          enddo
         enddo
        enddo
        do zatom = 1,3
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(Ta(1,zyx1,zyx2,zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,0,igen,
     &          .false.)
          enddo
         enddo
        enddo
      else
      call errquit(
     &        'intdd_1eh1P: cant do sphericals',
     &        911, INT_ERR)
      endif
c
c now shuffle transformed buffers to contiguous space
c
      nint_x = i_nbf_x*j_nbf_x
      nint_s = i_nbf_s*j_nbf_s
      if (nint_s.gt.nint_x) then
        call errquit
     &      ('intdd_1eh1: nint_s >.nint_x diff=',(nint_s-nint_x), 
     &      INT_ERR)
      elseif (nint_s.eq.nint_x) then
        return
      else
        call int_c2s_mv  ! do both H1a and Ta at the same time
     &      (H1a,nint_x,nint_s,(27*nat+27),scr,lscr,'intdd_1eh1')
      endif
c
      end
C> @}
