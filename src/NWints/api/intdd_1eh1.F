      subroutine intdd_1eh1(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lH1a,H1a)
C $Id: intdd_1eh1.F,v 1.8 2000-01-04 00:29:03 windus Exp $
      implicit none
#include "stdio.fh"
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
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] ``i'' contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] ``j'' contraction index
      integer lscr      ! [input] length of scratch space
      integer lH1a      ! [input] number of h1 integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      double precision scr(lscr) ! [input] scratch array
      double precision H1a(*)    ! [output] derivative integrals
c
c::local
      integer nint, offset, nat
c
      nat = ncenter(ibs_geom((i_basis + Basis_Handle_Offset)))
c
      nint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*3*nat*3.gt.lH1a) then
        write(luout,*) 'nint*3*3*nat*3 = ',nint*3*3*nat*3
        write(luout,*) 'lH1a       = ',lH1a
        call errquit('intd_1eh1: nint>lH1a error',911)
      endif
c
      offset = nint*3*3*nat*3 + 1  ! T is held seperately in H1a
c
      call intdd_1eh1P(i_basis,ish,j_basis,jsh,
     &       lscr,scr(offset),nint,H1a,H1a(offset),nat)
c
      end
      subroutine intdd_1eh1P(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,H1a,Ta,nat)
      implicit none
#include "stdio.fh"
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
      integer nat       ! [input] number of atoms
      integer nint      ! [input] number of integrals in shells ish and jsh
c                       ! NOTE: nint*3*3 integral derivatives returned per unique center
      double precision scr(lscr) ! [input] scratch array
      double precision H1a(nint,3,3,*)    ! [output] derivative integrals (nint,3,3,n_atoms)
      double precision Ta(nint,3,3,3)     ! [scratch] space for kinetic integrals
c::local
      logical doT
      integer ucont
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom
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
        call errquit('intdd_1eh1P: geom error ',911)
      endif
c
      if (iatom.eq.jatom) then
        doT = .false.
      else
        doT = .true.
      endif
c       
      if (cando_hnd_1e(i_basis,ish,0).and.cando_hnd_1e(j_basis,jsh,0))
     &    then
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
     &                0)
      endif
c
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
      call errquit(
     &        'intdd_1eh1P: cant do sphericals',
     &        911)
c

      end
      
