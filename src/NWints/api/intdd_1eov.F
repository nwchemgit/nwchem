C> \ingroup nwint
C> @{
C>
C> \brief Compute the overlap integral 2nd derivatives 
C>
C> Compute the overlap integral 2nd derivatives as given by
C> \f{eqnarray*}{
C>   \frac{\partial^2 (\mu|\nu)}{\partial X_x\partial X_y} 
C>     &=& \int \frac{\partial^2 [g_\mu(X_\mu,r_1)g_\nu(X_\nu,r_1)]}{\partial X_x\partial X_y}dr_1
C> \f}
C> The integrals are stored in `Ova` in an order consistent with the declaration
C> `Ova(nint,ncoord,ncoord)`, where `nint` is the number of integrals associated
C> with the shell pair, and `ncoord` is the number of Cartesian coordinates
C> associated with the two expansion centers. The variable `idatom` stores
C> the lexical indeces of the atoms associated with shells `ish` and `jsh`
C> respectively.
C>
      subroutine intdd_1eov(i_basis,ish,j_basis,jsh,lscr,scr,
     &       lOva,Ova,idatom)
c $Id$
      implicit none
#include "stdio.fh"
#include "errquit.fh"
c::functions
      integer int_nint_cart
      external int_nint_cart
c::passed
      integer i_basis   !< [Input] ish basis set handle
      integer ish       !< [Input] `i` contraction index
      integer j_basis   !< [Input] jsh basis set handle
      integer jsh       !< [Input] `j` contraction index
      integer lscr      !< [Input] length of scratch space `scr`
      integer lOva      !< [Input] number of overlap integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3 integral derivatives returned per unique center
      integer idatom(*) !< [Output] array identifying centers for derivatives
c                       ! e.g., the first nint*3  derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
c
c Order is...   nint*3*3*2 (3=> xyz, 2=atoms)
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
      double precision Ova(lOva) !< [Output] The overlap integral 2nd derivatives
*
      integer nint
*
      nint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (nint*3*3*4.gt.lOva) then
        write(luout,*) 'nint*3*3*4 = ',nint*3*3*4
        write(luout,*) 'lOva       = ',lOva
        call errquit('intdd_1eov: nint>lOva error',911, INT_ERR)
      endif
      call intdd_1eovP(i_basis,ish,j_basis,jsh,lscr,scr,nint,
     &    OVa,idatom)
      end
      
      subroutine intdd_1eovP(i_basis,ish,j_basis,jsh,lscr,scr,
     &       nint,Ova,idatom)
      implicit none
c
c returns 1e overlap second derivative integrals for shells/contractions ish,jsh
c
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
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
      logical cando_hnd_1edd
      logical cando_nw
      external cando_hnd_1edd
      external cando_nw
c::passed
      integer i_basis   ! [input] ish basis set handle
      integer ish       ! [input] "i" contraction index
      integer j_basis   ! [input] jsh basis set handle
      integer jsh       ! [input] "j" contraction index
      integer lscr      ! [input] length of scratch space "scr"
      integer nint      ! [input] number of overlap integral derivatives in shells ish and jsh
c                       ! NOTE: nint*3*3 integral derivatives returned per unique center
      integer idatom(*) ! [output] array identifying centers for derivatives
c                       ! e.g., the first nint*3  derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
c
c Order is...   nint*3*3*2 (3=> xyz, 2=atoms)
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
      double precision Ova(nint,3,3,4) ! [output] derivative integrals
c::local
      integer ucont
      integer ibas,iatom,inp,igen,iexp,icf,itype,igeom
      integer jbas,jatom,jnp,jgen,jexp,jcf,jtype,jgeom
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer nint_x, nint_s
      integer zatom, zyx1, zyx2
c
c  Temporary variable that needs to be taken out after testing!
c
c     integer itemp
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c
c  check if gencont, or spshell
c
      call int_nogencont_check(i_basis,'intdd_1eovP:i_basis')
      call int_nogencont_check(j_basis,'intdd_1eovP:j_basis')
      call int_nospshell_check(i_basis,'intdd_1eovP:i_basis')
      call int_nospshell_check(j_basis,'intdd_1eovP:j_basis')
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
        write(luout,*)'intdd_1eovP: two different geometries for',
     &         ' derivatives?'
        call errquit('intdd_1eovP: geom error ',911, GEOM_ERR)
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
     &       Ova,scr,scr,nint,
c............overlap, k-e,     pot-e,
     &       .true.,  .false., .false., 
     &       scr,lscr)
c     write(luout,*) 'S 1e- derivative integrals from hondo'
c     do itemp = 1,nint
c       write(luout,*) itemp,' 1,1,',Ova(itemp,1,1)
c       write(luout,*) itemp,' 2,1,',Ova(itemp,2,1)
c       write(luout,*) itemp,' 3,1,',Ova(itemp,3,1)
c       write(luout,*) itemp,' 1,2,',Ova(itemp,1,2)
c       write(luout,*) itemp,' 2,2,',Ova(itemp,2,2)
c       write(luout,*) itemp,' 3,2,',Ova(itemp,3,2)
c     enddo
      else
        call errquit('intdd_1eov: could not do hnd integrals',
     $        0, INT_ERR)
      endif
c
*     Ova now has the cartesian integral block  
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
        do zatom = 1,4
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(Ova(1,zyx1,zyx2,zatom),scr,
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
        do zatom = 1,4
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(Ova(1,zyx1,zyx2,zatom),scr,
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
        do zatom = 1,4
         do zyx2 = 1,3
          do zyx1 = 1,3
            call spcart_tran1e(Ova(1,zyx1,zyx2,zatom),scr,
     &          j_nbf_x,i_nbf_x,Jtype,jgen,
     &          j_nbf_s,i_nbf_s,0,igen,
     &          .false.)
          enddo
         enddo
        enddo
      else
        call errquit(
     &        'intdd_1eovP: cant do sphericals',
     &        911, INT_ERR)
      endif
c
c now shuffle transformed buffers to contiguous space
c
      nint_x = i_nbf_x*j_nbf_x
      nint_s = i_nbf_s*j_nbf_s
      if (nint_s.gt.nint_x) then
        call errquit
     &      ('intdd_1eov: nint_s >.nint_x diff=',(nint_s-nint_x),
     &       INT_ERR)
      elseif (nint_s.eq.nint_x) then
        return
      else
        call int_c2s_mv
     &      (Ova,nint_x,nint_s,(36),scr,lscr,'intdd_1eov')
      endif
c
      end
C> @}
