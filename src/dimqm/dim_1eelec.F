c $Id: int_1eelec.F 19696 2010-10-29 16:53:42Z d3y133 $
*
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1eelec}
c:tex-This routine computes the 1 electron electronic wave function (nder=-2) for shell ish, 
c:tex-electronic density (nder=-1), electrostatic potential (nder=0), 
c:tex-electric field (nder=1), and electric field gradient (nder=2).
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine dim_1eelec(i_basis,ish,j_basis,jsh,lscr,scr,lelec,elec,
     &                      nder,xyzpt,npts)
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e_prp
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e_prp
c::passed
c:tex-\begin{verbatim}
      integer i_basis ! [input] basis set handle for ish
      integer ish     ! [input] i shell/contraction
      integer j_basis ! [input] basis set handle for jsh
      integer jsh     ! [input] j shell/contraction
      integer lscr    ! [input] length of scratch array
      double precision scr(lscr)   ! [scratch] scratch array
      integer lelec                ! [input] length of elec buffer
      double precision elec(lelec) ! [output] elec integrals
      integer nder    ! [input] nder of integral
      double precision xyzpt(3,*) ! [input] points where integral is calculated
      integer npts    ! [input] number of points where integral is calculated
c:tex-\end{verbatim}
c::local
      integer igeom, jgeom, ibas, jbas, ucont
      integer itype, inp, igen, iexp, icent, icf, iatom
      integer jtype, jnp, jgen, jexp, jcent, jcf, jatom
c
      logical any_spherical, trani, tranj, shells_ok
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer ipts, ncartint ,i,j
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_1eelec'))
     &       call errquit('int_1eelec: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1eelec: invalid contraction/shell',0,
     &       INT_ERR)
c
c  check if gencont
c
      call int_nogencont_check(i_basis,'int_1eelec:i_basis')
      call int_nogencont_check(j_basis,'int_1eelec:j_basis')
c
      ibas = i_basis + basis_handle_offset
      jbas = j_basis + basis_handle_offset
c
      ucont = (sf_ibs_cn2ucn(ish,ibas))
      itype = infbs_cont(CONT_TYPE ,ucont,ibas)
      inp   = infbs_cont(CONT_NPRIM,ucont,ibas)
      igen  = infbs_cont(CONT_NGEN ,ucont,ibas)
      iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      icf   = infbs_cont(CONT_ICFP ,ucont,ibas)
      iatom = (sf_ibs_cn2ce(ish,ibas))
      igeom = ibs_geom(ibas)
c
      ucont = (sf_ibs_cn2ucn(jsh,jbas))
      jtype = infbs_cont(CONT_TYPE ,ucont,jbas)
      jnp   = infbs_cont(CONT_NPRIM,ucont,jbas)
      jgen  = infbs_cont(CONT_NGEN ,ucont,jbas)
      jexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      jcf   = infbs_cont(CONT_ICFP ,ucont,jbas)
      jatom = (sf_ibs_cn2ce(jsh,jbas))
      jgeom = ibs_geom(jbas)
c
      if (igeom.ne.jgeom) then
        write(luout,*)'int_1eelec: two different geometries for',
     &         ' properties?'
        call errquit('int_1eelec: geom error ',911, GEOM_ERR)
      endif
c
c     Determine # of cartesian integrals in block
c
      ncartint = int_nbf_x(itype)*int_nbf_x(jtype)
c
      call dim_elfder(
     &       coords(1,iatom,igeom),
     &       dbl_mb(mb_exndcf(iexp,ibas)),
     &       dbl_mb(mb_exndcf(icf,ibas)),
     &       inp,igen,itype,
c
     &       coords(1,jatom,jgeom),
     &       dbl_mb(mb_exndcf(jexp,jbas)),
     &       dbl_mb(mb_exndcf(jcf,jbas)),
     &       jnp,jgen,jtype,
c
     &       nder,ncartint,elec,scr,lscr,xyzpt,npts)
c
c     elec now has the cartesian integral block  
c     nder=-2 : (iblock,npts)
c     nder=-1 : (jblock,iblock,npts)
c     nder= 0 : (jblock,iblock,npts)
c     nder= 1 : (jblock,iblock,npts,3)
c     nder= 2 : (jblock,iblock,npts,6)
c      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
c      write(luout,*) "spherical?", any_spherical
      if (.not.any_spherical) return
c
      i_nbf_x = int_nbf_x(itype)
c
c     Make sure that in case of electronic wave function we only look at ish (jsh is not used)
c
      if (nder.eq.-2) then
         j_nbf_x  = 1
      else
         j_nbf_x = int_nbf_x(jtype)
      endif
c
c... assume we need to transform both i and j integrals
c
      trani = .true.
      tranj = .true. 
*.. do not tranform i component
      if (.not.bas_spherical(ibas)) trani = .false.
*.. do not tranform j component
      if (nder.eq.-2.or..not.bas_spherical(jbas)) tranj = .false.
c        
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4.
c
      if (itype.eq.-1) igen = 1
      if (jtype.eq.-1) jgen = 1
      ipts = npts*(max(3*nder,1))
      call spcart_2cBtran(elec,scr,lscr,
     &    j_nbf_x,int_nbf_s(jtype),jtype,jgen,tranj,
     &    i_nbf_x,int_nbf_s(itype),itype,igen,trani,
     &    ipts,.false.)
c
c     We now have the integrals in array (nsph_ints,npts,(max(3*nder,1)))
c
      return
      end
