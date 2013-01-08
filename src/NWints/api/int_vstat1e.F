C> \ingroup nwint
C> @{
C> 
C> \brief Compute the nuclear attaction integrals for pseudo nucleii
C>
C> This subroutine is modified from int_1epe. It evaluates the nuclear 
C> attraction integrals for each pseudo nucleus at positions specified by the
C> array qxyz, which are just the quardrature points. The effective charges
C> of these pseudo nuclei are determined by the transition densities and/or 
C> the quardrature weights.
C>
c $Id$
      subroutine int_vstat1e(i_basis,ish,j_basis,jsh,lscr,scr,nqpts,
     &     vint,qxyz,zq,dens)
c
c     This subroutine is modified from int_1epe. It evaluates the nuclear 
c     attraction integrals for each pseudo nucleus at positions specified by the
c     array qxyz, which are just the quardrature points. The effective charges
c     of these pseudo nuclei are determined by the transition densities and/or 
c     the quardrature weights.
c
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
c
c     external subroutines used
c
      logical cando_hnd_1e
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e
      external cando_nw_1e
      external cando_nw
c
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr)    !< [Scratch] scratch array
      integer nqpts                 !< [Input] length of potential buffer
      double precision vint(nqpts)  !< [Output] potential integrals
      double precision qxyz(3,nqpts) !< [Input] coordinates of q points
      double precision zq(nqpts)    !< [Input] effective charges on q points
      double precision dens(*)      !< [Input] elements of density matrix
c::local
      logical shells_ok, orel, oirel, ojrel, oNR
      integer i_geom, j_geom, ibas, jbas, ucont
      integer lbas, sbas, isbas, jsbas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_icfpS
c
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer rel_dbg, rel_typ
      integer ig,jg
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_vstat1e'))
     &       call errquit('int_vstat1e: int_init was not called' ,0,
     &            INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_vstat1e: invalid contraction/shell',0,
     &               BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      isbas = ibas
      jsbas = jbas
      oNR = .true.
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS  = infbs_cont(CONT_ICFP ,ucont,isbas)
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS  = infbs_cont(CONT_ICFP ,ucont,jsbas)
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_vstat1e: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
      if (cando_hnd_1e(i_basis,ish,0).and.cando_hnd_1e(j_basis,jsh,0))
     &     then
         call hnd_vstat(
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),
     &        j_prim, j_gen, Lj,
     &        qxyz, zq, nqpts,dens,
c..............................doS     doT     doV
     &        scr,scr,vint,nqpts,.false.,.false.,.true.,
     &        scr,lscr)
c
      else
         call errquit('int_vstat1e: could not do hnd, try sp or nw 
     &        integrals?', 0, INT_ERR)
      endif
c      write(*,*)"vints in int_vstat1e"
c      write(*,"(4f12.8)")((qxyz(ig,jg),ig=1,3),vint(jg),jg=1,nqpts)
c
*     vint now has the cartesian integral block  (jlo:jhi,ilo:ihi)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
      if (any_spherical)
     &     call errquit("int_vstat1e: spherical basis not supported 
     &     yet", 100, BASIS_ERR)
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
c
      if (li.eq.-1) i_gen = 1
      if (lj.eq.-1) j_gen = 1
c     
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*...  transform both i and j integrals
         i_nbf_x = int_nbf_x(Li)
         i_nbf_s = int_nbf_s(Li)
         j_nbf_x = int_nbf_x(Lj)
         j_nbf_s = int_nbf_s(Lj)
c     
         call spcart_tran1e(vint,scr,
     &        j_nbf_x,i_nbf_x,Lj,j_gen,
     &        j_nbf_s,i_nbf_s,Li,i_gen,
     &        .false.)
      else if (bas_spherical(ibas)) then
*..   transform on i component
         i_nbf_x = int_nbf_x(Li)
         i_nbf_s = int_nbf_s(Li)
         j_nbf_x = int_nbf_x(Lj)
         j_nbf_s = j_nbf_x
         call spcart_tran1e(vint,scr,
     &        j_nbf_x,i_nbf_x,0,j_gen,
     &        j_nbf_s,i_nbf_s,Li,i_gen,
     &        .false.)
      else if (bas_spherical(jbas)) then
*..   transform on j component
         i_nbf_x = int_nbf_x(Li)
         i_nbf_s = i_nbf_x
         j_nbf_x = int_nbf_x(Lj)
         j_nbf_s = int_nbf_s(Lj)
         call spcart_tran1e(vint,scr,
     &        j_nbf_x,i_nbf_x,Lj,j_gen,
     &        j_nbf_s,i_nbf_s,0,i_gen,
     &        .false.)
      else
         call errquit(
     &        'int_vstat1e: should never reach transform blocked else',
     &        911, INT_ERR)
      endif
c     
      end
C> @}
