      Subroutine grid_quadv0a(g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess, qwght, qxyz, xyz,iscr,nscr,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq, rdens_atom, cetobfr, 
     ,ibf_ao,rqbuf,rchi_atom,rdelchi_atom,
     ,Ex, Ec, Amat,Cmat,g_vxc,
     ,Amat2, Cmat2,
     ,nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r,do_xdm, tmat, l3d,
     &triplet,ttau,Mmat,Mmat2,
     & Lmat,
     & StericEnergyTotal,
cc AJL/Begin/FDE
c
     & do_gxc, Amat3, Cmat3, fde_option, g_dens_fde, ncenters_fde, 
     & rdens_atom_fde, Fmat_fde, Pmat_fde, Bmat_fde,
     & ttau_fde, rho_fde, delrho_fde, cetobfr_fde,
     & Amat_fde, Amat2_fde, Amat3_fde,
     & Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     & ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     & iniz_fde, ifin_fde, ibf_fde,
     & xyz_fde, rchi_atom_fde, rdelchi_atom_fde, rq_fde, rho_int_fde)
c
cc AJL/End
c
      implicit none
#include "errquit.fh"
c
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "util.fh"
#include "grid_quad.fh"
#include "steric.fh"
#include "msgids.fh"
c
      double precision rqbuf(*)
      double precision rad
      
      integer iqsh
      integer nqpts, ictr_buf
      integer ncube,istep,ntot_cube,ncontrset
      integer lbas_cent_info, ibas_cent_info,
     &        lbas_cset_info, ibas_cset_info,
     &        ldocset, idocset,i_iscratch,l_iscratch
      integer ncontrsetx,lbas_cent_xinfo, ibas_cent_xinfo,
     &        lbas_cset_xinfo, ibas_cset_xinfo,
     .     ldocsetx, idocsetx
      integer basis

      logical grid_file_rewind
      external grid_file_rewind
      logical xc_chkgrad, xc_chktau, xc_chklap, kske, kslap
      external xc_chkgrad, xc_chktau, xc_chklap
      logical l3d
      logical triplet
      integer k_becke
      integer l_becke
      logical oprint_qbecke
      double precision StericEnergyTotal, StericEnergy
c
c AJL/Begin/FDE
      integer ncontrset_fde
      integer lbas_cset_info_fde, ibas_cset_info_fde
      integer lbas_cent_info_fde, ibas_cent_info_fde
      integer ldocset_fde, idocset_fde
      character*32 pname

      pname = 'grid_quadv0a: '
      ibas_cset_info_fde = 0
      ibas_cent_info_fde = 0
      idocset_fde = 0
c AJL/End
c
c     see if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad()   ! checks if we have a gradient xc functional
      kske=xc_chktau()    ! checks if we have a tau xc functional
      kslap=xc_chklap()   ! checks if we have a lap xc functional
      basis = ao_bas_han
      if (adft) basis = xc_bas_han

      if(n_rec_in_file.eq.0) goto 100
c
      ibas_cent_xinfo = 0
      ibas_cset_xinfo = 0
      idocsetx = 0
c
c     becke charges
c
      oprint_qbecke = util_print('becke charges', print_high)
      if(oprint_qbecke) 
     R    call dft_initqbecke(natoms,k_becke,l_becke)
c
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit(pname//'rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
      if (.not.bas_numcont(basis, ncontrset))
     &     call errquit('Exiting in grid_quadv0a',0, BASIS_ERR)
c
c     Allocate and create info for new basis function evaluator
c
      if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit(pname//'cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit(pname//'cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(basis, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natoms)
      if (xcfit.and.(.not.adft)) then
         if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_xinfo',
     &        lbas_cent_xinfo, ibas_cent_xinfo))
     &        call errquit(pname//'cant alloca bas_cent_info',0,
     &       MA_ERR)
      if (.not.bas_numcont(xc_bas_han, ncontrsetx))
     &     call errquit('Exiting in grid_quadv0a',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrsetx, 'bas_cset_xinfo',
     &     lbas_cset_xinfo, ibas_cset_xinfo))
     &     call errquit(pname//'cannot allocate bas_cset_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrsetx, 'docsetx',
     &     ldocsetx, idocsetx))
     &     call errquit(pname//'cannot allocate docset',
     .     ncontrset, MA_ERR)
      call xc_make_basis_info(xc_bas_han, int_mb(ibas_cent_xinfo),
     &     int_mb(ibas_cset_xinfo), natoms)
      endif
      if (.not.MA_Push_Get(mt_int, 6*n_per_rec, 'iscratch',
     &     l_iscratch, i_iscratch))
     &     call errquit(pname//'cannot allocate bas_cset_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit(pname//'cannot allocate docset',
     .     ncontrset, MA_ERR)
c AJL/Begin/FDE
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then
       if (.not.bas_numcont(FD_bas_han, ncontrset_fde))
     &     call errquit(pname//'bas_numcont',0, BASIS_ERR)
       if (.not.MA_Push_Get(mt_int, 3*ncenters_fde, 'bas_cent_info_fde',
     &   lbas_cent_info_fde, ibas_cent_info_fde))
     &   call errquit(pname//'cannot allocate bas_cent_info_fde',0,
     &       MA_ERR)
       if (.not.MA_Push_Get(mt_int, 6*ncontrset_fde,'bas_cset_info_fde',
     &   lbas_cset_info_fde, ibas_cset_info_fde))
     &   call errquit(pname//'cannot allocate bas_cset_info_fde',0,
     &       MA_ERR)
       call xc_make_basis_info(FD_bas_han, int_mb(ibas_cent_info_fde),
     &     int_mb(ibas_cset_info_fde), ncenters_fde)
       if (.not.MA_Push_Get(mt_log, ncontrset_fde, 'docset_fde',
     &     ldocset_fde, idocset_fde))
     &     call errquit(pname//'cannot allocate docset_fde',
     .     ncontrset_fde, MA_ERR)
      endif
c AJL/End
c
      ntot_cube=0
      StericEnergyTotal = 0.d0 ! steric energy
c
c     write(*,*) "n_rec_in_file=",n_rec_in_file

      do 200 iqsh = 1, n_rec_in_file
c
c     Define the current range of radial shells and integration center.
c     
         call grid_file_read(n_per_rec, nqpts, ictr_buf, 
     &        rad,rqbuf,nsubb)
         if(nqpts.gt.buffer_size) call
     &        errquit(' buffersize exceed by qpts ',nqpts, DISK_ERR)
        if(nqpts.eq.0) goto 200
        istep=0

c       write(*,*) "nsubb=",nsubb

        do  ncube=1,nsubb
c     
c        temp put buf into currently used arrays qxyz and qwght
c     
         call grid_repack(rqbuf, qxyz, qwght, nqpts, rad,istep)
         if(nqpts.ne.0) then
            StericEnergy = 0.d0
            call grid_quadv0b(
     &       ictr_buf,qxyz, qwght, nqpts,rad,grad,
     &       g_dens, ncoef_max, natoms, 
     &       rho_int, nqtot, 
     &       iniz, ifin, rho, delrho, 
     &       hess,  xyz,iscr,nscr,
     &       expo, ccoef, Bmat, Fmat, 
     &       Pmat, rq, rdens_atom,  cetobfr, 
     &       ibf_ao,int_mb(ibas_cent_info),int_mb(ibas_cset_info),
     &       int_mb(i_iscratch),
     &       log_mb(idocset),rchi_atom,rdelchi_atom,
     &       Ex, Ec, Amat,Cmat,g_vxc,Amat2,Cmat2,
     &       nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm, 
     T       tmat, l3d,triplet,
     &       int_mb(ibas_cent_xinfo),int_mb(ibas_cset_xinfo),
     &       log_mb(idocsetx),
     &       ttau, kske, Mmat, Mmat2, 
     &       kslap, Lmat,
     &       StericEnergy,
     1       do_gxc, Amat3, Cmat3,
     Q       dbl_mb(k_becke),
cc AJL/Begin/FDE
     &       fde_option, g_dens_fde, ncenters_fde, rdens_atom_fde,
     &       Fmat_fde, Pmat_fde, Bmat_fde,ttau_fde, rho_fde, delrho_fde,
     &       cetobfr_fde, Amat_fde, Amat2_fde, Amat3_fde,
     &       Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     &       ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &       iniz_fde, ifin_fde, ibf_fde, xyz_fde, rchi_atom_fde, 
     &       rdelchi_atom_fde, rq_fde, int_mb(ibas_cent_info_fde), 
     &       int_mb(ibas_cset_info_fde),log_mb(idocset_fde),rho_int_fde)
cc AJL/End
c
             ntot_cube=ntot_cube+1
c
c            == Steric: Tally up the steric energy ==
             if (lsteric) then 
                  StericEnergyTotal = StericEnergyTotal + StericEnergy
             end if
c
         endif
      enddo
 200  continue
c
c     == Steric: Tally up over all nodes ==
      if (lsteric) then
         call ga_dgop(Msg_Excrho, StericEnergyTotal , 1, '+')  ! sum energy over all nodes
      end if
c
      if (.not.ma_chop_stack(lbas_cent_info))
     &     call errquit('grid_quadv0a: pop stack failed.',1, MA_ERR)
 100  continue
      if(oprint_qbecke) then
         call dft_printqbecke(geom,basis,ncenters,dbl_mb(k_becke))
         if (.not.ma_pop_stack(l_becke)) 
     Z        call errquit('grid_loop: pop stack failed.',0, MA_ERR)
      endif
c
      call ga_sync
      return
      end
c $Id$
      subroutine dft_printqbecke(geom,ao_bas_han,ncenters,
     C     q_becke)
      implicit none
#include "global.fh"
#include "sym.fh"
#include "stdio.fh"
#include "geom.fh"
#include "errquit.fh"
      integer geom
      integer ao_bas_han
      integer ncenters
      double precision q_becke(*)
c
      integer i_a
      double precision symfact
      character*2 symbol
      character*16 tag
      character*32 element
      integer atn
      double precision coord(3),qnuc
c
      call ga_dgop(2015, q_becke , ncenters, '+') 
      if(ga_nodeid().eq.0) then
         call dft_header
     &     (' Total Density - Becke Population Analysis')
        write(LuOut,1)
    1   format(/'    Atom         Charge   ')
        write(luout,11)
 11     format( ' -----------   -------- ')
         do i_a=1,ncenters
            if (.not.sym_atom(geom, i_a, symfact))  symfact=1d0
            q_becke(i_a)=q_becke(i_a)/symfact
            if(q_becke(i_a).gt.1d-12) then
               if (.not. geom_cent_get(geom, i_a, tag, coord, qnuc)) 
     $              call errquit('mull_pop: geom_cent_tag failed',0,
     &              GEOM_ERR)
               if (.not. geom_tag_to_element(tag, symbol, element, atn))
     $              symbol = 'X'
               write(LuOut,2) i_a,symbol,nint(qnuc),q_becke(i_a)
            endif
 2          format(1x,i4,1x,a2,i4,1x,f10.3)
         enddo
        write(LuOut,*)
      endif
      return
      end
      subroutine dft_initqbecke(natoms,k_becke,l_becke)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "mafdecls.fh"
      integer natoms
      integer k_becke
      integer l_becke
c
c     becke charges
c
      if (.not.MA_Push_Get(mt_dbl, natoms, 'becke charges',
     &     l_becke, k_becke))
     &     call errquit('dft_initqb: cant allocate becke q',0,MA_ERR)
      call dfill(natoms,0d0,dbl_mb(k_becke),1)
      return
      end
