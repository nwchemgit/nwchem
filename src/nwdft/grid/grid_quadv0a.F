      Subroutine grid_quadv0a(g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess, qwght, qxyz, xyz,iscr,nscr,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq, rdens_atom, cetobfr, 
     ,ibf_ao,rqbuf,rchi_atom,rdelchi_atom,
     ,Ex, Ec, Amat,Cmat,g_vxc,
     ,Amat2, Cmat2,
     ,nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, tmat, l3d,
     &triplet,tau,ttau,Mmat,Mmat2,
     & StericEnergyTotal)

c
c$Id: grid_quadv0a.F,v 1.46 2008-02-01 22:27:45 niri Exp $
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
#include "zora.fh"
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

      logical grid_file_rewind
      external grid_file_rewind
      logical xc_chkgrad, xc_chktau, kske
      external xc_chkgrad, xc_chktau
      logical l3d
      logical triplet
      double precision StericEnergyTotal, StericEnergy

c
c     see if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad()   ! checks if we have a gradient xc functional
      kske=xc_chktau()    ! checks if we have a tau xc functional

      if(n_rec_in_file.eq.0) goto 100
c
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit('grid_quadv0a: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in grid_quadv0a',0, BASIS_ERR)
c
c     Allocate and create info for new basis function evaluator
c
      if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natoms)
      if (xcfit) then
         if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_xinfo',
     &        lbas_cent_xinfo, ibas_cent_xinfo))
     &        call errquit('grid_quadv0a: cant alloca bas_cent_info',0,
     &       MA_ERR)
      if (.not.bas_numcont(xc_bas_han, ncontrsetx))
     &     call errquit('Exiting in grid_quadv0a',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrsetx, 'bas_cset_xinfo',
     &     lbas_cset_xinfo, ibas_cset_xinfo))
     &     call errquit('grid_quadv0a: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrsetx, 'docsetx',
     &     ldocsetx, idocsetx))
     &     call errquit('grid_quadv0a: cannot allocate docset',
     .     ncontrset, MA_ERR)
      call xc_make_basis_info(xc_bas_han, int_mb(ibas_cent_xinfo),
     &     int_mb(ibas_cset_xinfo), natoms)
      endif
      if (.not.MA_Push_Get(mt_int, 6*n_per_rec, 'iscratch',
     &     l_iscratch, i_iscratch))
     &     call errquit('grid_quadv0a: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('grid_quadv0a: cannot allocate docset',
     .     ncontrset, MA_ERR)
c
      ntot_cube=0
      StericEnergyTotal = 0.d0 ! steric energy
c 
      do 200 iqsh = 1, n_rec_in_file
c
c     Define the current range of radial shells and integration center.
c     
         call grid_file_read(n_per_rec, nqpts, ictr_buf, 
     &        rad,rqbuf,nsubb)
         if(nqpts.gt.buffer_size) call
     '        errquit(' buffersize exceed by qpts ',nqpts, DISK_ERR)
        if(nqpts.eq.0) goto 200
        istep=0

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
     &       nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, tmat,
     &       l3d,triplet,
     &       int_mb(ibas_cent_xinfo),int_mb(ibas_cset_xinfo),
     &       log_mb(idocsetx),
     &       tau, ttau, kske, Mmat, Mmat2, StericEnergy)
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
c
      call ga_sync
      return
      end
