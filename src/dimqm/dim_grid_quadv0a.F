      Subroutine dim_grid_quadv0a(ncoef_max, natoms, iniz, ifin,
     $                            rho, qwght, qxyz, xyz, expo, ccoef,
     $                            Bmat, Fmat, Pmat, rq, cetobfr,
     $                            ibf_ao, rqbuf, rchi_atom, g_vdim,
     $                            nmat, do_scf, do_cpks_l, l3d,
     $                            ipm, imag, dimxyz, muind)

      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "util.fh"
#include "steric.fh"
#include "msgids.fh"
#include "dimqm.fh"
c
c   Input Variables
      double precision rqbuf(*)
      double precision rad
      integer g_vdim
      double precision qxyz(3,n_per_rec)
      double precision qwght(n_per_rec)
      integer ncoef_max
      integer natoms
      integer iniz(natoms)
      integer ifin(natoms)
      double precision rho(*)
      double precision xyz(3,natoms)
      double precision expo(nbf_ao_mxprim)
      double precision Bmat(nbf_ao_mxnbf_ce*n_per_rec)
      double precision Pmat(*)
      double precision Fmat(*)
      double precision rq(n_per_rec, natoms)
      integer cetobfr(2,natoms)
      integer ibf_ao(nbf_ao)
      double precision rchi_atom(natoms)
      integer nmat
      logical do_scf
      logical do_cpks_l
      double precision ccoef(ncoef_max)
      
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
      integer ipm
      integer imag
      double precision dimxyz(3, nDIM)
      double precision muind(3, nDIM, nmat)
c
c      if(ldebug) then
c        write(luout,*) "dim_grid_quadv0a start"
c      end if
      
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

      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('grid_quadv0a: cannot allocate docset',
     .     ncontrset, MA_ERR)
c
      ntot_cube=0
c      write(luout,*) "n_recs_in_file:", n_rec_in_file
      do 200 iqsh = 1, n_rec_in_file
c      write(luout,*) "Top of loop", iqsh
c
c     Define the current range of radial shells and integration center.
c     
         call grid_file_read(n_per_rec, nqpts, ictr_buf, 
     &        rad,rqbuf,nsubb)
         if(nqpts.gt.buffer_size) call
     '        errquit(' buffersize exceed by qpts ',nqpts, DISK_ERR)
        if(nqpts.eq.0) goto 200
        istep=0
c        write(luout,*) iqsh, "nsubb:", nsubb
        do  ncube=1,nsubb
c     
c        temp put buf into currently used arrays qxyz and qwght
c     
         call grid_repack(rqbuf, qxyz, qwght, nqpts, rad,istep)
         if(nqpts.ne.0) then
            call dim_grid_quadv0b(
     $              ictr_buf, qxyz, qwght, nqpts, rad,
     $              ncoef_max, natoms, iniz, ifin, rho,
     $              xyz, expo, Bmat, Fmat, Pmat, rq,
     $              cetobfr, ibf_ao, int_mb(ibas_cent_info),
     $              int_mb(ibas_cset_info), log_mb(idocset), rchi_atom,
     $              g_vdim, nmat, do_scf, do_cpks_l, l3d, ipm, imag,
     $              dimxyz, muind)
c
             ntot_cube=ntot_cube+1
c
c
         endif
      enddo
c      write(luout,*) "Bottom of loop", iqsh
 200  continue
      call ga_igop(Msg_Excrho, ntot_cube , 1, '+')
c
      if (.not.ma_chop_stack(lbas_cent_info))
     &     call errquit('grid_quadv0a: pop stack failed.',1, MA_ERR)
 100  continue
c
c      if(ldebug) then
c        write(luout,*) "dim_grid_quadv0a end"
c      end if
      call ga_sync
      return
      end
