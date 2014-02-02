      subroutine grid_quadv0(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,tmat)
c
C$Id$
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(*)
      double precision rho_n
      double precision Exc(2)
      double precision tmat(*)  ! XC fit

      call grid_quadv0_gen(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,
     ,     1,1,tmat,.false.,.false.)

      return
      end
c
      subroutine grid_quadv0_gen(rtdb, g_dens, g_vxc,
     &      nexc, rho_n, Exc,
     &      nmat,  calc_type ,tmat,l3d,triplet)

      implicit none

#include "errquit.fh"
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(4)

      Logical oprint,grad,xc_chkgrad,oprintw
      Logical xc_chktau, kske
      external xc_chkgrad,xc_rep_close, xc_chktau
      logical xc_rep_close
      double precision Exc(2)
      double precision rho_n
      integer nmat              ! [input] number of XC matrices to make;
      double precision tmat(*)  ! XC fit
c                               !         only used for lhs cpks
      integer calc_type         ! what to calculate [input]
c                               !   1 = XC matrices for SCF
c                               !   2 = XC matrices for CPKS LHS
c                               !   3 = XC matrices for CPKS RHS
c                               !   4 = XC matrices for GIAO NMR RHS
c                               !   5 = XC matrices for third derivatives
      logical l3d               ! true if 3d GA is used
      logical triplet           ! true if triplet excitation energy calc
      double precision StericEnergyTotal
c
#include "cgridfile.fh"
#include "geom.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "dft2drv.fh"
#include "steric.fh"
#include "dft3drv.fh"
c
      logical rtdb_status       
      logical grid_reopen
      external grid_reopen
      logical lsa ! Qin
c
      integer me, nproc, nxyz, nctrs, 
     &        nbf_ao_mxnbf_ce2, nqtot, 
     &        ncoef_max,
     &        npol
c
      integer lxyz, ixyz, lcharge, icharge, 
     &        lqxyz, iqxyz, lqwght, iqwght, 
     &        lniz, iniz, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, ldelrho, idelrho, lBmat, iBmat,
     &        lhess, ihess, lrq, irq, 
     &        lfin, ifin, ltags, itags,
     &        lmmat, immat, lttau, ittau
c
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr,lamat,iamat,icmat,lcmat

      double precision nelec
c
      integer libf_ao,iibf_ao,irqbuf,lrqbuf,lbp,ibp,
     .     irchi_atom,irdelchi_atom,lrchi_atom,lrdelchi_atom
      integer iqs,lqs,nqlen,mxnrad,
     , laqbuf, iaqbuf,lRij,iRij,ndbl,l_qat,k_qat,irqlen
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2, l_Mmat2, k_Mmat2
      integer len_rho, len_delrho, len_Amat2, len_Cmat2, len_Mmat2
      integer len_mmat, len_ttau
      logical do_scf, do_cpks_l, do_cpks_r, do_nmrks_r
      integer max_elem_ang_scr, max_component
      integer nexpo,nibf
      double precision rho_int, ex, ec,dummy(3),job_grid_acc,rho_err
      integer ii,isp
      equivalence (dummy(1),Ex), (dummy(2),Ec), (dummy(3),rho_int)
c
      logical fon ! fractional occupancy
      double precision pstrace  ! a measure of the total charge
c
      logical do_gxc
      integer l_Amat3, k_Amat3, l_Cmat3, k_Cmat3
      integer len_Amat3, len_Cmat3
      logical lcpksl_gxc

#ifdef SOLARIS
c     needed for -stackvar option
      save ex, dummy, ec, rho_int
#endif      
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('quadrature', print_high)
      oprintw = util_print('quadrature warning', print_medium)
      grad=xc_chkgrad()
      kske = xc_chktau()
c
      irdelchi_atom = 0
c
      npol = (ipol*(ipol+1))/2
      if (.not.rtdb_get(rtdb,'dft:largenode', mt_log, 1, largenode))
     &     largenode=.false.
c
c     Open grid pts file
c     
cng
      if (grid_written .and. n_rec_in_file .le. 0 ) then
c
c       Somebody (fock_xc?) closed the grid file, better reopen it
c       before continuing. Note, fock_xc tends to open the grid file
c       at the beginning and close it at the end.
c
        if (.not.grid_reopen(geom)) call errquit(
     &    'grid_quadv0: could not open grid file',0,DISK_ERR)
      endif
cng
      if (.not. grid_written .or. n_per_rec .le. 0 ) then
         call grid_packing_info
      end if
c
      if (oprint) then      
        write(*,*) " *** grid_quadv0: calc_type:",calc_type 
        write(*,*) " *** grid_quadv0: grid_written:",grid_written 
        write(*,*) " *** grid_quadv0: grad:",grad 
        write(*,*) " *** grid_quadv0: kske:",kske 
      end if
c
c     Set up logicals denoting calculation type
c
      do_scf = calc_type .eq. 1
      do_cpks_l = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      do_nmrks_r = calc_type .eq. 4
      do_gxc = calc_type .eq. 5 ! for 3rd derivatives
      if (.not. (do_scf .or. do_cpks_l .or. do_cpks_r .or. 
     &    do_nmrks_r .or. do_gxc)) then
         call errquit('grid_quadv0_gen: illegal calculation type',1,
     &       INPUT_ERR)
      endif
c
      if (oprint) then
        write(*,*) " *** grid_quadv0: do_scf:",do_scf 
        write(*,*) " *** grid_quadv0: n_per_rec:",n_per_rec
      end if
c
      len_mmat = 1
      len_ttau = 1
      len_Amat2 = 1
      len_Cmat2 = 1
      len_Mmat2 = 1
      len_Amat3 = 1
      len_Cmat3 = 1 
c
      if (do_scf .or. do_nmrks_r) then
         len_rho = npol*n_per_rec
         len_delrho = 3*ipol*n_per_rec
         if (kske) len_ttau = npol*n_per_rec
         if (kske) len_mmat =ipol*n_per_rec  
      elseif (do_cpks_l) then
         len_rho = ipol*n_per_rec*(nmat+1)+n_per_rec
         len_delrho = 3*ipol*n_per_rec*(nmat+1)+3*n_per_rec
         if (kske) len_ttau = ipol*n_per_rec*(nmat+1)+n_per_rec
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
         if (kske) len_Mmat2 = n_per_rec*NCOL_MMAT2
      elseif (do_cpks_r) then
         call errquit('grid_quadv0: what for do_cpks_r?',0,
     &       UNKNOWN_ERR)
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
c
      elseif (do_gxc) then
c Here, len_rho and len_delrho reflect the fact that we have ipol ground
c state density matrices and ipol*nmat perturbed density matrices, where
c nmat depends on the number of roots requested.  These should be the
c correct dimensions.  The extra n_per_rec or 3*n_per_rec accounts for
c the total ground state density (or gradient) for unrestricted 
c calculations.
         len_rho = ipol*n_per_rec*(nmat+1) + n_per_rec
         len_delrho = 3*ipol*n_per_rec*(nmat+1) + 3*n_per_rec
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
         len_Amat3 = n_per_rec*NCOL_AMAT3
         len_Cmat3 = n_per_rec*NCOL_CMAT3
      endif
      if (.not. geom_ncent(geom, nctrs)) 
     &   call errquit('grid_den: geom_ncent failed',geom, GEOM_ERR)
      nxyz = 3*nctrs
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('grid_den: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'charge',lcharge,icharge))
     &   call errquit('grid_den: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,nctrs*16,'tags',ltags,itags))
     &   call errquit('grid_den: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, nctrs, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('grid_den: geom_cart_get failed',74, GEOM_ERR)
c
c      write(*,*)' Coordinates from grid_den.F '
c      call output(dbl_mb(ixyz), 1, 3, 1, nctrs, 3, nctrs, 1)
c     &                    
c
      if (.not.ma_pop_stack(ltags))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
c
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (do_nmrks_r) nbf_ao_mxnbf_ce2=nbf_ao_mxnbf_ce2*3
c
c     Divide the multi-center integrand into a set of single center integrands.
c
      nqtot = 0
      rho_int = 0.D0
c
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('grid_den: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('grid_den: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom list',lniz,iniz))
     &   call errquit('grid_den: cannot allocate atom list',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom nbf',lfin,ifin))
     &   call errquit('grid_den: cannot allocate atom nbf',0, MA_ERR)
      nexpo=nbf_ao_mxprim
      if(xcfit)nexpo=max(nbf_ao_mxprim,nbf_xc_mxprim)
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_den: cannot allocate expo',0, MA_ERR)
c
      ncoef_max = nbf_ao_mxnbf_cn * nbf_ao_mxprim
      if(xcfit) ncoef_max=max(ncoef_max,nbf_xc_mxnbf_cn * nbf_xc_mxprim)
c
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit('grid_den: cannot allocate coef',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',lFmat,iFmat))
     &   call errquit('grid_den: cannot allocate DENSF',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',lPmat,iPmat))
     &   call errquit('grid_den: cannot allocate DENS',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho',lrho,irho))
     &   call errquit('grid_den: cannot allocate rho',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,len_delrho,'delrho',ldelrho, idelrho))
     &   call errquit('grid_den: cannot allocate delrho',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl,len_ttau,'ttau',lttau, ittau))
     &   call errquit('grid_den: cannot allocate ttau',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*n_per_rec,'Bmat',
     &   lBmat,iBmat))call errquit('grid_den: cannot allocate Bmat',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,1,'hesdum',lhess,ihess))
     &   call errquit('grid_den: cannot allocate hesdum',0, MA_ERR)
c
c     Assume this numerical density routine can be called from anywhere.
c     -compute mapping vectors
c 
      if (.not.ma_push_get
     &   (mt_int,nctrs*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), nctrs, nshells_ao)
      if (.not.ma_chop_stack(lcntoce))
     &   call errquit('grid_den: cannot pop stack',1, MA_ERR)
c
c     -compute reduced total density matrices over atoms
c
      if (.not.MA_Push_Get(MT_Dbl,ipol*nctrs*nctrs,'rdens_atom',
     &   lrdens_atom,irdens_atom))
     &   call errquit('grid_den: cannot allocate rdens_atom',0, MA_ERR)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce*(nmat+1)
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit('grid_den: cannot allocate scr',0, MA_ERR)
c Daniel (1-11-13): I think we can stay similar to calc_type = 2 for
c XC-third derivatives in this portion of the code, since this part
c doesn't depend on having 3-dimensional GAs.
      if (calc_type.eq.2.or.calc_type.eq.5) then
         do isp=1,ipol
            call util_ga_mat_reduce(nbf_ao, nctrs, int_mb(icetobfr), 
     A           g_dens(1+(isp-1)*(nmat+1)),  1, 
     A           dbl_mb(irdens_atom+(isp-1)*nctrs*nctrs), 'rms', 
     &           dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
            do ii=2,(nmat+1)
               call util_ga_mat_reduce(nbf_ao,nctrs,int_mb(icetobfr), 
     A              g_dens(ii+(isp-1)*(nmat+1)),  1, 
     A              dbl_mb(irdens_atom+(isp-1)*nctrs*nctrs), 'rms', 
     &              dbl_mb(iscr), nbf_ao_mxnbf_ce,.false.)
            enddo
         enddo
      else
      call util_ga_mat_reduce(nbf_ao, nctrs, int_mb(icetobfr), g_dens, 
     &                        ipol, dbl_mb(irdens_atom), 'rms', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
      endif
c
      if (.not.ma_pop_stack(lscr))
     &     call errquit('grid_den: cannot pop stack',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,ipol*n_per_rec,'Amat',lamat,iamat))
     &   call errquit('grid_quadv0: cannot allocate Amat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Amat2,'Amat2',l_Amat2,k_Amat2))
     &   call errquit('grid_quadv0: cannot allocate Amat2',0, MA_ERR)
      call dfill(len_Amat2, 0.0d0, dbl_mb(k_Amat2), 1)
      if (.not.MA_Push_Get(mt_dbl,len_Amat3,'Amat3',l_Amat3,k_Amat3))
     1  call errquit('grid_quadv0: cannot allocate Amat3', 0, MA_ERR)
      call dfill(len_Amat3, 0.0d0, dbl_mb(k_Amat3), 1)
c
      if (.not.MA_Push_Get(
     .     mt_dbl,3*n_per_rec*ipol,'Cmat',lcmat,icmat))
     &   call errquit('grid_quadv0: cannot allocate Cmat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Cmat2,'Cmat2',l_Cmat2,k_Cmat2))
     &   call errquit('grid_quadv0: cannot allocate Cmat2',0, MA_ERR)
      call dfill(len_Cmat2, 0.0d0, dbl_mb(k_Cmat2), 1)
      if (.not.MA_Push_Get(mt_dbl,len_Cmat3,'Cmat3',l_Cmat3,k_Cmat3))
     1  call errquit('grid_quadv0: cannot allocate Cmat3', 0, MA_ERR)
      call dfill(len_Cmat3, 0.0d0, dbl_mb(k_Cmat3), 1)
c
      if (.not.MA_Push_Get(
     .     mt_dbl,len_mmat,'Mmat',lmmat,immat))
     &   call errquit('grid_quadv0: cannot allocate Mmat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Mmat2,'Mmat2',l_Mmat2,k_Mmat2))
     &   call errquit('grid_quadv0: cannot allocate Mmat2',0, MA_ERR)
      call dfill(len_Mmat2, 0.0d0, dbl_mb(k_Mmat2), 1)
c Daniel (1-11-13): Eventually we'll need to put stuff here for third order
c functional derivative arrays of the kinetic energy density, but I'm
c avoiding these because we don't yet have second derivatives for
c meta-GGAs (as far as I'm aware).
c
c     allocate scratch array for grid compression
c
      nibf=nbf_ao
      if(xcfit) nibf=max(nibf,nbf_xc)
      if (.not.MA_Push_Get(mt_int, nibf, 'ibf_ao', libf_ao,
     &     iibf_ao))
     &     call errquit('grid_quadv0: cannot allocate ibf_ao',2,
     &       MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 3,
     &       MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit('dft_scf: cannot allocate rchi_atom',0, MA_ERR)
      if(grad) then
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rdelchi_atom',
     &     lrdelchi_atom,irdelchi_atom))call errquit
     &        ('dft_scf: cannot allocate rdelchi_atom',0, MA_ERR)
      endif
      irqlen=nctrs*n_per_rec
      if(.not.grid_written) irqlen=max(irqlen,nctrs*nqmax)
      if (.not.MA_Push_get(mt_dbl,irqlen,'rq',lrq,irq))
     &   call errquit('grid_den: cannot allocate rq',0, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,irqlen,'bp',lbp,ibp))
     &     call errquit('grid_quadv0: cannot allocate bp',4, MA_ERR)
c
c     compute variables for scratch space in function evaluation
c     assumes cartesian; spherical set would only require 2*MXang + 1
c     
      max_component = (nbf_ao_mxang+1)*(nbf_ao_mxang+2)/2 
c     
c     Allocate scratch space.
c     
      if (nbf_ao_mxang .eq. 0)then
         max_elem_ang_scr = max_component * 3
      elseif (nbf_ao_mxang .le. 3)then
         max_elem_ang_scr = max_component * 9
      else                      ! general case
         max_elem_ang_scr = max_component * 28
      endif 
      nscr = 3*n_per_rec*nbf_ao_mxcont + 
     &     max(3*n_per_rec*nbf_ao_mxprim,
     &     n_per_rec*max_elem_ang_scr) + 1
      if (xcfit) nscr=max(nscr, 3*n_per_rec*nbf_xc_mxcont + 
     &     max(3*n_per_rec*nbf_xc_mxprim,
     &     n_per_rec*max_elem_ang_scr) + 1)
c     
      if (.not.MA_Push_Get(mt_dbl,nscr,'func scr',lscr,iscr))
     &     call errquit('grid_quadv0: cannot allocate bp',4, MA_ERR)
c
      Ex=0d0
      Ec=0d0
      call ga_sync()
c
      if(.not.do_cpks_l.or..not.do_gxc)
     *  call xc_rep_init(rtdb, geom, ao_bas_han,g_dens,g_vxc,
     &     nbf_ao,ipol,do_scf.and.(.not.l3d),.false.)
c      
c     Grid quadrature
c
      if(grid_written) then
c
         if (iAOacc .le. 0) then
          if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc))then
             iAOacc=-nint(log(e_conv))
          else
             iAOacc=max(iAOacc,-nint(log(e_conv)))
          endif
         end if
c
         call grid_quadv0a(g_dens,ncoef_max,nctrs,rho_int,nqtot, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     .        dbl_mb(iscr),nscr,
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        dbl_mb(irdens_atom), int_mb(icetobfr),
     ,        int_mb(iibf_ao), dbl_mb(irqbuf),
     ,        dbl_mb(irchi_atom), dbl_mb(irdelchi_atom),
     ,        ex,ec,dbl_mb(iamat), dbl_mb(icmat),g_vxc,
     ,        dbl_mb(k_amat2), dbl_mb(k_cmat2),
     ,        nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, tmat,
     ,        l3d,triplet,
     ,        dbl_mb(ittau), 
     &        dbl_mb(immat), dbl_mb(k_mmat2),StericEnergyTotal,
     1        do_gxc, dbl_mb(k_Amat3), dbl_mb(k_Cmat3))
      else
c
         if (.not. rtdb_get(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &        call errquit('gridv0: rtdb_get failed', 119, RTDB_ERR)
         nqlen = ncenters*mxnrad
         if (.not.MA_Push_get(MT_int,4*nqlen,'shell list',
     &        lqs,iqs))
     &        call errquit('gridv0: cannot allocate shell list',0,
     &       MA_ERR)
         if (.not.MA_Push_get(MT_dbl, 4*nqmax,
     &        'ang quad pts buffer', laqbuf, iaqbuf))
     &        call errquit('gridv0: cant allocate quad pt', 0, MA_ERR)
         ndbl = (ncenters*(ncenters+1))/2
         if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lRij,iRij))
     &        call errquit('grid_v0: cannot allocate Rij',0, MA_ERR)
         if (.not.MA_Push_get(mt_int,nqshells,'qshat',l_qat,k_qat))
     &        call errquit('gridv0: cannot allocate qat',5, MA_ERR)
c
c Daniel (1-11-13): Added XC-third derivative stuff here also, although
c it should never be used when we need to construct the grid.  This is
c needed for preventing odd behavior in calls to the grid_quadv0b and
c xc_eval_fnl routine later on, when someone implements new functionality
c down the line for the DFT code.
         call grid_quadvw(rtdb,
     ,        int_mb(iqs),dbl_mb(iaqbuf),
     ,        dbl_mb(irij),int_mb(k_qat),
     ,        g_dens,ncoef_max,nctrs,rho_int,nqtot, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     .        dbl_mb(iscr),nscr,
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        dbl_mb(irdens_atom), int_mb(icetobfr),
     ,        int_mb(iibf_ao),dbl_mb(ibp),
     ,        dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),
     ,        ex,ec,dbl_mb(iamat),dbl_mb(icmat),g_vxc,
     ,        dbl_mb(k_amat2),dbl_mb(k_cmat2),
     ,        nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, tmat,
     ,        l3d,triplet,
     ,        dbl_mb(ittau), 
     &        dbl_mb(immat), dbl_mb(k_mmat2),
     1        do_gxc, dbl_mb(k_Amat3), dbl_mb(k_Cmat3),
     2        0, 0.0d0, 0.0d0)
c
         if (.not. MA_chop_Stack(lqs))
     &        call errquit('gridv0: pop stack failed.',1, MA_ERR)
c Qin
         if (.not. rtdb_get(rtdb, 'dft:grid_make_sa', mt_log, 1, lsa))
     &     lsa = .false.
         if(lsa) then
           call grid_make_sa(rtdb)
           lsa = .false.
           if (.not. rtdb_put(rtdb, 'dft:grid_make_sa', mt_log, 1, lsa))
     &       call errquit('rtdb_put lsa failed', 0, RTDB_ERR)
         endif
c end Qin

      endif  ! grid_written
c     
      if(.not.xc_rep_close(rtdb, nbf_ao,ipol,g_vxc,.false.)) call
     .     errquit(' gridquadv0: xcrepclose failed ',0, DISK_ERR)
!      call ga_print(g_vxc(1))
      
      if (.not.ma_chop_stack(lqxyz))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
      call ga_sync()
      if(do_scf.and.(.not.l3d)) then
c
c         write(0,*) ' symmetr scf l r',do_scf,do_cpks_l,do_cpks_r,
c     A        calc_type
c         call ga_print(g_vxc(1))
c
          call ga_symmetrize(g_vxc(1))
          if(ipol.gt.1) call ga_symmetrize(g_vxc(2))
c
c         call ga_print(g_vxc(1))
      endif
c
      call ga_igop(Msg_nqtot,nqtot,1,'+')
      call ga_dgop(Msg_Excrho,dummy,3,'+')
      if(xcfit) call ga_dgop(Msg_Tmat,Tmat,nbf_xc*ipol,'+')
      if (nExc.eq.1)then
         Exc(1) = Exc(1) + Ec + Ex
      else
         Exc(1) = Exc(1) + Ex
         Exc(2) = Exc(2) + Ec
      endif
c
c     == Print total steric energy ==
c
      if (lsteric) then
         if (me.eq.0) write(LuOut,2221) StericEnergyTotal
      end if
 2221  format('         Total Steric Energy =', f22.12/)
c
c     get grid accuracy for job
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('grid_acc_def: rtdb_put failed', 911, RTDB_ERR)
c
      fon = .false.
c     fractional occupations?
c
      if (rtdb_get(rtdb,'dft:fon',mt_log,1,fon)) then 
        if (.not.fon) call errquit(
     &     'grid_quadv0: fon stored in RTDB but not .true.', 1,
     &     RTDB_ERR)
        if (.not. rtdb_get(rtdb, 'dft:pstrace', mt_dbl, 1, pstrace))
     &     call errquit('grid_quadv0: rtdb_get pstrace failed',
     &     1, RTDB_ERR)
      else
        fon = .false.
      end if
c
c     get number of electrons
c
      if(ipol.eq.1) then
         nelec=dble(2*noc(1))
      else
         nelec=dble(noc(1)+noc(2))
      endif
c
c     for fractional occupations the correct electron number comes from tr[PS]
      if (fon) nelec = pstrace
c         
      if(oprintw) then
         if(sic_orb_occ.eq.0) then
            rho_err=abs(rho_int-nelec)/nelec
         else
            rho_err=abs(rho_int-1d0)
         endif
c Daniel (1-11-13): I changed the logic here to stop the warning from
c being printed for either the TDDFT solver or TDDFT gradient solver.
         lcpksl_gxc = (do_cpks_l.or.do_gxc)
         if(.not.lcpksl_gxc.and.rho_err.gt.job_grid_acc.and.me.eq.0)
     &      write(LuOut,123) rho_int,job_grid_acc
c         if(.not.do_cpks_l.and.rho_err.gt.job_grid_acc.and.me.eq.0)
c     &      write(LuOut,123) rho_int,job_grid_acc
      endif ! oprintw

c     some debug output to check the density integration error:
c$$$      if (me.eq.0) then
c$$$        write(luout,'(1x,a,t15,i10)') "noc(1):",noc(1)
c$$$        write(luout,'(1x,a,t15,i10)') "noc(2):",noc(2)
c$$$        write(luout,'(1x,a,t15,i10)') "nelec:",nelec
c$$$        write(luout,'(1x,a,t15,f15.8)') "rho_int:",rho_int
c$$$        write(luout,'(1x,a,t15,f15.8)')
c$$$     &     "rhodiff: ",abs(rho_int-dble(nelec))
c$$$      end if
c
      rho_n = rho_int  ! integrated density
c
      if (me.eq.0.and.oprint)then
         write(LuOut,*)
         write(LuOut,'(a,f14.6,a,i8,a,i8)')
     &        ' XC quadr. integrated density : ',
     &        rho_int,' npts: ',nqtot,' nctrs: ',nctrs
        write(LuOut,*)
      endif
c
      if (.not.MA_Free_Heap(lxyz))
     &   call errquit('grid_den: cannot free heap',0, MA_ERR)
c
 123  format(
     &   1x,'Grid integrated density: ',f20.12,/
     &   1x,'Requested integration accuracy: ',e10.2)
c
      return
      end
