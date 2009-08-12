      subroutine grid_quadv0(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,tmat)
c
C$Id: grid_quadv0.F,v 1.56 2008-02-01 22:27:45 niri Exp $
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
#include "zora.fh"
#include "steric.fh"

      logical rtdb_status       
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
     &        ltau, itau,  lmmat, immat, lttau, ittau
c
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr,lamat,iamat,icmat,lcmat,
     ,        noc(2),nelec
c
      integer libf_ao,iibf_ao,irqbuf,lrqbuf,lbp,ibp,
     .     irchi_atom,irdelchi_atom,lrchi_atom,lrdelchi_atom
      integer iqs,lqs,nqlen,mxnrad,
     , laqbuf, iaqbuf,lRij,iRij,ndbl,l_qat,k_qat,irqlen
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2, l_Mmat2, k_Mmat2
      integer len_rho, len_delrho, len_Amat2, len_Cmat2, len_Mmat2
      integer len_tau, len_mmat, len_ttau
      logical do_scf, do_cpks_l, do_cpks_r, do_nmrks_r
      integer max_elem_ang_scr, max_component
      integer nexpo,nibf
      double precision rho_int, ex, ec,dummy(3),job_grid_acc,rho_err
      integer ii,isp
      equivalence (dummy(1),Ex), (dummy(2),Ec), (dummy(3),rho_int)

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
      npol = (ipol*(ipol+1))/2
      if (.not.rtdb_get(rtdb,'dft:largenode', mt_log, 1, largenode))
     &     largenode=.false.
c
c     Open grid pts file
c     
      if (.not. grid_written .or. n_per_rec .le. 0) then
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
      if (.not. (do_scf .or. do_cpks_l .or. do_cpks_r .or. 
     &    do_nmrks_r)) then
         call errquit('grid_quadv0_gen: illegal calculation type',1,
     &       INPUT_ERR)
      endif
c
      if (oprint) then
        write(*,*) " *** grid_quadv0: do_scf:",do_scf 
        write(*,*) " *** grid_quadv0: n_per_rec:",n_per_rec
      end if
c
      len_tau = 1
      len_mmat = 1
      len_ttau = 1
      len_Amat2 = 1
      len_Cmat2 = 1
      len_Mmat2 = 1
c
      if (do_scf .or. do_nmrks_r) then
         len_rho = npol*n_per_rec
         len_delrho = 3*ipol*n_per_rec
         if (kske) len_tau = 3*ipol*n_per_rec
         if (kske) len_ttau = npol*n_per_rec
         if (kske) len_mmat =ipol*n_per_rec  
      elseif (do_cpks_l) then
         len_rho = ipol*n_per_rec*(nmat+1)+n_per_rec
         len_delrho = 3*ipol*n_per_rec*(nmat+1)+3*n_per_rec
         if (kske) len_tau = 3*ipol*n_per_rec*(nmat+1)+3*n_per_rec
         if (kske) len_ttau = ipol*n_per_rec*(nmat+1)+n_per_rec
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
         if (kske) len_Mmat2 = n_per_rec*NCOL_MMAT2
      elseif (do_cpks_r) then
         call errquit('grid_quadv0: what for do_cpks_r?',0,
     &       UNKNOWN_ERR)
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
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
      if (.not.MA_Push_Get(MT_Dbl,len_tau,'tau',ltau, itau))
     &    call errquit('grid_den: cannot allocate tau',0, MA_ERR)
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
      if(calc_type.eq.2) then
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

      if (.not.ma_pop_stack(lscr))
     &     call errquit('grid_den: cannot pop stack',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,ipol*n_per_rec,'Amat',lamat,iamat))
     &   call errquit('grid_quadv0: cannot allocate Amat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Amat2,'Amat2',l_Amat2,k_Amat2))
     &   call errquit('grid_quadv0: cannot allocate Amat2',0, MA_ERR)
      call dfill(len_Amat2, 0.0d0, dbl_mb(k_Amat2), 1)
c
      if (.not.MA_Push_Get(
     .     mt_dbl,3*n_per_rec*ipol,'Cmat',lcmat,icmat))
     &   call errquit('grid_quadv0: cannot allocate Cmat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Cmat2,'Cmat2',l_Cmat2,k_Cmat2))
     &   call errquit('grid_quadv0: cannot allocate Cmat2',0, MA_ERR)
      call dfill(len_Cmat2, 0.0d0, dbl_mb(k_Cmat2), 1)
c
      if (.not.MA_Push_Get(
     .     mt_dbl,len_mmat,'Mmat',lmmat,immat))
     &   call errquit('grid_quadv0: cannot allocate Mmat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Mmat2,'Mmat2',l_Mmat2,k_Mmat2))
     &   call errquit('grid_quadv0: cannot allocate Mmat2',0, MA_ERR)
      call dfill(len_Mmat2, 0.0d0, dbl_mb(k_Mmat2), 1)
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
     
      if (.not.MA_Push_Get(mt_dbl,nscr,'func scr',lscr,iscr))
     &     call errquit('grid_quadv0: cannot allocate bp',4, MA_ERR)

      Ex=0d0
      Ec=0d0
      call ga_sync()
      if(.not.do_cpks_l)
     *  call xc_rep_init(rtdb, geom, ao_bas_han,g_dens,g_vxc,
     &     nbf_ao,ipol,do_scf.and.(.not.l3d),.false.)
c      
      if(grid_written) then

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
     ,        dbl_mb(itau), dbl_mb(ittau), 
     &        dbl_mb(immat), dbl_mb(k_mmat2),StericEnergyTotal)
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
     ,        dbl_mb(itau), dbl_mb(ittau), 
     &        dbl_mb(immat), dbl_mb(k_mmat2))
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

      endif
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
czora...start
czora...only symmetrize the zora spin-free term 
          if (do_zora) then
            if (so_term.eq.0) then
              call ga_symmetrize(g_vxc(1))
              if(ipol.gt.1) call ga_symmetrize(g_vxc(2))
            end if
            call ga_sync()
            call util_flush(LuOut)
          else 
            call ga_symmetrize(g_vxc(1))
            if(ipol.gt.1) call ga_symmetrize(g_vxc(2))
          end if
czora...end
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
c
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('grid_acc_def: rtdb_put failed', 911, RTDB_ERR)
c
c     get N electr
c
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('dft_scf: rtdb_get of noc failed', 0, RTDB_ERR)
      if(ipol.eq.1) then
         nelec=2*noc(1)
      else
         nelec=noc(1)+noc(2)
      endif
         
      if(oprintw) then
         if(sic_orb_occ.eq.0) then
            rho_err=abs(rho_int-dble(nelec))/nelec
         else
            rho_err=abs(rho_int-1d0)
         endif
         if (me.eq.0) then
           if(.not.do_cpks_l.and.rho_err.gt.job_grid_acc.and.me.eq.0)
     &        write(LuOut,123) rho_err,job_grid_acc
         end if
      endif

      rho_n = rho_int
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
 123  format(' WARNING: error on integrated density =',
     , d10.2,/' greater than required accuracy of',d10.2)

      return
      end
