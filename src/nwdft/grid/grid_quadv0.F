      Subroutine grid_quadv0(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,
     ,     tmat)
c
C$Id: grid_quadv0.F,v 1.20 2001-01-16 01:59:47 edo Exp $
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(4)
      double precision rho_n
      double precision Exc(2)
      double precision tmat(*)  ! XC fit
      call grid_quadv0_gen(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,
     ,     1,1,tmat)
      return
      end
c
      Subroutine grid_quadv0_gen(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,
     ,      nmat,  calc_type ,tmat)
      implicit none
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(4)
      Logical oprint,grad,xc_chkgrad,oprintw
      external xc_chkgrad
      double precision Exc(2)
      double precision rho_n
      integer nmat              ! [input] number of XC matrices to make;
      double precision tmat(*)  ! XC fit
c                               !         only used for lhs cpks
      integer calc_type         ! what to calculate [input]
c                               !   1 = XC matrices for SCF
c                               !   2 = XC matrices for CPKS LHS
c                               !   3 = XC matrices for CPKS RHS
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
c !!! BGJ test
#include "bgj.fh"
c !!! BGJ test
c
      integer me, nproc, nxyz, nctrs, 
     &        nbf_ao_mxnbf_ce2, nqtot, 
     &        ncoef_max,
     &        npol
      integer lxyz, ixyz, lcharge, icharge, 
     &        lqxyz, iqxyz, lqwght, iqwght, 
     &        lniz, iniz, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, ldelrho, idelrho, lBmat, iBmat, 
     &        lhess, ihess, lrq, irq, 
     &        lfin, ifin, ltags, itags
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr,lamat,iamat,icmat,lcmat,
     ,     noc(2),nelec
      integer libf_ao,iibf_ao,irqbuf,lrqbuf,lbp,ibp,lrscr,irscr,
     .     irchi_atom,irdelchi_atom,lrchi_atom,lrdelchi_atom
      integer iqs,lqs,nqlen,mxnrad,nq_list,
     , laqbuf, iaqbuf,lRij,iRij,ndbl,l_qat,k_qat,irqlen
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2
      integer len_rho, len_delrho, len_Amat2, len_Cmat2
      logical do_scf, do_cpks_l, do_cpks_r
      logical oldway
      logical donewbf           !!! !!!! *** temp
      integer max_elem_ang_scr, max_component
      integer nexpo,nibf
      double precision rho_int, ex, ec,dummy(3),job_grid_acc,rho_err
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
c
      npol = (ipol*(ipol+1))/2
      if (.not. rtdb_get(rtdb, 'dft:oldway', mt_log, 1, oldway))
     &     oldway=.false.
c
c     Open grid pts file
c     
      if(.not.grid_written) call grid_packing_info
c
c
c     Set up logicals denoting calculation type
c
      do_scf = calc_type .eq. 1
      do_cpks_l = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      if (.not. (do_scf .or. do_cpks_l .or. do_cpks_r)) then
         call errquit('grid_quadv0_gen: illegal calculation type',1)
      endif
      if (bgj_print() .gt. 0)
     &     write(*,*)'*** grid_quadv0: nmat logicals', nmat, 
     &     do_scf, do_cpks_l, do_cpks_r
      if (do_scf) then
         len_rho = npol*n_per_rec
         len_delrho = 3*ipol*n_per_rec
         len_Amat2 = 1
         len_Cmat2 = 1
      elseif (do_cpks_l) then
         len_rho = ipol*n_per_rec*(nmat+1)
         len_delrho = 3*ipol*n_per_rec*(nmat+1)
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
      elseif (do_cpks_r) then
         call errquit('grid_quadv0: what for do_cpks_r?',0)
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
      endif
      if (.not. geom_ncent(geom, nctrs))
     &   call errquit('grid_den: geom_ncent failed',73)     
      nxyz = 3*nctrs
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('grid_den: cannot allocate xyz',0)
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'charge',lcharge,icharge))
     &   call errquit('grid_den: cannot allocate charge',0)
      if (.not.MA_Push_Get(MT_Byte,nctrs*16,'tags',ltags,itags))
     &   call errquit('grid_den: cannot allocate tags',0)
      if (.not. geom_cart_get(geom, nctrs, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('grid_den: geom_cart_get failed',74)
c
c      write(*,*)' Coordinates from grid_den.F '
c      call output(dbl_mb(ixyz), 1, 3, 1, nctrs, 3, nctrs, 1)
c     &                    
c
      if (.not.ma_pop_stack(ltags))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('grid_den: cannot pop stack',0)
c
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
c
c     Divide the multi-center integrand into a set of single center integrands.
c
      nqtot = 0
      rho_int = 0.D0
c
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('grid_den: cannot allocate qxyz',0)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('grid_den: cannot allocate qwght',0)
      if (.not.MA_Push_get(MT_int,nctrs,'atom list',lniz,iniz))
     &   call errquit('grid_den: cannot allocate atom list',0)
      if (.not.MA_Push_get(MT_int,nctrs,'atom nbf',lfin,ifin))
     &   call errquit('grid_den: cannot allocate atom nbf',0)
      nexpo=nbf_ao_mxprim
      if(xcfit)nexpo=max(nbf_ao_mxprim,nbf_xc_mxprim)
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_den: cannot allocate expo',0)
c
      ncoef_max = nbf_ao_mxnbf_cn * nbf_ao_mxprim
      if(xcfit) ncoef_max=max(ncoef_max,
     ,     nbf_xc_mxnbf_cn * nbf_xc_mxprim)
c
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit('grid_den: cannot allocate coef',0)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',lFmat,iFmat))
     &   call errquit('grid_den: cannot allocate DENSF',0)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',lPmat,iPmat))
     &   call errquit('grid_den: cannot allocate DENS',0)
      if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho',lrho,irho))
     &   call errquit('grid_den: cannot allocate rho',0)
      if (.not.MA_Push_Get(MT_Dbl,len_delrho,'delrho',ldelrho,
     &   idelrho))call errquit('grid_den: cannot allocate delrho',0)
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*n_per_rec,'Bmat',
     &   lBmat,iBmat))call errquit('grid_den: cannot allocate Bmat',0)
      if (.not.MA_Push_Get(MT_dbl,1,'hesdum',lhess,ihess))
     &   call errquit('grid_den: cannot allocate hesdum',0)
c
c     Assume this numerical density routine can be called from anywhere.
c     -compute mapping vectors
c 
      if (.not.ma_push_get
     &   (mt_int,nctrs*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit('grid_den:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('grid_den:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('grid_den:push_get failed', 13)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), nctrs, nshells_ao)
      if (.not.ma_chop_stack(lcntoce))
     &   call errquit('grid_den: cannot pop stack',1)
c
c     -compute reduced total density matrices over atoms
c
      if (.not.MA_Push_Get(MT_Dbl,ipol*nctrs*nctrs,'rdens_atom',
     &   lrdens_atom,irdens_atom))
     &   call errquit('grid_den: cannot allocate rdens_atom',0)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit('grid_den: cannot allocate scr',0)
      call util_ga_mat_reduce(nbf_ao, nctrs, int_mb(icetobfr), g_dens, 
     &                        ipol, dbl_mb(irdens_atom), 'rms', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce)
      if (.not.ma_pop_stack(lscr))
     &   call errquit('grid_den: cannot pop stack',0)
c
      if (.not.MA_Push_Get(mt_dbl,ipol*n_per_rec,'Amat',lamat,iamat))
     &   call errquit('grid_quadv0: cannot allocate Amat',0)
      if (.not.MA_Push_Get(mt_dbl,len_Amat2,'Amat2',l_Amat2,k_Amat2))
     &   call errquit('grid_quadv0: cannot allocate Amat2',0)
      call dfill(len_Amat2, 0.0d0, dbl_mb(k_Amat2), 1)
      if (.not.MA_Push_Get(
     .     mt_dbl,3*n_per_rec*ipol,'Cmat',lcmat,icmat))
     &   call errquit('grid_quadv0: cannot allocate Cmat',0)
      if (.not.MA_Push_Get(mt_dbl,len_Cmat2,'Cmat2',l_Cmat2,k_Cmat2))
     &   call errquit('grid_quadv0: cannot allocate Cmat2',0)
      call dfill(len_Cmat2, 0.0d0, dbl_mb(k_Cmat2), 1)
c
c     allocate scratch array for grid compression
c
      nibf=nbf_ao
      if(xcfit) nibf=max(nibf,nbf_xc)
      if (.not.MA_Push_Get(mt_int, nibf, 'ibf_ao', libf_ao,
     &     iibf_ao))
     &     call errquit('grid_quadv0a: cannot allocate ibf_ao',2)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 3)
      if (.not.MA_Push_Get(mt_dbl,n_per_rec,
     ,     'rho scr',lrscr,irscr))
     &     call errquit('grid_quadv0a: cannot allocate rho scr',5)

      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit('dft_scf: cannot allocate rchi_atom',0)
      if(grad) then
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rdelchi_atom',
     &     lrdelchi_atom,irdelchi_atom))call errquit
     &        ('dft_scf: cannot allocate rdelchi_atom',0)
      endif
      irqlen=nctrs*n_per_rec
      if(.not.grid_written) irqlen=max(irqlen,nctrs*nqmax)
      if (.not.MA_Push_get(mt_dbl,irqlen,'rq',lrq,irq))
     &   call errquit('grid_den: cannot allocate rq',0)
      if (.not.MA_Push_Get(mt_dbl,irqlen,'bp',lbp,ibp))
     &     call errquit('grid_quadv0a: cannot allocate bp',4)
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
     &     call errquit('grid_quadv0a: cannot allocate bp',4)

      Ex=0d0
      Ec=0d0
      call ga_sync()
      if (.not. rtdb_get(rtdb, 'bgj:new_bf', mt_log,
     &     1, donewbf)) donewbf = .false.
      if(grid_written) then
         call grid_quadv0a(g_dens,ncoef_max,nctrs,rho_int,nqtot, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     .        dbl_mb(iscr),nscr,
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        dbl_mb(irdens_atom), int_mb(icetobfr),
     ,        int_mb(iibf_ao),dbl_mb(irqbuf),dbl_mb(ibp),
     ,        dbl_mb(irscr),
     ,        dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),
     ,        ex,ec,dbl_mb(iamat),dbl_mb(icmat),g_vxc,
     ,        dbl_mb(k_amat2),dbl_mb(k_cmat2),
     ,     nmat, do_scf, do_cpks_l, do_cpks_r, tmat,donewbf)
      else
         if (.not. rtdb_get(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &        call errquit('gridv0: rtdb_get failed', 119)
         nqlen = ncenters*mxnrad
         if (.not.MA_Push_get(MT_int,4*nqlen,'shell list',
     &        lqs,iqs))
     &        call errquit('gridv0: cannot allocate shell list',0)
c
c     get list of grid pts
c     
         call grid_list(rtdb, int_mb(iqs), nqlen, nq_list)
         if (.not.MA_Push_get(MT_dbl, 4*nqmax,
     &        'ang quad pts buffer', laqbuf, iaqbuf))
     &        call errquit('gridv0: cant allocate quad pt', 0)
         ndbl = (ncenters*(ncenters+1))/2
         if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lRij,iRij))
     &        call errquit('grid_v0: cannot allocate Rij',0)
         if (.not.MA_Push_get(mt_int,nqshells,'qshat',l_qat,k_qat))
     &        call errquit('gridv0: cannot allocate qat',5)
         if(oldway) then
         call grid_quadvold(rtdb,
     ,     int_mb(iqs),nq_list,dbl_mb(iaqbuf),
     ,        dbl_mb(irij),int_mb(k_qat),
     ,     g_dens,ncoef_max,nctrs,rho_int,nqtot, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     .        dbl_mb(iscr),nscr,
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        dbl_mb(irdens_atom), int_mb(icetobfr),
     ,        int_mb(iibf_ao),dbl_mb(ibp),
     ,        dbl_mb(irscr),
     ,        dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),
     ,        ex,ec,dbl_mb(iamat),dbl_mb(icmat),g_vxc,
     ,        dbl_mb(k_amat2),dbl_mb(k_cmat2),
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat)
         else
         call grid_quadvw(rtdb,
     ,     int_mb(iqs),nq_list,dbl_mb(iaqbuf),
     ,        dbl_mb(irij),int_mb(k_qat),
     ,     g_dens,ncoef_max,nctrs,rho_int,nqtot, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     .        dbl_mb(iscr),nscr,
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        dbl_mb(irdens_atom), int_mb(icetobfr),
     ,        int_mb(iibf_ao),dbl_mb(ibp),
     ,        dbl_mb(irscr),
     ,        dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),
     ,        ex,ec,dbl_mb(iamat),dbl_mb(icmat),g_vxc,
     ,        dbl_mb(k_amat2),dbl_mb(k_cmat2),
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat,donewbf)
         endif
         if (.not. MA_chop_Stack(lqs))
     &        call errquit('gridv0: pop stack failed.',1)
      endif
c     
      if (.not.ma_chop_stack(lqxyz))
     &   call errquit('grid_den: cannot pop stack',0)
      call ga_sync()
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
      if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &   job_grid_acc))
     &   call errquit('grid_acc_def: rtdb_put failed', 911)
c
c     get N electr
c
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('dft_scf: rtdb_get of noc failed', 0)
      if(ipol.eq.1) then
         nelec=2*noc(1)
      else
         nelec=noc(1)+noc(2)
      endif
         
      if(oprintw) then
         rho_err=abs(rho_int-dble(nelec))/nelec
         if(.not.do_cpks_l.and.rho_err.gt.job_grid_acc.and.me.eq.0)
     *        write(LuOut,123) rho_err,job_grid_acc
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
     &   call errquit('grid_den: cannot free heap',0)
c
 123  format(' WARNING: error on integrated density =',
     , d10.2,/' greater than required accuracy of',d10.2)
      return
      end
