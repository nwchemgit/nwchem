      Subroutine grid_quadv0(rtdb, g_dens, g_vxc, nexc, rho_n, Exc)
c
C$Id: grid_quadv0.F,v 1.1 1999-06-25 21:48:24 d3h449 Exp $
c
      implicit none
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(4)
      Logical oprint
      double precision Exc(2)
c
#include "cgridfile.fh"
#include "geom.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "util.fh"
c
      double precision rho_n
      integer me, nproc, icount, nxyz, nctrs, 
     &        nbf_ao_mxnbf_ce2, nqtot, 
     &        ncoef_max,
     &        npol
      integer lxyz, ixyz, lcharge, icharge, 
     &        lqxyz, iqxyz, lqwght, iqwght, lrq0, irq0, 
     &        lniz, iniz, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, ldelrho, idelrho, lBmat, iBmat, 
     &        lhess, ihess, lrq, irq, 
     &        lfin, ifin, ltags, itags
      integer len_rho, len_delrho
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr,lamat,iamat,icmat,lcmat
      double precision rho_int, ex, ec
#ifdef SOLARIS
c     needed for -stackvar option
      save rho_int
#endif      
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('quadrature', print_high)
c
      icount=1
c
      npol = (ipol*(ipol+1))/2
c
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
      len_rho = npol*n_per_rec
      len_delrho = 3*ipol*n_per_rec
c
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('grid_den: cannot allocate qxyz',0)
      call dfill(3*n_per_rec, 0.0d0, dbl_mb(iqxyz), 1)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('grid_den: cannot allocate qwght',0)
      call dfill(n_per_rec, 0.0d0, dbl_mb(iqwght), 1)
      if (.not.MA_Push_get(mt_dbl,nctrs,'rq0',lrq0,irq0))
     &   call errquit('grid_den: cannot allocate rq0',0)
      call dfill(nctrs, 0.0d0, dbl_mb(irq0), 1)
      if (.not.MA_Push_get(mt_dbl,nctrs*n_per_rec,'rq',lrq,irq))
     &   call errquit('grid_den: cannot allocate rq',0)
      call dfill(nctrs*n_per_rec, 0.0d0, dbl_mb(irq), 1)
      if (.not.MA_Push_get(MT_int,nctrs,'atom list',lniz,iniz))
     &   call errquit('grid_den: cannot allocate atom list',0)
      call ifill(nctrs, 0.0d0, int_mb(iniz), 1)
      if (.not.MA_Push_get(MT_int,nctrs,'atom nbf',lfin,ifin))
     &   call errquit('grid_den: cannot allocate atom nbf',0)
      call ifill(nctrs, 0.0d0, int_mb(ifin), 1)
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_den: cannot allocate expo',0)
      call dfill(nbf_ao_mxprim, 0.0d0, dbl_mb(iexpo), 1)
c
      ncoef_max = nbf_ao_mxnbf_cn * nbf_ao_mxprim
c
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit('grid_den: cannot allocate coef',0)
      call dfill(ncoef_max, 0.0d0, dbl_mb(icoef), 1)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',lFmat,iFmat))
     &   call errquit('grid_den: cannot allocate DENSF',0)
      call dfill(nbf_ao_mxnbf_ce2, 0.0d0, dbl_mb(iFmat), 1)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',lPmat,iPmat))
     &   call errquit('grid_den: cannot allocate DENS',0)
      call dfill(nbf_ao_mxnbf_ce2, 0.0d0, dbl_mb(iPmat), 1)
      if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho',lrho,irho))
     &   call errquit('grid_den: cannot allocate rho',0)
      call dfill(len_rho, 0.0d0, dbl_mb(irho), 1)
      if (.not.MA_Push_Get(MT_Dbl,len_delrho,'delrho',ldelrho,
     &   idelrho))call errquit('grid_den: cannot allocate delrho',0)
      call dfill(len_delrho, 0.0d0, dbl_mb(idelrho), 1)
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*n_per_rec,'Bmat',
     &   lBmat,iBmat))call errquit('grid_den: cannot allocate Bmat',0)
      call dfill(nbf_ao_mxnbf_ce*n_per_rec, 0.0d0, dbl_mb(iBmat), 1)
      if (.not.MA_Push_Get(MT_dbl,1,'hesdum',lhess,ihess))
     &   call errquit('grid_den: cannot allocate hesdum',0)
c
c     Assume this numerical density routine can be called from anywhere.
c
c     -compute mapping vectors
c 
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('grid_den:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('grid_den:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,nctrs*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit('grid_den:push_get failed', 13)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), nctrs, nshells_ao)
c
c     -compute reduced total density matrices over atoms
c
      if (.not.MA_Push_Get(MT_Dbl,ipol*nctrs*nctrs,'rdens_atom',
     &   lrdens_atom,irdens_atom))
     &   call errquit('grid_den: cannot allocate rdens_atom',0)
      call dfill(ipol*nctrs*nctrs, 0.0d0, dbl_mb(irdens_atom), 1)
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
     &   call errquit('xc_quadv0: cannot allocate Amat',0)
      call dfill(ipol*n_per_rec, 0.0d0, dbl_mb(iAmat), 1)
      if (.not.MA_Push_Get(
     .     mt_dbl,3*n_per_rec*ipol,'Cmat',lcmat,icmat))
     &   call errquit('xc_quadv0: cannot allocate Cmat',0)
      call dfill(3*ipol*n_per_rec, 0.0d0, dbl_mb(iCmat), 1)
      Ex=0d0
      Ec=0d0
      call grid_quadv0a(rtdb, g_dens, ncoef_max, nctrs, rho_int, nqtot, 
     &                int_mb(iniz), int_mb(ifin), 
     &                dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &                dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     &                dbl_mb(iexpo), dbl_mb(icoef), 
     &                dbl_mb(iBmat), dbl_mb(iFmat), 
     &                dbl_mb(iPmat), dbl_mb(irq), dbl_mb(irq0), 
     &                dbl_mb(irdens_atom), int_mb(icntoce),
     &                int_mb(icntobfr), int_mb(icetobfr),
     ,     ex,ec,dbl_mb(iamat),dbl_mb(icmat),g_vxc)
      if (nExc.eq.1)then
         Exc(1) = Exc(1) + Ec + Ex
      else
         Exc(1) = Exc(1) + Ex
         Exc(2) = Exc(2) + Ec
      endif
c
      if (.not.ma_pop_stack(lcmat))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lamat))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lrdens_atom))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lcetobfr))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lcntobfr))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lcntoce))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lhess))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lBmat))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(ldelrho))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lrho))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lPmat))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lFmat))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lcoef))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lexpo))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lfin))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lniz))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lrq))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lrq0))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lqwght))
     &   call errquit('grid_den: cannot pop stack',0)
      if (.not.ma_pop_stack(lqxyz))
     &   call errquit('grid_den: cannot pop stack',0)
c
      call ga_igop(Msg_nqtot,nqtot,1,'+')
      call ga_dgop(Msg_Excrho,rho_int,1,'+')
c
      rho_n = rho_int
c
      if (me.eq.0.and.oprint)then
         write(LuOut,*)
         write(LuOut,'(a,f14.6,a,i8,a,i8)')
     &        ' XC quadr. integrated density (from den): ',
     &        rho_int,' npts: ',nqtot,' nctrs: ',nctrs
        write(LuOut,*)
C        write(LuOut,'(/a)')'         Num. Quad. of the Electron Density'
c        write(LuOut,'(a,f20.10,2x,f20.10)')
c     &        ' integrated electron density:  ',rho_int, rho_n
C        write(LuOut,'(a,i20)')   ' no. of quad. points/center :  ',
C    &        nqtot/nctrs
C        write(LuOut,*)
      endif
c
      if (.not.MA_Free_Heap(lxyz))
     &   call errquit('grid_den: cannot free heap',0)
c
      return
      end
