      Subroutine grid_quadv0a(rtdb, g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess, qwght, qxyz, xyz,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq, rq0, rdens_atom, cntoce,
     &                       cntobfr, cetobfr, Ex, Ec, amat,cmat,g_vxc)
c
c$Id: grid_quadv0a.F,v 1.12 2000-03-21 23:11:59 edo Exp $
c
      implicit none
c
      integer rtdb
      integer g_dens(2),g_vxc(4)
      integer ncoef_max, nqtot
      integer natoms
      double precision Ex, Ec,Amat(*),Cmat(*)
c

#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer cntoce(nshells_ao), cntobfr(2,nshells_ao), 
     &        cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      double precision rho(*), delrho(*), hess(1), 
     &                 Bmat(nbf_ao_mxnbf_ce*nqmax)
      double precision Pmat(nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce), 
     &                 Fmat(nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce)
      double precision qwght(nqmax), qxyz(3,nqmax)
      double precision xyz(3,natoms)
      integer iniz(natoms), ifin(natoms)
      integer max_elem_ang_scr, max_component
      double precision ccoef(ncoef_max), expo(nbf_ao_mxprim)
      double precision rq(nqmax, natoms)
      double precision rq0(natoms)
      double precision  rho_int
      double precision acc_AO_gauss, rad, dummy(1)
      integer nbf_ao_mxnbf_ce2, 
     &        mbf_ao, 
     &        npol, iqsh, nscr, me
      integer ibp, lbp, iibf_ao, libf_ao, lchi_ao, ichi_ao, 
     &        ldelchi_ao, idelchi_ao, iscr, lscr, 
     &        lrscr, irscr
      integer lrchi_atom, irchi_atom, 
     &        lrdelchi_atom, irdelchi_atom
      integer nat_tot,nat_scr
      integer nqpts_per_buf, nqpts, ictr_buf
      integer lrqbuf, irqbuf,ldmat,idmat,i,
     .     grid_nbfm,ncube,istep,ntot_cube
      logical grid_file_rewind
      external grid_file_rewind,grid_nbfm
      Logical GRAD, oprint
      logical xc_chkgrad
      external xc_chkgrad
c
      oprint = util_print('quadrature', print_high)
c      write(*,*)'grid_quadv0a: pmat 1'
c      call ga_print(g_dens(1))
c
c     see if any Gradient Corrected Functional is not null
c
      nat_tot=0
      nat_scr=0
      grad=.true.
      grad=xc_chkgrad(xfac,cfac)
      me=ga_nodeid()
c
c     allocate scratch array for grid compression
c
      if (.not.MA_Push_Get(mt_int, nbf_ao, 'ibf_ao', libf_ao,
     &     iibf_ao))
     &     call errquit('grid_quadv0a: cannot allocate ibf_ao',2)
c
c     Determine the XC energy functionals.
c
      npol = (ipol*(ipol+1))/2
      acc_AO_gauss = dble(iAOacc)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
c
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit('grid_setspac: rewinding gridpts?', 0)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 3)
      nqpts_per_buf = n_per_rec
      if (.not.MA_Push_Get(mt_dbl,nqpts_per_buf*natoms,
     ,     'bp',lbp,ibp))
     &     call errquit('grid_quadv0a: cannot allocate bp',4)
      if (.not.MA_Push_Get(mt_dbl,nqpts_per_buf,
     ,     'rho scr',lrscr,irscr))
     &     call errquit('grid_quadv0a: cannot allocate rho scr',5)

      if (.not.MA_Push_Get(MT_Dbl,natoms,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit('dft_scf: cannot allocate rchi_atom',0)
      if(grad) then
      if (.not.MA_Push_Get(MT_Dbl,natoms,'rdelchi_atom',
     &     lrdelchi_atom,irdelchi_atom))call errquit
     &        ('dft_scf: cannot allocate rdelchi_atom',0)
      endif
C
      ntot_cube=0
!      do 200 iqsh = 1, n_tot_tasks
      do 200 iqsh = 1, n_rec_in_file
c
c     Define the current range of radial shells and integration center.
c     
         call grid_file_read(nqpts_per_buf, nqpts, ictr_buf, 
     &        rad,dbl_mb(irqbuf),nsubb)
         if(nqpts.gt.buffer_size) call
     '        errquit(' buffersize exceed by qpts ',nqpts)
        if(nqpts.eq.0) goto 200
        istep=0
        do 110 ncube=1,nsubb
c     
c        temp put buf into currently used arrays qxyz and qwght
c     
         call grid_repack(dbl_mb(irqbuf), qxyz, qwght, nqpts,
     ,          rad,istep)
!         do iii=1,nqpts
!            write(0,'(A,4F16.9)') 
!     .           ' RR ',
!     ,           qxyz(1,iii),qxyz(2,iii),qxyz(3,iii),
!     .           qwght(iii)
!         enddo
         if(nqpts.eq.0) goto 110
c     
c
c           Evaluate the AO basis set at each of the quad. points.
c           allocate arrays for exponents and contraction coefficients
c           within int_eval_gbsets
c           Evaluate the distances (squared) between each center and the qpoints.
c
         mbf_ao=nbf_ao
         if(natoms.gt.1) then
         call ifill(mbf_ao, 0, int_mb(iibf_ao), 1)
         mbf_ao=grid_nbfm(AO_bas_han,  natoms, 
     &        ictr_buf,rad,xyz,qxyz,nqpts,
     &        int_mb(iibf_ao),  iniz, ifin, expo, 
     &        minexp,ldiff,acc_ao_gauss)
         if (mbf_ao.eq.0) goto 110
         else
            iniz(1)=1
            ifin(1)=nbf_ao
            do i=0,nbf_ao-1
               int_mb(iibf_ao+i)=i+1
            enddo
         endif
c
         if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_ao, 'chi_ao', 
     &        lchi_ao, ichi_ao))
     &        call errquit('grid_quadv0a: cannot allocate chi_ao',0)
c     
         if(GRAD) then
            if (.not.MA_Push_Get(mt_dbl, 3*nqpts*mbf_ao,
     .           'delchi_ao', ldelchi_ao, idelchi_ao))
     &           call errquit('grid_quadv0a: cant alloc delchi_ao',0)
         endif
c     
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
         else                   ! general case
            max_elem_ang_scr = max_component * 28
         endif 
c     
         nscr = 3*nqpts*nbf_ao_mxcont + 
     &        max(3*nqpts*nbf_ao_mxprim,
     &        nqpts*max_elem_ang_scr) + 1
c     
         if (.not.MA_Push_Get(mt_dbl,nscr,'scrCD',lscr,iscr))
     &        call errquit('grid_quadv0a: cannot allocate scrCD',12)
c     
         call qdist(rq0, rq, qxyz, xyz, nqpts, natoms)
         call new_eval_gbsets(ao_bas_han, grad, .false., 
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &        hess, int_mb(iibf_ao), 0, rq0, 
     &        rq, qxyz, xyz, nqpts, nqpts, 
     &        mbf_ao, natoms, dbl_mb(iscr), nscr, 
     &        expo, nbf_ao_mxprim, ccoef, ncoef_max, 
     &        iniz,ifin)
c     
         if (.not.ma_pop_stack(lscr))
     &        call errquit('grid_quadv0a: cannot pop stack',0)
         
c     
c        get reduced Xi(r) and dXi(r) over atoms
c     
c
         call grid_rchi3(ntot_cube+1,
     ,        irchi, irdelchi,
     ,        nqpts, natoms,iniz,ifin,
     ,        dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
     ,        dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),
     ,        grad)
c     
         call xc_rhogen(tol_rho, AO_bas_han, g_dens,
     &        nbf_ao_mxnbf_ce, natoms, ipol, nqpts, nbf_ao,
     &        mbf_ao, ldelley, GRAD, Fmat, Pmat, Bmat, 
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &        int_mb(iibf_ao), iniz, ifin, rho, delrho, 
     &        dbl_mb(ibp), dbl_mb(irscr), 
     &        dbl_mb(irchi_atom), dbl_mb(irdelchi_atom),
     &        rdens_atom, cetobfr)
c
c
c           Evaluate the integrated charge density.
c     
         nqtot = nqtot + nqpts
c
         rho_int = rho_int + 
     &        ddot(nqpts, qwght, 1, rho, 1)
c      write(6,*)' iqsh,  nq, rho_int: ',
c     &                  iqsh,  nqpts, rho_int
c
            call xc_eval_fnl(rho, delrho, Amat,
!     &           Amat2, Cmat, Cmat2, nqpts, Ex, Ec, qwght,
     &           dummy, Cmat, dummy, nqpts, Ex, Ec, qwght,
     &           GRAD, .false., dummy, .false.)

            if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &         'Dmat', lDmat,iDmat))
     &         call errquit('grdi_quadv0_a: cannot allocate Dmat',0)
c     
            call xc_tabcd(tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, 
     &           Dbl_MB(iDmat), Dbl_MB(ichi_ao), 
     &           Dbl_MB(idelchi_ao), ipol, nqpts, mbf_ao, 
     &           nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD, 
     &           AO_bas_han, 
     &           natoms, iniz, ifin, g_vxc, INT_MB(iibf_ao),
     &           dbl_mb(irchi_atom), dbl_mb(irdelchi_atom),
     &           rdens_atom, cetobfr)


            if (.not.ma_pop_stack(lDmat))
     &           call errquit('grid_quadv0_a: cannot pop stack', 1)
            if(grad) then
               if (.not.ma_pop_stack(ldelchi_ao)) call
     &              errquit('grid_quadv0a: cannot pop stack', 2)
            endif
         if (.not.ma_pop_stack(lchi_ao))
     &        call errquit('grid_quadv0a: cannot pop stack', 3)
 110      ntot_cube=ntot_cube+1
  200 continue
      if(grad) then
         if (.not.ma_pop_stack(lrdelchi_atom))
     &        call errquit('grid_quadv0a: cannot pop stack',5)
      endif
      if (.not.ma_pop_stack(lrchi_atom))
     &     call errquit('grid_quadv0a: cannot pop stack',6)
      if (.not.ma_pop_stack(lrscr))
     &     call errquit('grid_quadv0a: cannot pop stack',7)
      if (.not.ma_pop_stack(lbp))
     &     call errquit('grid_quadv0a: cannot pop stack',8)
c
      if (.not.ma_pop_stack(lrqbuf))
     &     call errquit('grid_quadv0a: cannot pop stack', 9)
c

      if (.not.ma_pop_stack(libf_ao))
     &     call errquit('grid_quadv0a: cannot pop stack', 10)
c
      call ga_sync
c
c      write(6,*)' num_neg(1), neg_d_max(1): ',num_neg(1), neg_d_max(1)
c      write(6,*)' num_neg(2), neg_d_max(2): ',num_neg(2), neg_d_max(2)
c      write(6,*)' num_neg(3), neg_d_max(3): ',num_neg(2), neg_d_max(2)
C      write(0,*) ' SSCCRR ',nat_scr,nat_tot,100d0-((nat_scr*100d0)/
C     /nat_tot)
      return
      end
