      Subroutine dftg_gridv0a(rtdb, iga_dens,ldew,
     &                         ncoef_max, NXprim, max_at_bf,
     &                         rho_int, nqtot, iniz, ifin, rho, delrho,
     &                         func, d_qwght, qwght, qxyz, xyz, 
     &                         iandex, expo, ccoef, Amat, Amat2, Cmat,
     &                         Cmat2, Bmat, Dmat, Fmat, Pmat,  rq,
     &                         rq0, hh,  
     .     ibf_ao,curatoms,rchi_atom,rdelchi_atom,ppp,ff,zzz,
     ,     force, hess, g_rhs, do_force,
     &                         do_hess, do_cpks_r, natoms, oactive,
     &                         nactive, rdens_atom, cntobfr,
     &                         cetobfr)
c      
c     $Id: dftg_gridv0a.F,v 1.4 2000-04-14 01:47:45 edo Exp $
c     
      implicit none
#include "bas.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "cgridfile.fh"
c
#include "dft2drv.fh"
c !!! BGJ test !!!
#include "bgj.fh"
c !!! BGJ test !!!
      integer rtdb
      integer natoms
      integer ncoef_max, NXprim, max_at_bf
      integer nqtot
      integer nactive
      logical oactive(natoms)
      logical grad, gradchi, hesschi, der3chi
      logical xc_chkgrad
      external xc_chkgrad
      logical ldew
      logical do_force, do_hess, do_cpks_r
c
      integer cntobfr(2,nshells_ao), 
     &        cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      double precision Force(3*natoms), hess(3*natoms,3*natoms)
      integer g_rhs(*)
      double precision hh(3,n_per_rec)
      double precision rho(n_per_rec,ipol*(ipol+1)/2)
      double precision delrho(n_per_rec,3,ipol)
      double precision Amat(n_per_rec,ipol), Cmat(n_per_rec,3,ipol),
     &                 Bmat(n_per_rec,*), Dmat(n_per_rec,3,max_at_bf)
      double precision Amat2(n_per_rec,NCOL_AMAT2)
      double precision Cmat2(n_per_rec,NCOL_CMAT2)
      double precision Pmat(max_at_bf*max_at_bf),
     &                 Fmat(max_at_bf*max_at_bf)
      double precision func(n_per_rec)
      double precision d_qwght(3,n_per_rec,*)
      double precision qwght(n_per_rec), qxyz(3,n_per_rec)
      double precision xyz(3,natoms)
      integer iandex(natoms)
      integer iniz(natoms), ifin(natoms)
      integer ibf_ao(nbf_ao),curatoms(natoms)
      double precision rchi_atom(natoms),rdelchi_atom(3,natoms)
      double precision ppp(*),ff(*),zzz(*)
      double precision ccoef(ncoef_max), expo(NXprim)
      double precision rq(n_per_rec,natoms)
      double precision rq0(natoms)
c
c     external
c
      double precision ddot
      integer nxtask
      External ddot, nxtask
      double precision eps
      parameter(eps=1.d-8)
c
c---> Matrix Elements of the XC Potential & Energy
c
      integer iga_dens(2)
c
c     local crap
c
      integer npol
      integer iqsh
      integer me
      integer  nbf_max
      integer NXang, llmax, MXcont, max_elem_ang_scr
      integer nbf_ao_mxnbf_ce2
      double precision acc_AO_gauss 
      double precision Ex, Ec, rho_int,rad
      integer nqpts_per_buf, nqpts, ictr_buf
      integer lrqbuf, irqbuf,
     .     grid_nbfm,ncube,istep,ntot_cube
      logical grid_file_rewind
      external grid_file_rewind,grid_nbfm
c
c     
      me = ga_nodeid()
      Ec = 0.d0
      Ex = 0.d0
      if (bgj_print() .gt. 0)
     &    write(*,*)'*** dftg_quadv0a: logicals ',
     &    do_force,do_hess,do_cpks_r
c
c     Determine the XC energy functionals.
c
c
c     Define the number of basis functions in the AO Gaussian basis set.
c
      if (.not. bas_ncontr_cn_max(ao_bas_han, MXcont))
     &   call errquit('Exiting in dftg_quadv0.',1)
      if (.not. bas_high_angular(AO_bas_han, NXang))
     &   call errquit('exiting in dftg_quadv0.',7)
      llMAX=((NXang+2)*(NXang+1))/2
c scratch for eval_func
      if (NXang .ge. 3)then
         max_elem_ang_scr = 28
      elseif (NXang .ge. 1)then
         max_elem_ang_scr = 9
      else 
         max_elem_ang_scr = 3
      endif 
      max_elem_ang_scr=max_elem_ang_scr*llMAX
c
      if (.not.bas_nbf_cn_max(AO_bas_han, nbf_max))
     &   call errquit('exiting in dftg_quadv0.',7)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
c
c     Define the max number of basis functions, contractions,
c     and primitives comprising any given shell of this basis set.
c
      acc_AO_gauss= iAOacc
c
c
c     weight derivatives switched off
c
      if(ldew.and.(me.eq.0)) write(0,*) ' derivtive weights on'
      if (ldelley.and.ldew)then
         if (me.eq.0)write(LuOut,*) 
     &      ' derivative weights not available for Delley weights '
         ldew=.false.
      endif
c
c     see if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad(xfac,cfac)
c
c     allocate scratch array for grid compression
c
      if (grad)then
         gradchi = .true.
         hesschi = .true.
         der3chi = do_hess
      else
         gradchi = .true.
         hesschi = do_hess
         der3chi = .false.
      endif
c      write(6,*)'*** dftg_quadv0a: hesschi der3chi',hesschi,der3chi
c
      npol = (ipol*(ipol+1))/2
c     
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit('grid_setspac: rewinding gridpts?', 0)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 3)
      nqpts_per_buf = n_per_rec
c
      ntot_cube=0
      do 200 iqsh = 1, n_rec_in_file
c         
c           Define the current range of radial shells and integration center.
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
         if(nqpts.eq.0) goto 110

         call dftg_gridv0b(nqpts,rad,ictr_buf,iga_dens,
     &                         ncoef_max, NXprim, max_at_bf,
     &                         rho_int, nqtot, iniz, ifin, rho, delrho,
     &                         func, d_qwght, qwght, qxyz, xyz, 
     &                         iandex, expo, ccoef, Amat, Amat2, Cmat,
     &                         Cmat2, Bmat, Dmat, Fmat, Pmat,  rq,
     &                         rq0, hh, 
     .     ibf_ao,curatoms,rchi_atom,rdelchi_atom,ppp,ff,zzz,
     ,     force, hess, g_rhs, do_force,
     ,     llmax,max_elem_ang_scr,nbf_ao_mxnbf_ce2,mxcont,
     ,     acc_ao_gauss, ldew,grad,gradchi,hesschi,der3chi,npol,
     &                         do_hess, do_cpks_r, natoms, oactive,
     &                         nactive, rdens_atom, 
     &                         cetobfr)

 110      ntot_cube=ntot_cube+1
  200 continue  ! bottom of shell loop
      if (.not.ma_pop_stack(lrqbuf))
     &     call errquit('dftg_quadv0a: cannot pop stack', 0)
      return 
c     
      end
