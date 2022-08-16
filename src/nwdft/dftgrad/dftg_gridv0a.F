      Subroutine dftg_gridv0a(rtdb,iga_dens,ldew,
     &            max_at_bf,
     &            rho_int, nqtot, iniz, ifin, rho, delrho,
     &            func, d_qwght, qwght, qxyz, xyz, xyzm,
     &            iandex, expo, Amat, Amat2, Cmat,
     &            Cmat2, Bmat, Fmat, Pmat, rq,
     &            rq0, hh,  
     .            ibf_ao,curatoms,rchi_atom,rdelchi_atom,ppp,ff,zzz,
     .            drhonuc,ddelrhonuc,d2scr,nscr,
     ,            force, hess, g_rhs, do_force,
     &            do_hess, do_cpks_r, natoms, oactive,
     &            nactive, rdens_atom, cetobfr, ttau, Mmat, Mmat2,
     &            Lmat,
     1            do_gxc, Amat3, Cmat3, ngxc,
     2            rhosave, delrhosave, triplet, 
     3            do_tddftvxc,cd_coef)
c      
c     $Id$
c     
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "cgridfile.fh"
c
#include "dft2drv.fh"
c Daniel (1-14-13): Added 3rd derivative matrix header
#include "dft3drv.fh"
c !!! BGJ test !!!
#include "bgj.fh"
c !!! BGJ test !!!
c Qin
      integer rtdb
      logical lcdft, ldbw
      integer ipop, nconstr, i
      integer iatms1(max_constr), iatms2(max_constr)
      integer fatms1(max_constr), fatms2(max_constr)
      integer dtype(max_constr)
      double precision multipliers(max_constr), fact(max_constr)  ! pass to dftg_gridv0b
      double precision cd_coef(nbf_cd,*)

      integer natoms
      integer max_at_bf
      integer nqtot
      integer nactive
      logical oactive(natoms)
      logical grad, gradchi, hesschi, der3chi
      logical xc_chkgrad, xc_chktau, xc_chklap
      external xc_chkgrad, xc_chktau, xc_chklap
      logical ldew, kske, kslap
      logical do_force, do_hess, do_cpks_r
c
      integer cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      double precision Force(3*natoms), hess(3*natoms,3*natoms)
      integer g_rhs(*)
      double precision hh(3,n_per_rec)
      double precision rho(n_per_rec,ipol*(ipol+1)/2)
      double precision delrho(n_per_rec,3,ipol)

      double precision Mmat(*), ttau(*)
      double precision Lmat(*)

      double precision Amat(n_per_rec,ipol), Cmat(n_per_rec,3,ipol),
     &                 Bmat(n_per_rec,*)
c
      double precision Amat2(n_per_rec,NCOL_AMAT2)
      double precision Cmat2(n_per_rec,NCOL_CMAT2)
      double precision Mmat2(n_per_rec,NCOL_MMAT2)
c
      double precision Pmat(max_at_bf*max_at_bf),
     &                 Fmat(max_at_bf*max_at_bf)
      double precision func(n_per_rec)
      double precision d_qwght(3,n_per_rec,*)
      double precision qwght(n_per_rec), qxyz(3,n_per_rec)
      double precision xyz(3,natoms),xyzm(3,natoms)
      integer iandex(natoms)
      integer iniz(natoms), ifin(natoms)
      integer ibf_ao(nbf_ao),curatoms(natoms)
      double precision rchi_atom(natoms),rdelchi_atom(natoms)
      double precision ppp(*),ff(*),zzz(*)
      double precision expo(*)
      double precision rq(n_per_rec,natoms)
      double precision rq0(natoms)
      double precision drhonuc(*),ddelrhonuc(*),d2scr(*)
      integer nscr ! [in]
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
      integer iqsh
      integer me
      integer  nbf_max
      integer nbf_ao_mxnbf_ce2, ncontrset
      double precision acc_AO_gauss 
      double precision Ex, Ec, rho_int,rad
      integer nqpts_per_buf, nqpts, ictr_buf
      integer lrqbuf, irqbuf,
     .     ncube,istep,ntot_cube
      integer lbas_cent_info, ibas_cent_info,
     &        lbas_cset_info, ibas_cset_info,
     &        ldocset, idocset
      logical grid_file_rewind,grid_file_read
      external grid_file_rewind,grid_file_read
c Daniel (1-11-13): Define XC-3rd derivatives
      logical do_gxc
      double precision Amat3(n_per_rec,NCOL_AMAT3)
      double precision Cmat3(n_per_rec,NCOL_CMAT3)
      integer ngxc
c Daniel (1-17-13): Matrices for storing perturbed densities.
      double precision rhosave(*)
      double precision delrhosave(*)
c Daniel (2-6-13): Variable for doing restricted triplet TDDFT
c gradients from the XC-kernel
      logical triplet
c Daniel (2-7-13): Logical variable for doing the dVxc*P contribution
c to the TDDFT gradients.
      logical do_tddftvxc
      integer basis

      if (adft) then
        basis = cd_bas_han
      else
        basis = ao_bas_han
      endif
c
c Qin
      ldbw = .false.
      if (.not. rtdb_get(rtdb, 'dft:lcdft', mt_log, 1, lcdft))
     &     lcdft = .false.
      if(lcdft) then
        if (.not. rtdb_get(rtdb, 'dft:ipop', mt_int, 1, ipop))
     &     ipop = 3
        if(ipop.eq.1) then
         ldbw = .true.
         if (.not. rtdb_get(rtdb, 'dft:nconstr', mt_int, 1, nconstr))
     &     call errquit('dft_gridv0a: get nconstr failed', 0, RTDB_ERR)
        if (.not. (
     $            rtdb_get(rtdb,'dft:iatms1', mt_int, nconstr,iatms1(1))
     $      .and. rtdb_get(rtdb,'dft:iatms2', mt_int, nconstr,iatms2(1))
     $      .and. rtdb_get(rtdb,'dft:fatms1', mt_int, nconstr,fatms1(1))
     $      .and. rtdb_get(rtdb,'dft:fatms2', mt_int, nconstr,fatms2(1))
     $      ))
     &   call errquit('dft_gridv0a: atms failed', 2700, RTDB_ERR)
        if (.not.rtdb_get(rtdb, 'dft:multipliers', mt_dbl, nconstr,
     $       multipliers(1)))
     &    call errquit('dft_gridv0a: get multipliers failed', 0)
        if (.not. rtdb_get(rtdb,'dft:dtype', mt_int, nconstr, dtype(1)))
     &   call errquit('dft_gridv0a: get dtype failed', 3000, RTDB_ERR)
        endif
      endif
c Qin / end
c     
      me = ga_nodeid()
      Ec = 0.d0
      Ex = 0.d0
      if (bgj_print() .gt. 0)
     &    write(*,*)'*** dftg_gridv0a: logicals ',
     &    do_force,do_hess,do_cpks_r
c
c     Determine the XC energy functionals.
c
c
      if (.not.bas_nbf_cn_max(basis, nbf_max))
     &   call errquit('exiting in dftg_gridv0.',7, BASIS_ERR)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
c
c     Define the max number of basis functions, contractions,
c     and primitives comprising any given shell of this basis set.
c
      acc_AO_gauss= iAOacc
c
      if (ldelley.and.ldew)then
         if (me.eq.0)write(LuOut,*) 
     &      ' derivative weights not available for Delley weights '
         ldew=.false.
      endif
c
c     See if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad()
      kske=xc_chktau()
      kslap=xc_chklap()
c
c     allocate scratch array for grid compression
c
      if (kslap)then
         gradchi = .true.
         hesschi = .true.
         der3chi = .true.
      elseif (grad.or.kske)then
         gradchi = .true.
         hesschi = .true.
         der3chi = do_hess
      else
         gradchi = .true.
         hesschi = do_hess
         der3chi = .false.
      endif
c      write(6,*)'*** dftg_gridv0a: hesschi der3chi',hesschi,der3chi
c
c     Allocate and create info for new basis function evaluator
c
      if (.not.bas_numcont(basis, ncontrset))
     &     call errquit('Exiting in dftg_gridv0a',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('dftg_gridv0a: cannot allocate docset',
     ,     ncontrset,0)
      if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('dftg_gridv0a: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('dftg_gridv0a: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(basis, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natoms)
c     
c     rewind grid pts file
c
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'grid pts buffer', lrqbuf, irqbuf))
     &   call errquit('dftg_gen: cannot allocate grid pt buffer', 3,
     &       MA_ERR)
      if(n_rec_in_file.eq.0) goto 2001
      if (.not. grid_file_rewind())
     $   call errquit('grid_setspac: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
      nqpts_per_buf = n_per_rec
c
      ntot_cube=0
      do 200 iqsh = 1, n_rec_in_file
c         
c       Define the current range of radial shells and integration center.
c          
        if(.not.grid_file_read(nqpts_per_buf, nqpts, ictr_buf, 
     &      rad,dbl_mb(irqbuf),nsubb)) goto 2001
        if(nqpts.gt.buffer_size) call
     *      errquit(' buffersize exceed by qpts',nqpts,UNKNOWN_ERR)
        if(nqpts.eq.0) goto 200
        istep=0
        do 110 ncube=1,nsubb
c     
c        temp put buf into currently used arrays qxyz and qwght
c     
         call grid_repack(dbl_mb(irqbuf), qxyz, qwght, nqpts,
     ,          rad,istep)
         if(nqpts.eq.0) goto 110
c Qin
         if(ldbw) then
           do i = 1, nconstr
            if(ictr_buf.ge.iatms1(i) .and. ictr_buf.le.fatms1(i)) then
               fact(i) = multipliers(i)
            elseif(ictr_buf.ge.iatms2(i).and.ictr_buf.le.fatms2(i))then
               fact(i) = -multipliers(i)
            else
               fact(i) = 0d0
            endif
          enddo
         endif
c Qin/end
c Daniel (1-16-13): Thread a variable through here that allows XC-third
c derivativeis to happen.
         call dftg_gridv0b(nqpts,rad,ictr_buf,iga_dens,
     &     max_at_bf,
     &     rho_int, nqtot, iniz, ifin, rho, delrho,
     &     func, d_qwght, qwght, qxyz, xyz, xyzm,
     &     iandex, expo, Amat, Amat2, Cmat,
     &     Cmat2, Bmat, Fmat, Pmat,  rq,
     &     rq0, hh, 
     .     ibf_ao,int_mb(ibas_cent_info),int_mb(ibas_cset_info),
     .     log_mb(idocset),curatoms,rchi_atom,rdelchi_atom,ppp,
     .     ff,zzz,drhonuc,ddelrhonuc,d2scr,
     ,     force, hess, g_rhs, do_force,
     ,     nbf_ao_mxnbf_ce2,
     ,     acc_ao_gauss, ldew,grad,gradchi,hesschi,der3chi,
     &     do_hess, do_cpks_r, natoms, oactive,
     &     nactive, rdens_atom, 
     &     cetobfr,
     &     ttau, kske, Mmat, Mmat2, 
     &     kslap, Lmat, ldbw,
     $     fact, dtype, nconstr,  ! Qin
     1     do_gxc, Amat3, Cmat3, ngxc,
     2     rhosave, delrhosave, triplet, 
     3     do_tddftvxc,cd_coef)

 110      ntot_cube=ntot_cube+1
  200 continue  ! bottom of shell loop
 2001 continue
      if (.not.ma_pop_stack(lrqbuf))
     &     call errquit('dftg_gridv0a: cannot pop stack', 0, MA_ERR)
c
         if (.not.MA_Pop_Stack(lbas_cset_info))
     &        call errquit('dftg_gridv0a: pop stack failed.',2, MA_ERR)
         if (.not.MA_Pop_Stack(lbas_cent_info))
     &        call errquit('dftg_gridv0a: pop stack failed.',3, MA_ERR)
      if (.not.MA_Pop_Stack(ldocset))
     &     call errquit('dftg_gridv0a: pop stack failed.',4, MA_ERR)
c
      return
      end
