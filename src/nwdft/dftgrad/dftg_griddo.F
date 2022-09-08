      Subroutine dftg_griddo(rtdb,
     &     g_dens,ldew,
     &     ncoef_max, max_at_bf,
     &     rho_int, nqtot, iniz, ifin, rho, delrho,
     &     func, d_qwght, qwght, qxyz, xyz,
     &     iandex, expo, ccoef, Amat, Amat2, Cmat,
     &     Cmat2, Bmat,  Fmat, Pmat, rq,
     &     rq0, hh,  
     .     ibf_ao,curatoms,rchi_atom,rdelchi_atom,ppp,ff,zzz,
     .     drhonuc,ddelrhonuc,d2scr,iscr,nscr,
     ,     force, hess, g_rhs, do_force,
     &     do_hess, do_cpks_r, natoms, oactive,
     &     nactive, rdens_atom, cetobfr, ttau, 
     &     Mmat, Mmat2,
     &     Lmat,
     1     do_3rd, Amat3, Cmat3, ngxc,
     2     rhosave, delrhosave,
     3     do_tddftvxc,cd_coef)
c      
c     $Id$
c     
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "util.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "cgridfile.fh"
#include "grid_cube.fh"
c
#include "dft2drv.fh"
c !!! BGJ test !!!
#include "bgj.fh"
c !!! BGJ test !!!
      integer rtdb
      integer natoms
      integer ncoef_max, max_at_bf
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
      double precision xyz(3,natoms)
      integer iandex(natoms)
      integer iniz(natoms), ifin(natoms)
      integer ibf_ao(nbf_ao),curatoms(natoms)
      double precision rchi_atom(natoms),rdelchi_atom(natoms)
      double precision ppp(*),ff(*),zzz(*)
      double precision ccoef(*), expo(*)
      double precision rq(n_per_rec,*)
      double precision rq0(*)
      double precision drhonuc(*),ddelrhonuc(*),d2scr(*),
     ,     iscr(*)
      integer nscr ! [in]
      double precision cd_coef(*)
c     external
c
      double precision eps
      parameter(eps=1.d-8)
c
c---> Matrix Elements of the XC Potential & Energy
c
      integer g_dens(*)
c
c     local crap
c
      integer me
      integer  nbf_max
      integer nbf_ao_mxnbf_ce2
      double precision acc_AO_gauss 
      double precision Ex, Ec, rho_int
      integer lrqbuf, irqbuf
c
      integer k_qsh,l_qsh
      logical dosympts
      integer nops
      logical oprint
      integer nqlen,mxnrad,ndbl,
     ,     lrij,irij,l_qat,k_qat,ibp,lbp
      integer idum1,irqlen
      logical ldum1
      double precision ddum1
c
      logical grid_file_rewind,grid_file_read
      external grid_file_rewind,grid_file_read
c Daniel (1-11-13): These variables prevent XC-third derivatives from
c happening and prevent future issues when people add to the code.
      logical do_3rd
      double precision Amat3(*), Cmat3(*)
      double precision ngxc
      double precision rhosave(*)
      double precision delrhosave(*)
c Daniel (2-7-13): Variable for the alternative approach to doing
c the dVxc*P contribution in TDDFT gradients.
      logical do_tddftvxc
      integer basis
c
c     
      if (adft) then
        basis = cd_bas_han
      else
        basis = ao_bas_han
      endif
      me = ga_nodeid()
      Ec = 0.d0
      Ex = 0.d0
      oprint = util_print('quadrature', print_high)
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
c     See if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad()
      kske=xc_chktau()
      kslap=xc_chklap()
c
c     allocate scratch array for grid compression
c
      if (kslap) then
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
      if (.not. rtdb_get(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &     call errquit('dftgriddo: rtdb_get failed', 1, RTDB_ERR)
      nqlen = natoms*mxnrad
      if (.not.MA_Push_get(MT_int,4*nqlen,'shell list',
     &     l_qsh,k_qsh))
     &     call errquit('dftgriddo: cannot allocate shell list',0,
     &     MA_ERR)

      call grid_init(rtdb,geom,
     ,     natoms,
     ,     oskel,leb,dosympts,int_mb(k_qsh))

      if (.not.ma_push_get(MT_dbl, 4*buffer_size+4, 
     &                     'grid pts buffer', lrqbuf, irqbuf))
     &   call errquit('dftgriddo: cannot allocate grid pt buffer', 3,
     &       MA_ERR)
         ndbl = (natoms*(natoms+1))/2
         if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lrij,irij))
     &        call errquit('griddo: cannot allocate Rij',0, MA_ERR)
         if (.not.MA_Push_get(mt_int,nqshells,'qshat',l_qat,k_qat))
     &        call errquit('griddo: cannot allocate qat',5, MA_ERR)
      irqlen=natoms*n_per_rec
      irqlen=max(irqlen,natoms*nqmax)
      if (.not.MA_Push_Get(mt_dbl,irqlen,'bp',lbp,ibp))
     &     call errquit('grid_quadv0a: cannot allocate bp',4, MA_ERR)
      
c Daniel (1-11-13): Added XC-third derivative stuff to prevent problems
c with subroutine calls after this one.
cdebug      write(*,*) 'cd_coef', cd_coef(1:nbf_cd)

      call grid_loop(.true.,
     ,     nops,
     ,     oprint,grad,g_dens,
     ,     dosympts,sub4,largenode,grid_written,
     ,     ncoef_max,natoms,
     ,     nqtot,iniz, ifin, rho, delrho, 
     &     hess, qwght, qxyz, xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,dbl_mb(ibp),rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,idum1,
     ,     Amat2, Cmat2,
     ,     idum1, ldum1, ldum1, do_cpks_r,ldum1,ldum1,
     ,     cd_coef,ldum1,ldum1,
     &     ttau, kske, Mmat, Mmat2,
     &     kslap, Lmat,
     ,     int_mb(k_qat),
     ,     int_mb(k_qsh),dbl_mb(irqbuf),dbl_mb(irij),rho_int,
     D     func,d_qwght,rq0,hh,ppp,ff,zzz,
     D     drhonuc,ddelrhonuc,d2scr,force,
     D     acc_ao_gauss,
     L     do_force,do_hess,ldew,gradchi,hesschi,der3chi,oactive,
     I     iandex,max_at_bf,curatoms,g_rhs,
     I     nbf_ao_mxnbf_ce2,nactive, 
     1     do_3rd, Amat3, Cmat3, ngxc,
     2     rhosave, delrhosave,
     3     do_tddftvxc)
c
      if (.not.ma_chop_stack(l_qsh))
     &   call errquit('dftgriddo: cannot pop stack', 3, MA_ERR)
c
      return
      end
