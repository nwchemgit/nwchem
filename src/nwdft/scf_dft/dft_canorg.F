      subroutine dft_canorg( nmo, svals, g_svecs, g_xij)
c
c     Canonical Orthogonalization
c     Szabo & Ostlund page 144 (Dover Ed.)
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
c
      double precision svals(*) !  S evals [input]
      integer g_svecs ! GA handle for S evecs [input]
      integer g_xij  ! GA [ouput]
      integer nmo ! [out]
c
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
c     
      integer me, nproc, i
      integer ltmpm, itmpm,ndep
      double precision aaz
      integer ilo,ihi,jlo,jhi,numi,numj,icount
      integer jfirst,jlast,ptrfirst
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      ndep=nbf-nmo
c
c     
c        Build Xij
c     
      call ga_zero(g_xij)
      call ga_distribution(g_svecs,
     .           ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. ilo.le.ihi) then
         numi =  ihi-ilo+1
         numj =  jhi-jlo+1
         if (numi.gt.0 .and. numj.gt.0) then
            if (.not. MA_Push_Get(MT_Dbl, numi*numj, 'tmpm', 
     T           ltmpm, itmpm)) call errquit(
     T           'dftcano: failed to alloc tmpm',0, MA_ERR)
            call ga_get(g_svecs,ilo,ihi,jlo,jhi,
     .           dbl_mb(itmpm),numi)
            icount=itmpm
            jfirst=99999
            jlast=-1
            ptrfirst=0
            do i=jlo,jhi
               if(i.gt.ndep) then
                  if(ptrfirst.eq.0) ptrfirst=icount
                  jfirst=min(i-ndep,jfirst)
                  jlast=max(i-ndep,jlast)
                  aaz = svals(i)
                  aaz = 1.d0/sqrt(aaz)
                  call dscal(numi, aaz, DBL_MB(icount), 1)
               endif
               icount=icount+numi
            enddo
            if(jlast.ne.-1)
     C           call ga_put(g_xij,ilo,ihi,jfirst,jlast,
     .           dbl_mb(ptrfirst),numi)
            if(.not.MA_Pop_Stack(ltmpm)) call errquit(
     .           'failed popstack in dftcanorg ',0, 0)
         endif
      endif
c     
      return
      end
c
      integer function dft_lindep(nbf,toll_s,svals,oprint,luout)
      implicit none
      integer nbf ! [in]
      double precision toll_s ! [in]
      double precision svals(*) ! [in]
      logical oprint
      integer luout
c
      double precision smax,smin
      integer i,nmo
c
c     get linear dep
c
      smax=-1d30
      smin=1d30
      nmo = nbf
      do i=1,nbf
         if(svals(i).lt.toll_s) then
            nmo=nmo-1
            smax=max(smax,svals(i))
            smin=min(smin,svals(i))
         endif
      enddo
      if(oprint.and.nbf.ne.nmo) then
         write(luout,*)
         write(luout,'(A,i5,A)') ' WARNING : Found ',nbf-nmo,
     /        ' linear dependencies '
         write(luout,11) ' S eigenvalue threshold:',toll_s
 11      format(A,1pe16.5)
         write(luout,11) ' Smallest S eigenvalue :',smin
         write(luout,11) ' Largest  S eigenvalue :',smax
         write(luout,*)
      endif
      dft_lindep=nmo
      return
      end
c
      subroutine dft_guessin(movecs_in,ldmix,ncanorg,fon,
     ,     vecs_or_dens,
     ,     ipol,nbf_ao,g_movecs,g_gmovecs,
     ,     toll_s,svals)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "dft_scaleig.fh"
#include "global.fh"
#include "stdio.fh"
      character*(*) movecs_in
      logical ldmix ! [out]
      integer ncanorg ! [out]
      logical fon ! [in]
      character*(*) vecs_or_dens ! [out]
      integer ipol,g_movecs(2),g_gmovecs(2),nmo(2)
      integer nbf_ao
      double precision svals(*),toll_s
c
      external dft_lindep
      integer dft_lindep,nmo_out,me
c
      me=ga_nodeid()
      vecs_or_dens='density'
      ldmix=.false.
      ncanorg=0
      if((movecs_in.ne.'atomic'.and.
     *     movecs_in.ne.'fragment'.and.movecs_in.ne.'hcore')
     .     .or.fon) then
         vecs_or_dens='vectors'
         ldmix=.true.
         ncanorg=1
      endif
cedo with partial eigeensolver, need to get starting canonic. orthog
      if(yev.eq.'m') ncanorg=0
      nmo_out=dft_lindep(nbf_ao,toll_s,
     .        svals,me.eq.0,luout)
      if(nmo_out.ne.nbf_ao) then
         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'gvecs', 
     &        nbf_ao, 0, g_gmovecs(1))) 
     .        call errquit(' dftscf: ga_create failed ',0, GA_ERR)
         if(ipol.eq.2) then
            if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'gvecs', 
     &           nbf_ao, 0, g_gmovecs(2))) 
     .           call errquit(' dftscf: ga_create failed ',0, GA_ERR)
         endif
         call ga_copy(g_movecs(1),g_gmovecs(1))
         if(ipol.eq.2)
     /        call ga_copy(g_movecs(2),g_gmovecs(2))
      else
         g_gmovecs(1)=g_movecs(1)
         g_gmovecs(2)=g_movecs(2)
      endif
      if(movecs_in.eq.'atomic') then
         nmo(1)=nmo_out
         nmo(2)=nmo(1)
      endif
      return
      end
      subroutine dft_guessout(nmo,nbf_ao,g_gmovecs,g_movecs,ipol)
      implicit none
#include "errquit.fh"
#include "global.fh"
      integer nmo(2),nbf_ao,g_gmovecs(2),g_movecs(2),ipol
c
      if(nmo(1).ne.nbf_ao) then
         call ga_zero(g_movecs(1))
         call ga_copy_patch('n', g_gmovecs(1), 1, nbf_ao, 
     .        1, nmo(1),
     ,        g_movecs(1), 1,nbf_ao,1,nmo(1))
         if (.not. ga_destroy(g_gmovecs(1))) call errquit
     &      ('dft_scf: could not destroy g_gmovecs', 0, GA_ERR)
         if(ipol.eq.2) then
            nmo(2)=nmo(1)
            call ga_zero(g_movecs(2))
            call ga_copy_patch('n', g_gmovecs(2), 1, nbf_ao, 
     .           1, nmo(2),
     ,           g_movecs(2), 1,nbf_ao,1,nmo(2))
            if (.not. ga_destroy(g_gmovecs(2))) call errquit
     &         ('dft_scf: could not destroy g_gmovecs', 0, GA_ERR)
         endif
      endif
      return
      end
      subroutine util_tolls(rtdb,oprint,toll_s,n_dep,svals,nbf)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      integer nbf
      double precision svals(*)
      logical oprint
c
      integer n_dep
      double precision toll_s
c
c     Get user input parameters ... only one allowed
c
c     lindep:n_dep  ... force number of linearly dep vectors vectors
c     lindep:toll_s    ... tolerance 
c
      toll_s=1d-5
       if (rtdb_get(rtdb, 'lindep:tol', mt_dbl, 1, toll_s).or.
     .    rtdb_get(rtdb, 'dft:toll_s', MT_DBL, 1, toll_s)) then
          if (oprint) then
             write(luout, 3) toll_s
 3           format(/' Setting tolerance for linear dependence to',
     $            1p,d9.2)
             call util_flush(luout)
          endif
       endif
       if (rtdb_get(rtdb, 'lindep:n_dep', mt_int, 1, n_dep)) then
          if (oprint) then
             write(luout, 2) n_dep
 2           format(/' Forcing the number of linearly dependent ',
     $            'vectors to be', i3)
             call util_flush(luout)
          endif
       else
          n_dep = -1
       endif
       if(n_dep.ne.-1) then
          if(n_dep.eq.0) then
             toll_s=-1.d9
          else
             toll_s=svals(n_dep+1)
          endif
       endif
c
c     this to keep sodft happy
c
       if (.not. rtdb_put(rtdb, 'dft:toll_s', MT_DBL, 1,
     .      toll_s)) call errquit(
     .      'utiltools: rtdbput failed',0, RTDB_ERR)

       if(oprint) then
          write(luout,*) ' lindep tolerance = ',toll_s
       endif
c
!       if ((toll_s.ne.-1.0d0) .and. (n_dep.ne.-1)) call errquit
!     $      ('scf_lindep: only one of n_dep and tol permitted',0)
       return
       end
      subroutine dft_convpam(oprint,iter,etnew,etold,enuc)
      implicit none
#include "cdft.fh"
#include "stdio.fh"
#include "global.fh"
#include "dft_conv.fh"
#include "util.fh"
c
      logical oprint
      logical use_nwxc
      integer iter
      double precision etnew,etold,enuc
c
      double precision current_cpu
      logical gotsmear
      double precision ss
      gotsmear(ss)=abs(ss).gt.1d-9
c
      delta = -etold+etnew
      if (abs(delta).lt.5d-2.or.ltight)then
c     
c        Set coulomb acc to max (e.g., input parameter).
c        (note, may also require re-initializing DIIS)
c     
         itol2e = itol_max
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         use_nwxc = util_module_avail("nwxc")
         if (use_nwxc) then
            call nwxc_getvals("nwxc_is_on",use_nwxc)
         endif
         if (use_nwxc) call nwxc_reset_rho_tolerance(tol_rho)
         if(iter.ge.nfock/2) iswitc = iswitc+1
         if(abs(delta).lt.1d-6) iswitc = iswitc+2
      else
         iswitc=0
      endif
c     
      iter = iter + 1
c     
c     Check convergence parameters.
c     
      if(oprint) then
         write(luout,*) ' delta = ',delta
      endif
      if ((abs(delta).lt.dampon.and.abs(delta).gt.dampoff).
     .     or.iter.le.ncydp.or.iter.le.1)then
         damping = .true.
         ndamp = ndamp_input
      else
         damping = .false.
         ndamp = ndamp_def
      endif
c     
      if ((abs(delta).lt.levlon.and.abs(delta).gt.levloff).or.
     &    (iter.le.ncysh))then
         if (homo_lumo_gap.lt.hl_tol.and.(.not.gotsmear(ssmear)))then
            levelshifting = .true.
            rlshift = rlshift_input
            if (check_shift)then
               if (lumo .lt. homo)then
                  levelshifting = .false.
                  if (ga_nodeid().eq.0 .and. oprint)
     &               write(LuOut,2224)homo, lumo
               endif
            endif
         else
            levelshifting = .false.
            rlshift = rlshift_def
         endif 
      else
         levelshifting = .false.
         rlshift = rlshift_def
      endif
c     
      if ((abs(delta).lt.diison.and.abs(delta).gt.diisoff).or.
     &     iter.le.ncyds.or.keep_diis_on)then
         diising = .true.
c     
c        Once started, keep DIIS on until diisoff threshold.
c     
         keep_diis_on = .true.
      else
         diising = .false.
      endif
      if (abs(delta).lt.diisoff.or.
     .     (ncyds.gt.0.and.iter.gt.ncyds))then
         diising = .false.
         keep_diis_on = .false.
      endif
c     
      if (nodamping)damping = .false.
      if (nolevelshifting) then 
        levelshifting = .false.
        rlshift=rlshift_def
      endif        
      if (nodiis)diising = .false.
c
      lmaxov = lmaxov_sv
c      don't turn off max_ovl
c      if ((lumo - homo).lt.-hl_tol.and.lmaxov)then
c         lmaxov = .false.
c         if (ga_nodeid().eq.0 .and. oprint)
c     &      write(LuOut,224)homo, lumo
c      endif
c     
      if (oprint.and.ga_nodeid().eq.0)then         
         current_cpu = util_cpusec()
         if (nexc.le.1)then
            write(LuOut,222)etnew+enuc, ecore, Ecoul, Exc(1), enuc, 
     &                      rho_n, current_cpu
         else
            write(LuOut,223)etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     &                      enuc, rho_n, current_cpu
         endif
      endif
      return
 222  format(//
     &     '         Total DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '    Exchange-Corr. energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12//
     &     ' Numeric. integr. density =', f22.12//
     &     '     Total iterative time =', f9.1,'s'//)
 223  format(//
     &     '         Total DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '          Exchange energy =', f22.12/
     &     '       Correlation energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12//
     &     ' Numeric. integr. density =', f22.12//
     &     '     Total iterative time =', f9.1,'s'//)
 224  format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     &     /,10x,'Significant orbital reordering with',
     &     ' maximum overlap',
     &     /,10x,'turned ON.  Turning max_ovl OFF.')
 2224 format(10x,'HOMO = ',f6.2,' LUMO (with shift) = ',f6.2,
     &     /,10x,'Unshifted LUMO is less than HOMO.',
     &     /,10x,'Turning levelshifting OFF this iteration.')

      end
      double precision function dft_rabuck(ipol,iter,n_rabuck,noc,
     ,     derr,evala,evalb)
      implicit none
      integer ipol,iter,n_rabuck,noc(2)
      double precision evala(*),evalb(*),derr(*),delta
c
      integer nn,nss
      double precision ss,s_rabuck,derrt
      double precision smear_in
      s_rabuck(nn,ss,nss,delta)=min(ss*exp((1d0/(nss-1d0))*
     *     log(abs(1d0-10d0/ss))*(nn-1d0)),delta*1d-4)
c
      derrt=derr(1)
      if(ipol.eq.2) derrt=max(derr(1),derr(2))
      
      smear_in=max(abs(evala(noc(1))-evala(noc(1)+1)),1d-6)
      if(ipol.eq.2) smear_in=max(smear_in,
     ,     abs(evalb(noc(2))-evalb(noc(2)+1)))
      dft_rabuck=s_rabuck(iter,smear_in,n_rabuck,derrt)
c     ssmear=s_rabuck(iter,smear_in,n_rabuck,derrt)
c     ssmear=s_rabuck(iter,smear_in,n_rabuck,delta)
      
      return
      end
      integer function dft_nonzvecs(g_vecs,nbf,noc)
      implicit none
#include "errquit.fh"
c
      integer g_vecs ! GA handle for  vecs [input]
      integer nbf
      integer noc
c
#include "mafdecls.fh"
#include "global.fh"
#include "dft_scaleig.fh"
c     
      integer me, nproc, i
      integer ltmpm, itmpm,l_nz,k_nz
      double precision dabsmax,aaz,toll
      parameter(toll=1d-9)
      external dabsmax
      integer ilo,ihi,jlo,jhi,numi,numj,icount
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c     
      if (.not. MA_Push_Get(MT_int, nbf , 'nonzero', 
     T     l_nz, k_nz)) call errquit(
     T           'dftnz: failed to alloc nz',0, MA_ERR)
      call ifill(nbf,0d0,int_mb(k_nz),1)
      call ga_distribution(g_vecs,
     .           ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. ilo.le.ihi) then
         numi =  ihi-ilo+1
         numj =  jhi-jlo+1
         if (numi.gt.0 .and. numj.gt.0) then
            if (.not. MA_Push_Get(MT_Dbl, numi*numj, 'tmpm', 
     T           ltmpm, itmpm)) call errquit(
     T           'dftcano: failed to alloc tmpm',0, MA_ERR)
            call ga_get(g_vecs,ilo,ihi,jlo,jhi,
     .           dbl_mb(itmpm),numi)
            icount=itmpm
            do i=jlo,jhi
               aaz=dabsmax(numi,dbl_mb(icount))
               if(aaz.gt.toll) int_mb(k_nz+i-1)=1d0
               icount=icount+numi
            enddo
            if(.not.MA_Pop_Stack(ltmpm)) call errquit(
     .           'failed popstack in dftcanorg ',0, 0)
         endif
      endif
      call ga_igop(375, int_mb(k_nz), nbf, '+')
      do i=nbf,1,-1
         if(int_mb(k_nz+i-1).ne.0) then
            dft_nonzvecs=i
            goto 1
         endif
      enddo
 1    continue
c     we should be using all the electrons .. right?
      dft_nonzvecs=max(noc,dft_nonzvecs)
#if 0
      write(6,*) ' nonzvecs ',dft_nonzvecs
      call util_flush(6)
#endif
c     
      if(.not.MA_Pop_Stack(l_nz)) call errquit(
     .     'failed popstack in dftnz ',0, 0)
      return
      end
c
