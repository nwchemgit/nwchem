      Subroutine dft_fitcd(nfit,CD_coef, i3c_ERI, Ecoul1, 
     &                     g_dens, nTotEl, n_batch, n_bmax, iwhat,
     &                     n3c_dbl, iwhat_max, n_semi_bufs, fd, IOLGC,
     &                     natoms,
     &                     derfit,  Tvec, have_Tvec)
*
* $Id: dft_fitcd.F,v 1.41 2004-09-07 06:11:45 windus Exp $
*
      implicit none
#include "errquit.fh"
c
#include "bas.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "schwarz.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "sym.fh"
#include "geom.fh"
#include "inp.fh"
c
      integer natoms
      integer nfit           ! number of different fits to perform
      integer g_dens(nfit)   ! GA DM handles
      integer nTotEl
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax         ! max no. batches
      integer iwhat(n_bmax)  ! what batches have been computed [out]
      integer n3c_dbl
      integer iwhat_max
      integer n_semi_bufs
      integer fd
c
      double precision CD_coef(nbf_cd,nfit)
      double precision i3c_ERI(n3c_dbl)
      double precision Ecoul1
      double precision Tvec(nbf_cd,nfit)
c
      Logical IOLGC, derfit, have_Tvec
c
      integer current_nbuf, n3c_b, n_batch_read
      double precision offset
      double precision length, q1
      double precision F(2)
      double precision fac, THRESHOLD, Pmat_max, ERI_est, Q, AZ, alpha,
     &                 rho_int
      integer nproc, me, 
     &        ischw1, icount, Nao2, n1, n2, ishc, ifirstc, ilastc,
     &        nshbfc, ishd, ifirstd,  ilastd, nshbfd, istart,
     &        ishp, ifirstp, ilastp, nshbfp, Nintegrals, nscr,
     &        max1e, nscr1e, j, next,  mscratch_2e3c,
     &        ifirst, ilast, nshbf, maxg, ish, lmiss, ifit
      integer iTvec, lTvec, lgvec, igvec, lERI, iERI, lcgtf, icgtf,
     &     lscr, iscr, iqatom, lqatom, iqshell_cd, lqshell_cd, 
     &     iqshell_lab, lqshell_lab, lqshtyp_cd, iqshtyp_cd,
     &     lTcgtf, iTcgtf,ninteg_api
      logical oprint, dowork, oprint_ioinfo
      integer nlast, ntasks
      integer ishc_be, ishd_be
      integer ndone, ipoint
      character*(nw_max_path_len) cdfit_filename
      double precision dabsmax
      integer nxtask
      integer          itype,nprimo,isphere
      external dabsmax,nxtask
      logical spherical
      integer LU,ierr
      character*255 errmsg
      
c     
c     Fit electronic charge density. The fitting coefficients are obtained by 
c     solving the normal equations for a standard linear least squares problem.
c      
      if (.not.derfit .and. nfit.ne.1) then
         call errquit('dft_fitcd: nfit must be 1 for regular fits',0,
     &       INPUT_ERR)
      endif
      oprint = util_print('coulomb fit', print_high)
      oprint_ioinfo = util_print('io info',print_debug)
      nproc  = ga_nnodes()
      me = ga_nodeid()
c      
      THRESHOLD = 10.d0**(-itol2e)
c      
c     Project the Coulomb potential onto each of the CD fitting functions.
c
c     Allocate scratch space.
c      
      spherical=bas_is_spherical(CD_bas_han)
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd*nfit,'Tvec',lTvec,iTvec))
     &  call errquit('dft_fitcd: cannot allocate Tvec',0, MA_ERR)
c
c     Initialize the T vector.
c
      if (have_Tvec) then
         call dcopy(nbf_cd*nfit, Tvec, 1, dbl_mb(iTvec), 1)
      else
         call dfill(nbf_cd*nfit, 0.D0, Dbl_MB(iTvec), 1)
      endif
c
c     If we already have the right side vectors for the fit, skip right
c     to the determination of the fit coefficients
c
      if (have_Tvec) goto 500
c      
c     Parallelize the loop over the products of AO shells.
c      
      ischw1 = 0
cedo      call ga_sync
cedo      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      icount = 1
      icount = max(1,icount)
c      
      Nao2 = nbf_ao_mxnbf_cn*nbf_ao_mxnbf_cn
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nfit,'gather vector',lgvec,igvec))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0,
     &       MA_ERR)
c
      call int_mem_2e3c(maxg, mscratch_2e3c)
      if(.not.MA_Push_Get(MT_Dbl,maxg*nshells_cd,'ERI',lERI,iERI))
     &  call errquit('dft_fitcd: cannot allocate ERI',0, MA_ERR)
      nscr =  mscratch_2e3c 
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0, MA_ERR)
c
      if(.not.incore)then
        iwhat_max = -1
        goto 25
      endif
c
c     Exhaust in-core 3c 2e- buffers first.
c
      ndone = 1
      ipoint = 1
      if(ndone.gt.n_batch)goto 25
c
c     check for semi-direct, if so, load 1st buffer
c
      n3c_b = 8*n3c_dbl
      current_nbuf = 0
c      write(6,*)' n_semi_bufs = ',n_semi_bufs
      if (n_semi_bufs.gt.0)then
c
c        check to make sure file is there
c
c         write(6,*)' fd = ',fd
         if (eaf_length(fd, length) .ne. 0) call errquit
     $      ('dft_fitcd: unable to determine 3c2e file length',0,
     &       DISK_ERR)
c         write(6,*)' length = ',length
         current_nbuf = current_nbuf + 1
         offset = (current_nbuf-1)*n3c_b
         if (eaf_read(fd, offset, i3c_eri, n3c_b) .ne. 0)call errquit
     &      ('dft_fitcd: 3c-2e integral buffer read failed',0, DISK_ERR)
         n_batch = nint(i3c_eri(n3c_dbl))
c         write(6,*)' n_batch = ',n_batch
c
c        temp stop
c
c         call errquit(' temp stop in dft_3cincor ',0)
      endif
c
      do 24 ishc = 1,nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han, ishc, ifirstc, ilastc)) 
     &       call errquit('Exiting in dft_fitcd',111, BASIS_ERR)
        nshbfc = ilastc-ifirstc+1
c
        do 23 ishd = 1,ishc
          dowork = iwhat(ndone).eq.((ishc*(ishc-1))/2 + ishd-1)
          if(dowork) then
c
c            write(LuOut,*)
c     &      'me,ishc,ishd,ndone,iwhat(ndone),ipoint', 
c     &       me,ishc,ishd,ndone,iwhat(ndone),ipoint
c            call util_flush(LuOut)
c
            if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, ilastd)) 
     &         call errquit('Exiting in dft_fitcd',1111, BASIS_ERR)
            nshbfd = ilastd-ifirstd+1
            Nao2 = nshbfd*nshbfc
c
c           Gather the needed density matrix elements.
c
            do ifit = 1, nfit
               call ga_get(g_dens(ifit), ifirstd, ilastd, ifirstc,
     &                     ilastc, dbl_mb(igvec+(ifit-1)*Nao2), nshbfd)
            enddo
c            
c           Determine the maximum density matrix element for this product of AO shells.
c
            Pmat_max = dabsmax(Nao2*nfit,DBL_MB(igvec))
c            
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c
            ERI_est = schwarz_shell(ishc,ishd)
c            
            if( Pmat_max*ERI_est.gt.THRESHOLD )then
c              
              ischw1 = ischw1 + nshells_cd
c              
              istart = iERI
              fac = 2.0d0
              if(ishc.eq.ishd) fac = 1.d0
c
c             Use dgemv for case of 1 fit (more efficient), else use dgemm
c
              if (nfit.eq.1) then
                 call dgemv('T',Nao2,nbf_cd,fac,i3c_ERI(ipoint),Nao2,
     &                Dbl_MB(igvec),1,1.d0,Dbl_MB(iTvec),1)
              else
                 call dgemm('T','N',nbf_cd,nfit,Nao2,fac,
     &                i3c_ERI(ipoint),Nao2,dbl_mb(igvec),nbf_cd,
     &                1d0,dbl_mb(iTvec),nbf_cd)
                 call errquit('not impl in dft_fitcd',0, UNKNOWN_ERR)
              endif
            end if
c
            ndone = ndone+1
c            write(6,*)' ndone, n_batch ',ndone,n_batch
            if(ndone.gt.n_batch)then
c
c             if not reading 3c2e integrals from disk - then done
c             start "direct"
c
              if(n_semi_bufs.eq.0)goto 25
c
c             Finished processing this batch of integrals, either get
c
              if (current_nbuf .eq. n_semi_bufs)goto 25
              current_nbuf = current_nbuf + 1
              offset = (current_nbuf-1)*n3c_b
              if (eaf_read(fd, offset, i3c_eri, n3c_b) .ne. 0)
     &           call errquit
     &           ('dft_fitcd: 3c-2e integral buffer read failed',0,
     &       DISK_ERR)
              n_batch_read = nint(i3c_eri(n3c_dbl))
c              write(6,*)' n_batch_read = ',n_batch_read
              n_batch = n_batch + n_batch_read
              ipoint = 1
            else
              ipoint = ipoint+Nao2*nbf_cd
            endif
c
c            write(LuOut,*)' Done with: ishc, ishd, n_batch, ndone: ',
c     &                                 ishc, ishd, n_batch, ndone
c            call util_flush(LuOut)
c
          end if
c
   23   continue
   24 continue
   25 continue
c
c     Finally, process any 3c 2e- integrals still needed in direct mode.
c
c     Find last ishc and ishd; add 1 to iwhat_max since starting at 0.
c
      nlast = iwhat_max + 1
c
c     Determine total number of shell pairs.
c
      ntasks = (nshells_ao*(nshells_ao+1))/2
c
      if (nlast.eq.ntasks)goto 241
c
      do n1 = 1, nshells_ao
        if ((n1*(n1-1))/2.gt.nlast)then
          ishc_be = n1 - 1
          goto 219
        endif
      enddo
      ishc_be = nshells_ao
  219 continue
      ishd_be = nlast - (ishc_be*(ishc_be-1))/2
c
      ishd_be = ishd_be + 1
      if(ishd_be.gt.ishc_be)then
        ishc_be = ishc_be + 1
        ishd_be = 1
      endif
c
c      call ga_sync
c      do n1 = 0, nproc
c        call ga_sync
c        if(me.eq.n1)then
c          write(LuOut,*)' me, ishc_be, ishd_be ', me, ishc_be, ishd_be
c          call util_flush(LuOut)
c        endif
c        call ga_sync
c      enddo
c      call ga_sync
c
      n1 = nxtask(nproc,icount)
      n2 = 0
c
      do 240 ishc = ishc_be, nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han, ishc, ifirstc, ilastc)) 
     &       call errquit('Exiting in dft_fitcd',111, BASIS_ERR)
        nshbfc = ilastc-ifirstc+1
c
        do 230 ishd = ishd_be, ishc
          dowork = n1.eq.n2
          if(dowork) then
c
c            write(LuOut,*)'me,ishc,ishd,n1,n2', me,ishc,ishd,n1,n2
c            call util_flush(LuOut)
c
            if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, ilastd)) 
     &         call errquit('Exiting in dft_fitcd',1111, BASIS_ERR)
            nshbfd = ilastd-ifirstd+1
            Nao2 = nshbfd*nshbfc
c
c           Gather the needed density matrix elements.
c
            do ifit = 1, nfit
               call ga_get(g_dens(ifit), ifirstd, ilastd, ifirstc,
     &              ilastc, dbl_mb(igvec+(ifit-1)*Nao2), nshbfd)
            enddo
c            
c           Determine the maximum density matrix element for this product of AO shells.
c
            Pmat_max = dabsmax(Nao2*nfit,DBL_MB(igvec))
c            
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c
            ERI_est = schwarz_shell(ishc,ishd)
c            
            if( Pmat_max*ERI_est.gt.THRESHOLD )then
c              
              ischw1 = ischw1 + nshells_cd
c              
              istart = iERI
              fac = 2.0d0
              if(ishc.eq.ishd) fac = 1.d0
              do 220 ishp = 1,nshells_cd
c                  
c               Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c                  
                if (.not. bas_cn2bfr( CD_bas_han,ishp,ifirstp,ilastp))
     &               call errquit('Exiting in fitcd.',5, BASIS_ERR)
                nshbfp = ilastp-ifirstp+1
c                  
                Nintegrals = nshbfp*Nao2
                Ninteg_api = maxg
c
c               check for use of symmetry
c
c     Must turn symmetry off for derivative fits since derivative
c     densities will not be totally symmetric in general
c
                if (oskel .and. .not.derfit) then
                   if (sym_shell(cd_bas_han, ishp, q1))then
                      call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,
     &                ishd,nscr,Dbl_MB(iscr),Ninteg_api,Dbl_MB(istart))
                      call dscal(Nintegrals, q1, Dbl_MB(istart), 1)
                   else
                      call dcopy(nintegrals,0d0,0,Dbl_MB(istart), 1)
                   endif
                else
                   call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,
     &             ishd,nscr,Dbl_MB(iscr),Ninteg_api,Dbl_MB(istart))
                endif                   
c
                istart = istart+Nintegrals
c
  220         continue
c               
c             Multiply the matrix of 3-ctr ERIs by the density matrix elements.
c               
              if (nfit.eq.1) then
                 call dgemv('T',Nao2,nbf_cd,fac,Dbl_MB(iERI),Nao2,
     &                Dbl_MB(igvec),1,1.d0,Dbl_MB(iTvec),1)
              else
                 call errquit('dft_fitcd: not impl',1, BASIS_ERR)
              endif
            end if
c
c            write(LuOut,*)' Done with: ishc,ishd; n1,n2 ',
c     &                                 ishc,ishd, n1,n2
c            call util_flush(LuOut)
c
            n1 = nxtask(nproc,icount)
c
          end if
c          
          n2 = n2 + 1
c          
  230   continue
c
c       Reset ishd_be to 1.
c
        ishd_be = 1
  240 continue
  241 continue
c
      n1 = nxtask(-nproc,icount)
c
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitcd: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(lERI))
     &  call errquit('dft_fitcd: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(lgvec))
     &  call errquit('dft_fitcd: cannot pop stack',0, MA_ERR)
c
      call ga_dgop(Msg_gop_Tvec,Dbl_MB(iTvec),nbf_cd*nfit,'+')
c
c
      if (oskel .and. .not.derfit) then
         do ifit = 1, nfit
            call sym_vec_symmetrize(geom, cd_bas_han,
     &           Dbl_MB(iTvec+(ifit-1)*nbf_cd))
         enddo
      endif
c
 500  continue
c      
c     Compute the CD expansion coefficients.
c     
c     Allocate scratch space.
c      
c
      call int_mem_1e(max1e,nscr1e) 
c
c     correct for last shell length
c     necessary for spherical when calling int_1cg
c
      max1e=nbf_cd
      if(spherical) then
        if(.not.bas_continfo(CD_bas_han,nshells_cd,
     &       itype,nprimo,nshbf,isphere))
     &       call errquit('Exiting in fitcd.',44, BASIS_ERR)
        nshbf=((itype+1)*(itype+2))/2
        max1e=nbf_cd+nshbf
      endif
      if(.not.MA_Push_Get(MT_Dbl,max1e,'cgtf',lcgtf,icgtf))
     &     call errquit('dft_fitcd: cannot allocate cgtf',0, MA_ERR)
      nscr = max(nscr1e, nbf_cd)
c
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0, MA_ERR)
c      
c     Integrate CD fitting functions in each shell
c
      call ga_sync
      do ish = 1,nshells_cd
        if (.not. bas_cn2bfr( CD_bas_han,ish,
     &                       ifirst,ilast))
     &      call errquit('Exiting in fitcd.',44, BASIS_ERR)
          if(.not.bas_continfo(CD_bas_han,ish,
     &         itype,nprimo,nshbf,isphere))
     &      call errquit('Exiting in fitcd.',44, BASIS_ERR)
          if(spherical) then
            nshbf=((itype+1)*(itype+2))/2
          else
            nshbf=ilast-ifirst+1
          endif
              
        call int_1cg( CD_bas_han, ish, nscr, Dbl_MB(iscr), 
     &              nshbf, Dbl_MB(icgtf+ifirst-1) )
      enddo 
c
c     Compute the orbital integrated electron density.
c
      if (derfit) then
         Q = 0d0
      else
         Q = Dble(nTotEl)
      endif
c      
      if (IOLGC)then
c         call gacrea(geom,CD_bas_han,
c     &       nbf_cd,nbf_cd,'CD 2cERInv',g_cdinv,'row')
*ga:1:0
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERInv', 
     &                       nbf_cd, 0, g_cdinv)) 
     &      call errquit('dft_fitcd: error creating g_cdinv',0, GA_ERR)
        if (me.eq.0)then
           lmiss = 0
           call dft_invio('CDI',g_cdinv,nbf_cd,'READ',lmiss)
        endif
        call ga_brdcst(2000, lmiss,mitob(1),0)
        if(lmiss.ne.0)call errquit
     &       (' dft_fitcd; dft_invio - abnormal read of CDI ',0, GA_ERR)
cedo        call ga_sync
        if(oprint_ioinfo)then
          call ga_print(g_cdinv)
        endif
      endif
#if 0
c
c     Compute the Lagrange Multiplier, alpha.
c     
      F(1) = 0.D0
      F(2) = 0.D0
      next = nxtask(nproc, 1)+1
      do j = 1,nbf_cd
        if (j.eq.next) then
          next = nxtask(nproc, 1)+1
          call get_col(g_cdinv,nbf_cd,j,Dbl_MB(iscr))
          AZ = ddot(nbf_cd,DBL_mb(iscr),1,DBL_MB(icgtf),1)
          F(1) = F(1)+AZ*DBL_MB(iTvec+j-1)
          F(2) = F(2)+AZ*DBL_MB(icgtf+j-1)
        endif
      enddo
      next = nxtask(-nproc, 1)
cedo      call ga_sync
      call ga_dgop(Msg_Lagrange,F,2,'+')
cedo      call ga_sync
      alpha = (Q - F(1))/F(2)
c
c     Qfrac = (Q - F(1))/Q
c
c     Multiply: S**(-1) * (alpha*N + T)
c      
c     Initialize the CD_coef vector.
      
      call dfill(nbf_cd,0.D0,CD_coef,1)
c
      next = nxtask(nproc, 1)+1
      do j = 1,nbf_cd
        if (j.eq.next) then
          next = nxtask(nproc, 1) + 1
          AZ = DBL_MB(iTvec+j-1) + alpha*DBL_MB(icgtf+j-1)
          call get_col(g_cdinv,nbf_cd,j,Dbl_MB(iscr))
          call daxpy(nbf_cd,AZ,DBL_MB(iscr),1,CD_Coef,1)
        endif
      enddo
      next = nxtask(-nproc, 1)
      call ga_dgop(Msg_CD_coef,CD_coef,nbf_cd,'+')
cedo      call ga_sync
#else
c
c     Multiply T and cgtf by S**(-1) and store
c
c     CD_coef is initially used to store S**(-1) * T, the values
c     the coefficients would have without the charge constraint
c
      call dfill(nbf_cd*nfit,0.D0,CD_coef,1)
c
c     Need another temp array to store S**(-1) * cgtf
c
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd,'Tcgtf',lTcgtf,iTcgtf))
     &     call errquit('dft_fitcd: cannot allocate Tcgtf',0, MA_ERR)
      call dfill(nbf_cd,0.D0,dbl_mb(iTcgtf),1)
c
      next = nxtask(nproc, 1)+1
      do j = 1,nbf_cd
        if (j.eq.next) then
          next = nxtask(nproc, 1)+1
          call get_col(g_cdinv,nbf_cd,j,Dbl_MB(iscr))
          AZ = ddot(nbf_cd,DBL_mb(iscr),1,DBL_MB(icgtf),1)
          dbl_mb(iTcgtf+j-1) = dbl_mb(iTcgtf+j-1) + AZ
          do ifit = 1, nfit
             AZ = ddot(nbf_cd,DBL_mb(iscr),1,
     &            DBL_MB(iTvec+(ifit-1)*nbf_cd),1)
             CD_coef(j,ifit) = CD_coef(j,ifit) + AZ
          enddo
        endif
      enddo
      next = nxtask(-nproc, 1)
cedo      call ga_sync
      call ga_dgop(Msg_CD_coef,CD_coef,nbf_cd*nfit,'+')
c!!! What Msg should be here? !!!
      call ga_dgop(Msg_Tcgtf,dbl_mb(iTcgtf),nbf_cd,'+')
cedo      call ga_sync
c
c     Apply constraint correction to fits
c
      F(2) = ddot(nbf_cd,DBL_mb(icgtf),1,DBL_MB(iTcgtf),1)
c
      do ifit = 1, nfit
c
c        Compute Lagrange multiplier, alpha
c
         F(1) = ddot(nbf_cd,DBL_mb(icgtf),1,CD_coef(1,ifit),1)
            alpha = (Q - F(1))/F(2)
c         write(*,*)'F alpha',f(1),f(2),alpha
c
c        Correct unconstrained coefficients by alpha*S**(-1)*cgtf
c
         do j = 1,nbf_cd
            CD_coef(j,ifit) = CD_coef(j,ifit)
     &                      + alpha*DBL_MB(iTcgtf+j-1)
         enddo
      enddo
c
      if (.not.ma_pop_stack(lTcgtf))
     &   call errquit('dft_fitcd: cannot pop stack',0, MA_ERR)
#endif
      if(me.eq.0) then
         if (derfit) then
            call util_file_name('cdfit_deriv', .true., .false.,
     &           cdfit_filename)
         else
            call util_file_name('cdfit', .true., .false.,
     &           cdfit_filename)
         endif
         LU=79
         ierr=eaf_open(cdfit_filename(1:inp_strlen(cdfit_filename)),
     ,        eaf_rw, LU)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('dft_fitcd: eaf_open failed', ierr, GA_ERR)
        endif
         do ifit = 1, nfit
            offset=(ifit-1)*8d0*nbf_cd
            ierr = eaf_write(LU, offset, cd_coef(1,ifit),
     .           8*nbf_cd)
            if (ierr .ne. 0) then
               call eaf_errmsg(ierr, errmsg)
               write(LuOut,*) ' IO offset ', offset
               write(LuOut,*) ' IO error message ',
     $              errmsg(1:inp_strlen(errmsg))
               call errquit('dft_fitcd: write failed',ifit, INPUT_ERR)
               endif
         enddo
         ierr=eaf_close(LU)
         if (ierr .ne. 0) then
            write(luout,*) ga_nodeid(),' closing FD =',LU,
     ,           cdfit_filename(1:inp_strlen(cdfit_filename)) 
            call eaf_errmsg(ierr, errmsg)
            write(LuOut,*) ' IO error message ',
     $           errmsg(1:inp_strlen(errmsg))
            call util_flush(luout)
            call errquit('dft_fitcd: eaf_close failed',0, DISK_ERR)
         endif
      endif
c      
c     Compute the integrated electron density from the fit.
c      
      do ifit = 1, nfit
         rho_int = ddot(nbf_cd,CD_coef(1,ifit),1,Dbl_MB(icgtf),1)
         if (dabs(rho_int-Q).gt.0.001d0)then
            write(*,*)' rho_int check 1, density',ifit,' = ',rho_int
            call errquit('dft_fitcd: problem w/fitted density',911,
     &       UNKNOWN_ERR)
c     write(*,*)' Dbl_MB(icgtf) = ',(Dbl_MB(icgtf+i),i=0,nbf_cd-1)
c     write(*,*)' CD_coef = ',(CD_coef(i,ifit),i=1,nbf_cd)
         endif
      enddo
cedo      call ga_sync
c      if (derfit) then
c         write(*,*)'Fit coefficients:'
c         call output(cd_coef, 1, nbf_cd, 1, nfit, nbf_cd, nfit, 1)
c      endif
c
c     Get mulliken population analysis of charge density fit
c
      if (.not.derfit) then
c
      if (.not. geom_ncent(geom, natoms))
     &     call errquit('dft_fitcd: geom_ncent failed',73,
     &       GEOM_ERR)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'qatom',lqatom,iqatom))
     &   call errquit('dft_fitcd: cannot allocate qatom',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nshells_cd,'qshell_cd',lqshell_cd, 
     &   iqshell_cd))
     &   call errquit('dft_fitcd: cannot allocate qshell_cd',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Int,nshells_cd*natoms,'qshtyp_cd',
     &   lqshtyp_cd,iqshtyp_cd))
     &   call errquit('dft_fitcd: cannot allocate qshtyp_cd',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,nshells_cd,'qshell_lab',lqshell_lab,
     &   iqshell_lab))
     &   call errquit('dft_fitcd: cannot allocate qshell_cd',0, MA_ERR)
c
      call mull_pop_fit(CD_bas_han, natoms, nshells_cd, nbf_cd, 
     &                  dbl_mb(iqatom), dbl_mb(iqshell_cd), cd_coef, 
     &                  DBL_MB(icgtf), int_mb(iqshtyp_cd), 
     &                  Byte_MB(iqshell_lab))
      if (.not.ma_pop_stack(lqshell_lab))
     &   call errquit('dft_fitcd: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lqshtyp_cd))
     &   call errquit('dft_fitcd: cannot pop stack',1, MA_ERR)
      if (.not.ma_pop_stack(lqshell_cd))
     &   call errquit('dft_fitcd: cannot pop stack',2, MA_ERR)
      if (.not.ma_pop_stack(lqatom))
     &   call errquit('dft_fitcd: cannot pop stack',3, MA_ERR)
c
      endif                     ! .not.derfit
c
c      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 7'
c      if(IOLGC) call gadest(g_cdinv)
      if (IOLGC)then
         if (.not. ga_destroy(g_cdinv)) call errquit
     &      ('dft_fitcd: could not destroy g_cdinv', 0, GA_ERR)
      endif         
c
c     calculation of Ecoul1
c
      if (.not.derfit) then
         Ecoul1 = ddot(nbf_cd,DBL_mb(iTvec),1,CD_coef,1)
      endif
c      
c     Compute the integrated electron density from the fit.
c      
      do ifit = 1, nfit
         rho_int = ddot(nbf_cd,CD_coef(1,ifit),1,Dbl_MB(icgtf),1)
         if (dabs(rho_int-Q).gt.0.001d0) then
            write(*,*)' rho_int check 2, density',ifit,' = ',rho_int
            call errquit('dft_fitcd: problem w/fitted density',911,
     &       UNKNOWN_ERR)
c     write(*,*)' Dbl_MB(icgtf) = ',(Dbl_MB(icgtf+i),i=0,nbf_cd-1)
c     write(*,*)' CD_coef = ',(CD_coef(i,ifit),i=1,nbf_cd)
         endif
      enddo
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitcd: cannot pop stack',4, MA_ERR)
      if(.not.ma_pop_stack(lcgtf))
     &  call errquit('dft_fitcd: cannot pop stack',5, MA_ERR)
      if(.not.ma_pop_stack(lTvec))
     &  call errquit('dft_fitcd: cannot pop stack',6, MA_ERR)
      
      
      if(oprint.and.me.eq.0.and..not.derfit )then
        write(LuOut,*)
        write(LuOut,'(a,f14.6,a,1pe13.6)')
     &       '   CD fit integrated density:',rho_int,
     &       '    Lagrange Mult.:  ',alpha
        write(LuOut,*)
        call util_flush(LuOut)
      end if
      
      return
      end
      Subroutine mull_pop_fit(basis, natoms, nshells, nbf_cd,
     &     qatom, qshell, cd_coef, cgtf, shell_types, shell_labels)
      Implicit none
#include "errquit.fh"
      integer basis
      integer natoms, nshells, shell_types(natoms,nshells)
      integer ifirst, ilast, nbf_cd
      integer ish1, ish2, ish, iat, i, iend
      integer me, nproc, it, itype, nprim, ncont, icart
      double precision qatom(natoms), qshell(nshells)
      double precision cd_coef(nbf_cd), cgtf(nbf_cd)
      character*1 shell_labels(nshells), ang_mom_label(11)
      logical oprint_mull_fit
c
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "util.fh"
c     
      data ang_mom_label /'S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 
     &                    'K', 'L', 'X'/
      oprint_mull_fit = util_print('mulliken of fitted density',
     &                              print_high)
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      call dfill(natoms,0.D0,qatom,1)
      call dfill(nshells,0.D0,qshell,1)
      call ifill(nshells*natoms,0,shell_types,1)
c
      do iat = me+1, natoms, nproc
         if (.not.bas_ce2cnr(basis, iat, ish1, ish2))
     &      call errquit(' exiting in mull_pop_fit',4, BASIS_ERR)
         do ish = ish1, ish2
            if (.not. bas_cn2bfr(basis, ish, ifirst, ilast))
     &         call errquit(' exiting in mull_pop_fit.',5, BASIS_ERR)
            if (.not. bas_continfo(basis, ish, itype, nprim,
     &                             ncont, icart))
     &         call errquit(' exiting in mull_pop_fit.',5, BASIS_ERR)
            it = ilast - ifirst + 1
c            if (itype.le.9)then
c               shell_labels(ish) = ang_mom_label(itype+1)
c            else
c               shell_labels(ish) = ang_mom_label(11)
c            endif
            shell_types(iat,ish) = itype
            qshell(ish) = ddot(it,cd_coef(ifirst),1,cgtf(ifirst),1)
            qatom(iat) = qatom(iat) + qshell(ish)
         enddo
      enddo
cedo      call ga_sync
      call ga_dgop(Msg_Mull1, qatom, natoms,'+')
      call ga_dgop(Msg_Mull2, qshell, nshells,'+')
      call ga_igop(Msg_Mull1, shell_types, nshells*natoms,'+')
      if (me.eq.0.and.oprint_mull_fit) then
         call dft_header(
     &       'Fitted Density - Mulliken Population Analysis')
         write(LuOut,10)
         do iat = 1, natoms
            if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &         call errquit(' exiting in mull_pop',4, BASIS_ERR)
            do i = 1, nshells
               itype = shell_types(iat,i)
               if (itype.le.9)then
                  shell_labels(i) = ang_mom_label(itype+1)
               else
                  shell_labels(i) = ang_mom_label(11)
               endif
            enddo
            do ish = ish1, ish2, 10
               iend = min(ish+9,ish2)
               write(LuOut,20)(shell_labels(i),i=ish,iend)
               if (ish.eq.ish1)then
                  write(LuOut,30)iat,qatom(iat),(qshell(i),i=ish,iend)
               else
                  write(LuOut,35)(qshell(i),i=ish,iend)
               endif
            enddo
         enddo
      endif
      call ga_sync
   10 format(/' Atom   Charge   Shell Charges')
   20 format(15x,10(3x,a1,2x))
   30 format(1x,i4,1x,f8.2,1x,10(f6.2))
   35 format(15x,10(f6.2))
c
      return
      end
