      Subroutine dft_fitcd(nfit,CD_coef, i3c_ERI, Ecoul1, 
     &                     g_dens, nTotEl, n_batch, n_bmax, iwhat,
     &                     n3c_dbl, IOLGC,
     &                     natoms,
     &                     derfit,  Tvec, have_Tvec)
*
* $Id$
*
      implicit none
#include "errquit.fh"
c
#include "bas.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
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
      integer n3c_dbl        ! max no. ints [in]
      integer fd
c
      double precision CD_coef(nbf_cd,nfit)
      double precision i3c_ERI(n3c_dbl)
      double precision Ecoul1
      double precision Tvec(nbf_cd,nfit)
c
      Logical IOLGC, derfit, have_Tvec
c
      double precision offset
      double precision F(2)
      double precision THRESHOLD,  Q(nfit), AZ, alpha,
     &                 rho_int
      integer nproc, me, 
     &        nscr,max1e, nscr1e, j,   
     &        ifirst, ilast, nshbf, ish, lmiss, ifit, jfit
      integer iTvec, lTvec, lcgtf, icgtf,
     &     lscr, iscr, iqatom, lqatom, iqshell_cd, lqshell_cd, 
     &     iqshell_lab, lqshell_lab, lqshtyp_cd, iqshtyp_cd,
     &     lTcgtf, iTcgtf
      logical oprint, oprint_ioinfo
      character*(nw_max_path_len) cdfit_filename
      double precision dabsmax
      integer nxtask
      integer          itype,nprimo,isphere
      external dabsmax,nxtask
      logical spherical
      integer LU,ierr,ilo,ihi,jlo,jhi,nnii
      integer adrc,ldc,iptr,intdum
      character*255 errmsg
      double precision ddot
      external ddot
      
c     
c     Fit electronic charge density. The fitting coefficients are obtained by 
c     solving the normal equations for a standard linear least squares problem.
c      
      if (.not.derfit .and. nfit.gt.2) then
         call errquit('dft_fitcd: nfit must be <3 for regular fits',0,
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
         call dcopy(nbf_cd*nfit, 0.D0, 0, Dbl_MB(iTvec), 1)
      endif
c
c     If we already have the right side vectors for the fit, skip right
c     to the determination of the fit coefficients
c
      if (.not.have_Tvec) 
     T     call dft_getcdt(.false.,intdum,
     N     natoms, nfit, derfit, threshold, 
     C     .false., ! do not compute 3c ints
     N     n_batch, n_bmax, iwhat,
     &     n3c_dbl, i3c_eri,
     T     dbl_mb(itvec),g_dens)

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

      ! Only perform this step if the fit is being
      ! constrained
      if (.not.freecdfit) then
        do ish = 1,nshells_cd
          if (.not. bas_cn2bfr( CD_bas_han,ish,
     &                         ifirst,ilast))
     &        call errquit('Exiting in fitcd.',44, BASIS_ERR)
            if(.not.bas_continfo(CD_bas_han,ish,
     &           itype,nprimo,nshbf,isphere))
     &        call errquit('Exiting in fitcd.',44, BASIS_ERR)
            if(spherical) then
              nshbf=((itype+1)*(itype+2))/2
            else
              nshbf=ilast-ifirst+1
            endif
              
          call int_1cg( CD_bas_han, ish, nscr, Dbl_MB(iscr), 
     &                nshbf, Dbl_MB(icgtf+ifirst-1) )
        enddo 
      endif
c
c     Compute the orbital integrated electron density.
c
      if (derfit) then
         Q = 0d0
      else
         Q(1) = Dble(nTotEl)
         if (nfit.gt.1) then
           do ifit=1,nfit
             Q(ifit) = dble(noc(ifit))
           enddo
         endif
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
        call ga_brdcst(2000, lmiss,MA_sizeof(MT_INT,1,MT_BYTE),0)
        if(lmiss.ne.0)call errquit
     &       (' dft_fitcd; dft_invio - abnormal read of CDI ',0, GA_ERR)
cedo        call ga_sync
        if(oprint_ioinfo)then
          call ga_print(g_cdinv)
        endif
      endif
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
       call dcopy(nbf_cd,0.D0,0,dbl_mb(iTcgtf),1)
c
      call ga_sync()
      call  ga_distribution(g_cdinv,me,ilo,ihi,jlo,jhi)
      if (ilo.gt.0) then
      call ga_access(g_cdinv, ilo, ihi, jlo, jhi, adrc, ldc)
         nnii=ihi-ilo+1
         do j = jlo,jhi
            iptr=adrc+(j-jlo)*ldc
            AZ = ddot(nnii,DBL_mb(iptr),1,DBL_MB(icgtf+ilo-1),1)
            dbl_mb(iTcgtf+j-1) = dbl_mb(iTcgtf+j-1) + AZ
            do ifit = 1, nfit
               AZ = ddot(nnii,DBL_mb(iptr),1,
     &              DBL_MB(iTvec+(ifit-1)*nbf_cd+ilo-1),1)
               CD_coef(j,ifit) = CD_coef(j,ifit) + AZ
          enddo
       enddo
      call ga_release(g_cdinv, ilo, ihi, jlo, jhi)
      endif! own patch
      call ga_mask_sync(.true.,.false.)
      call ga_dgop(Msg_CD_coef,CD_coef,nbf_cd*nfit,'+')
      call ga_mask_sync(.false.,.true.)
      call ga_dgop(Msg_Tcgtf,dbl_mb(iTcgtf),nbf_cd,'+')
c
c     Apply constraint correction to fits
c
      F(2) = ddot(nbf_cd,DBL_mb(icgtf),1,DBL_MB(iTcgtf),1)
c
      if (.not.freecdfit) then
        do ifit = 1, nfit
c
c        Compute Lagrange multiplier, alpha
c
          F(1) = ddot(nbf_cd,DBL_mb(icgtf),1,CD_coef(1,ifit),1)
          alpha = (Q(ifit) - F(1))/F(2)
c
c        Correct unconstrained coefficients by alpha*S**(-1)*cgtf
c
         do j = 1,nbf_cd
            CD_coef(j,ifit) = CD_coef(j,ifit)
     &                      + alpha*DBL_MB(iTcgtf+j-1)
         enddo
        enddo
      endif
c
      if (.not.ma_pop_stack(lTcgtf))
     &   call errquit('dft_fitcd: cannot pop stack',0, MA_ERR)

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
      if (.not.freecdfit) then
        do ifit = 1, nfit
           rho_int = ddot(nbf_cd,CD_coef(1,ifit),1,Dbl_MB(icgtf),1)
           if (dabs(rho_int-Q(ifit)).gt.0.001d0)then
             write(*,*)' rho_int check 1, density',ifit,' = ',rho_int
             write(*,*)' CD cgtf',dabsmax(nbf_cd,cd_coef),
     D           dabsmax(nbf_cd,Dbl_MB(icgtf))
             call errquit('dft_fitcd: problem w/fitted density',911,
     &       UNKNOWN_ERR)
c            write(*,*)' Dbl_MB(icgtf) = ',(Dbl_MB(icgtf+i),i=0,nbf_cd-1)
c            write(*,*)' CD_coef = ',(CD_coef(i,ifit),i=1,nbf_cd)
         endif
        enddo
      endif
cedo      call ga_sync
c      if (derfit) then
c         write(*,*)'Fit coefficients:'
c         call output(cd_coef, 1, nbf_cd, 1, nfit, nbf_cd, nfit, 1)
c      endif
c
c     Get mulliken population analysis of charge density fit
c
      if (.not.derfit.and..false.) then
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
      if (IOLGC)then
         if (.not. ga_destroy(g_cdinv)) call errquit
     &      ('dft_fitcd: could not destroy g_cdinv', 0, GA_ERR)
      endif         
c
c     calculation of Ecoul1
c
      Ecoul1 = 0.0d0
      if (.not.derfit) then
         do ifit=1,nfit
           do jfit=1,nfit
             Ecoul1 = Ecoul1+
     $                ddot(nbf_cd,DBL_mb(iTvec+(ifit-1)*nbf_cd),1,
     $                CD_coef(1,jfit),1)
           enddo
         enddo
      endif
c      
c     Compute the integrated electron density from the fit.
c
      if (.not.freecdfit) then
      do ifit = 1, nfit
         rho_int = ddot(nbf_cd,CD_coef(1,ifit),1,Dbl_MB(icgtf),1)
         if (dabs(rho_int-Q(ifit)).gt.0.001d0) then
            write(*,*)' rho_int check 2, density',ifit,' = ',rho_int
            call errquit('dft_fitcd: problem w/fitted density',911,
     &       UNKNOWN_ERR)
c     write(*,*)' Dbl_MB(icgtf) = ',(Dbl_MB(icgtf+i),i=0,nbf_cd-1)
c     write(*,*)' CD_coef = ',(CD_coef(i,ifit),i=1,nbf_cd)
         endif
      if(oprint.and.me.eq.0.and..not.derfit )then
        write(LuOut,*)
        write(LuOut,'(a,f14.6,a,1pe13.6)')
     &       '   CD fit integrated density:',rho_int,
     &       '    Lagrange Mult.:  ',alpha
        write(LuOut,*)
        call util_flush(LuOut)
      end if
      enddo
      endif
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitcd: cannot pop stack',4, MA_ERR)
      if(.not.ma_pop_stack(lcgtf))
     &  call errquit('dft_fitcd: cannot pop stack',5, MA_ERR)
      if(.not.ma_pop_stack(lTvec))
     &  call errquit('dft_fitcd: cannot pop stack',6, MA_ERR)
      
      
      
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
      double precision ddot
      external ddot
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
      if (oprint_mull_fit) then
         call ga_dgop(Msg_Mull1, qatom, natoms,'+')
         call ga_dgop(Msg_Mull2, qshell, nshells,'+')
         call ga_igop(Msg_Mull1, shell_types, nshells*natoms,'+')
      endif
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
      subroutine dft_copyfrom(
     I     ld_dm,
     O     row_ld, col_ld,
     I     dm, 
     O     dmsub,  
     P     row_ptr,col_ptr)
      implicit none
      integer ld_dm                      ! [in]
      integer row_ld, col_ld             ! [in]
      double precision dm(ld_dm,*)       ! [in]
      double precision dmsub(1:row_ld,1:col_ld) ! [out]
      integer row_ptr,col_ptr            ! [in]
c
      integer j
c
      do j=1, col_ld
            
        call dcopy(row_ld,dm(row_ptr+1,col_ptr+j),1,dmsub(1,j),1)
      enddo
      return
      end
      subroutine dft_copyto(
     I     ld_dm,
     O     row_ld, col_ld,
     I     dm, 
     O     dmsub,  
     P     row_ptr,col_ptr)
      implicit none
      integer ld_dm                      ! [in]
      integer row_ld, col_ld             ! [in]
      double precision dm(ld_dm,*)       ! [out]
      double precision dmsub(1:row_ld,1:col_ld) ! [in]
      integer row_ptr,col_ptr            ! [in]
c
      integer j
c
      do j=1, col_ld
        call dcopy(row_ld,dmsub(1,j),1,dm(row_ptr+1,col_ptr+j),1)
      enddo
      return
      end
