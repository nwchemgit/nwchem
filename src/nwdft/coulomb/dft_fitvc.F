      Subroutine dft_fitvc(CD_coef, i3c_ERI, Ecoul2,  g_vc, 
     &                     n_batch, n_bmax, iwhat, n3c_dbl,
     &                     IOLGC,
     &                     dermat, nmat, do_init)
c     $Id$
      implicit none
#include "errquit.fh"
c      
#include "bas.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "schwarz.fh"
#include "util.fh"
#include "eaf.fh"
#include "cdft.fh"
#include "sym.fh"
#include "inp.fh"
c
      integer nmat           ! number of J matrices to make [input]
      integer g_vc(nmat)     ! GA handle for V Coulomb [output]
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax         ! max no. batches
      integer iwhat(n_bmax)  ! what batches have been computed [out]
      integer n3c_dbl
      integer LU,ierr
      character*(nw_max_path_len) cdfit_filename
      character*255 errmsg
c
      double precision CD_coef(nbf_cd,2*nmat)
      double precision i3c_ERI(n3c_dbl)
      double precision Ecoul2 ! [output]
c
      Logical IOLGC, dermat, do_init
c     
c     Compute the fitted Coulomb potential. The fitting coefficients are
c     obtained by solving the normal equations for a standard linear least
c     squares problem.
c
      integer me,nproc,i1,ishp,ishbf,lshbf,lmiss,ischw1,
     ,     icount,Nao2,i,nERI,maxg,mscratch_2e3c,nscr,
     ,     n1,ishc,ifirstc,ilastc,nshbfc,ishd,ifirstd,ilastd,
     ,     nshbfd,nfunz,ifirstp,ilastp,
     ,     nshbf_cd,imat,ninteg_api
      integer lcd_shell,icd_shell,lsvec,isvec,
     ,     iERI,lERI,lscr,iscr
      double precision THRESHOLD,CD_coef_max,shmax,
     ,     Vcoul_est,ERI_est
      double precision qsym,qscale,q1,kp_max
      integer nxtask
      logical dowork,doshell
      logical oprint_ioinfo,use_nxtask
      integer ipoint
      integer ilo, ihi, jlo, jhi, adrc, ldc,nnii,iptr,j
      integer atom_c,atom_d
      integer sh_lo_c,sh_hi_c, sh_lo_d, sh_hi_d, ishd_0, ishd_1
      integer bf_lo_c,bf_hi_c, bf_lo_d, bf_hi_d
      integer nbf_at_c,nbf_at_d
      integer k_p,l_p,icount_d,icount_c, n_batch_done
      integer k_at,l_at,atom_c_in,atom_d_in
      integer nbffac
      logical v_nonzero
      external nxtask
      double precision ddot
      external ddot
      double precision offset
c
      if(dermat) call errquit(
     C     'fitvc: dermat not coded yet',0,0)

      oprint_ioinfo = util_print('io info',print_debug)
c
      THRESHOLD=10.d0**(-itol2e)
c      
      me=ga_nodeid()
      nproc  = ga_nnodes()
      do imat=1,nmat
        call ga_zero(g_vc(imat))
      enddo
c      
c     Compute the fitted Vcoul.
c      
c     Allocate scratch space.
c      
      if(.not.MA_Push_Get(MT_Dbl,nshells_cd,'CD shell max',
     &  lcd_shell,icd_shell))
     &  call errquit('dft_fitvc: cannot allocate CD shell max',0,0)
c      
c      
c     Determine the largest magnitude of the CD expansion coefficients in
c     each shell.
c      
c     Note that for multiple matrices the max over all fitting coefficient
c     sets is taken
c
      CD_coef_max = 0.D0
      i1 = 1
      do  ishp = 1,nshells_cd
        if( .not. bas_cn2bfr(CD_bas_han,ishp,ishbf,lshbf) )then
          call errquit('Exiting in fitvcoul.',4,0)
        end if
        shmax = 0.D0
        do i = ishbf,lshbf
           do imat = 1, nmat
              shmax = max(abs(CD_coef(i,imat)), shmax)
           enddo
        enddo
        Dbl_MB((icd_shell-1)+ishp) = shmax
        CD_coef_max = max(shmax, CD_coef_max)
      enddo
      if (IOLGC)then
*ga:1:0
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERInv', 
     &                       nbf_cd, 0, g_2ceri)) 
     &      call errquit('dft_fitvc: error creating g_2ceri',0,0)
        if (me.eq.0)then
           lmiss = 0
           call dft_invio('CD',g_2ceri,nbf_cd,'READ',lmiss)
        endif
        call ga_brdcst(1999, lmiss,MA_sizeof(MT_INT,1,MT_BYTE),0)
        if(lmiss.ne.0)call errquit
     &     (' dft_fitvc; dft_invio - abnormal read of CD ',0,0)
        if(oprint_ioinfo)then
          call ga_print(g_2ceri)
        endif
      endif
c      
c     Parallelize the loop over the products of atoms
c      
c
      Nao2 = nbf_ao_mxnbf_cn**2
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nmat,'isvec',lsvec,isvec))
     &  call errquit('dft_fitvc: cannot allocate isvec',0,0)
      Nao2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nmat,'pp array',l_p,k_p))
     &  call errquit('dft_fitcd: cannot allocate p',0,  MA_ERR)
c
      call int_mem_2e3c(maxg, mscratch_2e3c)

      nERI    = maxg
      if(.not.MA_Push_Get(MT_Dbl,nERI,'ERI',lERI,iERI))
     &  call errquit('dft_fitvc: cannot allocate ERI',0,0)
      nscr    = mscratch_2e3c
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitvc: cannot allocate scr',0,0)
c
      if (.not.do_init) call errquit(
     C     ' dft_fitvc: do_init has to be true',0,0)
      do imat=1,nmat
        call ga_zero(g_vc(imat))
      enddo
      if(.not.MA_Push_Get(MT_int,ncenters,'at list',l_at,k_at))
     &  call errquit('dft_fitvc: cannot allocate at list',0, MA_ERR)
c     sort atoms list based on decreasing atomic weight
      call dft_sortatbyz(geom,ncenters,int_mb(k_at))
c
      n1=0
      use_nxtask=.not.incore
      if(use_nxtask) n1 = nxtask(nproc,1)
      icount = 0
      ipoint=1
c
      n_batch_done = 1
      do atom_c_in = 1,ncenters
         atom_c=int_mb(k_at+atom_c_in-1)
         if (.not. bas_ce2cnr( AO_bas_han, atom_c, sh_lo_c, sh_hi_c)) 
     &        call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
         if (.not. bas_ce2bfr( AO_bas_han, atom_c, bf_lo_c, bf_hi_c)) 
     &        call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
         nbf_at_c=bf_hi_c - bf_lo_c + 1
c     
         do atom_d_in= 1,atom_c_in
            atom_d=int_mb(k_at+atom_d_in-1)
            if(use_nxtask) then
               dowork=n1.eq.icount
            else
               dowork=iwhat(n_batch_done).eq.
     A              ((atom_c_in*(atom_c_in-1))/2 + atom_d_in-1)
            endif
            if(dowork) then
               v_nonzero=.false.
               if (.not.bas_ce2bfr(AO_bas_han,atom_d,bf_lo_d, bf_hi_d)) 
     &              call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
               nbf_at_d=bf_hi_d - bf_lo_d + 1
               if (.not.bas_ce2cnr(AO_bas_han,atom_d, sh_lo_d, sh_hi_d)) 
     &              call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
               icount_c=0
czero k_p
               call dcopy(nmat*nbf_at_d*nbf_at_c,0d0,0,dbl_mb(k_p),1)
               do ishc=sh_lo_c,sh_hi_c

                  if( .not. bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &                 call errquit('Exiting in fitvcoul.',3,0)
                  nshbfc=ilastc-ifirstc+1
                  
                  ishd_0 = sh_lo_d
                  ishd_1 = sh_hi_d
                  if (atom_c.eq.atom_d) ishd_1=ishc

                  icount_d=0
                  do ishd =  ishd_0,  ishd_1
                     
                     if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, 
     I                    ilastd))
     &                    call errquit('Exiting in fitvcoul.',4,0)
                     nshbfd=ilastd-ifirstd+1
                     Nao2 = nshbfc*nshbfd
c           
c     Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c           
                     ERI_est=schwarz_shell(ishc,ishd)
                     if(incore) CD_coef_max=1d0
c     schwartz with CDmax=1 for incore ints reading
                     if( ERI_est*CD_coef_max.gt.THRESHOLD )then

                        v_nonzero=.true.
                        qscale=1d0
                        if(ishc.ne.ishd) qscale=2d0
                        if(incore) then
                           if(ipoint+nao2*nbf_cd-1.gt.n3c_dbl) then
                              write(6,*) ' ip n3c ',
     A                             ipoint+nao2*nbf_cd-1,n3c_dbl
                              call errquit('dft_fitvc: out of i3c ',0,0)
                           endif
                           do imat=1,nmat
                             call dgemv('N',Nao2,nbf_cd,qscale,
     &                            i3c_eri(ipoint),Nao2,
     C                            CD_coef(1,imat),1,
     D                            0.D0,Dbl_MB(isvec+(imat-1)*nao2),1)
                           enddo
                           ipoint = ipoint + nao2*nbf_cd
                        else
                           call dcopy(Nao2*nmat,0.D0,0,dbl_mb(isvec),1)
                           nfunz=0
                           do ishp = 1,nshells_cd
c     
                              if( .not.bas_cn2bfr(CD_bas_han, ishp, 
     I                             ifirstp, ilastp))
     &                             call errquit('fitvc',0,BASIS_ERR)
                              nshbf_cd=ilastp-ifirstp+1
                              Vcoul_est=dbl_mb(icd_shell+ishp-1)*ERI_est
                              if( Vcoul_est*nshbf_cd.gt.THRESHOLD )then
                                 Ninteg_api=maxg
c                  
c     Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c     
c     check for use of symmetry
                                 qsym=1d0
                                 doshell=.true.
                                 if(oskel) doshell =
     Z                                sym_shell(cd_bas_han,ishp,qsym)

                                 if(doshell) then
                                    call int_2e3c(CD_bas_han,ishp,
     A                                   AO_bas_han,ishc,ishd,nscr,
     D                                   Dbl_MB(iscr),Ninteg_api,
     &                                   Dbl_MB(iERI))
c     
c     Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c     
                                    q1=qscale*qsym
                                    do imat=1,nmat
                                    call dgemv('N',Nao2,nshbf_cd,q1,
     &                                   Dbl_MB(iERI),Nao2,
     C                                   CD_coef(1+nfunz,imat),1,1.0D0,
     &                                   Dbl_MB(isvec+(imat-1)*nao2),1)
                                    enddo
                                 endif
                              endif ! schw on VCoul_est
c
                              nfunz=nfunz+nshbf_cd
c               
                           enddo ! ishp
                        endif   ! incore vs direct
c
c     Scatter the computed values of Vcoul.
c           
                          do imat=1,nmat
                           call dft_copyto(
     I                          nbf_at_d,nshbfd,nshbfc,
     I                          dbl_mb(k_p+(imat-1)*nbf_at_d*nbf_at_c),
     D                          dbl_mb(isvec+(imat-1)*nao2),
     A                          icount_d,icount_c)
                          enddo
                        endif   ! schwartz with cdmax

                     icount_d=icount_d+nshbfd
                  enddo         ! ishd
c     add extra bit to ipoint if ishd_1 is truncated
                  if (atom_c.eq.atom_d.and.ishd_1.ne.sh_hi_d) then

                     do ishd=ishd_1+1,sh_hi_d
                        if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, 
     I                       ilastd))
     &                       call errquit('Exiting in fitvcoul.',4,0)
                        nshbfd=ilastd-ifirstd+1
                        nao2 = nshbfc*nshbfd
                        ipoint=ipoint+nao2*nbf_cd
                     enddo
                  endif
                  icount_c=icount_c+nshbfc
               enddo            ! ishc
c     ga_put only if anything ever survived the schwarz scr
               
               if(v_nonzero) then
                    do imat=1,nmat
                      call ga_put(g_vc(imat), 
     B                     bf_lo_d, bf_hi_d,bf_lo_c,bf_hi_c,
     C                 dbl_mb(k_p+(imat-1)*nbf_at_c*nbf_at_d), nbf_at_d)
                    enddo
               endif

               if(use_nxtask) then
                  n1 = nxtask(nproc,1)
               else
                  n_batch_done = n_batch_done + 1
               endif
            endif
c     
          icount = icount + 1
c     
       enddo                    ! atom_d
      enddo                     ! atom_c

c
      if(.not.ma_chop_stack(lcd_shell))
     &     call errquit('dft_fitvc: cannot pop stack',0,0)
      if(use_nxtask) n1 = nxtask(-nproc,1)
      call ga_mask_sync(.not.use_nxtask,.false.)

      do imat=1,nmat
        call ga_symmetrize(g_vc(imat))
        if(oprint_ioinfo) call ga_print(g_vc(imat))
      enddo
c      
c
c     Compute the Coulombic self-energy of the fitted charge density.
c      
      if (dermat) then
         write(*,*)'*** Fix self-energy calculation in dft_fitvc'
      endif
c
c     Obtain the matrix of 2-ctr ERIs.
c      
c
      ! For ADFT, the cd_coef fields had CD_up + CD_dw + XC_up
      ! so we need to remove the XC part in order to compute
      ! Ecoul2
      if (ADFT) then
        call daxpy(nbf_cd,-1d0,cd_coef(1,nmat+1),1,cd_coef,1)
        if (ipol.gt.1) then
          call daxpy(nbf_cd,-1d0,cd_coef(1,nmat+2),1,cd_coef(1,2),1)
        endif
      endif
c
c     Compute Ecoul2 only for the first set of fit coefficients
c     !!! Should modify this so that nothing is done if nmat > 1
c     !!! Looks like g_2ceri not needed in that case
c
      Ecoul2=0.D0
      call ga_sync()
      call  ga_distribution(g_2ceri,me,ilo,ihi,jlo,jhi)
      if (ilo.gt.0) then
      call ga_access(g_2ceri, ilo, ihi, jlo, jhi, adrc, ldc)
         nnii=ihi-ilo+1
         do j = jlo,jhi
            iptr=adrc+(j-jlo)*ldc
            Ecoul2 = Ecoul2+ ddot(nnii,DBL_mb(iptr),1,
     C           cd_coef(ilo,1),1)*cd_coef(j,1)
       enddo
      call ga_release(g_2ceri, ilo, ihi, jlo, jhi)
      endif
c     

      Ecoul2=Ecoul2*0.5D0
c
c         write(*,*)'*** After destroy g_2ceri'
      call ga_mask_sync(.false.,.true.)
      call ga_dgop(Msg_Ecoul2,Ecoul2,1,'+')
      if (IOLGC)then
         if (.not. ga_destroy(g_2ceri)) call errquit
     &   ('dft_fitvc: could not destroy g_2ceri', 0,0)
      endif         
c
      ! For ADFT we need to save CD_up, CD_dw, XC_up, and XC_dw
      ! CD_up and CD_dw were already saved, so here we only add
      ! XC_up and XC_dw to the file
      if (adft) then
      if(me.eq.0) then
         call util_file_name('cdfit', .true., .false.,
     &         cdfit_filename)
         LU=79
         ierr=eaf_open(cdfit_filename(1:inp_strlen(cdfit_filename)),
     ,        eaf_rw, LU)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('dft_fitcd: eaf_open failed', ierr, GA_ERR)
        endif
         do imat = nmat+1, 2*nmat
            offset=(imat-1)*8d0*nbf_cd
            ierr = eaf_write(LU, offset, cd_coef(1,imat),
     .           8*nbf_cd)
            if (ierr .ne. 0) then
               call eaf_errmsg(ierr, errmsg)
               write(LuOut,*) ' IO offset ', offset
               write(LuOut,*) ' IO error message ',
     $              errmsg(1:inp_strlen(errmsg))
               call errquit('dft_fitcd: write failed',imat, INPUT_ERR)
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
      endif
c
      return
      end
