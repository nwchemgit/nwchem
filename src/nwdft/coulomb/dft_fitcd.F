      Subroutine dft_fitcd(CD_coef, i3c_ERI, Ecoul1, rtdb, g_cdinv,
     &                     g_dens, nTotEl, n_batch, n_bmax, iwhat,
     &                     n3c_dbl, iwhat_max, n_semi_bufs, fd, IOLGC)
c      
      implicit none
c
#include "bas.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "schwarz.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "sym.fh"
c
      double precision CD_coef(nbf_cd)
      double precision i3c_ERI(n3c_dbl)
      double precision Ecoul1
      integer rtdb
      integer g_cdinv
      integer g_dens         ! GA DM handle
      integer nTotEl
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax         ! max no. batches
      integer iwhat(n_bmax)  ! what batches have been computed [out]
      integer n3c_dbl
      integer iwhat_max
      integer n_semi_bufs
      integer fd
      Logical IOLGC
c
      integer current_nbuf, n3c_b, n_batch_read
      double precision offset
      double precision length, q1
      double precision F(2)
      double precision fac,THRESHOLD,Pmat_max,ERI_est,Q,AZ,alpha,
     ,     rho_int
      integer nproc,me,MXnshbf_ao,MXnshbf_cd,nshells_ao,
     ,     nshells_cd,ischw1,icount,Nao2,n1,n2,ishc,ifirstc,ilastc,
     ,     nshbfc,ishd,ifirstd, ilastd,nshbfd,istart,
     ,     ishp,ifirstp,ilastp,nshbfp,Nintegrals,nscr,
     ,     max1e,nscr1e,j,next,i,mscratch_2e3c,
     ,     ifirst,ilast,nshbf,maxg,ish,lmiss
      integer iTvec,lTvec,lgvec,igvec,lERI,iERI,lcgtf,icgtf,
     ,     lscr,iscr
      logical oprint,dowork,oprint_ioinfo
      integer nlast, ntasks
      integer ishc_be, ishd_be
      integer ndone,ipoint
      character*(nw_max_path_len) cdfit_filename
      double precision dabsmax
      integer nxtask
      external dabsmax,nxtask
      
#if defined(DFT_TIMER)
      integer nsh_blocks,niter
      double precision schw1
#include "timer.fh"
#endif
      
c     
c     Fit electronic charge density. The fitting coefficients are obtained by 
c     solving the normal equations for a standard linear least squares problem.
c     
c******************************************************************************
      

 
#if defined(DFT_TIMER)
      call ga_sync
      call t_stat(10,-1)
      call dft_second(ts1)
#endif
      oprint = util_print('coulomb fit', print_high)
      oprint_ioinfo = util_print('io info',print_debug)
      nproc  = ga_nnodes()
      me = ga_nodeid()
      
      THRESHOLD=10.d0**(-itol2e)
      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
      
      if( .not. bas_numbf(AO_bas_han,nbf_ao) )then
        call errquit('Exiting in fitcd.',0)
      end if
      if( .not. bas_nbf_cn_max(AO_bas_han,MXnshbf_ao) )then
        call errquit('Exiting in fitcd.',1)
      end if
      if( .not. bas_nbf_cn_max(CD_bas_han,MXnshbf_cd) )then
        call errquit('Exiting in fitcd.',2)
      end if
      
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in fitcd.',3)
      end if
      
      if( .not. bas_numcont(CD_bas_han,nshells_cd) )then
        call errquit('Exiting in fitcd.',4)
      end if
c      
c     Project the Coulomb potential onto each of the CD fitting functions.
c
c     Allocate scratch space.
      
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd,'Tvec',lTvec,iTvec))
     &  call errquit('dft_fitcd: cannot allocate Tvec',0)
c      call ga_sync
      
c     Initialize the T vector.
      
      call dfill(nbf_cd,0.D0,Dbl_MB(iTvec),1)
c      
c     Parallelize the loop over the products of AO shells.
c      
      ischw1 = 0
      call ga_sync
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      icount = max(1,icount)
      
      Nao2 = MXnshbf_ao*MXnshbf_ao
      if(.not.MA_Push_Get(MT_Dbl,Nao2,'gather vector',lgvec,igvec))
     &  call errquit('dft_fitcd: cannot allocate gather vector',0)

      call int_mem_2e3c(maxg, mscratch_2e3c)
      if(.not.MA_Push_Get(MT_Dbl,maxg*nshells_cd,'ERI',lERI,iERI))
     &  call errquit('dft_fitcd: cannot allocate ERI',0)
      nscr =  mscratch_2e3c 
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0)
c
      if(.not.incore)then
        iwhat_max = -1
        goto 25
      endif
c      call ga_sync
c      do n1 = 0, nproc
c        if(me.eq.n1)then
c          write(LuOut,*)' me, n_batch ', me, n_batch
c          call util_flush(LuOut)
c          call i_output(iwhat, 1, n_batch, 1, 1, n_batch, 1, 1)
c          call util_flush(LuOut)
c        endif
c        call ga_sync
c      enddo
c      call ga_sync
c
c     Exhaust in-core 3c 2e- buffers first.
c
      ndone=1
      ipoint=1
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
     $      ('dft_fitcd: unable to determine 3c2e file length',0)
c         write(6,*)' length = ',length
         current_nbuf = current_nbuf + 1
         offset = (current_nbuf-1)*n3c_b
         if (eaf_read(fd, offset, i3c_eri, n3c_b) .ne. 0)call errquit
     &      ('dft_fitcd: 3c-2e integral buffer read failed',0)
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
     &       call errquit('Exiting in dft_fitcd',111)
        nshbfc=ilastc-ifirstc+1
c
        do 23 ishd = 1,ishc
          dowork=iwhat(ndone).eq.((ishc*(ishc-1))/2 + ishd-1)
          if(dowork) then
c
c            write(LuOut,*)
c     &      'me,ishc,ishd,ndone,iwhat(ndone),ipoint', 
c     &       me,ishc,ishd,ndone,iwhat(ndone),ipoint
c            call util_flush(LuOut)
c
            if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, ilastd)) 
     &         call errquit('Exiting in dft_fitcd',1111)
            nshbfd=ilastd-ifirstd+1
            Nao2=nshbfd*nshbfc
c
c           Gather the needed density matrix elements.
c
            call ga_get(g_dens, ifirstd, ilastd, ifirstc, ilastc,
     &           dbl_mb(igvec), nshbfd)
c            
c           Determine the maximum density matrix element for this product of AO shells.
c
            Pmat_max = dabsmax(Nao2,DBL_MB(igvec))
c            
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c
            ERI_est=schwarz_shell(ishc,ishd)
c            
            if( Pmat_max*ERI_est.gt.THRESHOLD )then
c              
              ischw1 = ischw1 + nshells_cd
c              
              istart=iERI
              fac=2.0d0
              if(ishc.eq.ishd) fac=1.d0
              call dgemv('T',Nao2,nbf_cd,fac,i3c_ERI(ipoint),Nao2,
     &             Dbl_MB(igvec),1,1.d0,Dbl_MB(iTvec),1)
            end if
c
            ndone=ndone+1
c            write(6,*)' ndone, n_batch ',ndone,n_batch
            if(ndone.gt.n_batch)then
c
c             if not reading 3c2e integrals from disk - then done
c             start "direct"
c
              if(n_semi_bufs.eq.0)goto 25
c
c
c             Finished processing this batch of integrals, either get

c
              if (current_nbuf .eq. n_semi_bufs)goto 25
              current_nbuf = current_nbuf + 1
              offset = (current_nbuf-1)*n3c_b
              if (eaf_read(fd, offset, i3c_eri, n3c_b) .ne. 0)
     &           call errquit
     &           ('dft_fitcd: 3c-2e integral buffer read failed',0)
              n_batch_read = nint(i3c_eri(n3c_dbl))
c              write(6,*)' n_batch_read = ',n_batch_read
              n_batch = n_batch + n_batch_read
              ipoint=1
            else
              ipoint=ipoint+Nao2*nbf_cd
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
     &       call errquit('Exiting in dft_fitcd',111)
        nshbfc=ilastc-ifirstc+1
c
        do 230 ishd = ishd_be, ishc
          dowork=n1.eq.n2
          if(dowork) then
c
c            write(LuOut,*)'me,ishc,ishd,n1,n2', me,ishc,ishd,n1,n2
c            call util_flush(LuOut)
c
            if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, ilastd)) 
     &         call errquit('Exiting in dft_fitcd',1111)
            nshbfd=ilastd-ifirstd+1
            Nao2=nshbfd*nshbfc
c
c           Gather the needed density matrix elements.
c
            call ga_get(g_dens, ifirstd, ilastd, ifirstc, ilastc,
     &           dbl_mb(igvec), nshbfd)
c            
c           Determine the maximum density matrix element for this product of AO shells.
c
            Pmat_max = dabsmax(Nao2,DBL_MB(igvec))
c            
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c
            ERI_est=schwarz_shell(ishc,ishd)
c            
            if( Pmat_max*ERI_est.gt.THRESHOLD )then
c              
              ischw1 = ischw1 + nshells_cd
c              
              istart=iERI
              fac=2.0d0
              if(ishc.eq.ishd) fac=1.d0
              do 220 ishp = 1,nshells_cd
c                  
c               Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c                  
                if (.not. bas_cn2bfr( CD_bas_han,ishp,ifirstp,ilastp))
     &               call errquit('Exiting in fitcd.',5)
                nshbfp=ilastp-ifirstp+1
c                  
                Nintegrals = nshbfp*Nao2
                nscr=mscratch_2e3c
c
c               check for use of symmetry
c
                if (oskel) then
                   if (sym_shell(cd_bas_han, ishp, q1))then
                      call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,
     &                ishd,nscr,Dbl_MB(iscr),Nintegrals,Dbl_MB(istart))
                      call dscal(Nintegrals, q1, Dbl_MB(istart), 1)
                   else
                      call dzero(Dbl_MB(istart),Nintegrals)
                   endif
                else
                   call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,
     &             ishd,nscr,Dbl_MB(iscr),Nintegrals,Dbl_MB(istart))
                endif                   
c
                istart=istart+Nintegrals
c
  220         continue
c               
c             Multiply the matrix of 3-ctr ERIs by the density matrix elements.
c               
              call dgemv('T',Nao2,nbf_cd,fac,Dbl_MB(iERI),Nao2,
     &             Dbl_MB(igvec),1,1.d0,Dbl_MB(iTvec),1)
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
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(1)='  Tvec    '
      call tcllct(t_sub(1,1),ts2,ts1)
      call dft_second(ts1)
#endif
c
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitcd: cannot pop stack',0)
      if(.not.ma_pop_stack(lERI))
     &  call errquit('dft_fitcd: cannot pop stack',0)
      if(.not.ma_pop_stack(lgvec))
     &  call errquit('dft_fitcd: cannot pop stack',0)
      
      n1 = nxtask(-nproc,icount)
      call ga_dgop(Msg_gop_Tvec,Dbl_MB(iTvec),nbf_cd,'+')
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(2)='  DGOP    '
      call tcllct(t_sub(1,2),ts2,ts1)
      call dft_second(ts1)
#endif
c
      if (oskel)
     &   call sym_vec_symmetrize(geom, cd_bas_han,  Dbl_MB(iTvec))
c      
c     Compute the CD expansion coefficients.
c     
c     Allocate scratch space.
c      

      if(.not.MA_Push_Get(MT_Dbl,nbf_cd,'cgtf',lcgtf,icgtf))
     &  call errquit('dft_fitcd: cannot allocate cgtf',0)
      call int_mem_1e(max1e,nscr1e) 
      nscr = max(nscr1e, nbf_cd)

      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitcd: cannot allocate scr',0)
      
c      
c     Integrate CD fitting functions in each shell
c

      do ish = 1,nshells_cd
        if (.not. bas_cn2bfr( CD_bas_han,ish,
     &                       ifirst,ilast))
     &      call errquit('Exiting in fitcd.',44)
                nshbf=ilast-ifirst+1
        call int_1cg( CD_bas_han, ish, nscr, Dbl_MB(iscr), 
     &              nshbf, Dbl_MB(icgtf+ifirst-1) )
      enddo 
c
c     Compute the orbital integrated electron density.
c
      Q = Dble(nTotEl)
c      
#if defined(DFT_TIMER)
      call ga_sync
      call dft_second(ts2)
      subsect(3)='  Int_cgtf'
      call tcllct(t_sub(1,3),ts2,ts1)
      call dft_second(ts1)
#endif
c
c     Compute the Lagrange Multiplier, alpha.
c     
      if(IOLGC) then
        call gacrea(geom,CD_bas_han,
     &       nbf_cd,nbf_cd,'CD 2cERInv',g_cdinv,'row')
        if (me.eq.0)then
           lmiss = 0
           call dft_invio('CDI',g_cdinv,nbf_cd,'READ',lmiss)
           if(lmiss.ne.0)call errquit
     &     (' dft_fitcd; dft_invio - abnormal read of CDI ',0)
        endif
        call ga_sync
        if(me.eq.0.and.oprint_ioinfo)then
          write(6,*)' g_cdinv = ', g_cdinv
          call ga_print(g_cdinv)
        endif
      endif
#if defined(DFT_TIMER)
      call ga_sync
      call dft_second(ts2)
      subsect(4)='  RD CDinv'
      call tcllct(t_sub(1,4),ts2,ts1)
      call dft_second(ts1)
#endif
      
      F(1)=0.D0
      F(2)=0.D0
      next=nxtval(nproc)+1
      do j=1,nbf_cd
        if (j.eq.next) then
          next=nxtval(nproc)+1
          call get_col(g_cdinv,
     *         nbf_cd,j,Dbl_MB(iscr))
          AZ=ddot(nbf_cd,DBL_mb(iscr),1,DBL_MB(icgtf),1)
          F(1)=F(1)+AZ*DBL_MB(iTvec+j-1)
          F(2)=F(2)+AZ*DBL_MB(icgtf+j-1)
        endif
      enddo
      next=nxtval(-nproc)
      call ga_dgop(Msg_Lagrange,F,2,'+')
      alpha = (Q - F(1))/F(2)
      call ga_sync
c
c     Multiply: S**(-1) * (alpha*N + T)
c      
      do i=1,nbf_cd
        CD_coef(i)=0.D0
      enddo
      next=nxtval(nproc)+1
      do j=1,nbf_cd
        if (j.eq.next) then
          next=nxtval(nproc)+1
          AZ=DBL_MB(iTvec+j-1)+alpha*DBL_MB(icgtf+j-1)
          call get_col(g_cdinv,
     *         nbf_cd,j,Dbl_MB(iscr))
          call daxpy(nbf_cd,AZ,DBL_MB(iscr),1,CD_Coef,1)
        endif
      enddo
      next=nxtval(-nproc)
      call ga_dgop(Msg_CD_coef,CD_coef,nbf_cd,'+')
      call ga_sync 
      if(me.eq.0) then
         call util_file_name('cdfit', .true., .false., cdfit_filename)
         open(unit=79,file=cdfit_filename,access='sequential',
     &         form='unformatted',status='unknown')
         write(79) (CD_coef(i),i=1,nbf_cd)
         close(79)
      endif
      call ga_sync
      if(IOLGC) call gadest(g_cdinv)
c
c     calculation of Ecoul1
c
      Ecoul1=ddot(nbf_cd,DBL_mb(iTvec),1,CD_coef,1)
c      
c     Compute the integrated electron density from the fit.
c      
      rho_int = ddot(nbf_cd,CD_coef,1,Dbl_MB(icgtf),1)
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitcd: cannot pop stack',0)
      if(.not.ma_pop_stack(lcgtf))
     &  call errquit('dft_fitcd: cannot pop stack',0)
      if(.not.ma_pop_stack(lTvec))
     &  call errquit('dft_fitcd: cannot pop stack',0)
      
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(5)='  cd_coef '
      call tcllct(t_sub(1,5),ts2,ts1)
#endif
      
      if(oprint.and.me.eq.0 )then
        write(LuOut,*)
        write(LuOut,'(a,f14.6,a,1pe13.6)')
     &       '   CD fit integrated density:',rho_int,
     &       '    Lagrange Mult.:  ',alpha
        write(LuOut,*)
        call util_flush(LuOut)
      end if
      
#if defined(DFT_TIMER)
        
      call ga_igop(Msg_schw1,ischw1,1,'+')
      nsh_blocks = nshells_cd*((nshells_ao*(nshells_ao+1))/2)
        
      schw1 = (dble(nsh_blocks - ischw1)/dble(nsh_blocks))*100
        
      ltime = 7
        
      if( me.eq.0 )then
          
        niter  = (nshells_ao*(nshells_ao+1))/2
        ntasks = niter/icount
        write(ltime,*)
        write(ltime,9996) '      Analysis of FITCD'
        write(ltime,*)
        write(ltime,9995) ' number of nodes:  ',ga_nnodes()
        write(ltime,9995) ' number of tasks:  ',ntasks
        write(ltime,9995) ' chunk size     :  ',icount
        write(ltime,*)
        write(ltime,9996) 'Schwarz Screening - Density*(p|cd)'
        write(ltime,*)
        write(ltime,9997) 'Threshold :  ',THRESHOLD
        write(ltime,9998) '% Screened:  ',schw1
        write(ltime,*)
          
      end if
        
      subname=' Fitcd    '
      call t_stat(5,1)
      
 9995 format(10x,a,i10)
 9996 format(2x,a)
 9997 format(2x,a,e12.4)
 9998 format(2x,a,2f6.2)
      
#endif
c      call errquit(' temp stop semi-direct ',999)
      return
      end
