      Subroutine dft_fitcd( rtdb,AO_bas_han, CD_bas_han,geom,
     &     CD_coef,  nbf_cd,
     &     g_cdinv, g_dens, nTotEl, IOLGC,itol2e,
     &     Ecoul1,incore,
     &     n_batch,n_bmax,iwhat,i3c_ERI,n3c_dbl)
      
      implicit none
      
c****
c**** nwchem handles
c****
      integer AO_bas_han, CD_bas_han,geom,rtdb
      integer g_dens       ! GA DM handle
      integer g_cdinv
      integer nTotEl
      integer itol2e
      integer nbf_cd
      double precision Ecoul1
      logical incore ! true if 3c int are incore [input]
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl
      double precision i3c_ERI(n3c_dbl)
c-->  Charge Density Expansion Coefficients
      
      double precision CD_coef(nbf_cd)
c**** 
      Logical LResult,IOLGC, doesit, isit
      
#include "bas.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "context.fh"

#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "schwarz.fh"
#include "msgids.fh"

      double precision F(2)
      double precision fac,THRESHOLD,Pmat_max,ERI_est,Q,AZ,alpha,
     ,     rho_int
      integer nproc,me,nbf_ao,MXnshbf_ao,MXnshbf_cd,nshells_ao,
     ,     nshells_cd,ischw1,icount,Nao2,n1,n2,ishc,ifirstc,ilastc,
     ,     nshbfc,ishd,ifirstd, ilastd,nshbfd,istart,
     ,     ishp,ifirstp,ilastp,nshbfp,Nint,nscr,
     ,     max1e,nscr1e,j,next,i,mscratch_2e3c,
     ,     ifirst,ilast,nshbf,maxg,ish,lmiss
      integer iTvec,lTvec,lgvec,igvec,lERI,iERI,lcgtf,icgtf,
     ,     lscr,iscr
      logical oprint,dowork
      integer ndone,ipoint
      double precision dabsmax
      integer nxtask
      external dabsmax,nxtask
      
#if defined(DFT_TIMER)
      integer nsh_blocks,niter,ntasks
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
      
      LResult = MA_Push_Get(MT_Dbl,nbf_cd,'Tvec',lTvec,iTvec)
      call ga_sync
      
c     Initialize the T vector.
      
      call dfill(nbf_cd,0.D0,Dbl_MB(iTvec),1)
c      
c     Parallelize the loop over the products of AO shells.
c      
      ischw1 = 0
C      call ga_sync
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      icount = max(1,icount)
      
      Nao2 = MXnshbf_ao*MXnshbf_ao
      LResult = MA_Push_Get(MT_Dbl,Nao2,'gather vector',lgvec,igvec)

      call int_mem_2e3c(maxg, mscratch_2e3c)
      LResult = MA_Push_Get(MT_Dbl,maxg*nshells_cd,'ERI',lERI,iERI)
      nscr =  mscratch_2e3c 
      LResult = MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr)
      if(.not.incore) n1 = nxtask(nproc,icount)
      n2 = 0
      ndone=1
      ipoint=1
      do 240 ishc = 1,nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han, ishc, ifirstc, ilastc)) 
     *       call errquit('Exiting in dft_fitcd',111)
        nshbfc=ilastc-ifirstc+1

        do 230 ishd = 1,ishc
          dowork=.false.
          if(incore) then
            dowork=iwhat(ndone).eq.(ishc*(ishc-1)+ishd-1)
            if(ndone.gt.n_batch) then
              write(0,*) ' ndone nbatch ',ndone,n_batch
              call errquit(' fitcd: ndone gt nbatch ',ndone)
            endif
          else
            dowork=n1.eq.n2
          endif
          if(dowork) then
            if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, ilastd)) 
     *         call errquit('Exiting in dft_fitcd',1111)
            nshbfd=ilastd-ifirstd+1
            Nao2=nshbfd*nshbfc
C
c           Gather the needed density matrix elements.
c
            call ga_get(g_dens,
     *           ifirstd,ilastd,
     *           ifirstc,ilastc,
     *           dbl_mb(igvec),nshbfd)
            
            
c           Determine the maximum density matrix element for this product of AO shells.
            Pmat_max = dabsmax(Nao2,DBL_MB(igvec))
            
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
              ERI_est=schwarz_shell(ishc,ishd)
            
            if( Pmat_max*ERI_est.gt.THRESHOLD )then
              
              ischw1 = ischw1 + nshells_cd
              
              istart=iERI
              fac=2.0d0
              if(ishc.eq.ishd) fac=1.d0
              if(incore) then
                call dgemv('T',Nao2,nbf_cd,fac,i3c_ERI(ipoint),Nao2,
     &               Dbl_MB(igvec),1,1.d0,Dbl_MB(iTvec),1)
              else
                do 220 ishp = 1,nshells_cd
                  
c                 Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
                  
                  if (.not. bas_cn2bfr( CD_bas_han,ishp,ifirstp,ilastp))
     &                 call errquit('Exiting in fitcd.',5)
                  nshbfp=ilastp-ifirstp+1
                  
                  Nint = nshbfp*Nao2
                  nscr=mscratch_2e3c
                  
                  call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,
     &                 ishd,nscr,Dbl_MB(iscr),Nint,Dbl_MB(istart))
                  istart=istart+Nint
  220           continue
c               
c               Multiply the matrix of 3-ctr ERIs by the density matrix elements.
c               
                call dgemv('T',Nao2,nbf_cd,fac,Dbl_MB(iERI),Nao2,
     &               Dbl_MB(igvec),1,1.d0,Dbl_MB(iTvec),1)
              endif
              
            end if
            if(incore) then
              ndone=ndone+1
              if(ndone.gt.n_batch) goto 9999
              ipoint=ipoint+Nao2*nbf_cd
            else
              n1 = nxtask(nproc,icount)
            endif
          end if
          
          n2 = n2 + 1
          
  230   continue
  240 continue
 9999 continue
c
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(1)='  Tvec    '
      call tcllct(t_sub(1,1),ts2,ts1)
      call dft_second(ts1)
#endif
c
      LResult = MA_Pop_Stack(lscr)
      LResult = MA_Pop_Stack(lERI)
      LResult = MA_Pop_Stack(lgvec)
      
      n1 = nxtask(-nproc,icount)
      call ga_dgop(Msg_gop_Tvec,Dbl_MB(iTvec),nbf_cd,'+')
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(2)='  DGOP    '
      call tcllct(t_sub(1,2),ts2,ts1)
      call dft_second(ts1)
#endif
c      
c     Compute the CD expansion coefficients.
c     
c     Allocate scratch space.
c      

      LResult = MA_Push_Get(MT_Dbl,nbf_cd,'cgtf',lcgtf,icgtf)
      call int_mem_1e(max1e,nscr1e) 
      nscr = max(nscr1e, nbf_cd)

      LResult = MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr)
      
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
        if(me.eq.0)
     &       call dft_invio('CDI',g_cdinv,nbf_cd,'READ',lmiss)
        call ga_sync
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
      inquire(file='cdfit.dat',exist=doesit,opened=isit)
      if(.not.isit)then
        if(.not.doesit)then
          open(unit=79,file='cdfit.dat',access='sequential',
     &         form='unformatted',status='new')
        else
          open(unit=79,file='cdfit.dat',access='sequential',
     &         form='unformatted',status='old')
          endif
        endif
          write(79) (CD_coef(i),i=1,nbf_cd)
      rewind 79
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
      LResult = MA_Pop_Stack(lscr)
      LResult = MA_Pop_Stack(lcgtf)
      LResult = MA_Pop_Stack(lTvec)
      
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(5)='  cd_coef '
      call tcllct(t_sub(1,5),ts2,ts1)
#endif
      
      if(oprint.and.me.eq.0 )then
        write(LuOut,*)
        write(LuOut,'(a,f14.6,a,1pe13.6)')
     &       '   CD fit] integrated density:',rho_int,
     &       '    Lagrange Mult.:  ',alpha
        write(LuOut,*)
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
      return
      end
