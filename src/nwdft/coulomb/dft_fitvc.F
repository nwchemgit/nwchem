      Subroutine dft_fitvc(rtdb, geom, AO_bas_han, CD_bas_han, CD_coef, 
     &     g_vc, Ecoul2, nbf_cd, g_2ceri,IOLGC,itol2e,incore,
     &     n_batch,n_bmax,iwhat,i3c_ERI,n3c_dbl)
C     $Id: dft_fitvc.F,v 1.10 1996-10-18 00:52:58 d3h449 Exp $
      implicit none
      
c**** 
c**** nwchem handles
c**** 
      integer rtdb, geom, AO_bas_han, CD_bas_han
      integer g_vc    ! GA handle for V Coulomb [output]
      integer g_2ceri ! GA handle for repulsion int. [input]
      integer nbf_cd    ! no. CD functions [input]
      double precision Ecoul2 ! [output]
      integer itol2e ! bi-el int tolerance [input]
      logical incore ! true if 3c int are incore [input]
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl
      double precision i3c_ERI(n3c_dbl)
c**** 
c      
c-->  Charge Density Expansion Coefficients
c      
      double precision CD_coef(nbf_cd)

      Logical IOLGC
      
#include "bas.fh"
#include "rtdb.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "schwarz.fh"

#if defined(DFT_TIMER)
#include "timer.fh"
#endif
c     
c     Compute the fitted Coulomb potential. The fitting coefficients are
c     obtained by solving the normal equations for a standard linear least
c     squares problem.
c
      integer MXshbf_ao,nshells_ao,nshells_cd,MXshbf_cd
      integer me,nproc,i1,ishp,ishbf,lshbf,lmiss,ischw1,ischw2,
     ,     icount,Nao2,i,nERI,maxg,mscratch_2e3c,nscr,
     ,     n1,n2,ishc,ifirstc,ilastc,nshbfc,ishd,ifirstd,ilastd,
     ,     nshbfd,nfunz,ifirstp,ilastp,Nint,niter,next,nsh_blocks,
     ,     ntasks,nshbf_cd
      integer lcd_shell,icd_shell,lgvec,igvec,lsvec,isvec,
     ,     lsvec2,isvec2,iERI,lERI,lscr,iscr
      double precision THRESHOLD,CD_coef_max,shmax,
     ,     Vcoul_est,ERI_est,schw1,schw2
      integer nxtask
      logical dowork
      integer ndone,ipoint,ipoint2
      double precision ddot
      external ddot,nxtask
c
#if defined(DFT_TIMER)
      call ga_sync
      call t_stat(10,-1)
      call dft_second(ts1)
#endif
      THRESHOLD=10.d0**(-itol2e)
c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
      if( .not. bas_nbf_cn_max(AO_bas_han,MXshbf_ao) )then
        call errquit('Exiting in fitvcoul.',1)
      end if
      
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in fitvcoul.',2)
      end if
      
      if( .not. bas_numcont(CD_bas_han,nshells_cd) )then
        call errquit('Exiting in fitvcoul.',3)
      end if
      
      if( .not. bas_nbf_cn_max(CD_bas_han,MXshbf_cd) )then
        call errquit('Exiting in fitvcoul.',4)
      end if
      
      me=ga_nodeid()
      nproc  = ga_nnodes()
c      
c     Compute the fitted Vcoul.
c      
c     Allocate scratch space.
c      
      if(.not.MA_Push_Get(MT_Dbl,nshells_cd,'CD shell max',
     &  lcd_shell,icd_shell))
     &  call errquit('dft_fitvc: cannot allocate CD shell max',0)
c      
      if(.not.MA_Push_Get(MT_Dbl,MXshbf_cd,'gvec',lgvec,igvec))
     &  call errquit('dft_fitvc: cannot allocate gvec',0)
c      
c     Determine the largest magnitude of the CD expansion coefficients in
c     each shell.
c      
      CD_coef_max = 0.D0
      i1 = 1
      do 205 ishp = 1,nshells_cd
        if( .not. bas_cn2bfr(CD_bas_han,ishp,ishbf,lshbf) )then
          call errquit('Exiting in fitvcoul.',4)
        end if
        shmax = 0.D0
        do 200 i = ishbf,lshbf
          shmax = max(abs(CD_coef(i)), shmax)
  200   continue
        Dbl_MB((icd_shell-1)+ishp) = shmax
        CD_coef_max = max(shmax, CD_coef_max)
  205 continue
      if(IOLGC) then
        call gacrea(geom,CD_bas_han,
     &       nbf_cd,nbf_cd,'CD 2cERI',g_2ceri,'row')
        if(me.eq.0)
     &       call dft_invio('CD',g_2ceri,nbf_cd,'READ',lmiss)
        call ga_sync
      endif
c      
c     Parallelize the loop over the products of AO shells.
c      
      ischw1 = 0
      ischw2 = 0
      
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
c      
c     check for zero ... must be at least 1
c
      icount = max(icount,1)
      if(nproc.eq.1)icount = 1
c      
      Nao2 = MXshbf_ao**2
      if(.not.MA_Push_Get(MT_Dbl,Nao2,'isvec',lsvec,isvec))
     &  call errquit('dft_fitvc: cannot allocate isvec',0)
      if(.not.MA_Push_Get(MT_Dbl,Nao2,'isvec2',lsvec2,isvec2))
     &  call errquit('dft_fitvc: cannot allocate isvec2',0)
      nERI    = Nao2*MXshbf_cd
      if(.not.MA_Push_Get(MT_Dbl,nERI,'ERI',lERI,iERI))
     &  call errquit('dft_fitvc: cannot allocate ERI',0)
c
      call int_mem_2e3c(maxg, mscratch_2e3c)
      nscr    = mscratch_2e3c
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitvc: cannot allocate scr',0)

      call ga_zero(g_vc)
      if(.not.incore) n1 = nxtask(nproc,icount)
      n2 = 0
      ndone=1
      ipoint=1
c     
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(1)='  Pre-fit '
      call tcllct(t_sub(1,1),ts2,ts1)
      call dft_second(ts1)
#endif
c     
      do 240 ishc = 1,nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han,ishc,
     &       ifirstc,ilastc))
     &       call errquit('Exiting in fitvcoul.',3)
        nshbfc=ilastc-ifirstc+1
        
        do 230 ishd = 1,ishc
          

          if(incore) then
            dowork=iwhat(ndone).eq.(ishc*(ishc-1)+ishd-1)
            if(ndone.gt.n_batch) then
              write(0,*) ' ndone nbatch ',ndone,n_batch
              call errquit(' fitvc:ndone gt nbatch ',ndone)
            endif
          else
            dowork=n1.eq.n2
          endif
          if(dowork) then
            
            if (.not. bas_cn2bfr( AO_bas_han,ishd,
     &           ifirstd,ilastd))
     &           call errquit('Exiting in fitvcoul.',4)
            
            nshbfd=ilastd-ifirstd+1
            Nao2 = nshbfc*nshbfd
c           
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c           
            ERI_est=schwarz_shell(ishc,ishd)

            if( CD_coef_max*ERI_est.gt.THRESHOLD )then
              ischw1 = ischw1 + nshells_cd
              
              call dfill(Nao2,0.D0,dbl_mb(isvec),1)
              nfunz=0
              if(incore) ipoint2=ipoint
                  if(incore) then
c                 
c                   Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c                 
                    call dgemv('N',Nao2,nbf_cd,1.D0,
     &                   i3c_ERI(ipoint),Nao2,
     &                   CD_coef,1,1.D0,Dbl_MB(isvec),1)
                  else
              do 220 ishp = 1,nshells_cd
                
                if (.not. bas_cn2bfr( CD_bas_han,ishp,
     &               ifirstp,ilastp))
     &               call errquit('Exiting in fitvcoul.',5)
                nshbf_cd=ilastp-ifirstp+1
                Nint=nshbf_cd*Nao2
                Vcoul_est = Dbl_MB((icd_shell-1)+ishp)*ERI_est
                if( Vcoul_est.gt.THRESHOLD )then
                  
                  ischw2 = ischw2 + 1
                  
c                   
c                   Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c                   
                    call int_2e3c(CD_bas_han,ishp,AO_bas_han,
     &                   ishc,ishd,nscr,Dbl_MB(iscr),Nint,Dbl_MB(iERI))
c                   
c                   Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c                   
                    call dgemv('N',Nao2,nshbf_cd,1.D0,
     &                   Dbl_MB(iERI),Nao2,
     &                   CD_coef(1+nfunz),1,1.D0,Dbl_MB(isvec),1)
                  endif


                if(incore) ipoint2=ipoint2+Nint
                nfunz=nfunz+nshbf_cd
                
  220         continue
                endif
c           
c             Scatter the computed values of Vcoul.
c           
              call sca_lab_mat(AO_bas_han,ishc,ishd,
     &             nshbfc,nshbfd,g_vc,DBL_MB(isvec),dbl_mb(isvec2))

            endif
            
            if(incore) then
              ndone=ndone+1
              if(ndone.gt.n_batch) goto 9999
              ipoint=ipoint+Nao2*nbf_cd
            else
              n1 = nxtask(nproc,icount)
            endif
            
            
          end if
          n2 = n2 + 1
        
  230 continue
  240 continue
 9999 continue
      
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitvc: cannot pop stack',0)
      if(.not.ma_pop_stack(lERI))
     &  call errquit('dft_fitvc: cannot pop stack',0)
      if(.not.ma_pop_stack(lsvec2))
     &  call errquit('dft_fitvc: cannot pop stack',0)
      if(.not.ma_pop_stack(lsvec))
     &  call errquit('dft_fitvc: cannot pop stack',0)
c
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(2)='  Vcoul   '
      call tcllct(t_sub(1,2),ts2,ts1)
      call dft_second(ts1)
#endif
c
      n1 = nxtask(-nproc,icount)
c      
      if(.not.ma_pop_stack(lgvec))
     &  call errquit('dft_fitvc: cannot pop stack',0)
      if(.not.ma_pop_stack(lcd_shell))
     &  call errquit('dft_fitvc: cannot pop stack',0)
c      
c     Compute the Coulombic self-energy of the fitted charge density.
c      
c     Obtain the matrix of 2-ctr ERIs.
c      
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd,'scr',lscr,iscr))
     &  call errquit('dft_fitvc: cannot allocate scr',0)
      Ecoul2=0.D0
      next=nxtval(nproc)+1
      do i=1,nbf_cd
        if(i.eq.next) then
          next=nxtval(nproc)+1
          call get_col(g_2ceri,nbf_cd,i,DBL_MB(iscr))
          Ecoul2 = Ecoul2+ddot(
     *         nbf_cd,CD_coef,1,Dbl_MB(iscr),1)*CD_coef(i)
        endif
      enddo
      
      next=nxtval(-nproc)
      
      Ecoul2=Ecoul2*0.5D0
      if(IOLGC) call gadest(g_2ceri)
      call ga_dgop(Msg_Ecoul2,Ecoul2,1,'+')
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitvc: cannot pop stack',0)
c      
#if defined(DFT_TIMER)
      call dft_second(ts2)
      subsect(3)='  Ecoul2  '
      call tcllct(t_sub(1,3),ts2,ts1)
#endif
c      
#if defined(DFT_TIMER)
c      
      call ga_igop(Msg_schw2,ischw1,1,'+')
      call ga_igop(Msg_schw2,ischw2,1,'+')
      nsh_blocks = nshells_cd*((nshells_ao*(nshells_ao+1))/2)
      schw1 = (dble(nsh_blocks - ischw1)/dble(nsh_blocks))*100
      schw2 = (dble(nsh_blocks - ischw2)/dble(nsh_blocks))*100
      ltime = 7
      if( me.eq.0 )then
        niter  = (nshells_ao*(nshells_ao+1))/2
        ntasks = niter/icount
        write(ltime,*)
        write(ltime,9996) '      Analysis of FITVCOUL'
        write(ltime,*)
        write(ltime,9995) ' number of nodes:  ',nproc
        write(ltime,9995) ' number of tasks:  ',ntasks
        write(ltime,9995) ' chunk size     :  ',icount
        write(ltime,*)
        write(ltime,9996) 'Schwarz Screening - Density*(p|cd)'
        write(ltime,*)
        write(ltime,9997) 'Threshold :  ',THRESHOLD
        write(ltime,9998) '% Screened:  ',schw1,schw2
        write(ltime,*)
      end if
      subname=' Fitvcoul '
      call t_stat(3,1)
c      
 9995 format(10x,a,i10)
 9996 format(10x,a)
 9997 format(2x,a,e12.4)
 9998 format(2x,a,2f6.2)
c      
#endif
      return
      end
