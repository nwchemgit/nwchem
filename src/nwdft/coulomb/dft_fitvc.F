      Subroutine dft_fitvc(rtdb, geom, AO_bas_han, CD_bas_han, CD_coef, 
     &     iga_vc, Ecoul2, nbf_cd, iga_2ceri,IOLGC,itol2e)
C$Id: dft_fitvc.F,v 1.3 1995-06-19 18:24:15 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han
c**** 

      Logical LResult,IOLGC
      
#include <bas.fh>
#include "rtdb.fh"
#include "context.fh"

#include <mafdecls.h>
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"

#if defined(DFT_TIMER)
#include <timer.fh>
#endif
c      
c-->  Charge Density Expansion Coefficients
c      
      Dimension CD_coef(nbf_cd)
c     
c     Compute the fitted Coulomb potential. The fitting coefficients are
c     obtained by solving the normal equations for a standard linear least
c     squares problem.
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
      LResult = MA_Push_Get(MT_Dbl,nshells_cd,'CD shell max',
     &     lcd_shell,icd_shell)
      
      LResult = MA_Push_Get(MT_Dbl,MXshbf_cd,'gvec',lgvec,igvec)
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
     &       nbf_cd,nbf_cd,'CD 2cERI',iga_2ceri,'row')
        if(me.eq.0)
     &       call dft_invio('CD',iga_2ceri,nbf_cd,'READ',lmiss)
        call ga_sync
      endif
c      
c     Parallelize the loop over the products of AO shells.
c      
      ischw1 = 0
      ischw2 = 0
      
      nproc  = ga_nnodes()
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
c      
c     check for zero ... must be at least 1
c
      icount = max(icount,1)
      if(nproc.eq.1)icount = 1
c      
      Nao2 = MXshbf_ao**2
      LResult = MA_Push_Get(MT_Dbl,Nao2,'isvec',lsvec,isvec)
      LResult = MA_Push_Get(MT_Dbl,Nao2,'isvec2',lsvec2,isvec2)
      nERI    = Nao2*MXshbf_cd
      LResult = MA_Push_Get(MT_Dbl,nERI,'ERI',lERI,iERI)

      call int_mem_2e3c(maxg, mscratch_2e3c)
      nscr    = mscratch_2e3c
      LResult = MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr)

      call ga_zero(iga_vc)
      n1 = nxtask(nproc,icount)
      n2 = 0
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
     &                        ifirstc,ilastc))
     &          call errquit('Exiting in fitvcoul.',3)
        nshbfc=ilastc-ifirstc+1
 
        do 230 ishd = 1,ishc
          
c            
c         Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c         
          ERI_est=schwarz_shell(ishc,ishd)

          if( CD_coef_max*ERI_est.gt.THRESHOLD )then

            if( n1.eq.n2 )then
              
              ischw1 = ischw1 + nshells_cd
c              
c             Define the scatter index.
c              
                if (.not. bas_cn2bfr( AO_bas_han,ishd,
     &             ifirstd,ilastd))
     &             call errquit('Exiting in fitvcoul.',4)

                nshbfd=ilastd-ifirstd+1
                Nao2 = nshbfc*nshbfd
                           
                call dfill(Nao2,0.D0,dbl_mb(isvec),1)
                nfunz=0
              
              do 220 ishp = 1,nshells_cd
                
                Vcoul_est = Dbl_MB((icd_shell-1)+ishp)*ERI_est
                
                   if (.not. bas_cn2bfr( CD_bas_han,ishp,
     &                  ifirstp,ilastp))
     &                  call errquit('Exiting in fitvcoul.',5)
                   nshbf_cd=ilastp-ifirstp+1
                if( Vcoul_est.gt.THRESHOLD )then
                   
                   ischw2 = ischw2 + 1
                  
                   Nint=nshbf_cd*Nao2
c                 
c                 Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c
                  call int_2e3c(CD_bas_han,ishp,AO_bas_han,ishc,ishd,
     &                 nscr,Dbl_MB(iscr),Nint,Dbl_MB(iERI))
c              
c                 Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c              
                  call dgemv('N',Nao2,nshbf_cd,1.D0,Dbl_MB(iERI),Nao2,
     &                 CD_coef(1+nfunz),1,1.D0,Dbl_MB(isvec),1)
               endif

                nfunz=nfunz+nshbf_cd
                
 220         continue
c             
c             Scatter the computed values of Vcoul.
c          
              call sca_lab_mat(AO_bas_han,ishc,ishd,nshbfc,nshbfd,
     &             iga_vc,DBL_MB(isvec),dbl_mb(isvec2))

            n1 = nxtask(nproc,icount)

            end if
            
          n2 = n2 + 1
         
            
          end if
          
  230   continue
  240 continue
      
      LResult = MA_Pop_Stack(lscr)
      LResult = MA_Pop_Stack(lERI)
      LResult = MA_Pop_Stack(lsvec2)
      LResult = MA_Pop_Stack(lsvec)
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
      LResult = MA_Pop_Stack(lgvec)
      LResult = MA_Pop_Stack(lcd_shell)
c      
c     Compute the Coulombic self-energy of the fitted charge density.
c      
c     Obtain the matrix of 2-ctr ERIs.
c      
      LResult = MA_Push_Get(MT_Dbl,nbf_cd,'scr',lscr,iscr)
      Ecoul2=0.D0
      next=nxtval(nproc)+1
      do i=1,nbf_cd
        if(i.eq.next) then
          next=nxtval(nproc)+1
          call get_col(iga_2ceri,nbf_cd,i,DBL_MB(iscr))
          Ecoul2 = Ecoul2+ddot(
     *         nbf_cd,CD_coef,1,Dbl_MB(iscr),1)*CD_coef(i)
        endif
      enddo
      
      next=nxtval(-nproc)
      
      Ecoul2=Ecoul2*0.5D0
      if(IOLGC) call gadest(iga_2ceri)
      call ga_dgop(Msg_Ecoul2,Ecoul2,1,'+')
      LResult = MA_Pop_Stack(lscr)
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
