      Subroutine dft_fitvc(CD_coef, i3c_ERI, Ecoul2,  g_vc, 
     &                     n_batch, n_bmax, iwhat, n3c_dbl,
     &                     iwhat_max, n_semi_bufs, fd, IOLGC,
     &                     dermat, nmat, do_init)
c     $Id: dft_fitvc.F,v 1.40 2004-09-07 06:11:45 windus Exp $
      implicit none
#include "errquit.fh"
c      
#include "bas.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "schwarz.fh"
#include "util.fh"
#include "eaf.fh"
#include "cdft.fh"
#include "sym.fh"
c
      integer nmat    ! number of J matrices to make [input]
      integer g_vc     ! GA handle for V Coulomb [output]
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl
      integer iwhat_max
      integer n_semi_bufs
      integer fd
c
      double precision CD_coef(nbf_cd,nmat)
      double precision i3c_ERI(n3c_dbl)
      double precision Ecoul2 ! [output]
c
      Logical IOLGC, dermat, do_init
c     
c     Compute the fitted Coulomb potential. The fitting coefficients are
c     obtained by solving the normal equations for a standard linear least
c     squares problem.
c
      integer me,nproc,i1,ishp,ishbf,lshbf,lmiss,ischw1,ischw2,
     ,     icount,Nao2,i,nERI,maxg,mscratch_2e3c,nscr,
     ,     n1,n2,ishc,ifirstc,ilastc,nshbfc,ishd,ifirstd,ilastd,
     ,     nshbfd,nfunz,ifirstp,ilastp,Nintegrals,next,
     ,     ntasks,nshbf_cd,imat,ninteg_api
      integer lcd_shell,icd_shell,lsvec,isvec,
     ,     lsvec2,isvec2,iERI,lERI,lscr,iscr
      double precision THRESHOLD,CD_coef_max,shmax,
     ,     Vcoul_est,ERI_est
      integer current_nbuf, n3c_b, n_batch_read
      double precision offset
      double precision length, q1
      integer ishc_be, ishd_be, nlast
      integer nxtask
      logical dowork
      logical oprint_ioinfo
      integer ndone,ipoint,ipoint2
      external nxtask
c
c !!! BGJ test !!!
c      write(*,*)'*** dft_fitvc: dermat nmat =',dermat,nmat
      if (.not.dermat .and. nmat.ne.1) then
         call errquit('nmat must be 1 for regular fit in dft_fitvc',0,0)
      endif
c !!! BGJ test !!!
      oprint_ioinfo = util_print('io info',print_debug)
c
      THRESHOLD=10.d0**(-itol2e)
c      
      me=ga_nodeid()
      nproc  = ga_nnodes()
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
      do 5 ishp = 1,nshells_cd
        if( .not. bas_cn2bfr(CD_bas_han,ishp,ishbf,lshbf) )then
          call errquit('Exiting in fitvcoul.',4,0)
        end if
        shmax = 0.D0
        do 2 i = ishbf,lshbf
           do imat = 1, nmat
              shmax = max(abs(CD_coef(i,imat)), shmax)
           enddo
    2   continue
        Dbl_MB((icd_shell-1)+ishp) = shmax
        CD_coef_max = max(shmax, CD_coef_max)
    5 continue
      if (IOLGC)then
*ga:1:0
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERInv', 
     &                       nbf_cd, 0, g_2ceri)) 
     &      call errquit('dft_fitvc: error creating g_2ceri',0,0)
        if (me.eq.0)then
           lmiss = 0
           call dft_invio('CD',g_2ceri,nbf_cd,'READ',lmiss)
        endif
        call ga_brdcst(1999, lmiss,mitob(1),0)
        if(lmiss.ne.0)call errquit
     &     (' dft_fitvc; dft_invio - abnormal read of CD ',0,0)
        if(oprint_ioinfo)then
          call ga_print(g_2ceri)
        endif
      endif
c      
c     Parallelize the loop over the products of AO shells.
c      
      ischw1 = 0
      ischw2 = 0
      
c      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      icount = nshells_ao*(nshells_ao+1)/nproc/2
c      
c     check for zero ... must be at least 1
c
      icount = max(icount,1)
      if(nproc.eq.1)icount = 1
c      

c
      Nao2 = nbf_ao_mxnbf_cn**2
      if(.not.MA_Push_Get(MT_Dbl,Nao2*nmat,'isvec',lsvec,isvec))
     &  call errquit('dft_fitvc: cannot allocate isvec',0,0)
      if(.not.MA_Push_Get(MT_Dbl,nmat*Nao2,'isvec2',lsvec2,isvec2))
     &  call errquit('dft_fitvc: cannot allocate isvec2',0,0)
c
      call int_mem_2e3c(maxg, mscratch_2e3c)
c      nERI    = Nao2*nbf_cd_mxnbf_cn
      nERI    = maxg
      if(.not.MA_Push_Get(MT_Dbl,nERI,'ERI',lERI,iERI))
     &  call errquit('dft_fitvc: cannot allocate ERI',0,0)
      nscr    = mscratch_2e3c
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_fitvc: cannot allocate scr',0,0)
c
      if (do_init) then
            call ga_zero(g_vc)
      endif
c
      if(.not.incore)then
        iwhat_max = -1
        goto 25
      endif
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
      if (n_semi_bufs.gt.0)then
c
c        check to make sure file is there
c
         if (eaf_length(fd, length) .ne. 0) call errquit
     $      ('dft_fitcd: unable to determine 3c2e file length',0,0)
         current_nbuf = current_nbuf + 1
         offset = (current_nbuf-1)*n3c_b
         if (eaf_read(fd, offset, i3c_eri, n3c_b) .ne. 0)call errquit
     &      ('dft_fitcd: 3c-2e integral buffer read failed',0,0)
         n_batch = nint(i3c_eri(n3c_dbl))

      endif
c     


c     
      do 24 ishc = 1, nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han, ishc, ifirstc, ilastc))
     &       call errquit('Exiting in fitvcoul.',3,0)
        nshbfc=ilastc-ifirstc+1
c        
        do 23 ishd = 1, ishc
c          
          dowork=iwhat(ndone).eq.((ishc*(ishc-1))/2 + ishd-1)
          if(dowork) then
c
c            write(LuOut,*)
c     &      'me,ishc,ishd,ndone,iwhat(ndone),ipoint', 
c     &       me,ishc,ishd,ndone,iwhat(ndone),ipoint
c            call util_flush(LuOut)
c
            if (.not. bas_cn2bfr( AO_bas_han, ishd, ifirstd, ilastd))
     &           call errquit('Exiting in fitvcoul.',4,0)
            nshbfd=ilastd-ifirstd+1
            Nao2 = nshbfc*nshbfd
c           
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c           
            ERI_est=schwarz_shell(ishc,ishd)

            if( CD_coef_max*ERI_est.gt.THRESHOLD )then
              ischw1 = ischw1 + nshells_cd
              
c!!!              call dfill(Nao2*nmat,0.D0,dbl_mb(isvec),1)
              nfunz=0
              ipoint2=ipoint
c                 
c             Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c                 
              if (nmat.eq.1) then
c
c                Use dgemv when only 1 mat is being made - should be
c                more efficient than dgemm
c
                 call dgemv('N',Nao2,nbf_cd,1.D0,
     &                i3c_ERI(ipoint),Nao2,
     &                CD_coef,1,0.D0,Dbl_MB(isvec),1)
c!!!     &                CD_coef,1,1.D0,Dbl_MB(isvec),1)
c
              else
c
c                Case of more than 1 matrix
c
                 call dgemm('N','N',Nao2,nmat,nbf_cd,1d0,
     &                i3c_ERI(ipoint),Nao2,CD_coef,nbf_cd,
     &                0d0,dbl_mb(isvec),Nao2)
c
              endif
c           
c             Scatter the computed values of Vcoul.
c           
              if(dermat) then
                 call errquit('fitvc: not coded yet',0,0)
              else
                 call sca_lab_mat(AO_bas_han,ishc,ishd,
     &                nshbfc,nshbfd,g_vc,
!old     &                DBL_MB(isvec+(imat-1)*Nao2),dbl_mb(isvec2),
     &                DBL_MB(isvec),dbl_mb(isvec2),
     .                'put')
              endif
c
            endif
c            
            ndone=ndone+1
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
     &           ('dft_fitcd: 3c-2e integral buffer read failed',0,0)
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
      n1 = nxtask(nproc,icount)
      n2 = 0
c
      do 240 ishc = ishc_be, nshells_ao
        if( .not. bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &    call errquit('Exiting in fitvcoul.',3,0)
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
     &         call errquit('Exiting in fitvcoul.',4,0)
            nshbfd=ilastd-ifirstd+1
            Nao2 = nshbfc*nshbfd
c           
c           Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c           
            ERI_est=schwarz_shell(ishc,ishd)

            if( CD_coef_max*ERI_est.gt.THRESHOLD )then
              ischw1 = ischw1 + nshells_cd
c              
              call dfill(Nao2*nmat,0.D0,dbl_mb(isvec),1)
              nfunz=0
              do 220 ishp = 1,nshells_cd
c                
                if( .not.bas_cn2bfr(CD_bas_han, ishp, ifirstp, ilastp))
     &            call errquit('Exiting in fitvcoul.',5,0)
                nshbf_cd=ilastp-ifirstp+1
                Nintegrals=nshbf_cd*Nao2
                Ninteg_api=maxg
                Vcoul_est = Dbl_MB((icd_shell-1)+ishp)*ERI_est
                if( Vcoul_est.gt.THRESHOLD )then
c                  
                  ischw2 = ischw2 + 1
c                   
c                 Compute 3-ctr ERIs between a given pair of AO shells and a CD shell.
c                   
c                 check for use of symmetry
c
                  if (oskel .and. .not.dermat) then
                     if (sym_shell(cd_bas_han, ishp, q1))then
                        call int_2e3c(CD_bas_han,ishp,AO_bas_han,
     &                  ishc,ishd,nscr,Dbl_MB(iscr),Ninteg_api,
     &                  Dbl_MB(iERI))
                        call dscal(Nintegrals, q1, Dbl_MB(iERI), 1)
c                   
c                       Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c
                        if (nmat.eq.1) then
                           call dgemv('N',Nao2,nshbf_cd,1.D0,
     &                          Dbl_MB(iERI),Nao2,CD_coef(1+nfunz,1),
     &                          1,1.D0,Dbl_MB(isvec),1)
                        else
                           call dgemm('N','N',Nao2,nmat,nshbf_cd,1.D0,
     &                          Dbl_MB(iERI),Nao2,CD_coef(1+nfunz,1),
     &                          nbf_cd,1.D0,Dbl_MB(isvec),Nao2)
                        endif
                     else
                        call dcopy(nintegrals,0d0,0,Dbl_MB(ieri), 1)
                     endif
                  else
                     call int_2e3c(CD_bas_han,ishp,AO_bas_han,
     &               ishc,ishd,nscr,Dbl_MB(iscr),Ninteg_api,
     &               Dbl_MB(iERI))
c                   
c                    Multiply the matrix of 3-ctr ERIs by the CD expansion coefficients.
c
                     if (nmat.eq.1) then
                        call dgemv('N',Nao2,nshbf_cd,1.D0,
     &                       Dbl_MB(iERI),Nao2,
     &                       CD_coef(1+nfunz,1),1,1.D0,Dbl_MB(isvec),1)
                     else
                        call dgemm('N','N',Nao2,nmat,nshbf_cd,1.D0,
     &                       Dbl_MB(iERI),Nao2,CD_coef(1+nfunz,1),
     &                       nbf_cd,1.D0,Dbl_MB(isvec),Nao2)
                     endif
                  endif
                endif
c
                nfunz=nfunz+nshbf_cd
c               
  220         continue
c
c             Scatter the computed values of Vcoul.
c           
              if(dermat) then
                 call sca_lab_mat3d(AO_bas_han,ishc,ishd,
     &                nshbfc,nshbfd,g_vc,nmat,
     &                DBL_MB(isvec),dbl_mb(isvec2),
     .                'put')
              else
                    call sca_lab_mat(AO_bas_han,ishc,ishd,
     &                   nshbfc,nshbfd,g_vc,
     &                   DBL_MB(isvec),dbl_mb(isvec2),
     .                   'put')
              endif
            endif
c            
c            call util_flush(LuOut)
c
            n1 = nxtask(nproc,icount)
c
          endif
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
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitvc: cannot pop stack',0,0)
      if(.not.ma_pop_stack(lERI))
     &  call errquit('dft_fitvc: cannot pop stack',0,0)
      if(.not.ma_pop_stack(lsvec2))
     &  call errquit('dft_fitvc: cannot pop stack',0,0)
      if(.not.ma_pop_stack(lsvec))
     &  call errquit('dft_fitvc: cannot pop stack',0,0)
c
c
      n1 = nxtask(-nproc,icount)
c      
      if(.not.ma_pop_stack(lcd_shell))
     &  call errquit('dft_fitvc: cannot pop stack',0,0)
c
c     Compute the Coulombic self-energy of the fitted charge density.
c      
      if (dermat) then
         write(*,*)'*** Fix self-energy calculation in dft_fitvc'
      endif
c
c     Obtain the matrix of 2-ctr ERIs.
c      
c     !!! Not needed if not doing Ecoul2
      if(.not.MA_Push_Get(MT_Dbl,nbf_cd,'scr',lscr,iscr))
     &  call errquit('dft_fitvc: cannot allocate scr',0,0)
c
c
c     Compute Ecoul2 only for the first set of fit coefficients
c     !!! Should modify this so that nothing is done if nmat > 1
c     !!! Looks like g_2ceri not needed in that case
c
      Ecoul2=0.D0
      next=nxtask(nproc, 1)+1
      do i=1,nbf_cd
        if(i.eq.next) then
          next=nxtask(nproc, 1)+1
          call get_col(g_2ceri,nbf_cd,i,DBL_MB(iscr))
          Ecoul2 = Ecoul2+ddot(
     *         nbf_cd,CD_coef,1,Dbl_MB(iscr),1)*CD_coef(i,1)
        endif
      enddo
c     
      next=nxtask(-nproc, 1)
c     

      Ecoul2=Ecoul2*0.5D0
c
c         write(*,*)'*** After destroy g_2ceri'
c      call ga_sync
      call ga_dgop(Msg_Ecoul2,Ecoul2,1,'+')
      if (IOLGC)then
         if (.not. ga_destroy(g_2ceri)) call errquit
     &   ('dft_fitvc: could not destroy g_2ceri', 0,0)
      endif         
c
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_fitvc: cannot pop stack',0,0)
c
      return
      end
