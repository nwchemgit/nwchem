      Subroutine dft_3cincor(rtdb, geom, AO_bas_han, CD_bas_han, 
     &      nbf_cd, itol2e,n_batch,n_bmax,iwhat,i3c_ERI,n3c_dbl)
C$Id: dft_3cincor.F,v 1.2 1996-07-03 20:58:56 d3h449 Exp $
      implicit none
            
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han
      integer nbf_cd    ! no. CD functions [input]
      integer itol2e ! bi-el int tolerance [input]
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches [in]
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl ! max no. ints [in]
      double precision i3c_ERI(n3c_dbl) ! int values [out]
c**** 

      Logical LResult
      
#include "bas.fh"
#include "rtdb.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "schwarz.fh"

c     
c     compute the 3-centres integrals and store the in-core
c
      integer MXshbf_ao,nshells_ao,nshells_cd,MXshbf_cd
      integer me,nproc,i1,ishp,ishbf,lshbf,lmiss,
     ,     icount,Nao2,i,nERI,maxg,mscratch_2e3c,nscr,
     ,     n1,n2,ishc,ifirstc,ilastc,nshbfc,ishd,ifirstd,ilastd,
     ,     nshbfd,nfunz,ifirstp,ilastp,nint,niter,nsh_blocks,
     ,     ntasks,nshbf_cd,ndone,nb_done
      integer lscr,iscr
      double precision THRESHOLD,shmax,
     ,     Vcoul_est,ERI_est,schw1,schw2
      integer nxtask
      double precision fudge
      parameter (fudge=100.d0)
      external nxtask
C      write(0,*) ' in 3CINCO',itol2e
c
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
c     Parallelize the loop over the products of AO shells.
c      
      
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
c      
c     check for zero ... must be at least 1
c
      icount = max(icount,1)
      if(nproc.eq.1)icount = 1
c      
      Nao2 = MXshbf_ao**2
      nERI    = Nao2*MXshbf_cd
      call int_mem_2e3c(maxg, mscratch_2e3c)
      nscr    = mscratch_2e3c
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_3cincor: cannot allocate scr',0)
      call ga_sync
      n1 = nxtask(nproc,icount)
      n2 = 0
c
c      
      
      ndone=1
      n_batch=0
      do 240 ishc = 1,nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &          call errquit('exiting in fitvcoul.',3)
        nshbfc=ilastc-ifirstc+1
 
        do 230 ishd = 1,ishc
c            
c         screen integral blocks using the schwarz inequality:  (p|cd) .le. (cd|cd)
c         
          eri_est=schwarz_shell(ishc,ishd)/fudge

          if( eri_est.gt.threshold )then

            if( n1.eq.n2 )then
              
c              
c             define the scatter index.
c              
                if (.not. bas_cn2bfr( ao_bas_han,ishd,
     &             ifirstd,ilastd))
     &             call errquit('exiting in fitvcoul.',4)

                nshbfd=ilastd-ifirstd+1
                nao2 = nshbfc*nshbfd
                           
              do 220 ishp = 1,nshells_cd
                
                   if (.not. bas_cn2bfr( cd_bas_han,ishp,
     &                  ifirstp,ilastp))
     &                  call errquit('exiting in fitvcoul.',5)
                   nshbf_cd=ilastp-ifirstp+1
                  
                   nint=nshbf_cd*nao2
c                 
c              compute 3-ctr eris between a given pair of ao shells and a cd shell.
c
                  if((ndone+nint-1).gt.n3c_dbl) then
                    write(0,*) ' nlast n3cdbl ',ndone+nint-1,n3c_dbl
                  call errquit(' 3cincor: change factor in 3c ',999)
                    endif
                  call int_2e3c(cd_bas_han,ishp,ao_bas_han,ishc,ishd,
     &                 nscr,dbl_mb(iscr),nint,i3c_ERI(ndone))
C                  write(0,*) ' i3c ',ndone,i3c_ERI(ndone)
C                  call flush(0)
                  ndone=ndone+nint
 220         continue

             n_batch=n_batch+1
             if(n_batch.gt.n_bmax)
     .            call errquit(' 3cinc: n_bmax too small ',n_bmax)
C             iwhat(n_batch)=n2
             iwhat(n_batch)=ishc*(ishc-1)+ishd-1
            
             n1 = nxtask(nproc,icount)

            end if
            
          n2 = n2 + 1
         
            
          end if
  230   continue
  240 continue
      
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_3cincor: cannot pop stack',0)
c
c
      n1 = nxtask(-nproc,icount)
      
      call ga_sync
c      
      return
      end




