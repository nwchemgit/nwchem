      double precision function dft_n3cdbl( AO_bas_han, CD_bas_han,
     &     itol2e)
C$Id: dft_n3c.F,v 1.9 1999-06-22 01:20:43 edo Exp $
      implicit none
c****
c**** nwchem handles
c****
      integer  AO_bas_han, CD_bas_han
      integer itol2e ! bi-el int tolerance [input]
c**** 
#include "bas.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "schwarz.fh"
#include "stdio.fh"
#include "util.fh"
c     
c     compute no. of 3-ctr integrals non null 
c     according to Schwarz screening
c
      integer nshells_ao,me,nproc,ishc,ifirstc,ilastc,nshbfc,ishd,
     &        ifirstd,ilastd,nshbfd,nbf_cd,nbf_ao,mxshbf_ao,
     &        nao2_max,icount
      double precision THRESHOLD,ERI_est,n3,ischw1,batch_max
      logical oprint_3c2e
c
      oprint_3c2e = util_print('3c 2e integrals', print_default)
      THRESHOLD=10.d0**(-itol2e)
      nproc=ga_nnodes()
c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in dft_n3c.',2)
      end if
c      
      if ( .not. bas_numbf(AO_bas_han,nbf_ao) ) then
        call errquit('Exiting from dft_n3c',4)
      endif 
c
      if( .not. bas_nbf_cn_max(AO_bas_han,mxshbf_ao) )then
        call errquit('Exiting in dft_n3c.',1)
      end if
c      
      if ( .not. bas_numbf(CD_bas_han,nbf_cd) ) then
        call errquit('Exiting from dft_n3c',5)
      endif 
c      
      me=ga_nodeid()
c      
c     Determine number of 3 center 2e- integrals based on Schwarz screening.
c
      ischw1 = 0.0d0
c      
      do 240 ishc = 1,nshells_ao
c
        if( .not.bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &    call errquit('Exiting in dft_n3c.',3)
c
        nshbfc=ilastc-ifirstc+1
c        
        do 230 ishd = 1,ishc
c         
c         Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c         
          ERI_est=schwarz_shell(ishc,ishd)
c
          if( ERI_est.gt.THRESHOLD )then
c
            if (.not. bas_cn2bfr( AO_bas_han,ishd,
     &           ifirstd,ilastd))
     &           call errquit('Exiting in fitvcoul.',4)
c
            nshbfd=ilastd-ifirstd+1
            ischw1 = ischw1 + dble(nbf_cd*nshbfc*nshbfd)
c
          end if
c          
  230   continue
c
  240 continue
c
c     add max 3-center 2e- batch size
c      
c     Loops are parallelized over the products of AO shells
c     (check for zero ... must be at least 1).
c      
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
c      icount = (nshells_ao*(nshells_ao+1))/(2*nproc)
      icount = max(icount,1)
      if(nproc.eq.1)icount = 1
c      
      nao2_max = mxshbf_ao**2
      batch_max = dble(nao2_max)*dble(nbf_cd)*dble(icount)
c
c     Define n3 = maximum number of 3 center 2e- integrals. 
c
      n3=dble(nbf_ao**2)*dble(nbf_cd)
c
c     Define dft_n3cdbl=total number of non-zero 3 center 2e- ints.
c
      dft_n3cdbl=min(ischw1,n3)
c
c     Redefine dft_n3cdbl=buffer size of 3 center 2e- integrals per 
c     processor and allow for some load-unbalance, e.g., 1.5
c
      dft_n3cdbl = dft_n3cdbl/dble(nproc) + dble(batch_max) + 4096.0d0
c      dft_n3cdbl = dft_n3cdbl/dble(nproc) + dble(batch_max) + 1024.0d0
c
c     cannot be smaller than mxsh*mxsh*nbf_cd*icount
c
      dft_n3cdbl = max(dft_n3cdbl,batch_max)
c
      if (oprint_3c2e.and.me.eq.0)then
         write(LuOut,*)
         call util_print_centered
     &      (LuOut,'3 Center 2 Electron Integral Information',25,.true.)
         write(LuOut,1111)n3, ischw1, dft_n3cdbl, batch_max
      endif
 1111 format(10x,'Maximum number of 3-center 2e- integrals is: ',f15.0,
     &     /,10x,'  This is reduced with Schwarz screening to: ',f15.0,
     &     /,10x,'  Incore requires a per proc buffer size of: ',f15.0,
     &     /,10x,'        The minimum integral buffer size is: ',f15.0)
c      
      return
      end
      integer function dft_n3cint(AO_bas_han, CD_bas_han, 
     &      itol2e)
C$Id: dft_n3c.F,v 1.9 1999-06-22 01:20:43 edo Exp $
      implicit none
            
c****
c**** nwchem handles
c****
      integer AO_bas_han, CD_bas_han
      integer itol2e ! bi-el int tolerance [input]
c**** 
#include "bas.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
c
      integer nshells_ao
      integer nproc
c
      nproc=ga_nnodes()
c
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in dft_3cinc.',1)
      end if
c
      if(nproc.gt.1) then
c
c       Allow for some load-unbalance (factor of 1.2).
c
        dft_n3cint = nshells_ao*(nshells_ao+1)/nproc
        dft_n3cint=max(dft_n3cint,nshells_ao)
      else
        dft_n3cint = nshells_ao*nshells_ao
      endif
c
      return
      end
