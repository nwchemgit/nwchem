      integer function dft_n3cdbl( AO_bas_han, CD_bas_han,
     &     itol2e)
C$Id: dft_n3c.F,v 1.4 1996-11-23 00:50:43 d3h449 Exp $
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
c     
c     compute no. of 3-ctr integrals non null 
c     according to Schwarz screening
c
      integer nshells_ao
      integer me,nproc,ischw1,
     ,     ishc,ifirstc,ilastc,nshbfc,ishd,ifirstd,ilastd,
     ,     nshbfd,
     ,     nbf_cd,nbf_ao,n3,mxshbf_ao,nao2_max,batch_max
      double precision THRESHOLD,
     ,     ERI_est
c
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
      ischw1 = 0
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
            ischw1 = ischw1 + nbf_cd*nshbfc*nshbfd
c
          end if
c          
  230   continue
c
  240 continue
c
c     add max 3-center 2e- batch size
c
      nao2_max = mxshbf_ao**2
      batch_max = nao2_max*nbf_cd
c      
c     Define n3 = maximum number of 3 center 2e- integrals. 
c
      n3=nbf_ao*nbf_ao*nbf_cd
c
c     Define dft_n3cdbl=total number of non-zero 3 center 2e- ints.
c
      dft_n3cdbl=min(ischw1,n3)
c
c     Redefine dft_n3cdbl=buffer size of 3 center 2e- integrals per 
c     processor and allow for some load-unbalance, e.g., 1.5
c
      dft_n3cdbl = dft_n3cdbl/nproc + batch_max + 1024
c
c     cannot be smaller than mxsh*mxsh*nbf_cd
c
      dft_n3cdbl = max(dft_n3cdbl,4*batch_max)
c
      if(me.eq.0)write(LuOut,1111)n3, ischw1, dft_n3cdbl
 1111 format(/,10x,'Maximum number of 3-center 2e- integrals is: ',i15,
     &       /,10x,'  This is reduced with Schwarz screening to: ',i15,
     &       /,10x,'  Incore requires a per proc buffer size of: ',i15)
c      
      return
      end
      integer function dft_n3cint(AO_bas_han, CD_bas_han, 
     &      itol2e)
C$Id: dft_n3c.F,v 1.4 1996-11-23 00:50:43 d3h449 Exp $
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
        dft_n3cint = nshells_ao*nshells_ao*1.20d0/nproc
      else
        dft_n3cint = nshells_ao*nshells_ao
      endif
c
      return
      end
