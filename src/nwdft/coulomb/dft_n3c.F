      integer function dft_n3cdbl( AO_bas_han, CD_bas_han,
     &     itol2e)
C$Id: dft_n3c.F,v 1.1 1996-05-01 23:37:05 og845 Exp $
      implicit none
            
c****
c**** nwchem handles
c****
      integer  geom, AO_bas_han, CD_bas_han
      integer itol2e ! bi-el int tolerance [input]
c**** 

      Logical LResult
      
#include "bas.fh"
#include "tcgmsg.fh"
#include "global.fh"

#include "schwarz.fh"

c     
c     compute no. of 3-ctr integrals non null 
c     according to Schwarz screening
c
      integer nshells_ao
      integer me,nproc,ishp,ishbf,lshbf,lmiss,ischw1,
     ,     i,
     ,     ishc,ifirstc,ilastc,nshbfc,ishd,ifirstd,ilastd,
     ,     nshbfd,Nint,nsh_blocks,
     ,     nbf_cd,nbf_ao,n3,MXshbf
      double precision THRESHOLD,fudge,
     ,     ERI_est
      parameter(fudge=100.d0)
c
      THRESHOLD=10.d0**(-itol2e)/fudge
      nproc=ga_nnodes()
c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
      
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in fitvcoul.',2)
      end if
      
      if ( .not. bas_numbf(AO_bas_han,nbf_ao) ) then
        call errquit('Exiting from fitvcoul',4)
      endif 
      
      if ( .not. bas_numbf(CD_bas_han,nbf_cd) ) then
        call errquit('Exiting from fitvcoul',5)
      endif 
      
      me=ga_nodeid()
      
      ischw1 = 0
c
c      
      do 240 ishc = 1,nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han,ishc,
     &       ifirstc,ilastc))
     &       call errquit('Exiting in fitvcoul.',3)
        nshbfc=ilastc-ifirstc+1
        
        do 230 ishd = 1,ishc
          
c         
c         Screen integral blocks using the Schwarz inequality:  (p|cd) .le. (cd|cd)
c         
          ERI_est=schwarz_shell(ishc,ishd)

          if( ERI_est.gt.THRESHOLD )then

            if (.not. bas_cn2bfr( AO_bas_han,ishd,
     &           ifirstd,ilastd))
     &           call errquit('Exiting in fitvcoul.',4)

            nshbfd=ilastd-ifirstd+1
            ischw1 = ischw1 + nbf_cd*nshbfc*nshbfd

          end if
          
          
  230   continue
  240 continue
      
      n3=nbf_ao*nbf_ao*nbf_cd
      dft_n3cdbl=min(ischw1,n3)
      if(nproc.gt.1) then
c
c       allow for some load-unbalance
c
      if( .not. bas_nbf_cn_max(AO_bas_han,MXshbf) )then
        call errquit('Exiting in fitcd.',1)
      end if
        dft_n3cdbl = dft_n3cdbl*1.50d0/nproc
C
C       cannot be smaller than mxsh*mxsh*nbf_cd
C
        dft_n3cdbl = max(dft_n3cdbl,4*mxshbf*mxshbf*nbf_cd)
      endif
C      if(me.eq.0) then
C        write(6,*) ' schwarz screen and N^3'
C        write(6,*) ischw1,n3
C        write(6,*)
C      endif
c      
      return
      end
      integer function dft_n3cint(AO_bas_han, CD_bas_han, 
     &      itol2e)
C$Id: dft_n3c.F,v 1.1 1996-05-01 23:37:05 og845 Exp $
      implicit none
            
c****
c**** nwchem handles
c****
      integer AO_bas_han, CD_bas_han
      integer itol2e ! bi-el int tolerance [input]
c**** 

      Logical LResult
      
#include "bas.fh"
#include "rtdb.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"

      integer nshells_ao
      integer nproc

      nproc=ga_nnodes()

      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in dft_3cinc.',1)
      end if

      if(nproc.gt.1) then
c
c       allow for some load-unbalance
c
        dft_n3cint = nshells_ao*nshells_ao*1.20d0/nproc
      else
        dft_n3cint = nshells_ao*nshells_ao
      endif


      return
      end
