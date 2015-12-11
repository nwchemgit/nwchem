      double precision function dft_n3cdbl()
C$Id$
      implicit none
#include "errquit.fh"
c****
c**** nwchem handles
c****
!      integer  AO_bas_han, CD_bas_han
!      integer itol2e ! bi-el int tolerance [input]
c**** 
#include "cdft.fh"
#include "bas.fh"
#include "basP.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "schwarz.fh"
#include "stdio.fh"
#include "util.fh"
c     
c     compute no. of 3-ctr integrals non null 
c     according to Schwarz screening
c
      integer me,nproc,ishc,ifirstc,ilastc,nshbfc,ishd,
     &        ifirstd,ilastd,nshbfd,mxcebf_ao
      integer maxg, mscratch_2e3c
      integer itype,nprimo,isphere
      integer sh_lo_d, sh_hi_d, atom_d
      integer sh_lo_c, sh_hi_c, atom_c
      double precision THRESHOLD,ERI_est,n3,ischw1
      double precision imb_fact
      logical oprint_3c2e,spherical_ao,spherical_cd,
     C     oprint_3c2eh
      double precision dbldum
      integer intdum,ndone
      integer n_batch, n_bmax
      integer dft_nao2_max
      external dft_nao2_max
c
      oprint_3c2e = util_print('3c 2e integrals', print_default)
      oprint_3c2eh = util_print('3c 2e integrals details', print_high)
      THRESHOLD=10.d0**(-itol2e)
      nproc=ga_nnodes()
      spherical_ao=bas_spherical(ao_bas_han)
      spherical_cd=bas_spherical(cd_bas_han)
c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
c
      if( .not. bas_nbf_ce_max(AO_bas_han,mxcebf_ao) )then
        call errquit('Exiting in dft_n3c.',1, BASIS_ERR)
      end if
c      
c      
      me=ga_nodeid()
c      
c     Determine number of 3 center 2e- integrals based on Schwarz screening.
c
      ischw1 = 0.0d0
      call int_mem_2e3c(maxg, mscratch_2e3c)
c      
      do atom_c=1,ncenters
         if (.not. bas_ce2cnr( ao_bas_han, atom_c, sh_lo_c, sh_hi_c)) 
     &        call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
      do ishc = sh_lo_c,sh_hi_c
c
         if(spherical_ao) then
            if(.not.bas_continfo(ao_bas_han,ishc,
     &           itype,nprimo,nshbfc,isphere))
     &           call errquit('Exiting in fitcd.',44, BASIS_ERR)
            nshbfc=((itype+1)*(itype+2))/2*nshbfc
         else
            if( .not.bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &           call errquit('Exiting in dft_n3c.',3, BASIS_ERR)
            nshbfc=ilastc-ifirstc+1
         endif
c        
      do atom_d=1,atom_c
         if (.not. bas_ce2cnr( ao_bas_han, atom_d, sh_lo_d, sh_hi_d)) 
     &        call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
        do ishd = sh_lo_d, sh_hi_d
c         
c          Schwarz Screen integral blocks:  (p|cd) .le. (cd|cd)
c         
          ERI_est=schwarz_shell(ishc,ishd)
c
          if( ERI_est.gt.THRESHOLD )then
c
             if(spherical_ao) then
                if(.not.bas_continfo(ao_bas_han,ishd,
     &               itype,nprimo,nshbfd,isphere))
     &               call errquit('Exiting in fitcd.',44, BASIS_ERR)
                nshbfd=((itype+1)*(itype+2))/2*nshbfd
             else
                if (.not. bas_cn2bfr( AO_bas_han,ishd,
     &               ifirstd,ilastd))
     &               call errquit('Exiting in fitvcoul.',4, BASIS_ERR)
                nshbfd=ilastd-ifirstd+1
             endif
            ischw1 = ischw1 + dble(nshbfc*nshbfd)
c
          end if
c          
       enddo ! ishd
c
      enddo ! ishc
      enddo
      enddo
      ischw1 = ischw1*nbf_cd
c      
c
c     Define n3 = maximum number of 3 center 2e- integrals. 
c
      n3=dble(nbf_ao**2)*dble(nbf_cd)
c
c     Define dft_n3cdbl=total number of non-zero 3 center 2e- ints.
c
      dft_n3cdbl=min(ischw1,n3)
c
c     Redefine dft_n3cdbl=buffer size of 3 center 2e- integrals 
c
      call dft_getcdt(.true.,ndone,
     N     ncenters, 1, .false., threshold, .true.,
     N     intdum, 100000000, intdum,
     &     intdum, dbldum,
     T     dbldum, intdum)
      if(oprint_3c2eh) then
         write(6,*) ga_nodeid(),' ndone ',ndone
      endif
      dft_n3cdbl = ndone

c
      if (oprint_3c2e.and.me.eq.0)then
         write(LuOut,*)
         call util_print_centered
     &      (LuOut,'3 Center 2 Electron Integral Information',25,.true.)
         write(LuOut,1111)n3, ischw1, 
     D        dft_n3cdbl
      endif
 1111 format(10x,'Maximum number of 3-center 2e- integrals is:',f17.0,
     &     /,10x,'  This is reduced with Schwarz screening to:',f17.0,
     &     /,10x,'  Incore requires a per proc buffer size of:',f17.0)
c      
      return
      end
      integer function dft_n3cint()
C$Id$
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "cdft.fh"
c
      integer nproc
c
      nproc=ga_nnodes()
c
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in dft_3cinc.',1, BASIS_ERR)
      end if
c
      if(nproc.gt.1) then
c
c       Allow for some load-imbalance (factor of 2).
c
        dft_n3cint = ncenters*(ncenters+1)/nproc
        dft_n3cint = dft_n3cint*2
        dft_n3cint=max(dft_n3cint,ncenters)
#ifdef DEBUG
        if(ga_nodeid().eq.0) then
     C       write(6,*) ' dft_n3cint0 ',dft_n3cint, ' nctr0 ',ncenters
     C       write(6,*) ' dft_n3cint ',dft_n3cint, ' nctr ',ncenters
        endif
#endif
      else
         dft_n3cint = ncenters*(ncenters+1)
      endif
      return
      end
