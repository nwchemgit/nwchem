      Subroutine dft_rdinput(rtdb)
c
C$Id: dft_rdinput.F,v 1.16 1996-10-15 01:21:40 d3h449 Exp $
      implicit none
c      
      integer rtdb                               ! [input] database handle
c
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "sym.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer bases(3)
      integer iproc,ichg,itype,
     &     icenter,jcenter,iEx,iEc,i3
      integer nshells_cd,nshells_xc,nshells
      integer lcoord,icoord,lcharge,icharge,ltags,itags,
     &     lznuc,iznuc,lratom,iratom
      Logical LResult, lnewtype, oskel
      Logical int_normalize
      Logical int_norm_2c
c     
c     XC functionals coefficients
c     
      External int_normalize
      External int_norm_2c
c
c******************************************************************************

      iproc  = ga_nodeid()
c****
c**** recover charge from rtdb.
c****
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, rcharge)) then
        rcharge = 0.0d0
      endif 
      ichg = idint(1.000001d0 * rcharge) ! avoid bad truncation
c     
c>>>  Get logical parameters from rtdb.
c     
      if (.not. rtdb_get(rtdb, 'dft:nodamping', mt_log, 1, nodamping))
     $     call errquit('dft_input: rtdb_get failed', 1)
      if (.not. rtdb_get(rtdb, 'dft:nolevelshifting', mt_log, 1, 
     $  nolevelshifting))call errquit('dft_input: rtdb_get failed', 2)
      if (.not. rtdb_get(rtdb, 'dft:nodiis', mt_log, 1, nodiis))
     $     call errquit('dft_input: rtdb_get failed', 3)
      if (.not. rtdb_get(rtdb, 'dft:tight', mt_log, 1, ltight))
     $     call errquit('dft_input: rtdb_get failed', 4)
      if(.not.rtdb_get(rtdb,'dft:incore', mt_log, 1, incore)) 
     $     call errquit('dft_input: rtdb_get failed', 5)
      if(.not. rtdb_get(rtdb, 'dft:max_ovl', mt_log, 1, lmaxov ))
     $     call errquit('dft_input: rtdb_get failed', 6)
      if(.not. rtdb_get(rtdb, 'dft:store_wght', mt_log, 1, store_wght ))
     $     call errquit('dft_input: rtdb_get failed', 7)
      if(.not. rtdb_get(rtdb, 'dft:delley', mt_log, 1, ldelley ))
     $     call errquit('dft_input: rtdb_get failed', 8)
      if(.not. rtdb_get(rtdb, 'dft:lcfac', mt_log, 10, lcfac ))
     $     call errquit('dft_input: rtdb_get failed', 9)
      if(.not. rtdb_get(rtdb, 'dft:nlcfac', mt_log, 10, nlcfac ))
     $     call errquit('dft_input: rtdb_get failed', 10)
      if(.not. rtdb_get(rtdb, 'dft:lxfac', mt_log, 10, lxfac ))
     $     call errquit('dft_input: rtdb_get failed', 11)
      if(.not. rtdb_get(rtdb, 'dft:nlxfac', mt_log, 10, nlxfac ))
     $     call errquit('dft_input: rtdb_get failed', 12)
      if(.not. rtdb_get(rtdb, 'dft:xccomb', mt_log, 10, xccomb ))
     $     call errquit('dft_input: rtdb_get failed', 13)
c     
c>>>  Get integer parameters from rtdb.
c     
      if (.not. rtdb_get(rtdb, 'dft:itrscf', mt_int, 1, itrscf))
     $     call errquit('dft_input: rtdb_get failed', 101) 
      if (.not. rtdb_get(rtdb, 'dft:iscfcon', mt_int, 1, iscfcon))
     $     call errquit('dft_input: rtdb_get failed', 102)
      if (.not. rtdb_get(rtdb, 'dft:iguess', mt_int, 1, iguess))
     $     call errquit('dft_input: rtdb_get failed', 103)
      if (.not. rtdb_get(rtdb, 'dft:mult', mt_int, 1, mult))
     $     call errquit('dft_input: rtdb_get failed', 104)
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
     $     call errquit('dft_input: rtdb_get failed', 105)
      if (.not. rtdb_get(rtdb, 'dft:igcon', mt_int, 1, igcon))
     $     call errquit('dft_input: rtdb_get failed', 106)
      if (.not. rtdb_get(rtdb, 'dft:nfock', mt_int, 1, nfock))
     $     call errquit('dft_input: rtdb_get failed', 107)
      if (.not. rtdb_get(rtdb, 'dft:ncydp', mt_int, 1, ncydp))
     $     call errquit('dft_input: rtdb_get failed', 108)
      if (.not. rtdb_get(rtdb, 'dft:ncysh', mt_int, 1, ncysh))
     $     call errquit('dft_input: rtdb_get failed', 109)
      if (.not. rtdb_get(rtdb, 'dft:ncyds', mt_int, 1, ncyds))
     $     call errquit('dft_input: rtdb_get failed', 110)
      if (.not. rtdb_get(rtdb, 'dft:idirect', mt_int, 1, idirect))
     $     call errquit('dft_input: rtdb_get failed', 111)
      if (.not. rtdb_get(rtdb, 'dft:noio', mt_int, 1, noio))
     $     call errquit('dft_input: rtdb_get failed', 112)
      if (.not. rtdb_get(rtdb, 'dft:ndamp', mt_int, 1, ndamp))
     $     call errquit('dft_input: rtdb_get failed', 113)
      if (.not. rtdb_get(rtdb, 'dft:idecomp', mt_int, 1, idecomp))
     $     call errquit('dft_input: rtdb_get failed', 114)
      if (.not. rtdb_get(rtdb,'dft:nrad',mt_int,dft_ntags_bsmx,nrad))
     $     call errquit('dft_input: rtdb_get failed', 115)
      if (.not. rtdb_get(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     $     call errquit('dft_input: rtdb_get failed', 116)
      if (.not. rtdb_get(rtdb, 'dft:iqrot', mt_int, 1, iqrot))
     $     call errquit('dft_input: rtdb_get failed', 117)
      if (.not. rtdb_get(rtdb, 'dft:nrinc', mt_int, 1, nrinc))
     $     call errquit('dft_input: rtdb_get failed', 118)
      if (.not. rtdb_get(rtdb, 'dft:nrmax', mt_int, 1, nrmax))
     $     call errquit('dft_input: rtdb_get failed', 119)
      if (.not. rtdb_get(rtdb, 'dft:igeotype', mt_int, 1, igeotype))
     $     call errquit('dft_input: rtdb_get failed', 120)
      if (.not. rtdb_get(rtdb, 'dft:numgrp', mt_int, 1, numgrp))
     $     call errquit('dft_input: rtdb_get failed', 121)
      if (.not. rtdb_get(rtdb, 'dft:numset', mt_int, 1, numset))
     $     call errquit('dft_input: rtdb_get failed', 122)
      if (.not. rtdb_get(rtdb, 'dft:nquad_task', mt_int, 1,nquad_task))
     $     call errquit('dft_input: rtdb_get failed', 123)
      if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc))
     $     call errquit('dft_input: rtdb_get failed', 124)
      if (.not. rtdb_get(rtdb, 'dft:iXCacc', mt_int, 1, iXCacc))
     $     call errquit('dft_input: rtdb_get failed', 125)
      if (.not. rtdb_get(rtdb, 'dft:irqmax', mt_int, 1, irqmax))
     $     call errquit('dft_input: rtdb_get failed', 126)
      if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1, itol2e))
     $     call errquit('dft_input: rtdb_get failed', 127)
      if (.not. rtdb_get(rtdb, 'dft:imull', mt_int, 1, imull))
     $     call errquit('dft_input: rtdb_get failed', 128)
      if (.not. rtdb_get(rtdb, 'dft:iangquad', mt_int, 1, iangquad))
     $     call errquit('dft_input: rtdb_get failed', 129)
c     
c>>>  Get real parameters from rtdb.
c     
      if (.not. rtdb_get(rtdb, 'dft:dampon', mt_dbl, 1, dampon))
     $     call errquit('dft_input: rtdb_get failed', 201)
      if (.not. rtdb_get(rtdb, 'dft:dampoff', mt_dbl, 1, dampoff))
     $     call errquit('dft_input: rtdb_get failed', 202)
      if (.not. rtdb_get(rtdb, 'dft:levlon', mt_dbl, 1, levlon))
     $     call errquit('dft_input: rtdb_get failed', 203)
      if (.not. rtdb_get(rtdb, 'dft:levloff', mt_dbl, 1, levloff))
     $     call errquit('dft_input: rtdb_get failed', 204)
      if (.not. rtdb_get(rtdb, 'dft:diison', mt_dbl, 1, diison))
     $     call errquit('dft_input: rtdb_get failed', 205)
      if (.not. rtdb_get(rtdb, 'dft:diisoff', mt_dbl, 1, diisoff))
     $     call errquit('dft_input: rtdb_get failed', 206)
      if (.not. rtdb_get(rtdb, 'dft:rscale', mt_dbl, 1, rscale))
     $     call errquit('dft_input: rtdb_get failed', 207)
      if (.not. rtdb_get(rtdb, 'dft:cang', mt_dbl, 3, cang))
     $     call errquit('dft_input: rtdb_get failed', 208)
      if (.not. rtdb_get(rtdb, 'dft:cdist', mt_dbl, 3, cdist))
     $     call errquit('dft_input: rtdb_get failed', 209)
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, 10, cfac))
     $     call errquit('dft_input: rtdb_get failed', 210)
      if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, 10, xfac))
     $     call errquit('dft_input: rtdb_get failed', 211)
      if (.not. rtdb_get(rtdb, 'dft:lshift', mt_dbl, 1,rlshift))
     $     call errquit('dft_input: rtdb_get failed', 212)
      if (.not. rtdb_get(rtdb, 'dft:radius', mt_dbl, 1,radius))
     $     call errquit('dft_input: rtdb_get failed', 213)
      if (.not. rtdb_get(rtdb, 'dft:tol_rho', mt_dbl, 1,TOLL))
     $     call errquit('dft_input: rtdb_get failed', 214)
      if (.not. rtdb_get(rtdb, 'dft:rm', mt_dbl, 1,rm))
     $     call errquit('dft_input: rtdb_get failed', 215)
c     
c>>>  Get character parameters from rtdb.
c     
      if (.not. rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype))
     $     call errquit('dft_input: rtdb_cget failed', 301) 
      if (.not. rtdb_cget(rtdb, 'dft:cname', 10, cname))
     $     call errquit('dft_input: rtdb_cget failed', 302) 
      if (.not. rtdb_cget(rtdb, 'dft:xname', 10, xname))
     $     call errquit('dft_input: rtdb_cget failed', 303) 
      if (.not. rtdb_cget(rtdb, 'dft:xcname', 10, xcname))
     $     call errquit('dft_input: rtdb_cget failed', 304) 
c
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('rdinput: geom_create?',70)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('rdinput: no geometry ',72)
c
      itype = 1
      numgrp = 1
c****
c**** get number of atomic centers from geometry object
c****
      if (.not. geom_ncent(geom, ncenters))
     &     call errquit('rdinput: geom_ncent failed',73)      
c****
c**** allocate space for atomic coordinates and charges
c****
      if ( .not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &   icoord))call errquit('dft_rdinput: failed to alloc 
     &   coordinates',0)
      if ( .not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &   icharge))call errquit('dft_rdinput: failed to alloc 
     &   charges',0)
      if ( .not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit('dft_rdinput: failed to alloc 
     &   center tags',0)
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &     call errquit('rdinput: geom_cart_get failed',74)
c
c     Extract high level info from the data-base setting defaults
c
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     &     title = ' '
c
c     end lines from ddscf/scf
c     
      if (.not. bas_create(AO_bas_han, 'ao basis'))
     &     call errquit('rdinput: bas_create?',80)

      if (.not. bas_rtdb_load(rtdb, geom, AO_bas_han, 'ao basis'))
     &     call errquit('rdinput: no ao basis ',82)

      if (.not. int_normalize(AO_bas_han))
     &     call errquit('rdinput: int_normalize failed', 200)


      if (.not. bas_numbf(AO_bas_han, nbf))
     &     call errquit('rdinput:rdinput:',86)
      nbf_ao = nbf
c**** 
c       
c     Broadcast AO basis set to all nodes.
c****
c     Define the contracted Gaussian functions for the CD basis set.

      if (.not. bas_create(CD_bas_han, 'cd basis'))
     &   call errquit('rdinput: bas_create?',102)
      CDFIT=bas_rtdb_load(rtdb, geom, CD_bas_han, 'cd basis')
      if(CDFIT) then
        if (.not. int_norm_2c(CD_bas_han))
     &     call errquit('rdinput: int_norm_2c failed', 208)
        if (.not. bas_numbf(CD_bas_han, nbf_cd))
     &     call errquit('rdinput:rdinput:',86)
        if (.not. bas_numcont(CD_bas_han, nshells_cd))
     &     call errquit('rdinput:rdinput:',87)
      else
        if (.not. bas_destroy(CD_bas_han))
     &     call errquit('rdinput:not able to destroy CD_bas:',86)
        CD_bas_han=-99999
        nbf_cd=0
      endif
c       
      if (.not. bas_create(XC_bas_han, 'xc basis'))
     &   call errquit('rdinput: bas_create?',124)
      XCFIT=bas_rtdb_load(rtdb, geom, XC_bas_han, 'xc basis')
      if(XCFIT) then
        if (.not. int_normalize(XC_bas_han))
     &     call errquit('rdinput: int_normalize failed', 208)
        if (.not. bas_numbf(XC_bas_han, nbf_xc))
     &     call errquit('rdinput:rdinput:',86)
        if (.not. bas_numcont(XC_bas_han, nshells_xc))
     &     call errquit('rdinput:rdinput:',88)
      else
        if (.not. bas_destroy(XC_bas_han))
     &     call errquit('rdinput:not able to destroy XC_bas:',88)
        XC_bas_han=-99999
        nbf_xc=0
      endif
c****   
c****   Read the numerical quadrature input parameters.
c****   
      if(iproc.eq.0)then
        rscale = 0.0d0
        if(nrmax.gt.MXrad)then
          write(LuOut,*)'RDINPUT: Too many radial quadrature points.'
          call errquit(' rdinput: raise MXrad',1)
        end if
        if(nagrid.gt.MXagrid)then
          write(LuOut,*)'RDINPUT: Too many angular quadrature grids.'
          call errquit(' rdinput: raise MXagrid',3)
        end if
      end if
      call ga_sync
c****   
c****   generate number of atom types and atom type array iatype(icenter) 
c****   
      ntypes = 0
      do icenter = 1,ncenters
c****     
c****     is this a new type of atom?
c****     
        lnewtype = .true.
        do jcenter = 1,icenter-1
          if ( Dbl_MB(icharge+icenter-1) .eq.
     &      Dbl_MB(icharge+jcenter-1) ) then   ! same atom type
            lnewtype = .false.
            iatype(icenter)=iatype(jcenter)
            goto 8
          endif
        enddo 
    8   continue
        if ( lnewtype ) then
          ntypes = ntypes + 1
          iatype(icenter)=ntypes
        endif 
      enddo 
      if( ntypes.gt.dft_ntags_bsmx )then
        write(LuOut,*) 'RDINPUT:  Too many types of atoms.'
        call errquit(' rdinput: raise dft_ntags_bsmx',2)
      end if
c****
c**** set up type-indexed znuc array
c****
      if ( .not. Ma_Push_Get(MT_Dbl,ntypes,'nuc charge',lznuc,iznuc))
     &   call errquit('dft_rdinput: failed to alloc nuc charge',0)
c****
c**** for each type, place the nuclear charge in znuc array
c****
      do itype = 1,ntypes
        do icenter = 1,ncenters
          if ( iatype(icenter) .eq. itype ) then
c****       center icenter is of type itype; assign charge
            Dbl_MB(iznuc+itype-1) = dble( Dbl_MB(icharge+icenter-1) )
            goto 9 ! next type
          endif 
        enddo 
    9   continue 
      enddo 

      if ( .not. Ma_Push_Get(MT_Dbl,ntypes,'ratom',lratom,iratom))
     &     call errquit('dft_rdinput: failed to alloc ratom',0)
      nqshells = MXqshells
c
      call xc_setquad(rtdb,Dbl_MB(iznuc),Dbl_MB(iratom),dbl_mb(icoord))
c
      if ( .not. MA_Pop_Stack(lratom))
     &   call errquit('dft_rdinput: pop stack failed.',0)
      if ( .not. MA_Pop_Stack(lznuc))
     &   call errquit('dft_rdinput: pop stack failed.',0)
      if ( .not. MA_Pop_Stack(ltags))
     &   call errquit('dft_rdinput: pop stack failed.',0)
      if ( .not. MA_Pop_Stack(lcharge))
     &   call errquit('dft_rdinput: pop stack failed.',0)
      if ( .not. MA_Pop_Stack(lcoord))
     &   call errquit('dft_rdinput: pop stack failed.',0)
c
      if (.not. bas_numcont(AO_bas_han, nshells))
     &     call errquit('rdinput:rdinput:',86)
c
c     Check the user supplied input for any inconsistencies.
c
      call dft_inpana(rtdb)
c
c     inpana could have modified ipol
c
      if (.not. rtdb_put(rtdb, 'dft:ipol', MT_INT, 1, ipol ))
     $     call errquit('inpana: dft:ipol put failed', 0)
c
      if( iproc.eq.0 )then
        write(LuOut,*)
        write(LuOut,1) 
     &       '       AO basis - No. of functions:   ',nbf,
     &       '                  No. of shells   :   ',nshells
        if(nbf_cd.ne.0) then
          write(LuOut,1) 
     &         '       CD basis - No. of functions:   ',nbf_cd,
     &         '                  No. of shells   :   ',nshells_cd
        endif
        if(nbf_xc.ne.0) then
          write(LuOut,1) 
     &         '       XC basis - No. of functions:   ',nbf_xc,
     &         '                  No. of shells   :   ',nshells_xc
        endif
        write(LuOut,*)
        write(LuOut,1) 
     &       '       No. of quadrature shells   :   ',nqshells
    1   format(2(3x,A,1x,i6/))
        call flush(LuOut)
      endif
      call ga_sync
c****
c**** initialize the integral object
c****
      if(XCFIT.and.CDFIT) then
        i3 = 3
        bases(1) = AO_bas_han
        bases(2) = CD_bas_han
        bases(3) = XC_bas_han
      elseif((.not.XCFIT).and.CDFIT) then
        i3 = 2
        bases(1) = AO_bas_han
        bases(2) = CD_bas_han
      elseif((.not.CDFIT).and.XCFIT) then
        i3 = 2
        bases(1) = AO_bas_han
        bases(2) = XC_bas_han
      else
        i3 = 1
        bases(1) = AO_bas_han
      endif
c
      call int_init(rtdb, i3, bases)
c
      return
      end
