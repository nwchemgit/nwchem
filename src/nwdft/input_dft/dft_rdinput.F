      Subroutine dft_rdinput(rtdb)
c
C$Id: dft_rdinput.F,v 1.39 1997-03-29 01:15:58 d3h449 Exp $
c
c     The main purpose of this subroutine is to read all parameters
c     from the rtdb (which were read and placed by the input parser), 
c     set defaults for those which were not in the rtdb and load 
c     the whole lot into the rtdb and common (cdft.fh).
c
      implicit none
c      
      integer rtdb                               ! [input] database handle
c
#include "inp.fh"
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
#include "util.fh"
c
      integer bases(3)
      integer iproc,ichg,itype,
     &     icenter,jcenter,i3
      integer lcoord,icoord,lcharge,icharge,ltags,itags,
     &     lznuc,iznuc,lratom,iratom,lrot,irot
      integer iline, n
      integer three, six, fifteen
      Parameter ( three=3, six=6, fifteen=15 )
      character*40 blank_40c
      Logical lnewtype
      Logical int_normalize
      Logical int_norm_2c
      integer num_gridtypes, ind 
      parameter (num_gridtypes = 6)

      External int_normalize
      External int_norm_2c

      character*15 gridtypes(num_gridtypes)
      data gridtypes /'coarse', 'medium', 'fine', 'xfine', 
     &                'gausleg', 'lebedev'/
c     
c
c******************************************************************************
c     
      iproc  = ga_nodeid()
      do iline=1,15
          xfac(iline)=0.d0
          cfac(iline)=0.d0
         lcfac(iline)=.false.
        nlcfac(iline)=.false.
         lxfac(iline)=.false.
        nlxfac(iline)=.false.
        xccomb(iline)=.false.
      enddo
c****
c**** recover charge from rtdb.
c****
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, rcharge)) then
        rcharge = 0.0d0
      endif 
      ichg = idint(1.000001d0 * rcharge) ! avoid bad truncation
c     
c>>>  Get parameters from rtdb. If not available, set
c     default and put in the rtdb as well as loading in common.
c     
      if (.not. rtdb_get(rtdb, 'dft:nodamping', mt_log, 1,
     &   nodamping))then
         nodamping = .false.
         if (.not. rtdb_put(rtdb, 'dft:nodamping',
     &      mt_log, 1, nodamping))
     &      call errquit('dft_input: rtdb_put failed', 1)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nolevelshifting', mt_log, 1, 
     $   nolevelshifting))then
         nolevelshifting = .false.
         if (.not. rtdb_put(rtdb, 'dft:nolevelshifting',
     &      mt_log, 1, nolevelshifting))
     &      call errquit('dft_input: rtdb_put failed', 2)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:direct', mt_log, 1,
     &   direct))then
         direct = .false.
         if (.not. rtdb_put(rtdb, 'dft:direct',
     &      mt_log, 1, direct))
     &      call errquit('dft_input: rtdb_put failed', 3)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:diis', mt_log, 1,
     &   diis))then
         diis = .true.
         if (.not. rtdb_put(rtdb, 'dft:diis',
     &      mt_log, 1, diis))
     &      call errquit('dft_input: rtdb_put failed', 3)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nodiis', mt_log, 1,
     &   nodiis))then
         nodiis = .false.
         if (.not. rtdb_put(rtdb, 'dft:nodiis',
     &      mt_log, 1, nodiis))
     &      call errquit('dft_input: rtdb_put failed', 3)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ltight', mt_log, 1,
     &   ltight))then
         ltight=.false.
         if (.not. rtdb_put(rtdb, 'dft:ltight',
     &      mt_log, 1, ltight))
     &      call errquit('dft_input: rtdb_put failed', 4)
      endif
c
      if(.not.rtdb_get(rtdb,'dft:incore', mt_log, 1,
     &   incore)) then
         incore=.true.
         if(.not.rtdb_put(rtdb,'dft:incore',
     &      mt_log, 1, incore)) 
     &      call errquit('dft_input: rtdb_put failed', 5)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:max_ovl', mt_log, 1,
     &   lmaxov ))then
         lmaxov=.false.      
         if(.not. rtdb_put(rtdb, 'dft:max_ovl',
     &      mt_log, 1, lmaxov ))
     &      call errquit('dft_input: rtdb_put failed', 6)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:store_wght', mt_log, 1,
     &   store_wght ))then
         store_wght = .false.
         if(.not. rtdb_put(rtdb, 'dft:store_wght',
     &      mt_log, 1, store_wght ))
     &      call errquit('dft_input: rtdb_put failed', 7)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:delley', mt_log, 1,
     &   ldelley ))then
         ldelley = .false.
         if(.not. rtdb_put(rtdb, 'dft:delley',
     &      mt_log, 1, ldelley ))
     &      call errquit('dft_input: rtdb_put failed', 8)
      endif
c     
      if (.not. rtdb_get(rtdb, 'dft:iterations', mt_int, 1,
     &   iterations))then
         iterations=30
         if (.not. rtdb_put(rtdb, 'dft:iterations',
     &      mt_int, 1, iterations))
     &      call errquit('dft_input: rtdb_put failed', 101) 
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:mult', mt_int, 1,
     &   mult))then
         mult=1
         if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1,
     &      mult))
     &      call errquit('dft_input: rtdb_put failed', 104)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1,
     &   ipol))then
         ipol=1
         if (.not. rtdb_put(rtdb, 'dft:ipol',
     &      mt_int, 1, ipol))
     &      call errquit('dft_input: rtdb_put failed', 105)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))then
         e_conv=1.d-6
         if (.not. rtdb_put(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))
     &      call errquit('dft_input: rtdb_put failed', 102)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:d_conv', mt_dbl, 1, d_conv))then
         d_conv=1.d-4
         if (.not. rtdb_put(rtdb, 'dft:d_conv', mt_dbl, 1, d_conv))
     &      call errquit('dft_input: rtdb_put failed', 106)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:g_conv', mt_dbl, 1, g_conv))then
         g_conv=1.d-3
         if (.not. rtdb_put(rtdb, 'dft:g_conv', mt_dbl, 1, g_conv))
     &      call errquit('dft_input: rtdb_put failed', 106)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nfock', mt_int, 1,
     &   nfock))then
         nfock=10
         if (.not. rtdb_put(rtdb, 'dft:nfock',
     &      mt_int, 1, nfock))
     &      call errquit('dft_input: rtdb_put failed', 107)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:idirect', mt_int, 1,
     &   idirect))then
         idirect=0
         if (.not. rtdb_put(rtdb, 'dft:idirect',
     &      mt_int, 1, idirect))
     &      call errquit('dft_input: rtdb_put failed', 111)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:noio', mt_int, 1,
     &   noio))then
         noio=0
         if (.not. rtdb_put(rtdb, 'dft:noio',
     &      mt_int, 1, noio))
     &      call errquit('dft_input: rtdb_put failed', 112)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ndamp', mt_int, 1,
     &   ndamp))then
         damp = .true.
         if (.not. rtdb_put(rtdb, 'dft:damp',
     &      mt_log, 1, damp))
     &      call errquit('dft_input: rtdb_put failed', 2)
         ndamp=70
         if (.not. rtdb_put(rtdb, 'dft:ndamp',
     &      mt_int, 1, ndamp))
     &      call errquit('dft_input: rtdb_put failed', 113)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:idecomp', mt_int, 1,
     &   idecomp))then
         idecomp = 0
         if (.not. rtdb_put(rtdb, 'dft:idecomp',
     &      mt_int, 1, idecomp))
     &      call errquit('dft_input: rtdb_put failed', 114)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nquad_task', mt_int, 1,
     &   nquad_task))then
         nquad_task = 1
         if (.not. rtdb_put(rtdb, 'dft:nquad_task',
     &      mt_int, 1,nquad_task))
     &      call errquit('dft_input: rtdb_put failed', 123)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1,
     &   iAOacc))then
         iAOacc=20
         if (.not. rtdb_put(rtdb, 'dft:iAOacc',
     &      mt_int, 1, iAOacc))
     &      call errquit('dft_input: rtdb_put failed', 124)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:iXCacc', mt_int, 1,
     &   iXCacc))then
         iXCacc=20
         if (.not. rtdb_put(rtdb, 'dft:iXCacc',
     &      mt_int, 1, iXCacc))
     &      call errquit('dft_input: rtdb_put failed', 125)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:irqmax', mt_int, 1,
     &   irqmax))then
         irqmax=40
         if (.not. rtdb_put(rtdb, 'dft:irqmax',
     &      mt_int, 1, irqmax))
     &      call errquit('dft_input: rtdb_put failed', 126)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1,
     &   itol2e))then
         itol2e=12
         if (.not. rtdb_put(rtdb, 'dft:itol2e',
     &      mt_int, 1, itol2e))
     &      call errquit('dft_input: rtdb_put failed', 127)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:imull', mt_int, 1,
     &   imull))then
         imull=0
         if (.not. rtdb_put(rtdb, 'dft:imull',
     &      mt_int, 1, imull))
     &      call errquit('dft_input: rtdb_put failed', 128)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:iangquad', mt_int, 1,
     &   iangquad))then
         iangquad=0
         if (.not. rtdb_put(rtdb, 'dft:iangquad',
     &      mt_int, 1, iangquad))
     &      call errquit('dft_input: rtdb_put failed', 129)
      endif
c     
      if (.not. rtdb_get(rtdb, 'dft:hl_tol', mt_dbl, 1,
     &   hl_tol))then
         hl_tol = 2.d-1
         if (.not. rtdb_put(rtdb, 'dft:hl_tol',
     &      mt_dbl, 1, hl_tol))
     &      call errquit('dft_input: rtdb_put failed', 201)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:dampon', mt_dbl, 1,
     &   dampon))then
         dampon = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:dampon',
     &      mt_dbl, 1, dampon))
     &      call errquit('dft_input: rtdb_put failed', 201)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:dampoff', mt_dbl, 1,
     &   dampoff))then
         dampoff = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:dampoff',
     &      mt_dbl, 1, dampoff))
     &      call errquit('dft_input: rtdb_put failed', 202)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:levlon', mt_dbl, 1,
     &   levlon))then
         levlon = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:levlon',
     &      mt_dbl, 1, levlon))
     &      call errquit('dft_input: rtdb_put failed', 203)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:levloff', mt_dbl, 1,
     &   levloff))then
         levloff = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:levloff',
     &      mt_dbl, 1, levloff))
     &      call errquit('dft_input: rtdb_put failed', 204)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:diison', mt_dbl, 1,
     &   diison))then
         diison = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:diison',
     &      mt_dbl, 1, diison))
     &      call errquit('dft_input: rtdb_put failed', 205)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:diisoff', mt_dbl, 1,
     &   diisoff))then
         diisoff = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:diisoff',
     &      mt_dbl, 1, diisoff))
     &      call errquit('dft_input: rtdb_put failed', 206)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:rscale', mt_dbl, 1,
     &   rscale))then
         rscale = 0.0d0
         if (.not. rtdb_put(rtdb, 'dft:rscale',
     &      mt_dbl, 1, rscale))
     &      call errquit('dft_input: rtdb_put failed', 207)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:rlshift', mt_dbl, 1,
     &   rlshift))then
         levelshift = .true.
         if (.not. rtdb_put(rtdb, 'dft:levelshift',
     &      mt_log, 1, levelshift))
     &      call errquit('dft_input: rtdb_put failed', 2)
         rlshift=0.5d0
         if (.not. rtdb_put(rtdb, 'dft:rlshift',
     &      mt_dbl, 1,rlshift))
     &      call errquit('dft_input: rtdb_put failed', 212)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ncydp', mt_int, 1,
     &   ncydp))then
         ncydp = 2
         if (.not. rtdb_put(rtdb, 'dft:ncydp',
     &      mt_int, 1, ncydp))
     &      call errquit('dft_input: rtdb_put failed', 108)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ncysh', mt_int, 1,
     &   ncysh))then
         if(levelshift)then
           ncysh = iterations
         else
           ncysh = 0
         endif
         if (.not. rtdb_put(rtdb, 'dft:ncysh',
     &      mt_int, 1, ncysh))
     &      call errquit('dft_input: rtdb_put failed', 109)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ncyds', mt_int, 1,
     &   ncyds))then
         if(diis)then
           ncyds = iterations
         else
           ncyds = 0
         endif
         if (.not. rtdb_put(rtdb, 'dft:ncyds',
     &      mt_int, 1, ncyds))
     &      call errquit('dft_input: rtdb_put failed', 110)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:radius', mt_dbl, 1,
     &   radius))then
         radius=16.d0
         if (.not. rtdb_put(rtdb, 'dft:radius',
     &      mt_dbl, 1,radius))
     &      call errquit('dft_input: rtdb_put failed', 213)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:tol_rho', mt_dbl, 1,
     &   tol_rho))then
         tol_rho=1.D-15
         if (.not. rtdb_put(rtdb, 'dft:tol_rho',
     &      mt_dbl, 1,tol_rho))
     &      call errquit('dft_input: rtdb_put failed', 214)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:rm', mt_dbl, 1,
     &   rm))then
         rm = 2.0d0
         if (.not. rtdb_put(rtdb, 'dft:rm',
     &      mt_dbl, 1,rm))
     &      call errquit('dft_input: rtdb_put failed', 215)
      endif
c
c     get/set gridtype; set nagrid and nrad
c
      if (.not. rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype))then
c     
c        set default grid to  'MEDIUM'
c     
         gridtype='medium'
      endif
      if (.not. inp_match(num_gridtypes, .false., gridtype, gridtypes, 
     &   ind)) then
         write(6,*)' Unknown gridtype; defaulting to medium grid'
         gridtype='medium'
         ind = 2
      endif
      goto (10, 20, 30, 40, 50, 60, 99) ind
      call errquit('dft_rdinput: unknown gridtype', ind)
   10 continue
c
c     coarse grid;
c
      nagrid = 8
      nradpts = 35
      goto 99
   20 continue
c
c     medium grid;
c
      nagrid = 10
      nradpts = 50
      goto 99
   30 continue
c
c     fine grid;
c
      nagrid = 12
      nradpts = 75
      goto 99
   40 continue
c
c     xfine grid;
c
      nagrid = 14
      nradpts = 105
      goto 99
   50 continue
c
c     gausleg grid;
c
      if (.not. rtdb_get(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &   nagrid = 10
      if (.not. rtdb_get(rtdb,'dft:nradpts', mt_int, 1, nradpts))
     &   nradpts = 50
      goto 99
   60 continue
c
c     lebedev grid;
c
      if (.not. rtdb_get(rtdb, 'dft:iangquad', mt_int, 1, iangquad))then
         iangquad = 4
      else
         if(iangquad.lt.1)then
           write(6,*)' IANGQUAD must be in the range: 1-7; resetting. '
           iangquad = 1
         endif
         if(iangquad.gt.7)then
           write(6,*)' IANGQUAD must be in the range: 1-7; resetting. '
           iangquad = 7
         endif
      endif
      if (.not. rtdb_get(rtdb,'dft:nradpts', mt_int, 1, nradpts))
     &   nradpts = 50
   99 continue
      do n = 1, dft_ntags_bsmx
         nrad(n) = nradpts
      enddo 
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('dft_input: rtdb_put failed', 301) 
c
      if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &   call errquit('dft_input: rtdb_put failed', 116)
c
      if (.not. rtdb_put(rtdb,'dft:nrad', mt_int, dft_ntags_bsmx, nrad))
     &   call errquit('dft_input: rtdb_put failed', 115)
c
c     Get/set/put XC info.
c
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, 15,
     &   cfac))then
         if (.not. rtdb_put(rtdb, 'dft:cfac',
     &      mt_dbl, 15, cfac))
     &      call errquit('dft_input: rtdb_put failed', 210)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, 15,
     &   xfac))then
         if (.not. rtdb_put(rtdb, 'dft:xfac',
     &      mt_dbl, 15, xfac))
     &      call errquit('dft_input: rtdb_put failed', 211)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:lcfac', mt_log, 15,
     &   lcfac ))then
         if(.not. rtdb_put(rtdb, 'dft:lcfac',
     &      mt_log, 15, lcfac ))
     &      call errquit('dft_input: rtdb_put failed', 9)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:nlcfac', mt_log, 15,
     &   nlcfac ))then
         if(.not. rtdb_put(rtdb, 'dft:nlcfac',
     &      mt_log, 15, nlcfac ))
     &      call errquit('dft_input: rtdb_put failed', 10)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:lxfac', mt_log, 15,
     &   lxfac ))then
         if(.not. rtdb_put(rtdb, 'dft:lxfac',
     &      mt_log, 15, lxfac ))
     &      call errquit('dft_input: rtdb_put failed', 11)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:nlxfac', mt_log, 15,
     &   nlxfac ))then
         if(.not. rtdb_put(rtdb, 'dft:nlxfac',
     &      mt_log, 15, nlxfac ))
     &      call errquit('dft_input: rtdb_put failed', 12)
      endif
c
      if(.not. rtdb_get(rtdb, 'dft:xccomb', mt_log, 15,
     &   xccomb ))then
         if(.not. rtdb_put(rtdb, 'dft:xccomb',
     &      mt_log, 15, xccomb ))
     &      call errquit('dft_input: rtdb_put failed', 13)
      endif
c
c     Put XC names in rtdb.
c
      blank_40c = '                                        '
      do iline=1,15
         cname(iline)=blank_40c
         xname(iline)=blank_40c
        xcname(iline)=blank_40c
      enddo
c
c     Correlation functional name defaults.
c
      cname(1)='            VWN V Correlation Functional'
      cname(2)='    Lee-Yang-Parr Correlation Functional'
      cname(3)='      Perdew 1981 Correlation Functional'
      cname(4)='      Perdew 1986 Correlation Functional'
      cname(5)='      Perdew 1991 Correlation Functional'
      cname(6)='  Perdew 1991 LDA Correlation Functional'
      cname(7)='        VWN I RPA Correlation Functional'
      cname(8)='            VWN I Correlation Functional'
      cname(9)='           VWN II Correlation Functional'
      cname(10)='          VWN III Correlation Functional'
      cname(11)='           VWN IV Correlation Functional'
c
c     Exchange functional name defaults.
c
      xname(1)='           Hartree-Fock (Exact) Exchange'
      xname(2)='              Slater Exchange Functional'
      xname(3)='          Becke 1988 Exchange Functional'
c
c     Exchange-Correlation combination functional name defaults.
c
      xcname(1)='Adiabatic Connection Method XC Potential'
      xcname(2)='               B3LYP Method XC Potential'
      xcname(3)=' Becke half-and-half Method XC Potential'
c     
c>>>  place character XC parameters in rtdb
c     
      if (.not. rtdb_cput(rtdb, 'dft:cname', 
     &   15, cname))
     &   call errquit('dft_input: rtdb_put failed', 302) 
      if (.not. rtdb_cput(rtdb, 'dft:xname', 
     &   15, xname))
     &   call errquit('dft_input: rtdb_put failed', 303) 
      if (.not. rtdb_cput(rtdb, 'dft:xcname', 
     &   15, xcname))
     &   call errquit('dft_input: rtdb_put failed', 304) 
c
c     load geometry and symmetry info
c
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('rdinput: geom_create?',70)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('rdinput: no geometry ',72)
      if (.not.rtdb_get(rtdb, 'dft:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif
      if (.not.rtdb_get(rtdb, 'dft:adapt',MT_LOG, 1, oadapt)) then
         oadapt = sym_number_ops(geom) .gt. 0
      endif
c
c****
c**** allocate space for M and T (rotational invariance; 
c***  see B.G. Johnson et al, CPL, 220, 377 (1994), and auxiliary 
c***  species R_a - T and (R_a -T)(R_a - T)^T
c****
      if ( .not. Ma_Push_Get(MT_Dbl,fifteen,'rotation',lrot,
     &   irot))call errquit('dft_rdinput: failed to alloc rotation',0)
c
      call dzero(Dbl_MB(irot),fifteen)
c
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
c     Define the contracted Gaussian functions for the AO basis set.
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
c
c     Define the contracted Gaussian functions for the CD basis set.
c
      if (.not. bas_create(CD_bas_han, 'cd basis'))
     &   call errquit('rdinput: bas_create?',102)
      CDFIT=bas_rtdb_load(rtdb, geom, CD_bas_han, 'cd basis')
      if(CDFIT) then
        if (.not. int_norm_2c(CD_bas_han))
     &     call errquit('rdinput: int_norm_2c failed', 208)
        if (.not. bas_numbf(CD_bas_han, nbf_cd))
     &     call errquit('rdinput:rdinput:',86)
      else
        if (.not. bas_destroy(CD_bas_han))
     &     call errquit('rdinput:not able to destroy CD_bas:',86)
        CD_bas_han=-99999
        nbf_cd=0
      endif
c
c     Define the contracted Gaussian functions for the XC basis set.
c
      if (.not. bas_create(XC_bas_han, 'xc basis'))
     &   call errquit('rdinput: bas_create?',124)
      XCFIT=bas_rtdb_load(rtdb, geom, XC_bas_han, 'xc basis')
      if(XCFIT) then
        if (.not. int_normalize(XC_bas_han))
     &     call errquit('rdinput: int_normalize failed', 208)
        if (.not. bas_numbf(XC_bas_han, nbf_xc))
     &     call errquit('rdinput:rdinput:',86)
      else
        if (.not. bas_destroy(XC_bas_han))
     &     call errquit('rdinput:not able to destroy XC_bas:',88)
        XC_bas_han=-99999
        nbf_xc=0
      endif
c****   
c**** Read/set/put the numerical quadrature input parameters.
c****   
      nrinc = 1
      nrmax = 0
      do n = 1,dft_ntags_bsmx
        nrmax = max(nrmax,nrad(n))
      enddo 
      if (.not. rtdb_put(rtdb, 'dft:nrinc',
     &   mt_int, 1, nrinc))
     &   call errquit('dft_input: rtdb_put failed', 118)
      if (.not. rtdb_put(rtdb, 'dft:nrmax',
     &   mt_int, 1, nrmax))
     &   call errquit('dft_input: rtdb_put failed', 119)
      if(nrmax.gt.MXrad)then
        write(LuOut,*)'RDINPUT: Too many radial quadrature points.'
        call errquit(' rdinput: raise MXrad',1)
      end if
      if(nagrid.gt.MXagrid)then
        write(LuOut,*)'RDINPUT: Too many angular quadrature grids.'
        call errquit(' rdinput: raise MXagrid',3)
      end if
c
c     Set the switch to rotate the quadrature points.
c     This switch controls whether the grid is
c     rotationally invariant. For now, it is hard-wired to 
c     1 (rotationally invariant).
c     
      iqrot = 1
      if (.not. rtdb_put(rtdb, 'dft:iqrot',
     &   mt_int, 1, iqrot))
     &   call errquit('dft_input: rtdb_put failed', 117)
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
      call xc_setquad(rtdb,Dbl_MB(iznuc),Dbl_MB(iratom),dbl_mb(icoord),
     &                dbl_mb(icharge),dbl_mb(irot),dbl_mb(irot+three),
     &                dbl_mb(irot+six))
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
      if ( .not. MA_Pop_Stack(lrot))
     &   call errquit('dft_rdinput: pop stack failed.',0)
c
c     Check the user supplied input for any inconsistencies.
c
      call dft_inpana(rtdb)
c
c     inpana could have modified ipol
c
      if (.not. rtdb_put(rtdb, 'dft:ipol', mt_int, 1, ipol ))
     $     call errquit('inpana: dft:ipol put failed', 0)
c
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
      call int_app_set_no_texas(rtdb)
      call int_init(rtdb, i3, bases)
c
      return
      end
