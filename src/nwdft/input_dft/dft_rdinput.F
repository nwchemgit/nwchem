      Subroutine dft_rdinput(rtdb)
c
C$Id: dft_rdinput.F,v 1.139 2008-12-06 23:09:50 niri Exp $ dft_rdinput.F,v 1.55 1997/12/09 23:34:54 d3h449 Exp $
c
c     The main purpose of this subroutine is to read all parameters
c     from the rtdb (which were read and placed by the input parser), 
c     set defaults for those which were not in the rtdb and load 
c     the whole lot into the rtdb and common (cdft.fh).
c
      implicit none
#include "errquit.fh"
c      
      integer rtdb                               ! [input] database handle
c
#include "inp.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "cdft.fh"
#include "oep.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "sym.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "steric.fh"
c
      integer itype, icenter, jcenter
      integer lcoord, icoord, lcharge, icharge, ltags, itags,
     &        lznuc, iznuc
      integer iline
      integer nopen
      character*32 theory 
      character*40 blank_40c
      logical same_atom, same_bq, isbq
      Logical lnewtype
      Logical int_normalize
      Logical int_norm_2c
      Logical atom_tag_check
      integer me
c
      External int_normalize
      External int_norm_2c
      External atom_tag_check
c
      me=ga_nodeid()

      if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
     $     call errquit('dft_rdinput: failed getting theory',0,
     &       RTDB_ERR)
c
      do iline=1,numfunc
          xfac(iline)=0.d0
          cfac(iline)=0.d0
         lcfac(iline)=.false.
        nlcfac(iline)=.false.
         lxfac(iline)=.false.
        nlxfac(iline)=.false.
        xccomb(iline)=.false.
      enddo
c
c     recover charge from rtdb.
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, rcharge))then
        rcharge = 0.0d0
      endif 
c     
c     Get parameters from rtdb. If not available, set
c     default and put in the rtdb as well as loading in common.
c     
      if (.not. rtdb_get(rtdb, 'dft:nodamping', mt_log, 1,
     &   nodamping))then
         nodamping = .false.
         if (.not. rtdb_put(rtdb, 'dft:nodamping',
     &      mt_log, 1, nodamping))
     &      call errquit('dft_rdinput: rtdb_put failed', 1, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nolevelshifting', mt_log, 1, 
     $   nolevelshifting))then
         nolevelshifting = .false.
         if (.not. rtdb_put(rtdb, 'dft:nolevelshifting',
     &      mt_log, 1, nolevelshifting))
     &      call errquit('dft_rdinput: rtdb_put failed', 2, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:direct', mt_log, 1,
     &   direct))then
         direct = .false.
         if (.not. rtdb_put(rtdb, 'dft:direct',
     &      mt_log, 1, direct))
     &      call errquit('dft_rdinput: rtdb_put failed', 3, RTDB_ERR)
      endif
c
c     load dlb info
c
      call nxtask_init(rtdb)
c
      if (.not. rtdb_get(rtdb, 'dft:diis', mt_log, 1,
     &   diis))then
         diis = .true.
         if (.not. rtdb_put(rtdb, 'dft:diis',
     &      mt_log, 1, diis))
     &      call errquit('dft_rdinput: rtdb_put failed', 3, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nodiis', mt_log, 1,
     &   nodiis))then
         nodiis = .false.
         if (.not. rtdb_put(rtdb, 'dft:nodiis',
     &      mt_log, 1, nodiis))
     &      call errquit('dft_rdinput: rtdb_put failed', 3, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ltight', mt_log, 1,
     &   ltight))then
         ltight=.false.
         if (.not. rtdb_put(rtdb, 'dft:ltight',
     &      mt_log, 1, ltight))
     &      call errquit('dft_rdinput: rtdb_put failed', 4, RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'dft:incore', mt_log, 1,
     &   incore))then
         incore=.true.
         if (.not.rtdb_put(rtdb,'dft:incore',
     &      mt_log, 1, incore)) 
     &      call errquit('dft_rdinput: rtdb_put failed', 5, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:max_ovl', mt_log, 1,
     &   lmaxov))then
c
c     try to see if we see it with vectors lock
c
      if (.not. rtdb_get(rtdb, 'dft:movecs lock', mt_log, 1,
     &   lmaxov))  lmaxov=.false.      
         if (.not. rtdb_put(rtdb, 'dft:max_ovl',
     &      mt_log, 1, lmaxov))
     &      call errquit('dft_rdinput: rtdb_put failed', 6, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:lb94', mt_log, 1,
     &   lb94))then
         lb94 = .false.
         if (.not. rtdb_put(rtdb, 'dft:lb94',
     &      mt_log, 1, lb94))
     &      call errquit('dft_rdinput: rtdb_put failed', 9, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:cs00', mt_log, 1,
     &   cs00))then
         cs00 = .false.
         if (.not. rtdb_put(rtdb, 'dft:cs00',
     &      mt_log, 1, cs00))
     &      call errquit('dft_rdinput: rtdb_put failed', 10, RTDB_ERR)
      endif
c     
      if (.not. rtdb_get(rtdb, 'dft:iterations', mt_int, 1,
     &   iterations))then
         iterations=30
         if (.not. rtdb_put(rtdb, 'dft:iterations',
     &      mt_int, 1, iterations))
     &      call errquit('dft_rdinput: rtdb_put failed', 101, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:mult', mt_int, 1,
     &   mult))then
         mult=1
         if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1,
     &      mult))
     &      call errquit('dft_rdinput: rtdb_put failed', 104, RTDB_ERR)
         nopen = mult - 1
         if (.not.rtdb_put(rtdb, 'dft:nopen',mt_int,1,nopen))
     &       call errquit('dft_rdinput: rtdb_put failed', 9104,
     &       RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))then
         if (theory .eq. 'sodft') then
            ipol = 2
         else
            ipol=1
         endif
         if (.not. rtdb_put(rtdb, 'dft:ipol',
     &      mt_int, 1, ipol))
     &      call errquit('dft_rdinput: rtdb_put failed', 105,
     &       RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))then
         e_conv=1.d-6
         if (.not. rtdb_put(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))
     &      call errquit('dft_rdinput: rtdb_put failed', 102, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:d_conv', mt_dbl, 1, d_conv))then
         d_conv=1.d-5
         if (.not. rtdb_put(rtdb, 'dft:d_conv', mt_dbl, 1, d_conv))
     &      call errquit('dft_rdinput: rtdb_put failed', 106, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:g_conv', mt_dbl, 1, g_conv))then
         g_conv=5.d-4
         if (.not. rtdb_put(rtdb, 'dft:g_conv', mt_dbl, 1, g_conv))
     &      call errquit('dft_rdinput: rtdb_put failed', 106, RTDB_ERR)
      endif
c
      if (.not.rtdb_get(rtdb,'dft:delta_ac',mt_dbl,1,delta_ac))then
         delta_ac = 1.0d99
         if (.not.rtdb_put(rtdb,'dft:delta_ac',mt_dbl,1,delta_ac))
     &      call errquit('dft_rdinput: rtdb_put failed', 107, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nfock', mt_int, 1,
     &   nfock))then
         nfock=10
         if (.not. rtdb_put(rtdb, 'dft:nfock',
     &      mt_int, 1, nfock))
     &      call errquit('dft_rdinput: rtdb_put failed', 107, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:noio', mt_int, 1,
     &   noio))then
         noio=0
         if (.not. rtdb_put(rtdb, 'dft:noio',
     &      mt_int, 1, noio))
     &      call errquit('dft_rdinput: rtdb_put failed', 112, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:damp', mt_log, 1, damp))then
         damp = .false.
         if (.not. rtdb_put(rtdb, 'dft:damp', mt_log, 1, damp))
     &      call errquit('dft_rdinput: rtdb_put failed', 2, RTDB_ERR)
      endif
      if (.not. rtdb_get(rtdb, 'dft:ndamp', mt_int, 1, ndamp))then
         ndamp = 0
         if (.not. rtdb_put(rtdb, 'dft:ndamp', mt_int, 1, ndamp))
     &      call errquit('dft_rdinput: rtdb_put failed', 113, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:idecomp', mt_int, 1,
     &   idecomp))then
         idecomp = 0
         if (.not. rtdb_put(rtdb, 'dft:idecomp',
     &      mt_int, 1, idecomp))
     &      call errquit('dft_rdinput: rtdb_put failed', 114, RTDB_ERR)
      endif
      if (.not. rtdb_get(rtdb, 'dft:nquad_task', mt_int, 1,
     &   nquad_task))then
         nquad_task = 1
         if (.not. rtdb_put(rtdb, 'dft:nquad_task',
     &      mt_int, 1,nquad_task))
     &      call errquit('dft_rdinput: rtdb_put failed', 123, RTDB_ERR)
      endif
c
c
      if (.not. rtdb_get(rtdb, 'dft:iCDacc', mt_int, 1,
     &   iCDacc))then
         iCDacc=20
         if (.not. rtdb_put(rtdb, 'dft:iCDacc',
     &      mt_int, 1, iCDacc))
     &      call errquit('dft_rdinput: rtdb_put failed', 125, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:iXCacc', mt_int, 1,
     &   iXCacc))then
         iXCacc=20
         if (.not. rtdb_put(rtdb, 'dft:iXCacc',
     &      mt_int, 1, iXCacc))
     &      call errquit('dft_rdinput: rtdb_put failed', 125, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1,
     &   itol2e))then
         if(theory.eq.'sodft') then
            itol2e = 10
         else
            itol2e=8
         endif
         if (.not. rtdb_put(rtdb, 'dft:itol2e',
     &      mt_int, 1, itol2e))
     &      call errquit('dft_rdinput: rtdb_put failed', 127, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:imull', mt_int, 1,
     &   imull))then
         imull=0
         if (.not. rtdb_put(rtdb, 'dft:imull',
     &      mt_int, 1, imull))
     &      call errquit('dft_rdinput: rtdb_put failed', 128, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:iangquad', mt_int, 1,
     &   iangquad))then
         iangquad=0
         if (.not. rtdb_put(rtdb, 'dft:iangquad',
     &      mt_int, 1, iangquad))
     &      call errquit('dft_rdinput: rtdb_put failed', 129, RTDB_ERR)
      endif
c     
      if (.not. rtdb_get(rtdb, 'dft:hl_tol', mt_dbl, 1,
     &   hl_tol))then
         hl_tol = 0.01d0
         if (.not. rtdb_put(rtdb, 'dft:hl_tol',
     &      mt_dbl, 1, hl_tol))
     &      call errquit('dft_rdinput: rtdb_put failed', 201, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:dampon', mt_dbl, 1,
     &   dampon))then
         dampon = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:dampon',
     &      mt_dbl, 1, dampon))
     &      call errquit('dft_rdinput: rtdb_put failed', 201, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:dampoff', mt_dbl, 1,
     &   dampoff))then
         dampoff = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:dampoff',
     &      mt_dbl, 1, dampoff))
     &      call errquit('dft_rdinput: rtdb_put failed', 202, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:levlon', mt_dbl, 1,
     &   levlon))then
         levlon = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:levlon',
     &      mt_dbl, 1, levlon))
     &      call errquit('dft_rdinput: rtdb_put failed', 203, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:levloff', mt_dbl, 1,
     &   levloff))then
         levloff = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:levloff',
     &      mt_dbl, 1, levloff))
     &      call errquit('dft_rdinput: rtdb_put failed', 204, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:diison', mt_dbl, 1,
     &   diison))then
         diison = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:diison',
     &      mt_dbl, 1, diison))
     &      call errquit('dft_rdinput: rtdb_put failed', 205, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:diisoff', mt_dbl, 1,
     &   diisoff))then
         diisoff = 0.d0
         if (.not. rtdb_put(rtdb, 'dft:diisoff',
     &      mt_dbl, 1, diisoff))
     &      call errquit('dft_rdinput: rtdb_put failed', 206, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:levelshift',
     &   mt_log, 1, levelshift))then
         levelshift = .true.
         if (.not. rtdb_put(rtdb, 'dft:levelshift',
     &      mt_log, 1, levelshift))
     &      call errquit('dft_rdinput: rtdb_put failed', 212, RTDB_ERR)
      endif
      if (.not. rtdb_get(rtdb, 'dft:rlshift',
     &   mt_dbl, 1,rlshift))then
         rlshift=0.5d0
         if (.not. rtdb_put(rtdb, 'dft:rlshift',
     &      mt_dbl, 1,rlshift))
     &      call errquit('dft_rdinput: rtdb_put failed', 212, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ncydp', mt_int, 1,
     &   ncydp))then
         ncydp = 2
         if (.not. rtdb_put(rtdb, 'dft:ncydp',
     &      mt_int, 1, ncydp))
     &      call errquit('dft_rdinput: rtdb_put failed', 108, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:ncysh', mt_int, 1,
     &   ncysh))then
         if (levelshift)then
           ncysh = iterations
         else
           ncysh = 0
         endif
c
c     case with multiple tasks and we changed iterations ... (increased)
c
      else  if(ncysh.eq.30.and.iterations.gt.30) then
         ncysh=iterations
      endif
      if (.not. rtdb_put(rtdb, 'dft:ncysh',
     &     mt_int, 1, ncysh))
     &     call errquit('dft_rdinput: rtdb_put failed', 109, RTDB_ERR)
c     
      if (.not. rtdb_get(rtdb, 'dft:ncyds', mt_int, 1,
     &   ncyds))then
         if (diis)then
           ncyds = iterations
         else
           ncyds = 0
         endif
c
c     case with multiple tasks and we changed iterations ... (increased)
c
      else  if(ncyds.eq.30.and.iterations.gt.30) then
         ncyds=iterations
      endif
      if (.not. rtdb_put(rtdb, 'dft:ncyds',
     &     mt_int, 1, ncyds))
     &     call errquit('dft_rdinput: rtdb_put failed', 110, RTDB_ERR)
c
      if (.not. rtdb_get(rtdb, 'dft:tol_rho', mt_dbl, 1,
     &   tol_rho))then
         tol_rho=1.D-10
         if (.not. rtdb_put(rtdb, 'dft:tol_rho',
     &      mt_dbl, 1,tol_rho))
     &      call errquit('dft_rdinput: rtdb_put failed', 214, RTDB_ERR)
      endif
c
c
c     Get/set/put XC info.
c
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, numfunc,
     &   cfac))then
         if (.not. rtdb_put(rtdb, 'dft:cfac',
     &      mt_dbl, numfunc, cfac))
     &      call errquit('dft_rdinput: rtdb_put failed', 210, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc,
     &   xfac))then
         if (.not. rtdb_put(rtdb, 'dft:xfac',
     &      mt_dbl, numfunc, xfac))
     &      call errquit('dft_rdinput: rtdb_put failed', 211, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:lcfac', mt_log, numfunc,
     &   lcfac))then
         if (.not. rtdb_put(rtdb, 'dft:lcfac',
     &      mt_log, numfunc, lcfac))
     &      call errquit('dft_rdinput: rtdb_put failed', 9, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nlcfac', mt_log, numfunc,
     &   nlcfac))then
         if (.not. rtdb_put(rtdb, 'dft:nlcfac',
     &      mt_log, numfunc, nlcfac))
     &      call errquit('dft_rdinput: rtdb_put failed', 10, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:lxfac', mt_log, numfunc,
     &   lxfac))then
         if (.not. rtdb_put(rtdb, 'dft:lxfac',
     &      mt_log, numfunc, lxfac))
     &      call errquit('dft_rdinput: rtdb_put failed', 11, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nlxfac', mt_log, numfunc,
     &   nlxfac))then
         if (.not. rtdb_put(rtdb, 'dft:nlxfac',
     &      mt_log, numfunc, nlxfac))
     &      call errquit('dft_rdinput: rtdb_put failed', 12, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:xccomb', mt_log, numfunc,
     &   xccomb))then
         if (.not. rtdb_put(rtdb, 'dft:xccomb',
     &      mt_log, numfunc, xccomb))
     &      call errquit('dft_rdinput: rtdb_put failed', 13, RTDB_ERR)
      endif
c
      lsteric = .false.
      if (rtdb_get(rtdb, 'dft:steric', mt_log, 1, lsteric))then
       lsteric = .true.
       if (.not. rtdb_put(rtdb, 'dft:steric', mt_log,
     &   1, lsteric))
     &   call errquit('dft_rdinput: rtdb_put failed', 0, RTDB_ERR)
      endif
c
c     Put XC names in rtdb.
c
      blank_40c = '                                        '
      do iline=1,numfunc
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
      cname(12)=' PerdewBurkeErnz. Correlation Functional'
      cname(13)='            HCTH  Correlation Functional'
      cname(14)='        Becke 1997 Correlation Potential'
      cname(15)='      Becke 1997-1 Correlation Potential'
      cname(16)='         HCTH120  Correlation Functional'
      cname(17)='         HCTH147  Correlation Functional'
      cname(18)='       Becke 1998 Correlation Functional'
      cname(19)='   Becke 1997GGA-1 Correlation Potential'
      cname(20)='         HCTH407  Correlation Functional'
      cname(21)='         HCTHP1/4 Correlation Functional'
      cname(22)='            OPTC     Exchange Functional'
      cname(23)='      Becke 1997-2 Correlation Potential'
      cname(24)='  Filatov Thiel 1997 Method C Functional'
      cname(25)='   PKZB99 metaGGA Correlation Functional'
      cname(26)='         HCTH407+ Correlation Functional'
      cname(27)='   TPSS03 metaGGA Correlation Functional'
      cname(28)='      Becke 1997-3 Correlation Potential'
      cname(31)='               B95 Correlation Potential'
      cname(32)='            PW6B95 Correlation Potential'
      cname(33)='             PWB6K Correlation Potential'
      cname(34)='               M05 Correlation Potential'
      cname(35)='            M05-2X Correlation Potential'
      cname(36)='        OP_Becke88 Correlation Potential'
      cname(37)='          OP_PBE96 Correlation Potential'
      cname(38)='            M06-HF Correlation Potential'
      cname(39)='               M06 Correlation Potential'
      cname(40)='            M06-2X Correlation Potential'
      cname(41)='      Becke 1997-D Correlation Potential'
      cname(46)='    simplified PBE Correlation Potential'
      cname(47)='        ------  NO Correlation ---------'
c
c     Exchange functional name defaults.
c
      xname(1)='           Hartree-Fock (Exact) Exchange'
      xname(2)='              Slater Exchange Functional'
      xname(3)='          Becke 1988 Exchange Functional'
      xname(4)='               HCTH  Exchange Functional'
      xname(5)='          Becke 1997 Exchange Functional'
      xname(6)='        Becke 1997-1 Exchange Functional'
      xname(7)='PerdewBurkeErnzerhof Exchange Functional'
      xname(8)='         Gill 1996   Exchange Functional'
      xname(9)='       Perdew 1991   Exchange Functional'
      xname(10)='            HCTH120  Exchange Functional'
      xname(11)='            HCTH147  Exchange Functional'
      xname(12)='         Becke 1998  Exchange Functional'
      xname(13)='     Becke 1997GGA-1 Exchange Functional'
      xname(14)='            HCTH407  Exchange Functional'
      xname(15)='           HCTHP1/4  Exchange Functional'
      xname(16)='            OPTX     Exchange Functional'
      xname(17)='       Becke 1997-2  Exchange Functional'
      xname(18)='        PKZB metaGGA Exchange Functional'
      xname(19)='              mPW91  Exchange Functional'
      xname(20)='           HCTH407+  Exchange Functional'
      xname(21)='        TPSS metaGGA Exchange Functional'
      xname(22)='       Becke 1997-3  Exchange Functional'
      xname(25)='  Filatov Thiel 1997 Method X Functional'
      xname(26)='          PW6B95 GGA Exchange Functional'
      xname(27)='           PWB6K GGA Exchange Functional'
      xname(28)='         M05 metaGGA Exchange Functional'
      xname(29)='      M05-2X metaGGA Exchange Functional'
      xname(30)='                RPBE Exchange Functional'
      xname(31)='              RevPBE Exchange Functional'
      xname(33)='                VS98 Exchange Functional'
      xname(34)='               M06-L Exchange Functional'
      xname(35)='              M06-HF Exchange Functional'
      xname(36)='                 M06 Exchange Functional'
      xname(37)='              M06-2X Exchange Functional'
      xname(38)='         LC-BNL 2007 Exchange Functional'
      xname(39)='        Becke 1997-D Exchange Functional'
      xname(40)='         CAM-Becke88 Exchange Functional'
      xname(41)='           CAM-PBE96 Exchange Functional'
      xname(46)='               SSB-1 Exchange Functional'
      xname(47)='                 KT1 gradient correction'
c
c     Exchange-Correlation combination functional name defaults.
c
      xcname(1)='Adiabatic Connection Method XC Potential'
      xcname(2)='               B3LYP Method XC Potential'
      xcname(3)=' Becke half-and-half Method XC Potential'
      xcname(4)='          Becke 1997 Method XC Potential'
      xcname(5)='        Becke 1997-1 Method XC Potential'
      xcname(6)='              HCTH  Method XC Functional'
      xcname(7)='           HCTH120  Method XC Functional'
      xcname(8)='           HCTH147  Method XC Functional'
      xcname(9)='         Becke 1998 Method XC Functional'
      xcname(10)='               PBE0 Method XC Functional'
      xcname(11)='     Becke 1997GGA-1 Method XC Potential'
      xcname(12)='            HCTH407 Method XC Functional'
      xcname(13)='           HCTHP1/4 Method XC Functional'
      xcname(14)='         Becke 97-2 Method XC Functional'
      xcname(15)=' Filatov Thiel 1997 Method XC Functional'
      xcname(16)='              MPW1K Method XC Functional'
      xcname(17)='           HCTH407+ Method XC Functional'
      xcname(18)='        TPSS hybrid Method XC Functional'
      xcname(19)='         Becke 97-3 Method XC Functional'
      xcname(21)='            MPW1B95 Method XC Functional'
      xcname(22)='             MPWB1K Method XC Functional'
      xcname(23)='              B1B95 Method XC Functional'
      xcname(24)='               BB1K Method XC Functional'
      xcname(25)='             PW6B95 Method XC Functional'
      xcname(26)='              PWB6K Method XC Functional'
      xcname(27)='                BOP Method XC Functional'
      xcname(28)='              PBEOP Method XC Functional'
      xcname(29)='                M05 Method XC Functional'
      xcname(30)='             M05-2X Method XC Functional'
      xcname(32)='               VS98 Method XC Functional'
      xcname(33)='              M06-L Method XC Functional'
      xcname(34)='             M06-HF Method XC Functional'
      xcname(35)='                M06 Method XC Functional'
      xcname(36)='             M06-2X Method XC Functional'
      xcname(37)='         Becke 97-D Method XC Functional'
      xcname(38)='             LC-BNL Method XC Functional'
      xcname(46)='              SSB-D Method XC Functional'
      xcname(47)='                KT2 Method XC Functional'
c     
c     place character XC parameters in rtdb
c     
      if (.not. rtdb_cput(rtdb, 'dft:cname', 
     &   numfunc, cname))
     &   call errquit('dft_rdinput: rtdb_put failed', 302, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'dft:xname', 
     &   numfunc, xname))
     &   call errquit('dft_rdinput: rtdb_put failed', 303, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'dft:xcname', 
     &   numfunc, xcname))
     &   call errquit('dft_rdinput: rtdb_put failed', 304, RTDB_ERR)
c
c     load geometry and symmetry info
c
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('rdinput: geom_create?',70, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('rdinput: no geometry ',72, RTDB_ERR)
      if (.not.rtdb_get(rtdb, 'dft:skeleton',MT_LOG, 1, oskel))then
         oskel = sym_number_ops(geom) .gt. 0
      endif
      if (.not.rtdb_get(rtdb, 'dft:adapt',MT_LOG, 1, oadapt))then
         oadapt = sym_number_ops(geom) .gt. 0
      endif
c
c     get number of atomic centers from geometry object
c
      if (.not. geom_ncent(geom, ncenters))
     &     call errquit('rdinput: geom_ncent failed',73, GEOM_ERR)
c
c     allocate space for atomic coordinates and charges
c
      if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &   icoord))call errquit(
     .     'dft_rdinput: failed to alloc coordinates',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &   icharge))call errquit(
     .     'dft_rdinput: failed to alloc charges',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit(
     .     'dft_rdinput: failed to alloc center tags',0, MA_ERR)
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &     call errquit('rdinput: geom_cart_get failed',74, GEOM_ERR)
c
c     Extract high level info from the data-base setting defaults
c
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     &     title = ' '
c
c     Define the contracted Gaussian functions for the AO basis set.
c
      if (.not. bas_create(AO_bas_han, 'ao basis'))
     &   call errquit('rdinput: bas_create?', 86, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, AO_bas_han, 'ao basis'))
     &   call errquit('rdinput: no ao basis ', 86, BASIS_ERR)
      if (.not. int_normalize(rtdb,AO_bas_han))
     &   call errquit('rdinput: int_normalize failed', 86, INT_ERR)
      if (.not. bas_numbf(AO_bas_han, nbf))
     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
      nbf_ao = nbf
      if (.not. bas_nprim_cn_max(AO_bas_han,nbf_ao_mxprim))
     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
      if (.not. bas_high_angular(AO_bas_han,nbf_ao_mxang))
     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
      if (.not. bas_ncontr_cn_max(AO_bas_han,nbf_ao_mxcont))
     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_cn_max(AO_bas_han,nbf_ao_mxnbf_cn))
     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
      if (.not. bas_numcont(AO_bas_han,nshells_ao))
     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c
c     Define the contracted Gaussian functions for the CD basis set.
c
      if (.not. bas_create(CD_bas_han, 'cd basis'))
     &   call errquit('rdinput: bas_create?',102, BASIS_ERR)
      CDFIT=bas_rtdb_load(rtdb, geom, CD_bas_han, 'cd basis')
      if(theory.eq.'sodft') then
         if (me.eq.0) then
           write(luout,*)
     &        ' WARNING: CD fitting not compatible with spinorbit'
           write(luout,*)
     &        ' WARNING: disabling CD fitting'
         end if
         CDFIT=.false.
      endif
      if (CDFIT)then
         if (.not. int_norm_2c(rtdb,CD_bas_han))
     &      call errquit('rdinput: int_norm_2c failed', 208, BASIS_ERR)
         if (.not. bas_numbf(CD_bas_han, nbf_cd))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nprim_cn_max(CD_bas_han,nbf_cd_mxprim))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_high_angular(CD_bas_han,nbf_cd_mxang))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_ncontr_cn_max(CD_bas_han,nbf_cd_mxcont))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nbf_cn_max(CD_bas_han,nbf_cd_mxnbf_cn))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nbf_ce_max(CD_bas_han,nbf_cd_mxnbf_ce))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_numcont(CD_bas_han,nshells_cd))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. rtdb_put(rtdb, 'dft:cdfit', mt_log, 1, .true.))
     &      call errquit('rdinput: failed put_rtdb:', 99, RTDB_ERR)
c
c     texas breaks when AO_spherical and CD_cartesian
c
         if(bas_is_spherical(AO_bas_han).and.
     N        (.not.bas_is_spherical(CD_bas_han))) 
     T        call int_app_set_no_texas(rtdb)

      else
         if (.not. bas_destroy(CD_bas_han))
     &      call errquit('rdinput:not able to destroy CD_bas:',86,
     &       BASIS_ERR)
         CD_bas_han=-99999
         nbf_cd = 0
         nbf_cd_mxprim = 0
         nbf_cd_mxang = 0
         nbf_cd_mxcont = 0
         nbf_cd_mxnbf_cn = 0
         nbf_cd_mxnbf_ce = 0
         nshells_cd = 0
      endif
c
c     Define the contracted Gaussian functions for the XC basis set.
c
      if (.not. bas_create(XC_bas_han, 'xc basis'))
     &   call errquit('rdinput: bas_create?',124, BASIS_ERR)
      XCFIT=bas_rtdb_load(rtdb, geom, XC_bas_han, 'xc basis')
      if (XCFIT)then
         if (.not. int_normalize(rtdb,XC_bas_han))
     &      call errquit('rdinput: int_normalize failed', 208, INT_ERR)
         if (.not. bas_numbf(XC_bas_han, nbf_xc))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nprim_cn_max(XC_bas_han,nbf_xc_mxprim))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_high_angular(XC_bas_han,nbf_xc_mxang))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_ncontr_cn_max(XC_bas_han,nbf_xc_mxcont))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nbf_cn_max(XC_bas_han,nbf_xc_mxnbf_cn))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nbf_ce_max(XC_bas_han,nbf_xc_mxnbf_ce))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_numcont(XC_bas_han,nshells_xc))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
      else
         if (.not. bas_destroy(XC_bas_han))
     &      call errquit('rdinput:not able to destroy XC_bas:',88,
     &       BASIS_ERR)
         XC_bas_han=-99999
         nbf_xc = 0
         nbf_xc_mxprim = 0
         nbf_xc_mxang = 0
         nbf_xc_mxcont = 0
         nbf_xc_mxnbf_cn = 0
         nbf_xc_mxnbf_ce = 0
         nshells_xc = 0
      endif
c
c     Set the switch to rotate the quadrature points.
c     This switch controls whether the grid is
c     rotationally invariant. For now, it is hard-wired to 
c     1 (rotationally invariant).
c     
      call ga_sync()
      if (.not. rtdb_get(rtdb, 'dft:iqrot', mt_int, 1, iqrot))
     &   then
        iqrot = 1
        if (.not. rtdb_put(rtdb, 'dft:iqrot', mt_int, 1, iqrot))
     &       call errquit('dft_rdinput: rtdb_put failed', 117, RTDB_ERR)
      endif
c
c     generate number of atom types and atom type array iatype(icenter) 
c
      ntypes = 0
      do icenter = 1, ncenters
c
c        is this a new type of atom?
c
         isbq = geom_isbq(geom,icenter)
         lnewtype = .true.
         do jcenter = 1, icenter - 1
            same_atom = Dbl_MB(icharge + icenter - 1) .eq.
     &         Dbl_MB(icharge + jcenter - 1)
            same_bq = geom_isbq(geom,jcenter) .and. isbq
            same_atom = same_atom .or. same_bq
            if (same_atom .and.
     &          atom_tag_check(Byte_MB(itags + (icenter - 1)*16), 
     &                         Byte_MB(itags + (jcenter - 1)*16))
     &         )then   ! same atom type
               lnewtype = .false.
               iatype(icenter) = iatype(jcenter)
               goto 100
            endif
         enddo 
  100    continue
         if (lnewtype)then
            ntypes = ntypes + 1
            iatype(icenter) = ntypes
         endif 
      enddo 
      if (ntypes.gt.dft_ntags_bsmx)then
        write(LuOut,*) 'RDINPUT:  Too many types of atoms.'
        call errquit(' rdinput: raise dft_ntags_bsmx',2, INPUT_ERR)
      end if
c
c     definition of grid
c
c
c     Leb on by default
c
       if (.not. rtdb_get(rtdb, 'dft:lebedev', mt_log, 1, leb)) then
          leb=.true.
          if (.not.rtdb_put(rtdb, 'dft:lebedev', mt_log, 1, leb))
     .         call errquit(' dft_rdinput: put of leb failed',0,
     &       RTDB_ERR)
       endif
      if (.not. rtdb_cget(rtdb,'dft:radgrid', 1,radgrid)) then
         if(leb) then
             radgrid='murakn'
          else
            radgrid='eumac'
          endif
          if (.not. rtdb_cput(rtdb,'dft:radgrid', 1,radgrid)) 
     .         call errquit('dftrdinp: put radgrid failed',7,
     &       RTDB_ERR)
       endif
      if (.not.rtdb_get(rtdb, 'dft:trunc_neigh', mt_log, 1,
     &   trunc_neigh))trunc_neigh=.true.
c
c     partitioning function
c     erf1 on by default
c
      if (.not.rtdb_get(rtdb, 'dft:ssw', mt_log, 1, lssw)) then
         lssw=.true.
         if (.not.rtdb_put(rtdb, 'dft:ssw', mt_log, 1, lssw)) 
     .        call errquit('dftrdi: rtdb_put failed', 1710, RTDB_ERR)
         whichssw='erf1'
         if (.not.rtdb_cput(rtdb, 'dft:ssweights',  1, whichssw))
     &        call errquit('dftrdi: rtdb_cput failed', 1720, RTDB_ERR)
      else
         if(lssw) then
         if (.not.rtdb_cget(rtdb, 'dft:ssweights',  1, whichssw))
     &        call errquit('dftrdi: rtdb_cget failed', 1700, RTDB_ERR)
         endif
      endif
c
c     set up type-indexed znuc array
c
      if (.not. Ma_Push_Get(MT_Dbl,ntypes,'nuc charge',lznuc,iznuc))
     &   call errquit('dft_rdinput: failed to alloc nuc charge',0,
     &       MA_ERR)
c
c     for each type, place the nuclear charge in znuc array and put default 
c     radial cutoff in rad_cutoff
c
      do itype = 1, ntypes
cedo         rad_cutoff(1,itype) = irqmax
         rad_cutoff(2,itype) = -99
         do icenter = 1, ncenters
            if (iatype(icenter) .eq. itype)then
c
c              center icenter is of type itype; assign charge
c
               Dbl_MB(iznuc + itype - 1) = 
     &         dble(Dbl_MB(icharge + icenter - 1))
               goto 110 ! next type
            endif 
         enddo 
  110    continue 
      enddo 
c
c     Figure input/output MO vectors
c     
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $     movecs_in = 'atomic'
      if (.not. rtdb_cget(rtdb, 'dft:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ')then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
     $        movecs_in.eq.'project' .or. movecs_in.eq.'fragment'.or.
     $        movecs_in.eq.'rotate')then
            call util_file_name('movecs', .false.,.false.,movecs_out)
         else
            movecs_out = movecs_in
         endif
      endif
c
c     For spin-orbit calculations, turn off oskel 
c
      if(theory .eq. 'sodft') then
         oadapt=.false.
         oskel=.false.
      endif
c
c     Resolve names of MO files to full paths defaulting to the
c     permanent directory
c
      if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
     $        movecs_in.eq.'project' .or. movecs_in.eq.'fragment'.or.
     $        movecs_in.eq.'rotate')then
         continue
      else
         call util_file_name_resolve(movecs_in, .false.)
      endif
      call util_file_name_resolve(movecs_out, .false.)
c
      if (.not. MA_Pop_Stack(lznuc))
     &   call errquit('dft_rdinput: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(ltags))
     &   call errquit('dft_rdinput: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcharge))
     &   call errquit('dft_rdinput: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcoord))
     &   call errquit('dft_rdinput: pop stack failed.',0, MA_ERR)
c
      call ga_sync
c
      return
      end
      integer function dft_nbases(bases)
      implicit none
#include "cdft.fh"
      integer bases(3)
c
      dft_nbases=1
      if(CDFIT) dft_nbases=dft_nbases+1
      if(XCFIT) dft_nbases=dft_nbases+1
c
      if(XCFIT.and.CDFIT) then
        bases(1) = AO_bas_han
        bases(2) = CD_bas_han
        bases(3) = XC_bas_han
      elseif((.not.XCFIT).and.CDFIT) then
        bases(1) = AO_bas_han
        bases(2) = CD_bas_han
      elseif((.not.CDFIT).and.XCFIT) then
        bases(1) = AO_bas_han
        bases(2) = XC_bas_han
      else
        bases(1) = AO_bas_han
      endif
      return
      end
