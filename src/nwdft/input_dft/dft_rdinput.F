      Subroutine dft_rdinput( rtdb, iconst, nconst )

C$Id: dft_rdinput.F,v 1.12 1996-06-17 16:15:55 og845 Exp $
      implicit none
      
      integer rtdb                               ! [input] database handle
      integer nconst
      integer  iconst(nconst)
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "context.fh"
#include "msgids.fh"
#include "sym.fh"
#include "tcgmsg.fh"

#include "mafdecls.fh"
#include "quad.fh"
#include "stdio.fh"
c
      integer bases(3)
      character*255 title
      integer ncenters,iproc,ichg,iscfcon,iguess,mult,ipol,
     ,     igcon,icvgtyp,nfock,ncysh,idirect,idisoff,noio,
     ,     ndamp,idecomp,nagrid,nrinc,nrmax,
     ,     itype,numgrp,numset,itrscf,
     ,     icenter,jcenter,iEx,iEc,i3
      integer nbf,nshells_cd,nbf_xc,nshells_xc,nshells
      integer lcoord,icoord,lcharge,icharge,ltags,itags,
     ,     lznuc,iznuc,lratom,iratom
      double precision rcharge,rm,rscale,rlshift
      Logical LResult, lnewtype, XCFIT, CDFIT, oskel
      Logical int_normalize
      External int_normalize
      Logical int_norm_2c
      External int_norm_2c
c

c******************************************************************************

      iproc  = ga_nodeid()
c****
c**** recover logical parameters from rtdb
c****
C      if (.not.rtdb_get(rtdb, 'dft:skeleton',MT_LOG, 1, oskel)) then
C        write (6,*) "dft_rdinput: did not find dft:skeleton"
C        oskel = sym_number_ops(geom) .gt. 0
C        if (.not.rtdb_put(rtdb, 'dft:skeleton',MT_LOG, 1, oskel))
C     &       call errquit('rdinput: rtdb_put failed',88) 
C      else 
C        write (6,*) "dft_rdinput: found dft:skeleton"
C        write (6,*) "oskel: ",oskel
C      endif 
c****
c**** recover charge from rtdb
c****
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, rcharge)) then
        rcharge = 0.0d0
      endif 
      ichg = idint(1.000001d0 * rcharge) ! avoid bad truncation
c     
c     recover rdcalc integer parameters from rtdb
c     
      if (.not. rtdb_get(rtdb, 'dft:itrscf', mt_int, 1, itrscf))
     &     call errquit('rdinput: rtdb_get failed',2) 
      if (.not. rtdb_get(rtdb, 'dft:iscfcon', mt_int, 1, iscfcon))
     &     call errquit('rdinput: rtdb_get failed',4) 
      if (.not. rtdb_get(rtdb, 'dft:iguess', mt_int, 1, iguess))
     &     call errquit('rdinput: rtdb_get failed',6) 
      if (.not. rtdb_get(rtdb, 'dft:mult', mt_int, 1, mult))
     &     call errquit('rdinput: rtdb_get failed',10)
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &     call errquit('rdinput: rtdb_get failed',14)
      if (.not. rtdb_get(rtdb, 'dft:igcon', mt_int, 1, igcon))
     &     call errquit('rdinput: rtdb_get failed',16)
      if (.not. rtdb_get(rtdb, 'dft:icvgtyp', mt_int, 1, icvgtyp))
     &     call errquit('rdinput: rtdb_get failed',18)
      if (.not. rtdb_get(rtdb, 'dft:nfock', mt_int, 1, nfock))
     &     call errquit('rdinput: rtdb_get failed',20)
      if (.not. rtdb_get(rtdb, 'dft:ncysh', mt_int, 1, ncysh))
     &     call errquit('rdinput: rtdb_get failed',20)
      if (.not. rtdb_get(rtdb, 'dft:idirect', mt_int, 1, idirect))
     &     call errquit('rdinput: rtdb_get failed',22)
      if (.not. rtdb_get(rtdb, 'dft:idisoff', mt_int, 1, idisoff))
     &     call errquit('rdinput: rtdb_get failed',24)
      if (.not. rtdb_get(rtdb, 'dft:noio', mt_int, 1, noio))
     &     call errquit('rdinput: rtdb_get failed',32)
      if (.not. rtdb_get(rtdb, 'dft:lshift', mt_dbl, 1,rlshift))
     &     call errquit('rdinput: rtdb_get failed',34)
      if (.not. rtdb_get(rtdb, 'dft:ndamp', mt_int, 1, ndamp))
     &     call errquit('rdinput: rtdb_get failed',36)
      if (.not. rtdb_get(rtdb, 'dft:idecomp', mt_int, 1, idecomp))
     &     call errquit('rdinput: rtdb_get failed',38)
      if (.not. rtdb_get(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &     call errquit('rdinput: rtdb_get failed',40)
      if (.not. rtdb_get(rtdb, 'dft:iqrot', mt_int, 1, iqrot))
     &     call errquit('rdinput: rtdb_get failed',42)
      if (.not. rtdb_get(rtdb, 'dft:nrinc', mt_int, 1, nrinc))
     &     call errquit('rdinput: rtdb_get failed',44)
      if (.not. rtdb_get(rtdb, 'dft:nrmax', mt_int, 1, nrmax))
     &     call errquit('rdinput: rtdb_get failed',46)
c     
c     recover real parameters from rtdb
c     
      if (.not. rtdb_get(rtdb, 'dft:rm', mt_dbl, 1, rm))
     &     rm=2.d0
      if (.not. rtdb_get(rtdb, 'dft:rscale', mt_dbl, 1, rscale))
     &     call errquit('rdinput: rtdb_get failed',56)
c     
c     Read property evaluation options.

        iconst(39) = idecomp  ! only action in rdprop other than read in
        iconst(40) = idisoff

c     
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('rdinput: geom_create?',70)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('rdinput: no geometry ',72)

      itype = 1
      numgrp = 1
      iconst(44) = itype      
      iconst(45) = numgrp      
      iconst(46) = numset     
c****
c**** get number of atomic centers from geometry object
c****
      if (.not. geom_ncent(geom, ncenters))
     &     call errquit('rdinput: geom_ncent failed',73)      
c****
c**** allocate space for atomic coordinates and charges
c****
      LResult = Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &     icoord)
      LResult = Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &     icharge)
      LResult = Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &                       ltags, itags)
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &     call errquit('rdinput: geom_cart_get failed',74)
c
c     
c
      if (.not. rtdb_get(rtdb, 'dft:nrad', mt_int, ntags_bsmx, nrad))
     &     call errquit('rdinput: rtdb_get failed',76)
c
c     Extract high level info from the data-base setting defaults
c
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     &     title = ' '
c
c     end lines from ddscf/scf
c
c     
      if (.not. bas_create(AO_bas_han, 'ao basis'))
     &     call errquit('rdinput: bas_create?',80)

      if (.not. bas_rtdb_load(rtdb, geom, AO_bas_han, 'ao basis'))
     &     call errquit('rdinput: no ao basis ',82)

      if (.not. int_normalize(AO_bas_han))
     &     call errquit('rdinput: int_normalize failed', 200)


      if (.not. bas_numbf(AO_bas_han, nbf))
     &     call errquit('rdinput:rdinput:',86)
c**** 


c
c       
c     Broadcast AO basis set to all nodes.
c****
c     Define the contracted Gaussian functions for the CD basis set.


      if (.not. bas_create(CD_bas_han, 'cd basis'))
     &     call errquit('rdinput: bas_create?',102)
      CDFIT=bas_rtdb_load(rtdb, geom, CD_bas_han, 'cd basis')
      if(CDFIT) then
        if (.not. int_norm_2c(CD_bas_han))
     &       call errquit('rdinput: int_norm_2c failed', 208)

        
      if (.not. bas_numbf(CD_bas_han, nbf_cd))
     &     call errquit('rdinput:rdinput:',86)
          if (.not. bas_numcont(CD_bas_han, nshells_cd))
     &         call errquit('rdinput:rdinput:',87)
      else
        if (.not. bas_destroy(CD_bas_han))
     &       call errquit('rdinput:not able to destroy CD_bas:',86)
        CD_bas_han=-99999
        nbf_cd=0
      endif

c       
        
        if (.not. bas_create(XC_bas_han, 'xc basis'))
     &       call errquit('rdinput: bas_create?',124)
        XCFIT=bas_rtdb_load(rtdb, geom, XC_bas_han, 'xc basis')
        if(XCFIT) then
          if (.not. int_normalize(XC_bas_han))
     &         call errquit('rdinput: int_normalize failed', 208)


      if (.not. bas_numbf(XC_bas_han, nbf_xc))
     &     call errquit('rdinput:rdinput:',86)

          if (.not. bas_numcont(XC_bas_han, nshells_xc))
     &         call errquit('rdinput:rdinput:',88)
        else
          if (.not. bas_destroy(XC_bas_han))
     &         call errquit('rdinput:not able to destroy XC_bas:',88)
          XC_bas_han=-99999
          nbf_xc=0
        endif
c****   
c****   Read the numerical quadrature input parameters.
c****   
        if( iproc.eq.0 )then

          rscale = 0.0d0

          if( nrmax.gt.MXrad )then
            write(LuOut,*) 
     *           'RDINPUT:  Too many radial quadrature points.'
            call errquit(' rdinput: raise MXrad',1)
          end if

          if( nagrid.gt.MXagrid )then
            write(LuOut,*) 
     *           'RDINPUT:  Too many angular quadrature grids.'
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
     &           Dbl_MB(icharge+jcenter-1) ) then   ! same atom type
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
          if( ntypes.gt.ntags_bsmx )then
            write(LuOut,*) 'RDINPUT:  Too many types of atoms.'
            call errquit(' rdinput: raise ntags_bsmx',2)
          end if
c****
c**** set up type-indexed znuc array
c****
      LResult = Ma_Push_Get(MT_Dbl,ntypes,'nuc charge',lznuc,iznuc)
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
     &     call errquit('rdinput: failed to alloc ratom',101)


      nqshells = MXqshells

      call xc_setquad(rtdb,Dbl_MB(iznuc),Dbl_MB(iratom),dbl_mb(icoord),
     &     ncenters,
     &     rm,nrinc,nrmax,
     &     nagrid)


      LResult = MA_Pop_Stack(lratom)
      LResult = MA_Pop_Stack(lznuc)
      LResult = MA_Pop_Stack(ltags)
      LResult = MA_Pop_Stack(lcharge)
      LResult = MA_Pop_Stack(lcoord)


c
c     Write constant array to checkpoint file.
c     NEED TO LOSE THIS  (HT)

c     (These constants should be accessed through an integral API routine
c     when they are needed. I define these elements of the constant array
c     only to keep "old" routines working. Use of this constant array
c     should be eliminated at least for these Gaussian basis set parameters.)

      iconst( 8) = itrscf 
      iconst( 9) = iscfcon
      iconst(11) = ipol   
      iconst(18) = igcon  
      iconst(36) = icvgtyp
      iconst(37) = nfock
      iconst(38) = idirect
      iconst(49) = noio
C      iconst(50) = lshift
      iconst(51) = ndamp
      iconst(52) = ncysh

      iconst(12)  =  ichg  
      iconst(13)  =  mult  

      if (.not. bas_numcont(AO_bas_han, nshells))
     &     call errquit('rdinput:rdinput:',86)


c     Check the user supplied input for any inconsistencies.

        call dft_inpana(geom,rtdb,iconst,nconst,ncenters,
     &       itrscf,iscfcon,igcon,ipol,icvgtyp,idirect,
     &       nbf_xc,nbf_cd,iEx,iEc) 
c
c       inpana could have modified ipol
c
        ipol=iconst(11)

      if( iproc.eq.0 )then
        if(nfock.ne.0) then
          write(LuOut,'(3x,A,i2,A)') 
     &         '       ',nfock,' Fock matrices are used in DIIS'
          write(LuOut,*)
        endif
        if(abs(rlshift).gt.1.d-5) then
          write(LuOut,'(3x,A,f4.2,A,i4,A)') 
     &         '       Fock Matrix Level Shifted by ',rlshift,
     &         ' hartree for ',ncysh,' cycles'
          write(LuOut,*)
        endif
        if(ndamp.ne.0) then
          write(LuOut,'(3x,A,i2,A)') 
     &         '       Density Matrix Damping by ',ndamp,'%'
          write(LuOut,*)
        endif
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



      call ga_brdcst(Msg_iconst,iconst,mitob(nconst),0)


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

      call int_init(rtdb, i3, bases)


      return
      end

