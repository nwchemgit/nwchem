c***********************************************************************
c
c     inpana (input analysis)
c     Analyze input to deduce nature of system, and set key flags.
c
c***********************************************************************
      subroutine dft_inpana(rtdb)
C$Id: dft_inpana.F,v 1.14 1996-10-23 01:18:40 gutowski Exp $
      implicit none
      logical even, cksetd
      integer noc(2)
      Integer rtdb    !  runtime database handle
      integer me,ichg,nel
      integer noc1
      integer n
      double precision anucl_charg, ckfac
      double precision EPSILON
      character*9 local_c, nonlocal_c, nineb_c
      Parameter ( EPSILON=1.d-13 )
c
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "geom.fh"
#include "cdft.fh"
c
      me=ga_nodeid()
c
      if(me.eq.0) call dft_header('Input analysis')
c
c-->  Determine number of electrons in system.
c
      if (.not. geom_nuc_charge(geom, anucl_charg))
     $     call errquit('scf: geom_nuc_charge failed', 0)
c
      ichg = idint(1.000001d0 * rcharge) ! avoid bad truncation
      nel = int(anucl_charg+EPSILON) - ichg
c
c-->  Check to see if calculation type is allowed.
c
c-->  Even number of electrons required for RHF.
c     
      even=mod(nel,2).eq.0
c
c     odd # of electrons or not a singlet state --> LSD
c
      if((.not.even).or.(mult.ne.1)) ipol=2
      
      noc(2)=0
c
c-->  Calculate number of occupied orbitals.
c
      if(ipol.eq.1)then
        noc1 = nel/2
        noc(2)= 0
        noc(1)= noc1
      else
C
C       check consistency of no. elec and multiplicity
C
        even=mod((nel+mult-1),2).eq.0
        if(.not.even) call errquit(
     &    ' no. of electrons and multiplicity not compatible',nel)
        noc(2) = (nel - mult + 1)/2
        noc(1) = nel - noc(2)
        noc1  = noc(1) + noc(2)
      endif
c
c-->  Check to see if there are enough electrons for this
c-->  value of the multiplicity.
c
      if(noc(2).lt.0)then
        if(me.eq.0) write(LuOut,8100)mult
        stop
      endif
c****
c**** write noc (consistent with definition in ddscf) to rtdb
c****

      if (.not. rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     $     call errquit('inpana: rtdb_put of noc failed', 0) 
c
c-->  Write new data to checkpoint file.
c
c     Analyze any user specified XC functionals ... set if none
c
      cksetd = .true.
c
c     Check if user has specified some type of functional,
c     if so, do not set defaults.
c
      do n = 1, 15
        if(lcfac(n).or.nlcfac(n).or.lxfac(n).or.nlxfac(n))
     &    cksetd = .false.
      enddo
      if(cksetd)then
c
c       Set functional defaults.
c
         cfac(1) = 1.0d0
        lcfac(1) = .true.
         xfac(2) = 1.0d0
        lxfac(2) = .true.
c
c       Update rtdb.
c
        if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, 15, cfac))
     $     call errquit('dft_input: rtdb_put failed', 210)
        if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, 15, xfac))
     $     call errquit('dft_input: rtdb_put failed', 211)
        if (.not. rtdb_put(rtdb, 'dft:lcfac', mt_log, 15, lcfac ))
     $     call errquit('dft_input: rtdb_put failed', 9)
        if (.not. rtdb_put(rtdb, 'dft:lxfac', mt_log, 15, lxfac ))
     $     call errquit('dft_input: rtdb_put failed', 11)
      endif
      if(me.eq.0) then
c
c-->    Write to output.
c
        write(LuOut,8150)ncenters
c
        write(LuOut,8200)nel,ichg,mult
c
        write(LuOut,9000)
c     
        write(LuOut,9020)
c     
        if(ipol.eq.1)then
          write(LuOut,9050)
        elseif(ipol.eq.2)then
          write(LuOut,9055)
        endif
c
        write(LuOut,9030)itrscf
c     
        if(idirect.eq.1)write(LuOut,9110)
c     
        if(nbf_xc.gt.0)write(LuOut,9120)
c     
        if(nbf_cd.gt.0)write(LuOut,9130)
c     
        write(LuOut,9135)gridtype
        write(LuOut,9035)iscfcon
        if(igcon.gt.0)then
          write(LuOut,9040)igcon
        endif
c
c       Write out XC info. Combo info first, than X components,
c       than C components.
c
        write(LuOut,9222)
           local_c = 'local    '
        nonlocal_c = 'non-local'
           nineb_c = '         '
        do n = 1, 15
          if(xccomb(n))write(LuOut,9223) xcname(n)
        enddo
c
c       Do exact exchange differently.
c
        if(lxfac(1).or.nlxfac(1))
     &    write(LuOut,9224) xname(1), xfac(1), nineb_c
        do n = 2, 15
          if(lxfac(n).and.nlxfac(n))then
            write(LuOut,9224) xname(n), xfac(n), nineb_c
          elseif(lxfac(n).and.(.not.nlxfac(n)))then
            write(LuOut,9224) xname(n), xfac(n), local_c
          elseif((.not.lxfac(n)).and.nlxfac(n))then
            write(LuOut,9224) xname(n), xfac(n), nonlocal_c
          endif
        enddo
        do n = 1, 15
          if(lcfac(n).and.nlcfac(n))then
            write(LuOut,9224) cname(n), cfac(n), nineb_c
          elseif(lcfac(n).and.(.not.nlcfac(n)))then
            write(LuOut,9224) cname(n), cfac(n), local_c
          elseif((.not.lcfac(n)).and.nlcfac(n))then
            write(LuOut,9224) cname(n), cfac(n), nonlocal_c
          endif
        enddo
c
c       Check XC coefficients to make sure appropriate components
c       sum to 1.0
c        
        ckfac = 0.0d0
        do n = 1, 15
          if(lcfac(n))ckfac = ckfac + cfac(n)
        enddo
        if(abs(ckfac-1.0d0).gt.1.d-8)then
          write(LuOut,*)' WARNING: Sum of local correlation is ',ckfac
          write(LuOut,*)' Sum of components do not equal unity. '
        endif
        ckfac = 0.0d0
        do n = 1, 15
          if(nlcfac(n))ckfac = ckfac + cfac(n)
        enddo
        if(abs(ckfac-1.0d0).gt.1.d-8.and.abs(ckfac).gt.1.d-8)then
          write(LuOut,*)' WARNING: Sum of nonlocal correlation is ',
     &                  ckfac
          write(LuOut,*)' Sum of components do not equal unity or 0. '
        endif
        ckfac = 0.0d0
        do n = 1, 15
          if(lxfac(n))ckfac = ckfac + xfac(n)
        enddo
        if(abs(ckfac-1.0d0).gt.1.d-8)then
          write(LuOut,*)' WARNING: Sum of local exchange is ',ckfac
          write(LuOut,*)' Sum of components do not equal unity. '
        endif
        ckfac = 0.0d0
        do n = 1, 15
          if(nlxfac(n))ckfac = ckfac + xfac(n)
        enddo
        if(abs(ckfac-1.0d0).gt.1.d-8.and.abs(ckfac).gt.1.d-8)then
          write(LuOut,*)' WARNING: Sum of nonlocal exchange is ',ckfac
          write(LuOut,*)' Sum of components do not equal unity or 0. '
        endif
      endif
      return
 9000 format(//,5x,31('-'),/,5x,
     &     'Parameters for SCF calculation.',/,5x,31('-'),/)
 9020 format(10x,'SCF calculation type: DFT',/)
 9030 format(10x,'Maximum number of iterations = ',I3,/)
 9035 format(10x,'Convergence on energy requested = 10**(-',i2,')',/)
 9040 format(10x,'Convergence on density requested = 10**(-',i2,')',/)
 9050 format(10x,'Wavefunction type:  closed shell.',/)
 9055 format(10x,'Wavefunction type:  spin polarized.',/)
 9110 format(10x,'This is a Direct SCF calculation.',/)
 9120 format(10x,'An Exch-Corr fitting basis will be used.',/)
 9130 format(10x,'A Charge density fitting basis will be used.',/)
 9135 format(10x,'Grid used for XC integration:  ',a,/)
 9222 format(10x,'Exchange-Correlation Functional Information',/,
     &       10x,43('-'))
 9223 format(10x,a40)
 9224 format(10x,a40,1x,f6.2,1x,a9)
 8050 format(/,5x,'RHF not allowed for open-shell systems.',
     &       /,5x,'Stopping program in subroutine inpana.')
 8100 format(/,5x,'Not enough electrons for multiplicity ',i2,'.',
     &       /,5x,'Stopping program in subroutine inpana.')
 8150 format(/,10x,'No. of atoms      :',15x,i4)
 8200 format(/,10x,'No. of electrons  :',15x,i4,
     &      //,10x,'Charge            :',15x,i4,
     &      //,10x,'Spin multiplicity :',15x,i4)
      end
