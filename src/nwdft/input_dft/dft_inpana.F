c***********************************************************************
c
c     inpana (input analysis)
c     Analyze input to deduce nature of system, and set key flags.
c
c***********************************************************************
      subroutine dft_inpana(geom,rtdb,iconst,nconst,nat,
     ,     itrscf,iscfcon,igcon,ipol,
     &     idirect,nbf_xc,nbf_cd,iEx,iEc)
C$Id: dft_inpana.F,v 1.9 1996-10-01 21:22:18 d3h449 Exp $
      implicit none
      integer nat,nconst
      logical even
      integer iconst(nconst)
      integer noc(2)
      Integer rtdb    !  runtime database handle
      integer geom
      character*50 xname,cname
      integer ipol,itrscf,idirect,nbf_xc,nbf_cd,
     ,     iscfcon,igcon,iEx,iEc
c
c
C
      integer me,ichg,mult,nel
      integer noc1
      double precision anucl_charg
      character*10 gridtype
      double precision EPSILON
      Parameter ( EPSILON=1.d-13 )
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "geom.fh"
c
      me=ga_nodeid()

      if(me.eq.0) call dft_header('Input analysis')
c
c-->  Read parameters from checkpoint file.
c**** these should be changed to rtdb gets
c
      ipol    = iconst(11)
      ichg    = iconst(12)
      mult    = iconst(13)

c
c-->  Determine number of electrons in system.
c
      if (.not. geom_nuc_charge(geom, anucl_charg))
     $     call errquit('scf: geom_nuc_charge failed', 0)
c
      nel = int(anucl_charg+EPSILON) - ichg
c
c-->  Check to see if calculation type is allowed.
c
c
c-->  Even number of electrons required for RHF.
c     
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
C     check consistency of no. elec and multiplicity
C
        even=mod((nel+mult-1),2).eq.0
        if(.not.even) call errquit(
     (       ' no. of electrons and multiplicity not compatible',nel)
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

      if (.not. rtdb_put(rtdb, 'dft:ipol', MT_INT, 1, ipol ))
     $     call errquit('inpana: dft:ipol put failed', 0)
c
c-->  Write new data to checkpoint file.
c
      iconst(11) = ipol
      if(me.eq.0) then
c
c-->  Write to output.
c
         write(LuOut,8150)nat
c
         write(LuOut,8200)nel,ichg,mult

      endif
c
      if (.not. rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype))
     $     call errquit('dft_inpana: rtdb_get failed',0) 
      if (.not. rtdb_cget(rtdb, 'dft:cname', 1, cname))
     $     call errquit('dft_inpana: rtdb_get failed',1) 
      if (.not. rtdb_cget(rtdb, 'dft:xname', 1, xname))
     $     call errquit('dft_inpana: rtdb_get failed', 2) 
      if(ga_nodeid().eq.0) then
      write(LuOut,9000)
c     
      write(LuOut,9020)
c     
      if(ipol.eq.1)then
        write(LuOut,9050)
      elseif(ipol.eq.2)then
        write(LuOut,9055)
      endif
      write(LuOut,9030)itrscf
c     
c      if(icvgtyp.eq.0) then
c        write(LuOut,9080)
c      elseif(icvgtyp.eq.1) then
c        write(LuOut,9090)
c      elseif(icvgtyp.eq.2) then
c        write(LuOut,9100)
c      endif
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
      write(LuOut,9223) xname
      if(cname.ne.' ') write(LuOut,9224) cname
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
 9080 format(10x,'SCF convergence type: simple SCF.',/)
 9090 format(10x,'SCF convergence type: density matrix mixing.',/)
 9100 format(10x,'SCF convergence type: DIIS.',/)
 9110 format(10x,'This is a Direct SCF calculation.',/)
 9120 format(10x,'An Exch-Corr fitting basis will be used.',/)
 9130 format(10x,'A Charge density fitting basis will be used.',/)
 9135 format(10x,'Grid used for XC integration:  ',a,/)
 9210 format(10x,'Exact (Hartree-Fock)                 exchange',/)
 9223 format(10x,A,/)
 9224 format(10x,A,/)
8050  format(/,5x,'RHF not allowed for open-shell systems.',
     &       /,5x,'Stopping program in subroutine inpana.')
8100  format(/,5x,'Not enough electrons for multiplicity ',i2,'.',
     &       /,5x,'Stopping program in subroutine inpana.')
8150  format(/,10x,'No. of atoms      :',15x,i4)
8200  format(/,10x,'No. of electrons  :',15x,i4,
     &      //,10x,'Charge            :',15x,i4,
     &      //,10x,'Spin multiplicity :',15x,i4)
      end



