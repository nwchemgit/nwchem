c***********************************************************************
c
c     inpana (input analysis)
c     Analyze input to deduce nature of system, and set key flags.
c
c***********************************************************************
      subroutine dft_inpana(geom,rtdb,iconst,nconst,nat)
C$Id: dft_inpana.F,v 1.7 1996-01-09 17:54:15 og845 Exp $
      implicit none
      integer nat,nconst
      logical even
      integer iconst(nconst)
      integer noc(2)
      Integer rtdb    !  runtime database handle
      integer geom
C
      integer me,ipol,ichg,mult,nel
      integer noc1
      double precision anucl_charg
      double precision EPSILON
      Parameter ( EPSILON=1.d-13 )
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
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
         noc(2) = (nel - mult + 1)/2
         noc(1) = nel - noc(2)
         noc1  = noc(1) + noc(2)
      endif
c
c-->  Check to see if there are enough electrons for this
c-->  value of the multiplicity.
c
      if(noc(2).lt.0)then
         if(me.eq.0) write(6,9100)mult
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
         write(LuOut,9150)nat
c
         write(LuOut,9200)nel,ichg,mult

      endif
c
9000  format(/,5x,'There are ',i4,' electrons in this system.',
     &       /,5x,'RHF is not an appropriate SCF type.',
     &       /,5x,'Stopping program in subroutine inpana.')
9050  format(/,5x,'RHF not allowed for open-shell systems.',
     &       /,5x,'Stopping program in subroutine inpana.')
9100  format(/,5x,'Not enough electrons for multiplicity ',i2,'.',
     &       /,5x,'Stopping program in subroutine inpana.')
9150  format(/,10x,'No. of atoms      :',15x,i4)
9200  format(/,10x,'No. of electrons  :',15x,i4,
     &      //,10x,'Charge            :',15x,i4,
     &      //,10x,'Spin multiplicity :',15x,i4)
9250  format(/,10x,'Atomic center number ',i4,' is a ghost atom.')
      end
