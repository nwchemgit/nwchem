c***********************************************************************
c
c     inpana (input analysis)
c     Analyze input to deduce nature of system, and set key flags.
c
c***********************************************************************
      subroutine dft_inpana(rtdb,zan,c,iconst,nconst,nat)
C$Id: dft_inpana.F,v 1.1 1995-04-27 03:16:08 og845 Exp $
      implicit real*8 (a-h,o-z)
      logical even
      dimension zan(nat), c(3,nat), iconst(nconst)
      integer noc(2)
      Integer rtdb    !  runtime database handle
#include "stdio.fh"
#include <rtdb.fh>     ! for rtdb_get
#include <mafdecls.h>  ! for mt_int
#include "global.fh"
c
      me=ga_nodeid()

      if(me.eq.0) call dft_header('Input analysis')
c
c-->  Read parameters from checkpoint file.
c**** these should be changed to rtdb gets
c
      nbf     = iconst(5)
      ipol    = iconst(11)
      ichg    = iconst(12)
      mult    = iconst(13)
      igcon   = iconst(18)

c      if ( ipol.eq.2 .and. igcon.gt.0) then
c        igcon = -1
c        write (LuOut,*) "INPANA: found density convergence turned on"
c        write (LuOut,*) "This fails for open-shell; removing flag"
c      endif 
c
c-->  Scan for ghost atoms.
c
      nghost=0
      do 500 inat=1,nat
        if (zan(inat).lt.0.5)then
          nghost=nghost+1
          zan(inat)=0.0d0
          if(me.eq.0) write(LuOut,9250) inat
        endif 
  500 continue 
c
c-->  Determine number of electrons in system.
c
      izantot = 0
      do 1000 i=1,nat
         izantot = izantot + int(zan(i))
1000  continue
c
      nel = izantot - ichg
c
c-->  Check to see if calculation type is allowed.
c
c
c-->  Even number of electrons required for RHF.
c     
c
      even=mod(nel,2).eq.0
c
c     odd # of electrons --> LSD
c
      if(.not.even) ipol=2
c
c     not a singlet state --> LSD
c
      if(mult.ne.1) ipol=2

      nocb=0
c
c-->  Calculate number of occupied orbitals.
c
      if(ipol.eq.1)then
         noc1 = nel/2
         nocb= noc1
         noca= noc1
      else
         nocb = (nel - mult + 1)/2
         noca = nel - nocb
         noc1  = noca + nocb
      endif
c
c-->  Check to see if there are enough electrons for this
c-->  value of the multiplicity.
c
      if(nocb.lt.0)then
         if(me.eq.0) write(6,9100)mult
         stop
      endif
c****
c**** write noc (consistent with definition in ddscf) to rtdb
c****
      noc(1) = noca
      if ( ipol .eq. 1) then
        noc(2) = 0
      else 
        noc(2) = nocb
      endif 

      if (.not. rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     $     call errquit('inpana: rtdb_put of noc failed', 0) 

      if (.not. rtdb_put(rtdb, 'dft:ipol', MT_INT, 1, ipol ))
     $     call errquit('inpana: dft:ipol put failed', 0)
c
c-->  At later date insert geometry analysis here using c(3,nat).
c
c      call geom(c,nat)
c
c-->  Write new data to checkpoint file.
c
      iconst(15) = noc1
      iconst(16) = noca
      iconst(17) = nocb
      iconst(11) = ipol
      if(me.eq.0) then
c
c-->  Write to output.
c
         write(LuOut,9150)nat-nghost
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









