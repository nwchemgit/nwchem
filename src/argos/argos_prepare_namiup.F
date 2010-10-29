      logical function argos_prepare_namiup(lfnout,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm)
c
c $Id$
c
c     function to comply with IUPAC atom naming
c
c     in  : lfnout    = output file logical file number
c           kse
c           mseq      = dimension of the sequence list
c           nseq      = length of the sequence list
c
c     out : lseq(1,*) = segment numbers
c           lseq(2,*) = number of atoms in segment
c           lseq(3,*) = index to list of unique segments
c           cseq      = segment names
c           nseq      = actual length of the sequence list
c
      implicit none
c
#include "util.fh"
#include "argos_prepare_common.fh"
c
      integer argos_prepare_atnum
      external argos_prepare_atnum
      logical argos_prepare_protein
      external argos_prepare_protein
c
      integer lfnout
      integer mseq,matm
      integer nseq,natm
      integer lseq(6,mseq),latm(5,matm)
      character*6 catm(3,matm)
      character*5 cseq(2,mseq)
      character*2 elemnt
c
      integer i
      logical flag1
c
      if(util_print('sequence',print_debug)) then
      write(lfnout,2000)
 2000 format(/,'argos_prepare_namiup ')
      endif
c
      flag1=.false.
      do 1 i=1,natm
      if(argos_prepare_protein(ffield,cseq(1,latm(1,i)))) then
      if(catm(1,i)(1:2).eq.'1H') flag1=.true.
      endif
      catm(2,i)=catm(1,i)
    1 continue
c
c     IUPAC hydrogen naming
c
      if(flag1) then
      do 2 i=1,natm
      if(argos_prepare_protein(ffield,cseq(1,latm(1,i)))) then
      if(catm(2,i)(1:2).eq.'3H') catm(2,i)(1:2)='4H'
      if(catm(2,i)(1:2).eq.'2H') catm(2,i)(1:2)='3H'
      if(catm(2,i)(1:2).eq.'1H') catm(2,i)(1:2)='2H'
      endif
    2 continue
      if(util_print('sequence',print_high)) then
      write(lfnout,2001)
 2001 format(/,'Hydrogen atom prefix corrected ')
      endif
      endif
c
c     ILE CD1 -> CD
c
      do 3 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'ILE') then
      if(catm(2,i)(1:4).eq.' CD1') then
      catm(2,i)(1:4)=' CD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) latm(1,i)
 2002 format(/,'Renamed ILE',i3.3,' CD1 to CD ')
      endif
      endif
      if(catm(2,i)(2:4).eq.'HD1') then
      catm(2,i)(2:4)='HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2003) latm(1,i)
 2003 format(/,'Renamed ILE',i3.3,' HD1 to HD ')
      endif
      endif
      endif
    3 continue
c
      do 4 i=1,natm
      if(latm(2,i).le.0) latm(2,i)=argos_prepare_atnum(catm(2,i)(1:2))
      if(latm(2,i).eq.0) then
      elemnt(1:1)=' '
      elemnt(2:2)=catm(2,i)(2:2)
      latm(2,i)=argos_prepare_atnum(elemnt)
      if(util_print('sequence',print_high)) then
      write(lfnout,2004) latm(1,i)
 2004 format(/,'Unrecognized element in segment ',i4)
      endif
      endif
    4 continue
c
      argos_prepare_namiup=.true.
c
      return
      end
