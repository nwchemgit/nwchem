      logical function pre_namiup(lfnout,lprint,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm)
c
c     function to comply with IUPAC atom naming
c
c     in  : lfnout    = output file logical file number
c           lprint    = integer for print to output
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
      integer pre_atnum
      external pre_atnum
c
      integer lfnout
      integer mseq,matm
      integer nseq,natm
      integer lseq(4,mseq),latm(5,matm)
      character*6 catm(3,matm)
      character*5 cseq(2,mseq)
      character*255 filnam
      integer lprint
c
      integer isgm,nsgm,i,j,nchain,iatm
      logical flag1
c
      if(lprint.ge.5) then
      write(lfnout,2000)
 2000 format(/,'pre_namiup ')
      endif
c
      flag1=.false.
      do 1 i=1,natm
      if(catm(1,i)(1:2).eq.'1H') flag1=.true.
      catm(2,i)=catm(1,i)
    1 continue
c
c     IUPAC hydrogen naming
c
      if(flag1) then
      do 2 i=1,natm
      if(catm(2,i)(1:2).eq.'3H') catm(2,i)(1:2)='4H'
      if(catm(2,i)(1:2).eq.'2H') catm(2,i)(1:2)='3H'
      if(catm(2,i)(1:2).eq.'1H') catm(2,i)(1:2)='2H'
    2 continue
      endif
c
c     ILE CD1 -> CD
c
      do 3 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'ILE') then
      if(catm(2,i)(1:4).eq.' CD1') catm(2,i)(1:4)=' CD '
      if(catm(2,i)(2:4).eq.'HD1') catm(2,i)(2:4)='HD '
      endif
    3 continue
c
      do 4 i=1,natm
      latm(2,i)=pre_atnum(catm(2,i)(1:2))
    4 continue
c
      pre_namiup=.true.
c
      return
      end
