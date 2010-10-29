      logical function argos_prep_wrtxyz(lfnout,lfnxyz,
     + filxyz,lfnqqq,filqqq,
     + num,mat,isat,xs,qsa,msa,nsa,box)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      character*2 argos_prep_atnam
      external argos_prep_atnam
c
      integer lfnout,lfnxyz,lfnqqq,msa,nsa,mat
      character*255 filxyz,filqqq
      integer num(mat)
      integer isat(msa)
      real*8 xs(3,msa),qsa(msa),box(3)
c
      integer length
c
      integer i,k
c
c     open XYZ file
c
      length=index(filxyz,' ')-1
      open(unit=lfnxyz,file=filxyz(1:length),form='formatted',
     + status='unknown',err=9999)
c
      write(lfnxyz,1000) nsa
 1000 format('* Coordinates in XYZ format',/,i5)
      do 1 i=1,nsa
      write(lfnxyz,1001)
     +  argos_prep_atnam(num(isat(i))),(1.0d1*xs(k,i),k=1,3)
 1001 format(a2,3f12.6)
    1 continue
      close(unit=lfnxyz)
c
      if(util_print('files',print_default)) then
      write(lfnout,1002) filxyz(1:length)
 1002 format(' Created xyz',t40,a,/)
      endif
c
      filqqq=filxyz
      i=index(filxyz,'.xyz')
      if(i.gt.1) then
      filqqq(i:i+5)='.xyzq '
      else
      filqqq=filxyz(1:length)//'_q'
      endif
c
      length=index(filqqq,' ')-1
      open(unit=lfnqqq,file=filqqq(1:length),form='formatted',
     + status='unknown',err=9999)
c
      write(lfnqqq,2000) nsa,1,box
 2000 format(i7,i3,3f12.6)
      do 2 i=1,nsa
      write(lfnqqq,2001)
     +  argos_prep_atnam(num(isat(i))),(1.0d1*xs(k,i),k=1,3),
     + qsa(i)
 2001 format(a2,4f12.6)
    2 continue
      close(unit=lfnqqq)
c
      if(util_print('files',print_default)) then
      write(lfnout,2002) filqqq(1:length)
 2002 format(' Created q',t40,a,/)
      endif
c
      argos_prep_wrtxyz=.true.
      return
c
 9999 continue
      argos_prep_wrtxyz=.false.
      return
      end
