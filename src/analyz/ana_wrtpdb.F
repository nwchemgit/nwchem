      subroutine ana_wrtpdb(lfnpdb,filpdb,sgmnam,isel,xs)
c
c $Id: ana_wrtpdb.F,v 1.2 2000-05-04 06:26:37 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
c
      character*16 sgmnam(msa)
      integer isel(msa)
      integer lfnpdb
      real*8 xs(msa,3)
      character*255 filpdb
c
      integer i,k,length
c
      if(me.eq.0) then
c
      length=index(filpdb,' ')-1
      open(unit=lfnpdb,file=filpdb(1:length),form='formatted',
     + status='unknown')
c
      rewind(lfnpdb)
c
      write(lfnpdb,1000)
 1000 format(
     + 'HEADER',/,
     + 'TITLE     ',/,
     + 'TITLE    2',/,
     + 'TITLE    3',/,
     + 'REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.1, 25-OCT-1996')
c
      do 1 i=1,nsa
      if(isel(i).gt.0) then
      write(lfnpdb,1001) i,sgmnam(i)(6:9),sgmnam(i)(1:3),
     + sgmnam(i)(13:16),(1.0d1*xs(i,k),k=1,3)
 1001 format('ATOM',i7,1x,a4,1x,a3,2x,a4,4x,3f8.3)
      endif
    1 continue
c
      write(lfnpdb,1002)
 1002 format('END')
c
      close(unit=lfnpdb)
c
      endif
c
      return
      end
