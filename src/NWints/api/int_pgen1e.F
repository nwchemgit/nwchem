      subroutine int_pgen1e(msg, i_basis, j_basis, ish, jsh, 
     &       ilab, jlab, l1e, Gen1e)
c $Id: int_pgen1e.F,v 1.1 1994-04-08 15:59:46 d3e129 Exp $
#include "basP.fh"
#include "geobasmapP.fh"
      implicit none
c
c::external subroutines used
c dcopy
c errquit
c::function
      logical int_chk_init
      external int_chk_init
      logical int_chk_sh
      external int_chk_sh
c::passed
      character*(*) msg
      integer i_basis, j_basis, ish, jsh, l1e
      integer ilab(l1e),jlab(l1e)
      double precision Gen1e(l1e)
c::local
      integer ii
c
      write(6,*)'pgen1e:',msg
      write(6,10000)i_basis,j_basis,ish,jsh
10000 format(1x,'<ibasis:jbasis>-<',i4,':',i4,'>',5x,
     &       '<ish:jsh>-<',i3,':',i3,'>')
c
      do 00100 ii = 1,l1e
        write(6,10001)ilab(ii),jlab(ii),Gen1e(ii)
00100 continue
c
10001 format('i = ',i4,' j = ',i4,5x,'1eval = ',1pd20.10)
c
      end
