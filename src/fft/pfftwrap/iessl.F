*
* $Id$
*
c************************************************************************
      INTEGER FUNCTION initESSL(icontext)
      IMPLICIT NONE
#include "errquit.fh"
      INTEGER icontext
      integer nprow, npcol, myrow, mycol
      integer nnodes, nodeid, nproc, me
      integer icontext1, icontext2
      common /xparms/nprow,npcol,myrow,mycol

      nproc = nnodes()
      me = nodeid()
      nprow = 1
      npcol = nproc
c obtain system default and map to a row
      call blacs_get(0,0,icontext1)
      call blacs_gridinit(icontext1,'R',1,nproc)
c reconfigure (doesn't seem to be necessary)
c     call blacs_get(icontext1,10,icontext2)
c     call blacs_gridinit(icontext2,'R',nprow,npcol)
c     print*,me,' icontext2=',icontext2,nproc
      call blacs_gridinfo(icontext2,nprow,npcol,myrow,mycol)
      initESSL = 0
      if (nprow.NE.myrow.OR.npcol.NE.mycol) then
        initESSL = -1
        if (me.EQ.0) then
          write(6,*) 'proc grid incorrect: ',nprow,myrow,npcol,mycol
        endif 
        call errquit('initESSL error initializing',me, UNKNOWN_ERR)
      endif
      icontext = icontext1
      RETURN
      END 

c************************************************************************
