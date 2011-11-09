      program pippo
      implicit none

      character*255 stcr
      integer zcomp,zcomm
      integer valcomm,valcomp
      integer newdist
      integer ga_Nodeid
      external ga_Nodeid
      integer getdist
      external getdist

      stcr='comp'
      call initialize_dist()
      zcomp = getdist( stcr )

      valcomp=123456
      call updist(zcomp, valcomp-1)
      call updist(zcomp, valcomp+10)
      call updist(zcomp, valcomp-1)
      call updist(zcomp, valcomp-2)
      call updist(zcomp, valcomp)
      call updist(zcomp, valcomp+10)


      call finalize_dist( ga_Nodeid() )

      stop
      end
c $Id$
