c
c     $Id: MOTION_nframes.F,v 1.2 2003-02-12 20:12:19 edo Exp $
c

*     *********************************
*     *                               *
*     *        MOTION_nframes         *
*     *                               *
*     *********************************
*
*    This routine returns the number of frames in 
*  the MOTION file name filename.  This routine is
*  now needed because the format of the MOTION file
*  has changed.
*
      subroutine MOTION_nframes(filename,nframes,nion,omega)
      implicit none
      character*(*) filename
      integer nframes,nion
      real*8  omega

      real*8    time
      integer   i,nion1

      nframes = 0
      nion    = 0

      open(unit=19,file=filename,form='formatted',status='old')
      do while(.true.)
        read(19,*,ERR=30,END=30) time,nion1,omega
        if (nion1.gt.nion) nion = nion1
        do i=1,nion1
          read(19,*,ERR=30,END=30) 
        end do
        nframes = nframes + 1
      end do
 30   continue
      close(unit=19)

      return
      end

