c
c     $Id$
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
      subroutine MOTION_nframes(filename,nframes,nion,omega,unita)
      implicit none
      character*(*) filename
      integer nframes,nion
      real*8  omega,unita(3,3)

      real*8    time
      integer   i,nion1

      nframes = 0
      nion    = 0

      open(unit=19,file=filename,form='formatted',status='old')
      do while(.true.)
        read(19,*,ERR=30,END=30) time,nion1,omega,unita
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

*     *********************************
*     *                               *
*     *         xyz_nframes           *
*     *                               *
*     *********************************
*     
*    This routine returns the number of frames in
*  the xyz file name filename.  
*     
      subroutine xyz_nframes(filename,nframes,nion)
      implicit none
      character*(*) filename
      integer nframes,nion
      
      integer   i,nion1

      nframes = 0
      nion    = 0
      
      open(unit=19,file=filename,form='formatted',status='old')
      do while(.true.)
        read(19,*,ERR=30,END=30) nion1
        read(19,*,ERR=30,END=30)
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


