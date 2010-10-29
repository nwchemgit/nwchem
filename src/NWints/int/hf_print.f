      subroutine hf_print(msg,xyz,prims,coefs,npa,nca,la)
c $Id$
      implicit none
c
      character*(*) msg    ! info to print 
      integer npa          ! number of primitives
      integer nca          ! number of general contractions
      integer la           ! type of basis function
c
      double precision xyz(3)         ! coordinates
      double precision prims(npa)     ! exponents
      double precision coefs(npa,nca) ! coefs
c
      integer hfprintit
      common /hfprtc/ hfprintit
c
      integer i,j
c
      if (hfprintit.eq.0) return
c
      write(6,10000)
      write(6,'(a)')msg
      write(6,10000)
c
      write(6,10001)xyz
      write(6,'(a,i12,a,i12,a,i12,a)')
     &       ' <type:',la,'>  <npx:',npa,'>  <ncx:',nca,'>'
      write(6,*)' exponents | coefficients '
      do 00100 i=1,(min(npa,100))
        write(6,10002)i,prims(i),' | ',(coefs(i,j),j=1,nca)
00100 continue
c
      write(6,10000)
c
10000 format(80('='))
10001 format(' <x:',1pd12.4,'>  <y:',1pd12.4,'>  <z:',1pd12.4,'>')
10002 format(1x,i3,1pd14.6,a,1pd14.6,1pd14.6,1pd14.6,1pd14.6,
     &       1pd14.6,1pd14.6,1pd14.6,1pd14.6,1pd14.6,1pd14.6)
      end
      subroutine hf_print_set(value)
      implicit none
c
      integer hfprintit
      common /hfprtc/ hfprintit
c
      integer value
c
      hfprintit = value
      end



