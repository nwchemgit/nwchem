      subroutine argos_diana_anahdr(sgmnam)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      character*16 sgmnam(nsa)
      character*11 atom(4)
      integer i,j
c
      write(lfnana,1000) ndist,nangl,ntors,nimps
 1000 format(4i10)
c
      do 1 i=1,ndist
      do 11 j=1,2
      if(idist(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(idist(i,j,1))(12:16),':',
     + sgmnam(idist(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+idist(i,j,1),':',
     + wnam(idist(i,j,2))(6:10)
      endif
   11 continue
c
      write(lfnana,1001) (atom(j),j=1,2)
 1001 format(a11,'-',a11)
    1 continue
c
      do 2 i=1,nangl
      do 21 j=1,3
      if(iangl(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(iangl(i,j,1))(12:16),':',
     + sgmnam(iangl(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+iangl(i,j,1),':',
     + wnam(iangl(i,j,2))(6:10)
      endif
   21 continue
      write(lfnana,1002) (atom(j),j=1,3)
 1002 format(a11,'-',a11,'-',a11)
    2 continue
c
      do 3 i=1,ntors
      do 31 j=1,4
      if(itors(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(itors(i,j,1))(12:16),':',
     + sgmnam(itors(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+itors(i,j,1),':',
     + wnam(itors(i,j,2))(6:10)
      endif
   31 continue
      write(lfnana,1003) (atom(j),j=1,4)
 1003 format(a11,'-',a11,'-',a11,'-',a11)
    3 continue
c
      do 4 i=1,nimps
      do 41 j=1,4
      if(iimps(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(iimps(i,j,1))(12:16),':',
     + sgmnam(iimps(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+iimps(i,j,1),':',
     + wnam(iimps(i,j,2))(6:10)
      endif
   41 continue
      write(lfnana,1003) (atom(j),j=1,4)
    4 continue
c
      return
      end
c $Id$
