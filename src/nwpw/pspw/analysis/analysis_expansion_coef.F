c
c $Id$
c
      subroutine analysis_expansion_coef(filename,lorb0,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "btdb.fh"

      character*255 filename
      integer lorb0
      integer rtdb

*     **** local variables ****
      integer lmax,i,nrho,ind,lorb
      real*8  b0(0:5),u2(0:5),u3(0:5)
      real*8  x2,y2,x3,y3,dum
      character*2 atom
 
      character*255 id,full_filename


      lorb = lorb0

      call util_file_name_noprefix(filename,.false.,.false.,
     >                    full_filename)


      ind = index(full_filename,' ') - 1
      open(82,FILE=full_filename,STATUS='OLD',FORM='FORMATTED')
        read(82,'(A2)') atom
        read(82,*) dum,dum,lmax
        read(82,*)
        read(82,*) nrho,dum
        read(82,*)
        do i=1,nrho
          read(82,*)
        end do
        read(82,*) 
        read(82,*) x2,(u2(i),i=0,lmax)
        read(82,*) x3,(u3(i),i=0,lmax)
      close(82)
      if (lorb.lt.0) lorb = lmax-1


      do i=0,lorb
          Y2=U2(i)/X2**(i+1)
          Y3=U3(i)/X3**(i+1)
          b0(i)=Y2-X2*(Y2-Y3)/(X2-X3)
      end do

      WRITE(6,*)
      WRITE(6,1300) filename(1:ind)
      WRITE(6,1305) 'ATOM ','S','P','D','F'
      WRITE(6,1310) ATOM,(B0(i),i=0,lorb)
      WRITE(6,*)
      WRITE(6,*)
      call util_flush(6)
 1300 FORMAT('== ',A,' expansion coefficients =='/)
 1305 FORMAT(A5,6X,A,12X,A,12X,A,12X,A)
 1310 FORMAT(A2,' : ',4E13.5)

      id = 'analysis:lorb'//atom
      if (.not. btdb_put(rtdb,id,mt_int,1,lorb))
     >        call errquit(
     >        'analysis_input: btdb_put failed', 0, RTDB_ERR)

      id = 'analysis:expansion'//atom
      if (.not. btdb_put(rtdb,id,mt_dbl,(lorb+1),b0))
     >        call errquit(
     >        'analysis_input: btdb_put failed', 0, RTDB_ERR)

      return
      end
