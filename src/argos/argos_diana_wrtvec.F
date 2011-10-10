      logical function argos_diana_wrtvec(ivec,fil,ndex,xa)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ivec
      character*255 fil
      integer i,j,k,l
      real*8 value,xa(msa,3)
      integer ndex(nsa)
c
      if(me.eq.0) then
      open(unit=lfntcl,file=filtcl(1:index(filtcl,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnvec,file=filvec(1:index(filvec,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnecc,file=filecc(1:index(filecc,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnamp,file=filamp(1:index(filamp,' ')-1),
     + form='formatted',status='unknown',err=9999)
      k=0
      l=0
      do 15 i=1,nsel
      value=0.0d0
      do 16 j=1,3
      value=value+(dbl_mb(i_cov+k))**2
      k=k+1
   16 continue
      l=l+1
      write(lfntcl,5555)
     + (10.0*xa(ndex(i),j),j=1,3),
     + (dbl_mb(i_cov+k-4+j),j=1,3)
 5555 format('draw vector2 {',3f12.6,'} {',3f12.6,'} 1.0 10')
      write(lfnamp,'(i5,f12.6)') i,sqrt(value)
   15 continue
c
      write(lfnecc,3013)
 3013 format('%begin%input file%1%char',/,
     + 'geometry noautosym units nm')
      do 314 i=1,nsa
      if(int_mb(i_osel+i-1).eq.1) then
      if(byte_mb(i_snam+(i-1)*16+5).eq.'2'.or.
     + byte_mb(i_snam+(i-1)*16+5).eq.'3'.or.
     + byte_mb(i_snam+(i-1)*16+5).eq.'4') then
      write(lfnecc,3014) ' ',byte_mb(i_snam+(i-1)*16+6),
     + (dbl_mb(i_xadat+i-1+j*msa),j=0,2)
      else
      write(lfnecc,3014) (byte_mb(i_snam+(i-1)*16+j),j=5,6),
     + (dbl_mb(i_xadat+i-1+j*msa),j=0,2)
      endif
 3014 format(2a1,1x,3f12.6)
      endif
  314 continue
      write(lfnecc,3012) nselo
 3012 format('end',/,'basis "ao basis" spherical print',/,
     + ' C library "6-31+G*"',/,
     + ' N library "6-31+G*"',/,
     + ' O library "6-31+G*"',/,
     + ' H library "6-31+G*"',/,
     + ' S library "6-31+G*"',/,
     + ' P library "6-31+G*"',/,
     + ' Ca library "6-31+G*"',/,
     + ' Na library "6-31+G*"',/,
     + 'end',/,
     + 'task scf gradient',/,
     + '%end%input file%1%char',/,
     + 'task_gradient%begin%total gradient%3 ',i7,'%double')
c
      k=0
      do 315 i=1,nsa
      if(int_mb(i_osel+i-1).eq.1) then
      if(int_mb(i_isel+i-1).eq.0) then
      write(lfnecc,3015) 0.0d0,0.0d0,0.0d0
 3015 format(3f12.6)
      else
      write(lfnecc,3015) (dbl_mb(i_cov+k+j),j=0,2)
      k=k+3
      endif
      endif
  315 continue
      write(lfnecc,3011) nselo
 3011 format('task_gradient%end%total gradient%3 ',i7,'%double')
c
      close(unit=lfntcl)
      close(unit=lfnvec)
      close(unit=lfnecc)
      close(unit=lfnval)
      close(unit=lfnprj)
      close(unit=lfnamp)
      write(*,'(/,a,i5,a,a)') ' Projection of vector ',ivec,
     +' written to ',fil(1:index(fil,' ')-1)
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
      argos_diana_wrtvec=.true.
      return
 9999 continue
      argos_diana_wrtvec=.false.
      return
      end
c $Id$
