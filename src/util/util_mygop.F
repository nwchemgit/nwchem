      subroutine util_mygadgop(msgid,buff,size,op)
      implicit none
      integer msgid ! [in]
      double precision buff(*) ! [in/out]
      integer size ! [in]
      character*(*) op
#include "mafdecls.fh"      
#include "global.fh"
c#define DEBUG 1
c
      integer bigint,nsteps,i,istart,len
c
#ifdef EXT_INT
      bigint = (2**31)-1024
#else
      bigint = (2**30)-1024
#endif
      bigint=bigint/ma_sizeof(mt_int, 1, mt_byte)

      nsteps = ceiling(dble(size)/dble(bigint))
      istart=1
#ifdef DEBUG         
      if(ga_nodeid().eq.0) 
     w     write(6,'(a,i10,a,i10,a,i10)') ' size ',size,
     w     ' nsteps ',nsteps,
     w ' bigint ',bigint
#endif
      do i=1,nsteps
         len=bigint
         if (istart+len-1.gt.size) len=size - istart + 1
#ifdef DEBUG         
         if(ga_nodeid().eq.0) 
     w        write(6,'(i10,a,i10,a,i10)') i,' istart ',istart,
     w        ' len ',len
#endif
         call ga_dgop(msgid+(i-1)*1000,buff(istart),len,op)
         istart=istart+len
      enddo

      return
      end
      
      subroutine util_mygaigop(msgid,ibuff,size,op)
      implicit none
      integer msgid ! [in]
      integer ibuff(*) ! [in/out]
      integer size ! [in]
      character*(*) op
#include "mafdecls.fh"      
#include "global.fh"      
c
      integer bigint,nsteps,i,istart,len
c
#ifdef EXT_INT
      bigint = (2**31)-1024
#else
      bigint = (2**30)-1024
#endif
      bigint=bigint/ma_sizeof(mt_int, 1, mt_byte)

      nsteps = ceiling(dble(size)/dble(bigint))
      istart=1
      do i=1,nsteps
         len=bigint
         if (istart+len-1.gt.size) len=size - istart + 1
#ifdef DEBUG         
         if(ga_nodeid().eq.0) 
     w        write(6,'(i10,a,i10,a,i10)') i,' istart ',istart,
     w        ' len ',len
#endif
         call ga_igop(msgid+(i-1)*1000,ibuff(istart),len,op)
         istart=istart+len
      enddo

      return
      end
      
