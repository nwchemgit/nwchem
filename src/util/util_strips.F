      subroutine util_gacc_strips(g_a,ilo,ihi,jlo,jhi,
     C     buf,ld,alpha,lgth_in)
      implicit none
#include "global.fh"
      integer g_a
      integer ilo,ihi,jlo,jhi,ld ![in]
      double precision alpha     ![in]
      double precision buf(*)   ![in]
      integer lgth_in ![in]
c
      integer jj,acc_ilo,acc_ihi,acc_ilgt
      integer iptr
      integer npass
#define ACC_LENGTH 800
#define ACC_LENGTH_MAX 950
#define ACC_LENGTH_MAX_ 1950
c
      if(lgth_in.gt.0) then
c         npass=(ld*(jhi-jlo+1))/max(lgth_in+1
cold         npass=(ld*(jhi-jlo+1))/max(ACC_LENGTH,lgth_in)+1
         npass=(ld*(jhi-jlo+1))/min(ACC_LENGTH_MAX,lgth_in)+1
      else
         npass=(ld*(jhi-jlo+1))/ACC_LENGTH+1
      endif
      acc_ilgt=ld/npass
      do jj=1,npass
         iptr=(jj-1)*acc_ilgt
         acc_ilo=ilo+iptr
         acc_ihi=acc_ilo+acc_ilgt-1
         if(acc_ihi.gt.ihi) then
            write(6,*) ' acc_ihi = ',acc_ihi, ' gt ihi ',ihi
            call errquit(' booh ',0,0)
         endif
         call ga_acc(g_a,acc_ilo,acc_ihi,jlo,jhi,
     S        buf(iptr+1),
     &        ld,alpha)
      enddo
c     cleanup
      if(acc_ihi.lt.ihi) then
         iptr=npass*acc_ilgt
         call ga_acc(g_a,acc_ihi+1,ihi,jlo,jhi,
     S        buf(iptr+1),
     &        ld,alpha)
      endif
      return
      end
