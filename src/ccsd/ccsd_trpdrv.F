      subroutine ccsd_trpdrv(t1,buff,bufo,lbfo,bufv,lbfv,
     &                       f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &                       eccsdt,g_objo,g_objv,g_coul,g_exch,
     &                       ncor,nocc,nvir,ndel,iprt,emp4,emp5,
     &                       oseg_lo,oseg_hi)
      implicit none
      double precision t1(*),buff(*),bufo(*),
     &                 bufv(*),f1n(*),f1t(*),f2n(*),
     &                 f2t(*),f3n(*),f3t(*),f4n(*),f4t(*),eorb(*),
     &                 emp4,emp5
      integer lbfo,lbfv,g_objo,g_objv,ncor,nocc,nvir,ndel,iprt,g_coul,
     &        g_exch,oseg_lo,oseg_hi
c
      double precision gmoint2,eaijk,eccsdt
      integer a,b,c,i,j,k,aold,av,inode,t2j,eintc,eintx,
     &        dint,t2x,t2c,fintc,fintx,len,ad1,ad2,ad3,
     &        next,nxtval
c
      double precision zip,one,half,two
      data zip/0.0d00/,one/1.0d00/,half/0.5d00/,two/2.0d00/
C
      Integer Nodes, IAm
C
#include "ccsd_len.fh"
#include "global.fh"
c
capr      call ga_print(g_coul)
capr      call ga_print(g_exch)
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
C
      call ga_sync()
      call qenter('trpdrv',0)
      inode=-1
      aold=0
      next=nxtval(-nodes)
      next=nxtval(nodes)
      do a=oseg_lo,oseg_hi
       av=a-ncor-nocc
       do j=1,nocc
        inode=inode+1
        if (inode.eq.next)then
c         print *,'iam,a,j',iam,a,j
         t2j=1
         eintc=t2j+lnovv
         eintx=eintc+lnoov
         len=eintx+lnoov-1
         call ga_get(g_objo,1,len,j,j,bufo,len)
         if (iprt.gt.10)then
          print *,'t2j ',a,j,iam
          write(6,4859)(bufo(i),i=t2j,t2j+lnovv-1)
          print *,'eintc ',a,j,iam
          write(6,4859)(bufo(i),i=eintc,eintc+lnoov-1)
          print *,'eintx ',a,j,iam
          write(6,4859)(bufo(i),i=eintx,eintx+lnoov-1)
         endif
         if (aold.ne.a)then
          dint=1
          t2c=dint+lnoov
          t2x=t2c+lnoov
          len=t2x+lnoov-1
          call ga_get(g_objv,1,len,av,av,bufv,len)
          fintc=1
          fintx=fintc+lnovv
c - make fints
CRK note change - unsure
          ad1=(a-oseg_lo)*nocc
          call ga_get(g_coul,1,lnvv,ad1+1,ad1+nocc,buff(fintc),lnvv)
          call ga_get(g_exch,1,lnvv,ad1+1,ad1+nocc,buff(fintx),lnvv)
c          ad1=0
c          do i=1,nocc
c           do b=nocc+1,nocc+nvir
c            do c=nocc+1,nocc+nvir
c             print *,'coul ',a,i,b,c,buff(fintc+ad1),gmoint2(a,i,b,c)
c             print *,'exch ',a,i,b,c,buff(fintx+ad1),gmoint2(a,b,i,c)
c             buff(fintc+ad1)=gmoint2(a,i,b,c)
c             buff(fintx+ad1)=gmoint2(a,b,i,c)
c             ad1=ad1+1
c            enddo
c           enddo
c          enddo
          if (iprt.gt.10)then
           print *,'fintc ',a,iam
           write(6,4859)(buff(i),i=fintc,fintc+lnovv-1)
           print *,'fintx ',a,iam
           write(6,4859)(buff(i),i=fintx,fintx+lnovv-1)
 4859      format(1x,5e14.4)
           print *,'t2c ',a,iam
           write(6,4859)(bufv(i),i=t2c,t2c+lnoov-1)
           print *,'t2x ',a,iam
           write(6,4859)(bufv(i),i=t2x,t2x+lnoov-1)
           print *,'dint ',a,iam
           write(6,4859)(bufv(i),i=dint,dint+lnoov-1)
          endif
          aold=a
         endif
c
         do i=1,nocc
          do k=1,i
           call dfill(lnvv,zip,f1n,1)
           call dfill(lnvv,zip,f1t,1)
           call dfill(lnvv,zip,f2n,1)
           call dfill(lnvv,zip,f2t,1)
           call dfill(lnvv,zip,f3n,1)
           call dfill(lnvv,zip,f3t,1)
           call dfill(lnvv,zip,f4n,1)
           call dfill(lnvv,zip,f4t,1)
           ad1=(i-1)*lnvv
           ad2=(k-1)*lnvv
           call ccsd_dovvv(buff(ad1+fintc),buff(ad1+fintx),
     $                     bufo(ad2+t2j),f1n,f2n,f3n,f4n,nocc,nvir)
           call ccsd_dovvv(buff(ad2+fintc),buff(ad2+fintx),
     $                     bufo(ad1+t2j),f1t,f2t,f3t,f4t,nocc,nvir)
           ad1=(k-1)*lnov
           ad2=(i-1)*lnov
           call ccsd_doooo(bufo(ad1+eintc),bufo(ad1+eintx),
     $                     bufv(ad2+t2c),bufv(ad2+t2x),f1n,f2n,
     $                     f3n,f4n,nocc,nvir)
           call ccsd_doooo(bufo(ad2+eintc),bufo(ad2+eintx),
     $                     bufv(ad1+t2c),bufv(ad1+t2x),f1t,f2t,
     $                     f3t,f4t,nocc,nvir)
           if (iprt.gt.50)then
            call prtfmat(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,nvir)
           endif
           eaijk=eorb(ncor+i)+eorb(ncor+j)+eorb(ncor+k)-eorb(a)
           ad1=(j-1)*lnov+(i-1)*nvir
           ad2=(i-1)*lnov+(j-1)*nvir
           ad3=(k-1)*nvir+1
           call ccsd_tengy(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     &                     bufv(dint+ad1),bufv(dint+ad2),
     $                     t1(ad3),eorb,eaijk,emp4,emp5,
     $                     ncor,nocc,nvir,ndel)
c          write(6,2234)i,k,iam,a,j,emp4,emp5
c 2234     format(' te1 iam aijk',5i5,2e15.5)
           if (i.ne.k)then
            ad1=(j-1)*lnov+(k-1)*nvir
            ad2=(k-1)*lnov+(j-1)*nvir
            ad3=(i-1)*nvir+1
            call ccsd_tengy(f1t,f1n,f2t,f2n,f3t,f3n,f4t,f4n,
     &                      bufv(dint+ad1),bufv(dint+ad2),
     $                      t1(ad3),eorb,eaijk,emp4,emp5,
     $                      ncor,nocc,nvir,ndel)
c          write(6,3234)i,k,iam,a,j,emp4,emp5
c 3234     format(' te2 iam aijk',5i5,2e15.5)
           endif
          enddo
         enddo
         if (iprt.gt.50)then
          write(6,1234)iam,a,j,emp4,emp5
 1234     format(' iam aijk',3i5,2e15.5)
         endif
         next=nxtval(nodes)
        endif
       enddo
      enddo
c      print *,'end of loops'
c
      next=nxtval(-nodes)
      call qexit('trpdrv',0)
      return
      end







      subroutine ccsd_doooo(eintc,eintx,t2c,t2x,f1,f2,f3,f4,nocc,nvir)
      implicit none
      double precision eintc(*),eintx(*),t2c(*),t2x(*),
     &                 f1(*),f2(*),f3(*),f4(*)
      integer nocc,nvir
      double precision one
      data one/1.0d00/
      call qenter('doooo',0)
c
c      f1n(b,c)=f1n(b,c)-t2(b,l,a,i)*vooo(c,k,l,j)
      call dgemm('n','n',nvir,nvir,nocc,-one,t2c,nvir,
     &           eintx,nocc,one,f1,nvir)
c      f2n(b,c)=f2n(b,c)-t2(a,l,b,i)*vooo(c,k,l,j)
      call dgemm('n','n',nvir,nvir,nocc,-one,t2x,nvir,
     &           eintx,nocc,one,f2,nvir)
c      f3n(b,c)=f3n(b,c)-t2(b,l,a,i)*vooo(c,j,l,k)
      call dgemm('n','n',nvir,nvir,nocc,-one,t2c,nvir,
     &           eintc,nocc,one,f3,nvir)
c      f4n(b,c)=f4n(b,c)-t2(a,l,b,i)*vooo(c,j,l,k)
      call dgemm('n','n',nvir,nvir,nocc,-one,t2x,nvir,
     &           eintc,nocc,one,f4,nvir)
      call qexit('doooo',0)
      return
      end






      subroutine ccsd_dovvv(fintc,fintx,t2,f1,f2,f3,f4,nocc,nvir)
      implicit none
      real*8 fintc(*),fintx(*),t2(*),f1(*),f2(*),f3(*),f4(*)
      integer nocc,nvir
      double precision one
      data one/1.0d00/
      call qenter('dovvv',0)
c
c --  f1n(b,c)=f1n(b,c)+vvvo(b,d,a,i)*t2(d,j,c,k)
      call dgemm('n','t',nvir,nvir,nvir,one,fintc,nvir,
     &           t2,nvir,one,f1,nvir)
c --  f2n(b,c)=f2n(b,c)+vvvo(a,d,b,i)*t2(d,j,c,k)
      call dgemm('n','t',nvir,nvir,nvir,one,fintx,nvir,
     &           t2,nvir,one,f2,nvir)
c --  f3n(b,c)=f3n(b,c)+vvvo(b,d,a,i)*t2(d,k,c,j)
      call dgemm('n','n',nvir,nvir,nvir,one,fintc,nvir,
     &           t2,nvir,one,f3,nvir)
c --  f4n(b,c)=f4n(b,c)+vvvo(a,d,b,i)*t2(d,k,c,j)
      call dgemm('n','n',nvir,nvir,nvir,one,fintx,nvir,
     &           t2,nvir,one,f4,nvir)
      call qexit('dovvv',0)
      return
      end

