      subroutine ccsd_trpdrv(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     eccsdt,g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi,
     $     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     $     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
C     $Id: ccsd_trpdrv.F,v 2.8 1997-10-27 18:04:54 d3g681 Exp $
      implicit none
c
#include "global.fh"
#include "ccsd_len.fh"
c     
      double precision t1(*),
     &     f1n(*),f1t(*),f2n(*),
     &     f2t(*),f3n(*),f3t(*),f4n(*),f4t(*),eorb(*),
     &     emp4,emp5
      double precision Tij(*), Tkj(*), Tia(*), Tka(*), Xia(*), Xka(*), 
     $     Jia(*), Jka(*), Kia(*), Kka(*),
     $     Jij(*), Jkj(*), Kij(*), Kkj(*), Dja(*), Djka(*), Djia(*)

      integer g_objo,g_objv,ncor,nocc,nvir,iprt,g_coul,
     &     g_exch,oseg_lo,oseg_hi
c     
      double precision eaijk,eccsdt
      integer a,i,j,k,akold,av,inode,len,ad1,ad2,ad3,
     &     next,nxtval
      external nxtval
C     
      Integer Nodes, IAm
c     
*      double precision Tij(1000), Tkj(10000) ! lnvv and kchunk*lnvv
c     
*      double precision Tia(1000), Tka(10000) ! lnov and kchunk*lnov
*      double precision Xia(1000), Xka(10000) ! lnov and kchunk*lnov
c     
*      double precision Jia(1000), Jka(10000) ! lnvv and kchunk*lnvv
*      double precision Kia(1000), Kka(10000) ! lnvv and kchunk*lnvv
c     
*      double precision Jij(1000), Jkj(10000) ! lnov and kchunk*lnov
*      double precision Kij(1000), Kkj(10000) ! lnov and kchunk*lnov
c     
c     rjh: D is something to do with the (T) correction ?
*      double precision Dja(1000), Djka(1000), Djia(100) ! lnov, nvir*kchunk, nvir
c     
      integer klo, khi, start, end
      integer kchunk
C     
      double precision zip
      data zip/0.0d00/
c     
c     apr      call ga_print(g_coul)
c     apr      call ga_print(g_exch)
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
C     
      call ga_sync()
      call qenter('trpdrv',0)
      inode=-1
      next=nxtval(nodes)
c     
      do klo = 1, nocc, kchunk
         akold=0
         khi = min(nocc, klo+kchunk-1)
         do a=oseg_lo,oseg_hi
            av=a-ncor-nocc
            do j=1,nocc
               inode=inode+1
               if (inode.eq.next)then
c     
c     Get Dja = Dci,ja for given j, a, all ci
c     
                  start = 1 + (j-1)*lnov
                  len   = lnov
                  end   = start + len - 1
                  call ga_get(g_objv,start,end,av,av,Dja,len)
c     
c     Get Tkj = T(b,c,k,j) for given j, klo<=k<=khi, all bc
c     
                  start = (klo-1)*lnvv + 1
                  len   = (khi-klo+1)*lnvv
                  end   = start + len - 1
                  call ga_get(g_objo,start,end,j,j,Tkj,len)
c     
c     Get Jkj = J(c,l,k,j) for given j, klo<=k<=khi, all cl
c     
                  start = lnovv + (klo-1)*lnov + 1
                  len   = (khi-klo+1)*lnov
                  end   = start + len - 1
                  call ga_get(g_objo,start,end,j,j,Jkj,len)
c     
c     Get Kkj = K(c,l,k,j) for given j, klo<=k<=khi, all cl
c     
                  start = lnovv + lnoov + (klo-1)*lnov + 1
                  len   = (khi-klo+1)*lnov
                  end   = start + len - 1
                  call ga_get(g_objo,start,end,j,j,Kkj,len)
c     
                  if (akold .ne. a) then
                     akold = a
c     
c     Get Jka = J(b,c,k,a) for given a, klo<=k<=khi, all bc
c     
                     start = (a-oseg_lo)*nocc + klo
                     len   = (khi-klo+1)
                     end   = start + len - 1
                     call ga_get(g_coul,1,lnvv,start,end,Jka,lnvv)
c     
c     Get Kka = K(b,c,k,a) for given a, klo<=k<=khi, all bc
c     
                     start = (a-oseg_lo)*nocc + klo
                     len   = (khi-klo+1)
                     end   = start + len - 1
                     call ga_get(g_exch,1,lnvv,start,end,Kka,lnvv)
c     
c     Get Tka = Tbl,ka for given a, klo<=k<=khi, all bl
c     
                     start = 1 + lnoov + (klo-1)*lnov
                     len   = (khi-klo+1)*lnov
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,Tka,len)
c     
c     Get Xka = Tal,kb for given a, klo<=k<=khi, all bl
c     
                     start = 1 + lnoov + lnoov + (klo-1)*lnov
                     len   = (khi-klo+1)*lnov
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,Xka,len)
                  endif
c     
c     Get Djka = Dcj,ka for given j, a, klo<=k<=khi, all c
c     
                  do k = klo, khi
                     start = 1 + (j-1)*nvir + (k-1)*lnov
                     len   = nvir
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,
     $                    Djka(1+(k-klo)*nvir),len)
                  enddo
c     
                  do i=1,nocc
c     
c     Get Tij = T(b,c,i,j) for given j, i, all bc
c     
                     start = (i-1)*lnvv + 1
                     len   = lnvv
                     end   = start + len - 1
                     call ga_get(g_objo,start,end,j,j,Tij,len)
c     
c     Get Jij = J(c,l,i,j) for given j, i, all cl
c     
                     start = lnovv + (i-1)*lnov + 1
                     len   = lnov
                     end   = start + len - 1
                     call ga_get(g_objo,start,end,j,j,Jij,len)
c     
c     Get Kij = K(c,l,i,j) for given j, i, all cl
c     
                     start = lnovv + lnoov + (i-1)*lnov + 1
                     len   = lnov
                     end   = start + len - 1
                     call ga_get(g_objo,start,end,j,j,Kij,len)
c     
c     Get Jia = J(b,c,i,a) for given a, i, all bc
c     
                     start = (a-oseg_lo)*nocc + i
                     len   = 1
                     end   = start + len - 1
                     call ga_get(g_coul,1,lnvv,start,end,Jia,lnvv)
c     
c     Get Kia = K(b,c,i,a) for given a, i, all bc
c     
                     start = (a-oseg_lo)*nocc + i
                     len   = 1
                     end   = start + len - 1
                     call ga_get(g_exch,1,lnvv,start,end,Kia,lnvv)
c     
c     Get Dia = Dcj,ia for given j, i, a, all c
c     
                     start = 1 + (j-1)*nvir + (i-1)*lnov
                     len   = nvir
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,Djia,len)
c     
c     Get Tia = Tbl,ia for given a, i, all bl
c     
                     start = 1 + lnoov + (i-1)*lnov
                     len   = lnov
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,Tia,len)
c     
c     Get Xia = Tal,ib for given a, i, all bl
c     
                     start = 1 + lnoov + lnoov + (i-1)*lnov
                     len   = lnov
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,Xia,len)
c     
                     do k=klo,min(khi,i)
                        call dfill(lnvv,zip,f1n,1)
                        call dfill(lnvv,zip,f1t,1)
                        call dfill(lnvv,zip,f2n,1)
                        call dfill(lnvv,zip,f2t,1)
                        call dfill(lnvv,zip,f3n,1)
                        call dfill(lnvv,zip,f3t,1)
                        call dfill(lnvv,zip,f4n,1)
                        call dfill(lnvv,zip,f4t,1)
c     
c     sum(d) (Jia, Kia)bd * Tkj,cd -> Fbc
c     
                        call ccsd_dovvv(Jia, Kia,
     $                       Tkj(1+(k-klo)*lnvv),f1n,f2n,f3n,f4n,
     $                       nocc,nvir)
c     
c     sum(d) (Jka, Kka)bd * Tij,cd -> Fbc
c     
                        call ccsd_dovvv(Jka(1+(k-klo)*lnvv),
     $                       Kka(1+(k-klo)*lnvv),
     $                       Tij,f1t,f2t,f3t,f4t,nocc,nvir)
c     
                        ad1=(k-1)*lnov
                        ad2=(i-1)*lnov
c     
c     sum(l) (Jij, Kij)cl  * Tkl,ab -> Fbc
c     
                        call ccsd_doooo(Jkj(1+(k-klo)*lnov),
     $                       Kkj(1+(k-klo)*lnov),
     $                       Tia,Xia,
     $                       f1n,f2n,
     $                       f3n,f4n,nocc,nvir)
c     
c     sum(l) (Jkj, Kkj)cl  * Tli,ba -> Fbc
c     
                        call ccsd_doooo(Jij, Kij,
     $                       Tka(1+(k-klo)*lnov),Xka(1+(k-klo)*lnov),
     $                       f1t,f2t,
     $                       f3t,f4t,nocc,nvir)
                        if (iprt.gt.50)then
                           call prtfmat(f1n,f1t,f2n,f2t,f3n,f3t,f4n,
     $                          f4t, nvir)
                        end if
                        eaijk=eorb(ncor+i)+eorb(ncor+j)+eorb(ncor+k)-
     $                       eorb(a)
                        ad1=(j-1)*lnov+(i-1)*nvir
                        ad2=(i-1)*lnov+(j-1)*nvir
                        ad3=(k-1)*nvir+1
                        call ccsd_tengy(f1n,f1t,f2n,f2t,f3n,f3t,f4n,
     $                       f4t,
     &                       Dja(1+(i-1)*nvir),Djia,
     $                       t1(ad3),eorb,eaijk,emp4,emp5,
     $                       ncor,nocc,nvir)
c     
                        if (i.ne.k)then
                           ad1=(j-1)*lnov+(k-1)*nvir
                           ad2=(k-1)*lnov+(j-1)*nvir
                           ad3=(i-1)*nvir+1
                           call ccsd_tengy(f1t,f1n,f2t,f2n,
     $                          f3t,f3n,f4t,f4n,
     &                          Dja(1+(k-1)*nvir),Djka(1+(k-klo)*nvir), 
     $                          t1(ad3),eorb,eaijk,emp4,emp5,
     $                          ncor,nocc,nvir)
c     
                        end if
                     end do
                  end do
                  if (iprt.gt.50)then
                     write(6,1234)iam,a,j,emp4,emp5
 1234                format(' iam aijk',3i5,2e15.5)
                  end if
                  next=nxtval(nodes)
               end if
            end do
         end do
      end do
c     
      next=nxtval(-nodes)
      call ga_sync
      call qexit('trpdrv',0)
c     
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
      double precision fintc(*),fintx(*),t2(*),f1(*),f2(*),f3(*),f4(*)
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

