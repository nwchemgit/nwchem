      subroutine ccsd_trpdrv_nb(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     eccsdt,g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi,
     $     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     $     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
c
C     $Id: ccsd_trpdrv_nb.F,v 2.10 2004-12-16 23:02:58 edo Exp $
c
c     CCSD(T) non-blocking modifications written by 
c     Jeff Hammond, Argonne Leadership Computing Facility
c     Fall 2009
c
      implicit none
c
#include "global.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
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
     &     next,nxtask
      external nxtask
c     
      Integer Nodes, IAm
c     
      integer klo, khi, start, end
      integer kchunk
c
c==================================================
c
c  NON-BLOCKING stuff
c
c==================================================
c
c  Dependencies (global array, local array, handle): 
c
c      g_objv, Dja, nbh_objv1
c      g_objv, Tka, nbh_objv2
c      g_objv, Xka, nbh_objv3
c      g_objv, Djka(1+(k-klo)*nvir), nbh_objv4(k)
c      g_objv, Djia, nbh_objv5
c      g_objv, Tia, nbh_objv6
c      g_objv, Xia, nbh_objv7
c      g_objo, Tkj, nbh_objo1
c      g_objo, Jkj, nbh_objo2
c      g_objo, Kkj, nbh_objo3
c      g_objo, Tij, nbh_objo4
c      g_objo, Jij, nbh_objo5
c      g_objo, Kij, nbh_objo6
c      g_exch, Kka, nbh_exch1
c      g_exch, Kia, nbh_exch2
c      g_coul, Jka, nbh_coul1
c      g_coul, Jia, nbh_coul2
c
c  non-blocking handles
c
       integer nbh_objv1,nbh_objv2,nbh_objv3
       integer nbh_objv5,nbh_objv6,nbh_objv7
       integer nbh_objv4(nocc)
c
       integer nbh_objo1,nbh_objo2,nbh_objo3
       integer nbh_objo4,nbh_objo5,nbh_objo6
c
       integer nbh_exch1,nbh_exch2,nbh_coul1,nbh_coul2
c
c==================================================
c     
      double precision zip
      data zip/0.0d00/
c     
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
c     
      call ga_sync()
c
      if (occsdps) then
         call pstat_on(ps_trpdrv)
      else
         call qenter('trpdrv',0)
      endif
      inode=-1
      next=nxtask(nodes, 1)
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
                  call ga_nbget(g_objv,start,end,av,av,Dja,len,
     1                          nbh_objv1)
c     
c     Get Tkj = T(b,c,k,j) for given j, klo<=k<=khi, all bc
c     
                  start = (klo-1)*lnvv + 1
                  len   = (khi-klo+1)*lnvv
                  end   = start + len - 1
                  call ga_nbget(g_objo,start,end,j,j,Tkj,len,
     1                          nbh_objo1)
c     
c     Get Jkj = J(c,l,k,j) for given j, klo<=k<=khi, all cl
c     
                  start = lnovv + (klo-1)*lnov + 1
                  len   = (khi-klo+1)*lnov
                  end   = start + len - 1
                  call ga_nbget(g_objo,start,end,j,j,Jkj,len,
     1                          nbh_objo2)
c     
c     Get Kkj = K(c,l,k,j) for given j, klo<=k<=khi, all cl
c     
                  start = lnovv + lnoov + (klo-1)*lnov + 1
                  len   = (khi-klo+1)*lnov
                  end   = start + len - 1
                  call ga_nbget(g_objo,start,end,j,j,Kkj,len,
     1                          nbh_objo3)
c     
                  if (akold .ne. a) then
                     akold = a
c     
c     Get Jka = J(b,c,k,a) for given a, klo<=k<=khi, all bc
c     
                     start = (a-oseg_lo)*nocc + klo
                     len   = (khi-klo+1)
                     end   = start + len - 1
                     call ga_nbget(g_coul,1,lnvv,start,end,Jka,lnvv,
     1                             nbh_coul1)
c     
c     Get Kka = K(b,c,k,a) for given a, klo<=k<=khi, all bc
c     
                     start = (a-oseg_lo)*nocc + klo
                     len   = (khi-klo+1)
                     end   = start + len - 1
                     call ga_nbget(g_exch,1,lnvv,start,end,Kka,lnvv,
     1                             nbh_exch1)
c     
c     Get Tka = Tbl,ka for given a, klo<=k<=khi, all bl
c     
                     start = 1 + lnoov + (klo-1)*lnov
                     len   = (khi-klo+1)*lnov
                     end   = start + len - 1
                     call ga_nbget(g_objv,start,end,av,av,Tka,len,
     1                             nbh_objv2)
c     
c     Get Xka = Tal,kb for given a, klo<=k<=khi, all bl
c     
                     start = 1 + lnoov + lnoov + (klo-1)*lnov
                     len   = (khi-klo+1)*lnov
                     end   = start + len - 1
                     call ga_nbget(g_objv,start,end,av,av,Xka,len,
     1                             nbh_objv3)
                  endif
c     
c     Get Djka = Dcj,ka for given j, a, klo<=k<=khi, all c
c     
                  do k = klo, khi
                     start = 1 + (j-1)*nvir + (k-1)*lnov
                     len   = nvir
                     end   = start + len - 1
                     call ga_nbget(g_objv,start,end,av,av,
     1                    Djka(1+(k-klo)*nvir),len,nbh_objv4(k)) ! k <= nocc
                  enddo
c     
                  do i=1,nocc
c     
c     Get Tij = T(b,c,i,j) for given j, i, all bc
c     
                     start = (i-1)*lnvv + 1
                     len   = lnvv
                     end   = start + len - 1
                     call ga_nbget(g_objo,start,end,j,j,Tij,len,
     1                             nbh_objo4)
c     
c     Get Jij = J(c,l,i,j) for given j, i, all cl
c     
                     start = lnovv + (i-1)*lnov + 1
                     len   = lnov
                     end   = start + len - 1
                     call ga_nbget(g_objo,start,end,j,j,Jij,len,
     1                             nbh_objo5)
c     
c     Get Kij = K(c,l,i,j) for given j, i, all cl
c     
                     start = lnovv + lnoov + (i-1)*lnov + 1
                     len   = lnov
                     end   = start + len - 1
                     call ga_nbget(g_objo,start,end,j,j,Kij,len,
     1                             nbh_objo6)
c     
c     Get Jia = J(b,c,i,a) for given a, i, all bc
c     
                     start = (a-oseg_lo)*nocc + i
                     len   = 1
                     end   = start + len - 1
                     call ga_nbget(g_coul,1,lnvv,start,end,Jia,lnvv,
     1                             nbh_coul2)
c     
c     Get Kia = K(b,c,i,a) for given a, i, all bc
c     
                     start = (a-oseg_lo)*nocc + i
                     len   = 1
                     end   = start + len - 1
                     call ga_nbget(g_exch,1,lnvv,start,end,Kia,lnvv,
     1                             nbh_exch2)
c     
c     Get Dia = Dcj,ia for given j, i, a, all c
c     
                     start = 1 + (j-1)*nvir + (i-1)*lnov
                     len   = nvir
                     end   = start + len - 1
                     call ga_nbget(g_objv,start,end,av,av,Djia,len,
     1                             nbh_objv5)
c     
c     Get Tia = Tbl,ia for given a, i, all bl
c     
                     start = 1 + lnoov + (i-1)*lnov
                     len   = lnov
                     end   = start + len - 1
                     call ga_nbget(g_objv,start,end,av,av,Tia,len,
     1                             nbh_objv6)
c     
c     Get Xia = Tal,ib for given a, i, all bl
c     
                     start = 1 + lnoov + lnoov + (i-1)*lnov
                     len   = lnov
                     end   = start + len - 1
                     call ga_nbget(g_objv,start,end,av,av,Xia,len,
     1                             nbh_objv7)
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
                  next=nxtask(nodes, 1)
               end if
            end do
         end do
      end do
c     
      next=nxtask(-nodes, 1)
      call ga_sync
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
c     
      end

