      subroutine ccsd_trpdrv(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi,
     $     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     $     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
C     $Id$
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
      double precision eaijk
      integer a,b,c,i,j,k,akold,av,inode,len,
     &     next,nxtask
      external nxtask
C     
      Integer Nodes, IAm
c     
      integer klo, khi, start, end
      integer kchunk
c 
c
      double precision zip
      data zip/0.0d00/
c     
c     apr      call ga_print(g_coul)
c     apr      call ga_print(g_exch)
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
C     
      call ga_sync()
      if (occsdps) then
         call pstat_on(ps_trpdrv)
      else
         call qenter('trpdrv',0)
      endif
      inode=-1
      next=nxtask(nodes, 1)
c     
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
                  call ccsd_getdjka(djka,
     G                 g_objv,nvir,lnov,
     K                 j,av,klo,khi)
c     
c     Get Tia = Tbl,ia for given a,  all i,bl
c     
                     start = 1 + lnoov 
                     len   = lnov*nocc
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,Tia,len)
c     
c     Get Xia = Tal,ib for given a, all i,bl
c     
                     start = 1 + lnoov + lnoov 
                     len   = lnov*nocc
                     end   = start + len - 1
                     call ga_get(g_objv,start,end,av,av,Xia,len)
c     
c     Get Jij = J(c,l,i,j) for given j, all i,cl
c     
                     start = lnovv +  1
                     len   = lnov*nocc
                     end   = start + len - 1
                     call ga_get(g_objo,start,end,j,j,Jij,len)
c     
c     Get Kij = K(c,l,i,j) for given j, all,i cl
c     
                     start = lnovv + lnoov +  1
                     len   = lnov*nocc
                     end   = start + len - 1
                     call ga_get(g_objo,start,end,j,j,Kij,len)
c     
c     Get Dia = Dcj,ia for given j, i, a, all c
c     
                     call ccsd_getdjka(djia,
     G                    g_objv,nvir,lnov,
     K                    j,av,1,nocc)
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
                     do k=klo,min(khi,i)
                        call dcopy(lnvv,zip,0,f1n,1)
                        call dcopy(lnvv,zip,0,f1t,1)
                        call dcopy(lnvv,zip,0,f2n,1)
                        call dcopy(lnvv,zip,0,f2t,1)
                        call dcopy(lnvv,zip,0,f3n,1)
                        call dcopy(lnvv,zip,0,f3t,1)
                        call dcopy(lnvv,zip,0,f4n,1)
                        call dcopy(lnvv,zip,0,f4t,1)
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
c     sum(l) (Jij, Kij)cl  * Tkl,ab -> Fbc
c     
                        call ccsd_doooo(Jkj(1+(k-klo)*lnov),
     $                       Kkj(1+(k-klo)*lnov),
     $                       Tia(1+(i-1)*lnov),Xia(1+(i-1)*lnov),
     $                       f1n,f2n,
     $                       f3n,f4n,nocc,nvir)
c     
c     sum(l) (Jkj, Kkj)cl  * Tli,ba -> Fbc
c     
                        call ccsd_doooo(Jij(1+(i-1)*lnov), 
     K                       Kij(1+(i-1)*lnov),
     $                       Tka(1+(k-klo)*lnov),Xka(1+(k-klo)*lnov),
     $                       f1t,f2t,
     $                       f3t,f4t,nocc,nvir)
c
                        if (iprt.gt.50)then
                           call prtfmat(f1n,f1t,f2n,f2t,f3n,f3t,f4n,
     $                          f4t, nvir)
                        end if
c
                        eaijk=eorb(ncor+i)+eorb(ncor+j)+eorb(ncor+k)
     $                       -eorb(a)
c
                        call ccsd_tengy(f1n,f1t,f2n,f2t,f3n,f3t,f4n,
     $                       f4t,
     &                       Dja(1+(i-1)*nvir),Djia(1+(i-1)*nvir),
     $                       t1((k-1)*nvir+1),eorb,eaijk,emp4,emp5,
     $                       ncor,nocc,nvir)
c     
                        if (i.ne.k)then
                           call ccsd_tengy(f1t,f1n,f2t,f2n,
     $                          f3t,f3n,f4t,f4n,
     $                          Dja(1+(k-1)*nvir),Djka(1+(k-klo)*nvir), 
     $                          t1((i-1)*nvir+1),eorb,eaijk,emp4,emp5,
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
            if(ga_nodeid().eq.0) then
               write(6,4321) ' ccsd(t): done ',
     A              a-(ncor+nocc)+((klo-1)/kchunk)*nvir,
     O              ' out of ',(nocc/kchunk)*nvir,
     O              ' progress: ',
     O              ((a-(ncor+nocc)+((klo-1)/kchunk)*nvir)*100d0)/
     D              ((nocc/kchunk)*nvir),
     P              '%'
 4321          format(a,i8,a,i8,a,f6.1,a1)
                  endif
         end do
      end do
c     
c
      next=nxtask(-nodes, 1)
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
c     
      end

      subroutine ccsd_doooo(eintc,eintx,t2c,t2x,f1,f2,f3,f4,nocc,nvir)
      implicit none
      double precision eintc(*),eintx(*),t2c(*),t2x(*),
     &                 f1(*),f2(*),f3(*),f4(*)
      integer nocc,nvir
      double precision one
      data one/1.0d00/
#include "ccsdps.fh"
      if (occsdps) then
         call pstat_on(ps_doooo)
      else
         call qenter('doooo',0)
      endif
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
      if (occsdps) then
         call pstat_off(ps_doooo)
      else
         call qexit('doooo',0)
      endif
      return
      end






      subroutine ccsd_dovvv(fintc,fintx,t2,f1,f2,f3,f4,nocc,nvir)
      implicit none
      double precision fintc(*),fintx(*),t2(*),f1(*),f2(*),f3(*),f4(*)
      integer nocc,nvir
      double precision one,zero
      data one/1.0d00/
      data zero/0.0d00/
#include "ccsdps.fh"
      if (occsdps) then
         call pstat_on(ps_dovvv)
      else
         call qenter('dovvv',0)
      endif
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
      if (occsdps) then
         call pstat_off(ps_dovvv)
      else
         call qexit('dovvv',0)
      endif
      return
      end

      subroutine ccsd_getdjka(djka,
     G     g_objv,nvir,lnov,
     K     j,av,klo,khi)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      double precision djka(*)
      integer g_objv, nvir,lnov
      integer av,klo,khi,j
c
      integer k,l_kk,k_kk,lenk,start,end,iptr,i,startc
      logical failed
c
      start = 1 + (j-1)*nvir + (klo-1)*lnov
      end =  (j-1)*nvir + (khi)*lnov
      lenk=end-start+1
      lenk=min(lenk,(ma_inquire_avail(MT_DBL)*66)/100)
crestrict to 5M to avoid OOMs
      lenk=min(lenk,5000000)
 1    continue
      failed=.not.ma_push_get(MT_Dbl,lenk,
     Q     'tdjka',l_kk,k_kk)
      if(failed) then
         lenk=lenk/2
         goto 1
      endif
c     round it to the nearest multiple of lnov
      lenk=(lenk/lnov)*lnov
      iptr=1
      failed=.false.
 2    continue
      call ccsd_cpyloop(g_objv,start,end,lenk,
     V     av,nvir,lnov, iptr,
     A     dbl_mb(k_kk),djka)
      if(start-1.ne.end)then
         if(failed)
     A     call errquit('getdjka:error ', 0, 0)
         lenk=end-start+1
         goto 2
      endif
      if (.not.ma_pop_stack(l_kk))
     &     call errquit('ccsdgetdjka: cannot pop stack',11,0)

      return
      end
      subroutine ccsd_cpyloop(g_objv,start,end,lenk,
     V     av,nvir,lnov, iptr,
     A     kk,djka)
      implicit none
      integer g_objv,start,end,lenk,
     V     av,nvir,lnov
      double precision kk(*),djka(*)
c
      integer startc,iptr,i,k
c
      do i=1,(end-start+1)/lenk
         call ga_get(g_objv,start,start+lenk-1,av,av,
     $        kk,lenk)
         startc=1
         do k = 1,lenk/lnov
            call dcopy(nvir,kk(startc),1,Djka(iptr),1)
            iptr=iptr+nvir
            startc=startc+lnov
         enddo
         start=start+lenk
      enddo
      return
      end




