      subroutine selci_makeh(h,g,int12,int34,w1,w2,ioconf,indxci,hd,
     $     work1, work2, work3, f)
#include "implicit.fh"
#include "ciinfo.fh"
#include "cbuff.fh"
      dimension h(*),g(*),int12(*),int34(*),
     &     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1),
     &     ioconf(nintpo,noconf),indxci(noconf),hd(nci)
c
      parameter (lenhbf = 20000)
      dimension iocc(255),indbar(255), lists(255), listd(255),
     $     hbuf(lenhbf), jbuf(2,lenhbf)
c
      dimension work2(nfmax*nfmax),
     $     work1(nfmax*nfmax),work3(nfmax*nfmax),
     $     f(norbs*(norbs+1)/2),ind(4)
c
      integer selci_iodiff
c    
c     naively compute the hamiltonian matrix
c
c     For each orbital configuration I as many interactions
c     as possible are accumulated in hbuf. When this is full
c     or we start processing another I value this buffer is flushed
c     into the I/O buffer proper
c
c     no. of elements in hbuf
      ninhb = 0
      ninjb = 0
c     no. of elements in I/O buffer
      nrinb = 0
      niinb = 0
      nrec = 0
      nval = 0
*mdc*if parallel
      me = nodeid()
      nproc = nnodes()
      icount = -1
*mdc*endif
      call dfill(nci, 0.0d0, hd, 1)
c     
c     loop through I occupancies
c     
      do 10 iconf = 1,noconf
c     
c     get required information on I
c     
         call selci_upkcon(norbs, iocc, ioconf(1,iconf), nintpo, nbitpi)
         call selci_mkindb(norbs, iocc, indbar, listd, lists, ns, nd)
         call selci_makef(f, h, g, int12, int34, iocc, listd, lists,
     $        ns, nd)
         ibase = indxci(iconf)
*mdc*if debug
*        write(ifllog,*) ' makeh: configuration i ', iconf
*        call selci_wrtcon(ifllog,iocc,norbs)
*        write(ifllog,*) ' makeh: indbar'
*        write(ifflog,1101) (indbar(iqiq),iqiq=1,norbs)
*1101    format(20i4)
*        write(ifllog,*) ' makeh: f'
*        write(ifllog,1102) (f(iqiq),iqiq=1,norbs)
*1102    format(5d16.8)
*mdc* endif
c     
         do 20 jconf = iconf,noconf
*mdc*if parallel
            icount = icount + 1
            if (mod(icount, nproc) .ne. me) goto 20
*mdc*endif
c     
c     locate interacting J occupancies
c     all this can be vectorised at a later date
c     
            iexcit = selci_iodiff(ioconf(1,iconf), ioconf(1,jconf), 
     $           nintpo)
*mdc*if debug
*            write(ifllog,*) ' makeh: configuration j ',jconf
*            call selci_prtcon(ifllog, norbs, ioconf(1,jconf), 
*     &        nintpo, nbitpi)
*            write(6,*) ' makeh: iexcit = ',iexcit
*mdc*endif
            if (iexcit .gt. 4) goto 20
c     
c     jconf interacts with iconf
c     
            if (iexcit.eq.0) then
*               write(6,*)' icount,who,iconf ',icount,nodeid(),me,
*     &                   nproc,iconf
               call flush(6)
               call selci_makehd(work1,.false.,
     $              ns,nd,indbar,iocc,lists,listd,w1,w2,
     $              work2,work3,f,h,
     $              g,int12,int34,numf)
c    zero the upper half and diagonal of work1 to simplify putinb
               do 33 id = 1,numf
                  idid = (id-1)*numf
                  hd(ibase+id) = work1(idid+id)
                  do 34 jd = id,numf
                     work1(idid+jd) = 0.0d0
 34               continue
 33            continue
               numf2 = numf
            else  if (iexcit.eq.2) then
               call selci_getij(i,j,ioconf(1,iconf),ioconf(1,jconf),
     $              nintpo,nbitpi,iocc)
               call selci_makehs(work1,i,j,ns,indbar,iocc,lists,
     $              w1,w2,work2,work3,f,g,int12,int34,numf,numf2)
            else if(iexcit.eq.4) then
               call selci_gtijkl(ind,ioconf(1,iconf),ioconf(1,jconf),
     $              nintpo,nbitpi,iocc,icase)
               call selci_makhdb(work1,ind,icase,ns,indbar,iocc,
     $              w1,w2,work2,work3,g,int12,int34,numf,numf2)
            else
               call errquit('strange excitation value ',iexcit)
            endif
            if (numf*numf2+ninhb.gt.lenhbf)
     $           call selci_putinb(iflham,indxci(iconf),numf,hbuf,ninhb,
     $           jbuf,ninjb)
            call dcopy(numf*numf2,work1,1,hbuf(ninhb+1),1)
            ninhb = ninhb + numf*numf2
            ninjb = ninjb + 1
            jbuf(1,ninjb) = indxci(jconf)
            jbuf(2,ninjb) = numf2
 20      continue
         call selci_putinb(iflham,indxci(iconf),numf,hbuf,ninhb,jbuf,
     $        ninjb)
 10   continue
c     
      call selci_putinb(iflham,-1,0,hbuf,0,jbuf,0)
c     if running in parallel need to get all the diags together
*mdc*if parallel
c      call dgop (99, hd, nci, '+', work1)
*      write(6,*)' before hd(1)   = ',hd(1),nodeid()
*      write(6,*)' before hd(782) = ',hd(782),nodeid()
      call dgop (99, hd, nci, '+')
*      write(6,*)' after  hd(1)   = ',hd(1),nodeid()
*      write(6,*)' after  hd(782) = ',hd(782),nodeid()
      if (me .eq. 0) then
*mdc*endif
        rewind iflhdg
        call selci_swrite(iflhdg,hd,nci)
        close(iflhdg, status='keep')
*mdc*if parallel
      endif
      write(ifllog, 99) me, nrec, nval
c 99   format(/' process number              ',i3/
c     $        ' no. of hamiltonian records  ',i9/
c     $        ' no. of hamiltonian elements ',i9/)
 99   format(/' n_proc ',i3,
     $        '| n_recs ',i9,
     $        '| n_eles ',i9,'|')
*mdc*else
*        write(ifllog,99) nrec,nval
* 99     format(/' no. of hamiltonian records  ',i9/
*     $          ' no. of hamiltonian elements ',i9/)
*mdc*endif
      close(iflham, status='keep')
c
      end
