      subroutine ccsd_ht2pm(basis,nsh,ncor,nocc,nvir,nact,nbf,g_sht2,
     &                      g_nht2,cmo,scra,scrb,offsh)
C     $Id$
      implicit none
#include "errquit.fh"
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,g_sht2,g_nht2,
     &        offsh(nsh,nsh)
      double precision cmo(nbf,nbf),scra(nbf*nbf),scrb(nbf*nbf)

#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "ccsd_debug.fh"
#include "ccsdps.fh"
c
      integer g_sht2t,g_jlo,g_jhi,g_ilo,g_ihi,
     &        ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,
     &        x,y,xy,ymax,nfi,nfj,ipp,imm,ii,jj,iijj,
     &        jjii,i,j,k,l,ad1,ad2,ad3,iam,lnoo,a,b
      integer nxtask
      external nxtask
c
      if (occsdps) then
         call pstat_on(ps_ht2pm)
      else
         call qenter('ht2pm',0)
      endif
c
      iam=ga_nodeid()
      lnoo=nocc*nocc
c
      call ga_distribution(g_sht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
        nfi=ihi-ilo+1
        do jsh=1,ish
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
          nfj=jhi-jlo+1
          do x=1,nfi
            ymax=nfj
            if (ish.eq.jsh)ymax=x
            do y=1,ymax
              xy=offsh(ish,jsh)+(x-1)*nfj+y
              if (xy.ge.g_ilo.and.xy.le.g_ihi)then
                call ga_get(g_sht2,1,lnoo,xy,xy,scra,lnoo)
                ipp=0
                imm=nocc*(nocc+1)/2
                do ii=1,nocc
                  do jj=1,ii-1
                    ipp=ipp+1
                    imm=imm+1
                    iijj=(ii-1)*nocc+jj
                    jjii=(jj-1)*nocc+ii
                    scrb(iijj)=scra(ipp)+scra(imm)
                    scrb(jjii)=scra(ipp)-scra(imm)
                  enddo
                  iijj=(ii-1)*nocc+ii
                  ipp=ipp+1
                  scrb(iijj)=scra(ipp)
                enddo
                call ga_put(g_sht2,1,lnoo,xy,xy,scrb,lnoo)
                if (ish.ne.jsh.or.x.ne.y)then
                  xy=offsh(jsh,ish)+(y-1)*nfi+x
                  do ii=1,nocc
                    do jj=1,nocc
                     iijj=(ii-1)*nocc+jj
                     jjii=(jj-1)*nocc+ii
                     scra(iijj)=scrb(jjii)
                    enddo
                  enddo
                  call ga_put(g_sht2,1,lnoo,xy,xy,scra,lnoo)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
      call ga_sync()
c
c ------------------------------------------------------------
c - transform ao indices of ht2 array into the mo basis
c ------------------------------------------------------------

*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,lnoo,'sht2t',
     &                   nbf*nbf,0,g_sht2t))
     &     call errquit('ga_create g_sht2t failed',0, GA_ERR)

      call ga_transpose(g_sht2,g_sht2t)
c - redefine g_sht2
      if (.not.ga_destroy(g_sht2))
     &    call errquit('ga_dest g_sht2 fail',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nact*nact,lnoo,'sht2',
     &                   nact*nact,0,g_sht2))
     &     call errquit('ga_create g_sht2 failed',0, GA_ERR)

      call ga_distribution(g_nht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            ad3=(i-1)*nocc+j
              call ga_get(g_sht2t,1,nbf*nbf,ad3,ad3,
     &                    scrb,nbf*nbf)
              ad3=0
              do ksh=1,nsh
                if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     &              call errquit('vvvv: bas_cn2bfr',ksh, BASIS_ERR)
                do lsh=1,nsh
                  if (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     &                call errquit('vvvv: bas_cn2bfr',lsh, BASIS_ERR)
                  do k=klo,khi
                    do l=llo,lhi
                      ad3=ad3+1
                      scra((k-1)*nbf+l)=scrb(ad3)
                    enddo
                  enddo
                enddo
              enddo
              call dgemm('n','n',nbf,nact,nbf,1.0d00,scra,nbf,
     &                   cmo(1,ncor+1),nbf,0.0d00,scrb,nbf)
              call dgemm('t','n',nact,nact,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &                   scrb,nbf,0.0d00,scra,nact)
              ad3=(i-1)*nocc+j
              call ga_put(g_sht2,1,nact*nact,ad3,ad3,
     &                    scra,nact*nact)
c
              if (dob(1).eq.2)then
                ad3=0
                do a=nocc+1,nact
                  do b=nocc+1,nact
                    ad3=ad3+1
                    scrb(ad3)=scra((a-1)*nact+b)
                  enddo
                enddo
                call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                      scrb,nvir,1.0d00)
              endif
c
            endif
          enddo
        endif
      enddo
      if (.not.ga_destroy(g_sht2t))
     &    call errquit('ga_dest g_sht2t fail',0, GA_ERR)
c
      if (occsdps) then
         call pstat_off(ps_ht2pm)
      else
         call qexit('ht2pm',0)
      endif
c
      return
      end
      subroutine ccsd_t2pm(basis,nsh,ncor,nocc,nvir,nbf,g_st2,g_nt2,cmo,
     &                     t1,scra,scrb,nbfdim)
      implicit none
#include "errquit.fh"
      integer basis,nsh,ncor,nocc,nvir,nbf,g_st2,g_nt2,nbfdim
      double precision cmo(nbf,nbf),scra(nbf*nbf),scrb(nbf*nbf),
     &                 t1(nocc,nvir)
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "ccsdps.fh"
c
      integer g_st2t,lnoo,iam
      integer i,j,k,l,a,b,ad1,ad2,ad3,g_jlo,g_jhi,g_ilo,g_ihi,
     &        ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,
     &        xy,x,y,ipp,imm,ii,jj,iijj,jjii,ad4,nodes
c
      if (occsdps) then
         call pstat_on(ps_t2pm)
      else
         call qenter('t2pm',0)
      endif
c
      iam=ga_nodeid()
      nodes=ga_nnodes()
      lnoo=nocc*nocc
c
*ga:1:0
      do ad4=iam+1,lnoo,nodes
         i=(ad4-1)/nocc+1
         j=ad4-(i-1)*nocc
              ad1=(i-1)*nvir
              ad2=(j-1)*nvir
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              ad3=0
              do a=1,nvir
                do b=1,nvir
                  ad3=ad3+1
                  scra(ad3)=scra(ad3)+t1(i,a)*t1(j,b)
                enddo
              enddo
              call dgemm('n','t',nvir,nbf,nvir,1.0d00,scra,nvir,
     &                   cmo(1,ncor+nocc+1),nbf,0.0d00,scrb,nvir)
              call dgemm('n','n',nbf,nbf,nvir,1.0d00,cmo(1,ncor+nocc+1),
     &                   nbf,scrb,nvir,0.0d00,scra,nbf)
              ad3=0
              do ksh=1,nsh
                if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     &              call errquit('vvvv: bas_cn2bfr',ksh, BASIS_ERR)
                do lsh=1,ksh
                  if (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     &                call errquit('vvvv: bas_cn2bfr',lsh, BASIS_ERR)
                  do k=klo,khi
                    do l=llo,lhi
                      ad3=ad3+1
                      scrb(ad3)=scra((k-1)*nbf+l)
                    enddo
                  enddo
                enddo
              enddo
              ad3=(i-1)*nocc+j
              call ga_put(g_st2,ad3,ad3,1,nbfdim,
     &                    scrb,1)
      enddo
      call ga_sync()
c
c ------------------------------------------------------------
c - form t2+/- (see gustavos paper)
c - t2+ = t_ij^ab + t_ij^ba 
c - t2- = t_ij^ab - t_ij^ba 
c - some prefactors absorbed
c ------------------------------------------------------------
      call ga_distribution(g_st2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      xy=0
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
        do jsh=1,ish
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
          do x=ilo,ihi
            do y=jlo,jhi
              xy=xy+1
              if (xy.ge.g_ilo.and.xy.le.g_ihi)then
                call ga_get(g_st2,1,lnoo,xy,xy,scra,lnoo)
                ipp=0
                imm=nocc*(nocc+1)/2
                do ii=1,nocc
                  do jj=1,ii-1
                    iijj=(ii-1)*nocc+jj
                    jjii=(jj-1)*nocc+ii
                    ipp=ipp+1
                    imm=imm+1
                    scrb(ipp)=(scra(iijj)+scra(jjii))
                    scrb(imm)=(scra(iijj)-scra(jjii))
                  enddo
                  iijj=(ii-1)*nocc+ii
                  ipp=ipp+1
                  scrb(ipp)=scra(iijj)+scra(iijj)
                enddo
                if (x.eq.y)then
                  call dscal(lnoo,0.25d00,scrb,1)
                else
                  call dscal(lnoo,0.5d00,scrb,1)
                endif
                call ga_put(g_st2,1,lnoo,xy,xy,scrb,lnoo)
              endif
            enddo
          enddo
        enddo
      enddo
      call ga_sync()
c
      if (occsdps) then
         call pstat_off(ps_t2pm)
      else
         call qexit('t2pm',0)
      endif
c
      return
      end
      subroutine ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,
     &                    g_sht2,g_c,g_x,offsh,snsi,sisn,lssni,scre,
     &                    mem2,max2e,eri1,eri2,t1,cmo,t1ao,scra,
     &                    scrb,lscr,kscr,tol2e,iprt,tklst,maxints,
     &                    ish_idx,shinf,max_sht2_blk,nbfdim,a_st2,
     &                    use_ccsd_omp)
      implicit none
#include "errquit.fh"
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,g_sht2,
     &        g_c,g_x,lssni,
     &        mem2,max2e,lscr,kscr,iprt,max_sht2_blk,nbfdim
      integer offsh(nsh,nsh,2),tklst(nsh*(nsh+1)/2,2),maxints
      integer ish_idx(nsh*nsh+2),shinf(nsh,3)
      double precision tol2e,eri1(maxints),eri2(maxints),scre(mem2),
     &                 t1(nocc*nvir),cmo(nbf,nbf),scra(kscr),
     &                 scrb(lscr),snsi(lssni),sisn(lssni),
     &                 t1ao(nocc*nbf),tx(2),a_st2(nocc*nocc,nbfdim)
      logical schwarz1,schwarz2

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "eaf.fh"
*rak-s
#include "ccsd_time.fh"
*rak-e
#include "ccsdps.fh"
#include "ccsd_data.fh"
      integer ad1,ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,next,
     &        icnt,nfi,nfk,lnijkl,offjl,offlj,lnoo,nfj,nfl,nodes,iam
      logical flush, store, use_storage
      integer nxtask,blocksize,logsize,intsize,ijk(2),cnijkl
      integer icount(2),jcnt,btki,kcount,jcount,ishb,iblk,bnijkl,ival
      integer cart_2e4c,ad2
      double precision store_need,store_used
      double precision faddr,kaddr
      external nxtask, cart_2e4c
      integer trishel,iamoff,iii,jjj,istart
      logical, optional, intent(in) :: use_ccsd_omp
      if (.not.present(use_ccsd_omp)) then
          call errquit('ccsd_pzamp: use_ccsd_omp not present!',0,0)
      endif
c
      if (occsdps) then
         call pstat_on(ps_sxy)
      else
         call qenter('sxy',0)
      endif
c
      lnoo=nocc*nocc
      iam=ga_nodeid()
      nodes=ga_nnodes()
c
c ------------------------------------------------------------
c - work out an index array for shell offset
c ------------------------------------------------------------
      call ga_sync()
      ad1=0
      ad2=0
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
        nfi=ihi-ilo+1
        shinf(ish,1)=nfi
        shinf(ish,2)=ilo
        shinf(ish,3)=ihi
        do jsh=1,nsh
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
          nfj=jhi-jlo+1
          offsh(ish,jsh,1)=ad1
          ad1=ad1+nfi*nfj
          if (jsh.le.ish) then
             offsh(ish,jsh,2)=ad2
             ad2=ad2+nfi*nfj
          endif
        enddo
      enddo
c ------------------------------------------------------------
c create t1ao
c ------------------------------------------------------------
      call dgemm('n','t',nbf,nocc,nvir,1.0d00,cmo(1,ncor+nocc+1),
     &           nbf,t1,nocc,0.0d00,t1ao,nbf)
c
c ------------------------------------------------------------
c - loop over the integral generation this is what needs to be
c - efficient, so gather timings for this
c - at the moment we are shell blocked, we could go to atom
c - blocking if this is inefficient
c - note integrals computed 4 times minimal list
c ------------------------------------------------------------
c
c      call ga_sync()
      tx(1)=tcgtime()
      if (iam.eq.0.and.iprt.gt.5)
     &    print *,' begin parallel integral generation'
      call ga_zero(g_c)
c      call ga_sync()
      faddr=0.d0
      kaddr=0.d0
      logsize=ma_sizeof(MT_LOG,2,MT_BYTE)
      intsize=ma_sizeof(MT_INT,1,MT_BYTE)
c
c If this is the second call to ccsd_sxy we may have the integrals and other
c info on disk. We should read from there....
c
      if (repeat.and.use_disk) then
c
c Get the number of lsh,lsh blocks this node is doing
c
      if (eaf_read(sxy_hl,faddr,icount,intsize*2).ne.0)
     &   call errquit('ccsd_sxy: read failed',1,DISK_ERR)
      faddr=faddr+intsize*2
c
c Start the main loop, using the integrals that are stored on disk
c
      do jcnt=1,icount(1)
c
c Get the jsh and lsh shell info
c
         if (eaf_read(sxy_hl,faddr,ijk,intsize*2).ne.0)
     &      call errquit('ccsd_sxy: read failed',2,DISK_ERR)
         faddr=faddr+intsize*2
c
         jsh=(tklst(ijk(1)+1,2)-1)/nsh+1
         lsh=tklst(ijk(1)+1,2)-(jsh-1)*nsh
c
         nfj=shinf(jsh,1)
         jlo=shinf(jsh,2)
         jhi=shinf(jsh,3)
         nfl=shinf(lsh,1)
         llo=shinf(lsh,2)
         lhi=shinf(lsh,3)
c
         call dcopy(lnoo*nfj*nfl,0.0d00,0,scrb,1)
         call dcopy(nfj*nfl*nbf*nocc,0.0d00,0,snsi,1)
         call dcopy(nfj*nfl*nbf*nocc,0.0d00,0,sisn,1)
c
c ijk(2) contains number of ksh*(ish-block) blocks stored
c
         do iblk=1,ijk(2)
c
c Get data for ksh and ish-block. First get number of blocks
c for the ksh to be calculated, then get whole block of appropriate
c length
c
            call qenter('r_read',0)
            if (eaf_read(sxy_hl,faddr,ival,intsize).ne.0)
     &         call errquit('ccsd_sxy: read failed',
     &         3,DISK_ERR)
            blocksize=intsize*ival
            if (eaf_read(sxy_hl,faddr,ish_idx,blocksize).ne.0)
     &         call errquit('ccsd_sxy: read failed',
     &         4,DISK_ERR)
            faddr=faddr+blocksize
c
c Get integrals, number of integrals in ish_idx(1,2)
c
            blocksize=ma_sizeof(MT_DBL,ish_idx(2),MT_BYTE)
            if (eaf_read(sxy_hl,faddr,eri1,blocksize).ne.0)
     &         call errquit('ccsd_sxy: read failed',
     &         5,DISK_ERR)
            faddr=faddr+blocksize
            if (eaf_read(sxy_hl,faddr,eri2,blocksize).ne.0)
     &         call errquit('ccsd_sxy: read failed',
     &         6,DISK_ERR)
            faddr=faddr+blocksize
            call qexit('r_read',0)
            call t2eri(ish_idx(1),ish_idx,jlo,jhi,nfj,llo,lhi,nfl,nsh,
     &                 eri1,eri2,scra,scrb,lnoo,nocc,offsh,nbf,g_st2,
     &                 shinf,max_sht2_blk,snsi,sisn,t1ao,nbfdim,a_st2,
     &                 use_ccsd_omp)
         enddo
         ad1=offsh(jsh,lsh,1)
         call ga_put(g_sht2,1,lnoo,ad1+1,ad1+nfj*nfl,scrb,lnoo)
         offjl=offsh(jsh,lsh,1)
         offlj=offsh(lsh,jsh,1)
         if (use_ccsd_omp) then
           call ccsd_idx2_omp(snsi,sisn,cmo,lscr,
     &                        nfj,nfl,ncor,nocc,nact,nbf,
     &                        jlo,jhi,llo,lhi,offjl,offlj,g_x,g_c)
         else
           call ccsd_idx2(snsi,sisn,cmo,scra,scrb,lscr,
     &                    nfj,nfl,ncor,nocc,nact,nbf,
     &                    jlo,jhi,llo,lhi,offjl,offlj,g_x,g_c)
         endif
      enddo
c
c     End reading integrals from disk in repeat iteration
c     
      endif!!!!if (repeat.and.use_disk) then
c
      if (.not.repeat.or..not.use_disk.or.icount(2).gt.0) then
c
c We got here in one of two cases:
c
c 1. First time through, use dynamical (nxtask) distribution to 
c    calculate and store the integrals....
c    - Start from iam
c    - use_storage = true
c    - reset icount array
c 2. Repeat time through, needing to recalculate some integrals
c    as they could not be stored to disk...
c    If no disk was used, case is like 1, but with use_storage = false
c    - Start from icount(2)
c    - use_storage = false
c
      if (repeat.and.use_disk) then
         istart = icount(2)
         use_storage=.false.
      else
         istart = iam
         icount(1)=0
         icount(2)=0
         if (use_disk) then
            if (eaf_write(sxy_hl,faddr,icount,intsize*2).ne.0)
     &         call errquit('ccsd_sxy: write failed',1,DISK_ERR)
            faddr=faddr+intsize*2
         endif
         use_storage=use_disk
         store_used=0.d0
      endif
      do icnt=istart,nsh*(nsh+1)/2-1,nodes
c
c ijk contains jsh and lsh info plus number of ksh and ish blocks
c
            ijk(1)=icnt
            jsh=(tklst(icnt+1,2)-1)/nsh+1
            lsh=tklst(icnt+1,2)-(jsh-1)*nsh
c
            nfj=shinf(jsh,1)
            jlo=shinf(jsh,2)
            jhi=shinf(jsh,3)
            nfl=shinf(lsh,1)
            llo=shinf(lsh,2)
            lhi=shinf(lsh,3)
c
            call dcopy(lnoo*nfj*nfl,0.0d00,0,scrb,1)
            call dcopy(nfj*nfl*nbf*nocc,0.0d00,0,snsi,1)
            call dcopy(nfj*nfl*nbf*nocc,0.0d00,0,sisn,1)
c
            bnijkl=1
            kcount=2
            flush=.false.
            jcount=0
            trishel=nsh*(nsh+1)/2
            iamoff=iam*(trishel/nodes)
c
c Determine if we can still store the integrals to disk
c If not, calculate and set integral recalculation point icnt
c
            if (use_storage) then
              store_need=dble(ma_sizeof(MT_DBL,2*nfj*nfl*nbf*(nbf+1),
     $                                    MT_BYTE))
              store_need=store_need+dble(ma_sizeof(MT_INT,nsh*nsh+2,
     $                                               MT_BYTE))
              if ((store_need+store_used).gt.store_avail) then
                use_storage=.false.
                icount(2)=icnt 
              else
                icount(1)=icount(1)+1
                kaddr=faddr
                if (eaf_write(sxy_hl,faddr,ijk,intsize*2).ne.0) call
     &                 errquit('ccsd_sxy: write failed',2,DISK_ERR)
                faddr=faddr+intsize*2
              endif
            endif
c
            do iii=iamoff,trishel+iamoff-1 !begin loop
               jjj=iii
               if(iii.gt.trishel-1) jjj=iii-trishel
               ish=(tklst(jjj+1,2)-1)/nsh+1
               ksh=tklst(jjj+1,2)-(ish-1)*nsh
               nfk=shinf(ksh,1)
                schwarz1=schwarz_shell(ish,jsh)*
     &                   schwarz_shell(ksh,lsh).ge.tol2e
                schwarz2=schwarz_shell(ish,lsh)*
     &                   schwarz_shell(ksh,jsh).ge.tol2e
                nfi=shinf(ish,1)
                lnijkl=nfi*nfj*nfk*nfl
                cnijkl=cart_2e4c(basis,ish,jsh,ksh,lsh)
c
c Check if the integral buffers are full or we have done all ish and ksh
c If so store and process integrals, else add another block
c               
                store=(bnijkl+cnijkl.gt.maxints)
  111           continue
                if ((store.and.(schwarz1.or.schwarz2)).or.flush) then
c
c Store the integrals and process them in t2eri and ccsd_idx1
c
c First four spots in ish_idx are used to store additional data:
c  (1)=kcount  -> number of ish+ksh blocks plus 2 for this data (2*2)
c  (2)=bnijkl  -> number of integrals
c    
                   ish_idx(1)=kcount
                   ish_idx(2)=bnijkl
                   if (use_storage) then
                     jcount=jcount+1
                     blocksize=intsize*kcount
                     call qenter('f_write',0)
                     if (eaf_write(sxy_hl,faddr,ish_idx,blocksize).ne.0)
     &                  call errquit('ccsd_sxy: write failed',
     &                  3,DISK_ERR)
                     faddr=faddr+blocksize
                     store_used=store_used+blocksize
                     blocksize=ma_sizeof(MT_DBL,bnijkl,MT_BYTE)
                     if (eaf_write(sxy_hl,faddr,eri1,blocksize).ne.0) 
     &                  call errquit('ccsd_sxy: write failed',
     &                  4,DISK_ERR)
                     faddr=faddr+blocksize
                     if (eaf_write(sxy_hl,faddr,eri2,blocksize).ne.0) 
     &                  call errquit('ccsd_sxy: write failed',
     &                  5,DISK_ERR)
                     faddr=faddr+blocksize
                     call qexit('f_write',0)
                     store_used=store_used+2*blocksize
                   endif
c
c Process the integral block
c
                   call t2eri(kcount,ish_idx,jlo,jhi,nfj,llo,lhi,nfl,
     &                        nsh,eri1,
     &                        eri2,scra,scrb,lnoo,nocc,offsh,nbf,g_st2,
     &                        shinf,max_sht2_blk,snsi,sisn,t1ao,
     &                        nbfdim,a_st2,
     &                        use_ccsd_omp)
c     
c Reset some indices
c
                   bnijkl=1
                   kcount=2
                endif
c
c Add next block of integrals
c               
                if (.not.flush.and.(schwarz1.or.schwarz2)) then
                   call qenter('ints',0)
                   kcount=kcount+1
                   ish_idx(kcount)=(ksh-1)*nsh+ish
                   if (schwarz1) then
                      if(bnijkl+max2e.gt.maxints) then
                         write(6,'(I6,A,I20,A)') ga_nodeid(),
     W                        ' : increase stack memory ',
     W                     (trishel+iamoff-iii+1)*2*max2e*8/1024/1024,
     M                        ' MBytes '
                         call util_flush(6)
                         call errquit(' ccsdsxy: not enough MA',0,0)
                      endif
                     call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                             mem2, scre, max2e, eri1(bnijkl) )
                   else
                     call dcopy(lnijkl,0.0d00,0,eri1(bnijkl),1)
                   endif
                   if (schwarz2) then
                     if (jsh.eq.lsh) then
                        call dcopy(lnijkl,eri1(bnijkl),1,eri2(bnijkl),1)
                     else
                      if(bnijkl+max2e.gt.maxints) then
                         write(6,*) ' increase stack memory'
                         call util_flush(6)
                         call errquit(' ccsdsxy: not enough MA',0,0)
                      endif
                        call int_2e4c(basis,ish,lsh,basis,ksh,jsh,
     $                                mem2,scre,max2e,eri2(bnijkl))
                     endif
                   else
                     call dcopy(lnijkl,0.0d00,0,eri2(bnijkl),1)
                   endif
                   bnijkl=bnijkl+lnijkl
                   call qexit('ints',0)
                endif
                if (iii.eq.(trishel+iamoff-1).and..not.flush) then
                   flush=.true.
                   goto 111
                endif
                flush=.false.
            enddo!end loop
c
            if (use_storage) then
              ijk(2)=jcount
              call qenter('f_write',0)
              if (eaf_write(sxy_hl,kaddr,ijk,intsize*2).ne.0)
     &           call errquit('ccsd_sxy: write failed',6,DISK_ERR)
              call qexit('f_write',0)
            endif
c
            ad1=offsh(jsh,lsh,1)
            call ga_put(g_sht2,1,lnoo,ad1+1,ad1+nfj*nfl,scrb,lnoo)
            offjl=offsh(jsh,lsh,1)
            offlj=offsh(lsh,jsh,1)
            if (use_ccsd_omp) then
              call ccsd_idx2_omp(snsi,sisn,cmo,lscr,
     &                           nfj,nfl,ncor,nocc,nact,nbf,
     &                           jlo,jhi,llo,lhi,offjl,offlj,g_x,g_c)
            else
              call ccsd_idx2(snsi,sisn,cmo,scra,scrb,lscr,
     &                       nfj,nfl,ncor,nocc,nact,nbf,
     &                       jlo,jhi,llo,lhi,offjl,offlj,g_x,g_c)
            endif
      enddo
      faddr=0.d0
      if (.not.repeat.and.use_disk) then
         if (eaf_write(sxy_hl,faddr,icount,intsize*2).ne.0) call
     &       errquit('ccsd_sxy: write failed',7,DISK_ERR)
      endif
      repeat=.true.
      endif!!!if (.not.repeat.or..not.use_disk.or.icount(2).gt.0)
cedo
      call ga_sync()
c sync before 3 and 4 index transformation
      if (use_ccsd_omp) then
        call ccsd_idx34_omp(basis,cmo,
     &                      nsh,ncor,nocc,nact,nbf,
     &                      g_x,g_c)
      else
        call ccsd_idx34(basis,cmo,scra,scrb,
     &                  nsh,ncor,nocc,nact,nbf,
     &                  g_x,g_c)
      endif
c
      call ga_sync()
      tx(2)=tcgtime()
      if (iam.eq.0) then
*rak     write(6,*)'Time around main block',tx(2)-tx(1)
         main_block_time = tx(2)-tx(1)
      endif
c
      if (occsdps) then
         call pstat_off(ps_sxy)
      else
         call qexit('sxy',0)
      endif
c
      return
      end
