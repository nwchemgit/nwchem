      Subroutine dft_3cincor( n_batch, n_bmax, iwhat, i3c_ERI, 
     &                       n3c_dbl, iwhat_max, n_semi_bufs, fd)
c 
C$Id: dft_3cincor.F,v 1.20 2000-05-03 22:45:08 edo Exp $
c
      implicit none
c    
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches [in]
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl ! max no. ints [in]
      double precision i3c_ERI(n3c_dbl) ! int values [out]
c
#include "bas.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "schwarz.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "basP.fh"
#include "eaf.fh"
#include "util.fh"
#include "sym.fh"
c     
c     Compute the 3-center integrals and store as many of them 
c     as possible in-core.
c
      integer me, nproc, ishp, icount, Nao2, maxg, 
     &        mscratch_2e3c, nscr, n1, n2, ishc, ifirstc,
     &        ilastc, nshbfc, ishd, ifirstd, ilastd,
     &        nshbfd, ifirstp, ilastp, nintegrals,ninteg_api,
     &        nshbf_cd, ndone
      integer lscr, iscr, n_semi_bufs
      integer nlast, ntasks, nxtbuf, nbufs, n_chk
      integer iwhat_max, nao2_max, batch_max, neri_max
c
      integer availkb, current_nbuf, n3c_b
      character*20 fstype
      character*255 filename, actualname
      integer fd, other_stuff, n_batch_buf,errnob
      logical oprint_semi_file, file_opened
      logical oprint_3c2e
      double precision offset
c
      integer itype,nprimo,isphere,nbf_cd_cart
      double precision THRESHOLD, ERI_est, ratio
      double precision q1
c
      THRESHOLD=10.d0**(-itol2e)
      oprint_semi_file = util_print('semi-direct info', print_high)
      oprint_3c2e = util_print('3c 2e integrals', print_high)
      me = ga_nodeid()
      call dzero(i3c_eri, n3c_dbl)
      call int_mem_2e3c(maxg, mscratch_2e3c)
c
c     Check whether of not to enable semi-direct; check direct & noio
c
      if(noio.eq.1)n_semi_bufs = 0
c
c     temp set n_semi_bufs = 10
c
c      n_semi_bufs = 10
c      write(6,*)' n_semi_bufs = ',n_semi_bufs,me
c
      n3c_b = 8*n3c_dbl
      other_stuff = 8*(2*nbf_ao**2 + 2*nbf_cd**2 + 2*nbf_xc**2
     &                 + 1024)
      current_nbuf = 0
      if (n_semi_bufs.gt.0)then
c
c        inquire of available disk space
c
         call util_file_name('3ceri', .true., .false., filename)
         call util_pname(filename, actualname)
c
         if (eaf_stat(actualname, availkb, fstype) .ne. 0) 
     &      call errquit('dft_3cincor: unable to get disk info',0)
c
c        need global space available check here
c
         call ga_igop(msg_max_3c_incore, availkb, 1, 'min')
c
         if (availkb .gt. (2*n3c_b + other_stuff)/1024)then
            if (eaf_open(actualname, eaf_rw, fd) .ne. 0)
     &        call errquit('dft_3cincor; eaf_open failed', 0)
            file_opened = .true.
            if (me.eq.0.and.oprint_semi_file)then
               write(6,*) ' dft_3cincor: opening 3c-2e integral file: ',
     &         actualname(1:inp_strlen(actualname))
               write(6,*)' 1 availkb: ',availkb,' fstype: ',fstype,
     &         ' fd: ',fd
               call util_flush(6)
            endif
         endif
      endif
c      
      me=ga_nodeid()
      nproc  = ga_nnodes()
c      
c     Parallelize the loop over the products of AO shells.
c      
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
c      
c     check for zero ... must be at least 1
c
      icount = max(icount,1)
      if(nproc.eq.1)icount = 1
c      
      nao2_max=0
      do ishc=1,nshells_ao
         if(.not.bas_continfo(ao_bas_han,ishc,
     &        itype,nprimo,nshbfc,isphere))
     &        call errquit('Exiting in fitcd.',44)
         nao2_max=max(nao2_max,
     ,        ((itype+1)*(itype+2))/2*nshbfc)
      enddo
      nbf_cd_cart=0
      do ishc=1,nshells_cd
         if(.not.bas_continfo(cd_bas_han,ishc,
     &        itype,nprimo,nshbfc,isphere))
     &        call errquit('Exiting in fitcd.',44)
         nbf_cd_cart=nbf_cd_cart+
     ,        (itype+1)*(itype+2)/2*nshbfc
      enddo
      batch_max = nao2_max*nbf_cd_cart*icount
c
      nERI_max = Nao2_max*nbf_cd_mxnbf_cn
c
c
c     Determine/allocate scratch memory for 3 center 2e- integrals.
c
      nscr    = mscratch_2e3c
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_3cincor: cannot allocate scr',0)
c
      call ga_sync
c
      nxtbuf = batch_max
      nbufs = icount
      call nxtask_nleft (nproc, icount, n1, nbufs)
      n2 = 0
c      
      ndone=1
      n_batch=0
      n_batch_buf = 0
      do 240 ishc = 1,nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &     call errquit('exiting in fitvcoul.',3)
        nshbfc=ilastc-ifirstc+1
c 
        do 230 ishd = 1,ishc
c            
c         screen integral blocks using schwarz ineq:  (p|cd) .le. (cd|cd)
c         
          eri_est=schwarz_shell(ishc,ishd)
c
          if( eri_est.gt.threshold )then
c
            if( n1.eq.n2 )then
c              
c             Define the scatter index.
c              
              if( .not. bas_cn2bfr( ao_bas_han,ishd,
     &          ifirstd,ilastd))call errquit('exiting in fitvcoul.',4)
c
              nshbfd=ilastd-ifirstd+1
              nao2 = nshbfc*nshbfd
c                           
              do 220 ishp = 1,nshells_cd
c
                if ( .not. bas_cn2bfr( cd_bas_han,ishp,ifirstp,ilastp))
     &              call errquit('exiting in fitvcoul.',5)

                nshbf_cd=ilastp-ifirstp+1

                nintegrals=nshbf_cd*nao2
                Ninteg_api = maxg

c               compute the 3 center 2e- integrals between this given pair 
c               of ao shells and cd shell.
c
c                n_chk = ndone+nintegrals-1
                n_chk = ndone+maxg-1
                if( (n_chk) .gt. n3c_dbl) then
                  write(LuOut,*) ' ishc ishd ishp ',ishc,ishd,Ishp
                  write(LuOut,*) ' nshao nshcd ',
     ,                 nshells_ao,nshells_cd
                  write(LuOut,*) ' Nintegrals ',nintegrals
                  write(Luout,*) ' batchmx maxg',batch_max,maxg
                  write(LuOut,*) ' nlast n3cdbl ',n_chk, n3c_dbl
                  call errquit('dft_3cincor: blew out buf of dbles',999)
                endif
c
c               check for use of symmetry
c
                if (oskel) then
                   if (sym_shell(cd_bas_han, ishp, q1))then
                      call int_2e3c(cd_bas_han,ishp,ao_bas_han,ishc,
     &                ishd,nscr,dbl_mb(iscr),ninteg_api,i3c_ERI(ndone))
                      call dscal(nintegrals, q1, i3c_ERI(ndone), 1)
                    else
                      call dzero(i3c_ERI(ndone),Nintegrals)
                   endif
                else
                   call int_2e3c(cd_bas_han,ishp,ao_bas_han,ishc,
     &             ishd,nscr,dbl_mb(iscr),ninteg_api,i3c_ERI(ndone))
                endif                   

                ndone=ndone+nintegrals
c
 220          continue
c
              n_batch = n_batch + 1
              n_batch_buf = n_batch_buf + 1
              if( n_batch.gt.n_bmax) then
                 call errquit('dft_3cincor: blew out buf of ints',
     .                n_bmax)
              else
                 iwhat(n_batch)=(ishc*(ishc-1))/2 + ishd-1
c
c            
c             Check to see if there is room to store another batch of
c             integrals in the local buffer.  If so, do next task,
c             otherwise go wait for global sync.
c
              if(nbufs.eq.0)nxtbuf = ndone + batch_max
              endif
c
c              write(LuOut,1111)me,ishc,ishd,n1,n2,nxtbuf,n3c_dbl, 
c     &                         n_batch, iwhat(n_batch), ndone
c 1111         format(' Done: me,ishc,ishd,n1,n2,nxtbuf,n3c_dbl, ',
c     &                     ' n_batch, iwhat(n_batch), ndone ',/,10i7)
c              call util_flush(LuOut)
c
              if( nxtbuf.lt.n3c_dbl) then
                call nxtask_nleft (nproc, icount, n1, nbufs)
c                n1 = nxtask(nproc,icount)
              else
c
c               The most simplest of semi-direct implementations.
c               When this semi-incore memory buffer is full, dump it
c               and start over.
c
c               if I/O allowed, sync (so that all procs have similar size buffers
c               to process later), write out buffer and continue.
c
                if(n_semi_bufs.eq.0)go to 241
c
c               write buffer of 3c2e ints
c
                current_nbuf = current_nbuf + 1
                offset = (current_nbuf-1)*n3c_b
c
c               put number of batches as last element in buffer
c
                i3c_eri(n3c_dbl) = dble(n_batch_buf)
c  
                errnob=eaf_write(fd, offset, i3c_eri, n3c_b)
                if (errnob .ne. 0)
     &             call errquit(
     .               'dft_3cincor: 3c-2e integral buffer write failed',
     .               errnob)
c                
c               now clear buffer and continue
c
                call dzero(i3c_eri, n3c_dbl)
                nxtbuf = batch_max
                ndone = 1
                n_batch_buf = 0
c
c               check to see if there is still more disk space available:
c
                if (eaf_stat(actualname, availkb, fstype) .ne. 0) 
     &            call errquit('dft_3cincor: unable to get disk info',0)
c
c               need global space available check here
c
c                call ga_igop(msg_max_3c_incore, availkb, 1, 'min')
c
                if (availkb*1024 .gt. (2*n3c_b + other_stuff))then
                   if (me.eq.0.and.oprint_semi_file)then
                      write(6,*)' 2 availkb: ',availkb,' fstype: ',
     &                fstype,' fd: ',fd
                      call util_flush(6)
                   endif
c
c                  get next task
c
                   call nxtask_nleft (nproc, icount, n1, nbufs)
                else
                   go to 242
                endif
c
              endif
c
            end if
c            
            n2 = n2 + 1
c            
          end if
c
  230   continue
c
  240 continue
c
  241 continue
c
c     Either all integrals have been computed and are incore or all
c     local buffers are full. 
c
      if (n_semi_bufs.gt.0)then
c
c        dump last buffer     
c
         current_nbuf = current_nbuf + 1
         offset = (current_nbuf-1)*n3c_b
c
c        put number of batches as last element in buffer
c
         i3c_eri(n3c_dbl) = dble(n_batch_buf)
c  
         if (eaf_write(fd, offset, i3c_eri, n3c_b) .ne. 0)call errquit
     &      ('dft_3cincor: 3c-2e integral buffer write failed',0)
c                
c        now clear buffer and continue
c
         call dzero(i3c_eri, n3c_dbl)
      endif
c
  242 continue

      call nxtask_nleft (-nproc, icount, n1, nbufs)
c
c     Determine last shell-pair computed.
c
      iwhat_max = iwhat(n_batch)
      call ga_igop(msg_max_3c_incore, iwhat_max, 1, 'max')
c
c     Determine total number of shell pairs.
c
      ntasks = (nshells_ao*(nshells_ao+1))/2
c
      nlast = iwhat_max + 1

      ratio = dble(nlast)/dble(ntasks)
c
c     Find last ishc and ishd; add 1 to iwhat_max since starting at 0.
c
      do n1 = 1, nshells_ao
        if ((n1*(n1-1))/2.gt.nlast)then
          ishc = n1 - 1
          goto 243
        endif
      enddo
      ishc = nshells_ao
  243 continue
      ishd = nlast - (ishc*(ishc-1))/2
c
      if (oprint_3c2e)then
         if (me.eq.0)then
            if (n_semi_bufs .eq. 0)then
               write(LuOut,3241)ratio*100.d0
            else
               write(LuOut,3242)ratio*100.d0, current_nbuf
            endif
         endif
      endif
c
c     set and store number of 3c2e integral buffers
c
      if (noio .eq. 1)then
         n_semi_bufs = 0
      else
         n_semi_bufs = current_nbuf
      endif
c
c      write(LuOut,*)' ntasks, nlast, ishc, ishd: ', 
c     &                ntasks, nlast, ishc, ishd
c      call util_flush(LuOut)
c
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_3cincor: cannot pop stack',0)
c
      call nxtask_nleft (-nproc, icount, n1, nbufs)
 3241 format(10x,'The percent of 3c 2e- integrals held in-core is: ',
     &       f6.2,/)
 3242 format(10x,'The percent of 3c 2e- integrals held on-disk is: ',
     &       f6.2,/,
     &       10x,' in ',i6,' integral buffers.',/)
      return
      end
      subroutine nxtask_nleft (nproc, ichunk, nxtask_n, nstill)
      implicit none
#include "global.fh"      
c
      integer nproc, ichunk     ! [input]
c
      integer icount, nleft, nxtval, junk, nxtask_n, nstill
      external nxtval
      save icount, nleft
      data nleft, icount /0, 0/
c     
c     This is a copy of nxtask (from ddscf) but returning
c     copy of nleft as argument.
c     
      if(nproc.gt.0) then
         if(nleft.eq.0) then
            icount = nxtval(nproc) * ichunk
            nleft = ichunk
         endif
         nxtask_n = icount
         icount = icount + 1
         nleft = nleft -1
      else
         nleft = 0
         nxtask_n = -1
         junk = nxtval(nproc)
      endif
      nstill = nleft
      return
      end
