      Subroutine dft_3cincor(rtdb, geom, AO_bas_han, CD_bas_han, 
     &      nbf_cd, itol2e,n_batch,n_bmax,iwhat,i3c_ERI,n3c_dbl,
     &      iwhat_max)
C$Id: dft_3cincor.F,v 1.4 1996-11-22 19:08:23 d3h449 Exp $
      implicit none
            
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han
      integer nbf_cd    ! no. CD functions [input]
      integer itol2e ! bi-el int tolerance [input]
      integer n_batch        ! no. batches computed by this node [out]
      integer n_bmax  ! max no. batches [in]
      integer iwhat(n_bmax) ! what batches have been computed [out]
      integer n3c_dbl ! max no. ints [in]
      double precision i3c_ERI(n3c_dbl) ! int values [out]
c**** 
#include "bas.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "schwarz.fh"
#include "stdio.fh"
c     
c     Compute the 3-center integrals and store as many of them 
c     as possible in-core.
c
      integer MXshbf_ao,nshells_ao,nshells_cd,MXshbf_cd
      integer me,nproc,ishp,
     ,     icount,Nao2,nERI,maxg,mscratch_2e3c,nscr,
     ,     n1,n2,ishc,ifirstc,ilastc,nshbfc,ishd,ifirstd,ilastd,
     ,     nshbfd,ifirstp,ilastp,nint,
     ,     nshbf_cd,ndone
      integer lscr,iscr
      integer nlast, ntasks, nxtbuf
      integer iwhat_max, nao2_max, batch_max, neri_max
      double precision THRESHOLD, ERI_est, ratio
      integer nxtask
      external nxtask
C      write(0,*) ' in 3CINCO',itol2e
c
      THRESHOLD=10.d0**(-itol2e)
c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
      if( .not. bas_nbf_cn_max(AO_bas_han,MXshbf_ao) )then
        call errquit('Exiting in fitvcoul.',1)
      end if
c      
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in fitvcoul.',2)
      end if
c      
      if( .not. bas_numcont(CD_bas_han,nshells_cd) )then
        call errquit('Exiting in fitvcoul.',3)
      end if
c      
      if( .not. bas_nbf_cn_max(CD_bas_han,MXshbf_cd) )then
        call errquit('Exiting in fitvcoul.',4)
      end if
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
      nao2_max = MXshbf_ao**2
      batch_max = nao2_max*nbf_cd
      nERI_max = Nao2_max*MXshbf_cd
c
c      write(LuOut,*)' me, rtdb, geom, AO_bas_han, CD_bas_han: ',
c     &                me, rtdb, geom, AO_bas_han, CD_bas_han
c      call util_flush(LuOut)
c      write(LuOut,*)' me, mxshbf_ao, mxshbf_cd, nbf_cd: ',
c     &                me, mxshbf_ao, mxshbf_cd, nbf_cd
c      call util_flush(LuOut)
c      write(LuOut,*)' me, icount, nao2_max, batch_max, neri_max: ',
c     &                me, icount, nao2_max, batch_max, neri_max
c
c      call util_flush(LuOut)
c
c     Determine/allocate scratch memory for 3 center 2e- integrals.
c
      call int_mem_2e3c(maxg, mscratch_2e3c)
      nscr    = mscratch_2e3c
      if(.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &  call errquit('dft_3cincor: cannot allocate scr',0)
c
      call ga_sync
c
      n1 = nxtask(nproc,icount)
      n2 = 0
c      
      ndone=1
      n_batch=0
      do 240 ishc = 1,nshells_ao
        if (.not. bas_cn2bfr( AO_bas_han,ishc,ifirstc,ilastc))
     &     call errquit('exiting in fitvcoul.',3)
        nshbfc=ilastc-ifirstc+1
c 
        do 230 ishd = 1,ishc
c            
c         screen integral blocks using the schwarz inequality:  (p|cd) .le. (cd|cd)
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
                if( .not. bas_cn2bfr( cd_bas_han,ishp,
     &            ifirstp,ilastp))
     &            call errquit('exiting in fitvcoul.',5)
                nshbf_cd=ilastp-ifirstp+1
c                  
                nint=nshbf_cd*nao2
c                 
c               compute the 3 center 2e- integrals between this given pair 
c               of ao shells and cd shell.
c
                if( (ndone+nint-1).gt.n3c_dbl) then
                  write(LuOut,*) ' nlast n3cdbl ',ndone+nint-1,n3c_dbl
                  call errquit('dft_3cincor: blew out buf of dbles',999)
                endif
                call int_2e3c(cd_bas_han,ishp,ao_bas_han,ishc,ishd,
     &                        nscr,dbl_mb(iscr),nint,i3c_ERI(ndone))
C                  write(0,*) ' i3c ',ndone,i3c_ERI(ndone)
C                  call util_flush(LuOut)
c
                ndone=ndone+nint
c
 220          continue
c
              n_batch=n_batch+1
              if( n_batch.gt.n_bmax)
     &          call errquit('dft_3cincor: blew out buf of ints',n_bmax)
c
              iwhat(n_batch)=(ishc*(ishc-1))/2 + ishd-1
c            
c             Check to see if there is room to store another batch of
c             integrals in the local buffer.  If so, do next task,
c             otherwise go wait for global sync.
c
              nxtbuf = ndone + batch_max
c
c              write(LuOut,1111)me,ishc,ishd,n1,n2,nxtbuf,n3c_dbl, 
c     &                         n_batch, iwhat(n_batch), ndone
c 1111         format(' Done: me,ishc,ishd,n1,n2,nxtbuf,n3c_dbl, ',
c     &                     ' n_batch, iwhat(n_batch), ndone ',/,10i5)
c              call util_flush(LuOut)
c
              if( nxtbuf.lt.n3c_dbl) then
                n1 = nxtask(nproc,icount)
              else
                go to 241
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
      call ga_sync
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
          goto 242
        endif
      enddo
      ishc = nshells_ao
  242 continue
      ishd = nlast - (ishc*(ishc-1))/2
c
      if(me.eq.0)write(LuOut,3241)ratio*100.d0
c
c      write(LuOut,*)' ntasks, nlast, ishc, ishd: ', 
c     &                ntasks, nlast, ishc, ishd
c      call util_flush(LuOut)
c
      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_3cincor: cannot pop stack',0)
c
      n1 = nxtask(-nproc,icount)
c
c      call ga_sync
c      do n1 = 0, nproc
c        if(me.eq.n1)then
c          write(LuOut,*)' me, n_batch ', me, n_batch
c          call util_flush(LuOut)
c          call i_output(iwhat, 1, n_batch, 1, 1, n_batch, 1, 1)
c          call util_flush(LuOut)
c        endif
c        call ga_sync
c      enddo
c      call ga_sync
c
c      call errquit(' temp stop semi-direct ',999)
c      
 3241 format(/,10x,'The percent of 3c 2e- integrals held in-core is',
     &       /,10x,f6.2)
 3242 format(/,10x,'All 3-center 2e- integrals are in-core. ')
      return
      end
      subroutine i_output (z,rowlow,rowhi,collow,colhi,rowdim,coldim,
     $     nctl)
c.......................................................................
c output prints a real*8 matrix in formatted form with numbered rows
c and columns.  the input is as follows;
c        matrix(*,*).........matrix to be output
c        rowlow..............row number at which output is to begin
c        rowhi...............row number at which output is to end
c        collow..............column number at which output is to begin
c        colhi...............column number at which output is to end
c        rowdim..............row dimension of matrix(*,*)
c        coldim..............column dimension of matrix(*,*)
c        nctl................carriage control flag; 1 for single space
c                                                   2 for double space
c                                                   3 for triple space
c the parameters that follow matrix are all of type integer*4.  the
c program is set up to handle 5 columns/page with a 1p5d24.15 format for
c the columns.  if a different number of columns is required, change
c formats 1000 and 2000, and initialize kcol with the new number of
c columns.
c author;  nelson h.f. beebe, quantum theory project, university of
c          florida, gainesville
c.......................................................................
C$Id: dft_3cincor.F,v 1.4 1996-11-22 19:08:23 d3h449 Exp $
      implicit none
      integer rowlow,rowhi,collow,colhi,rowdim,coldim,begin,kcol
      integer nctl, i, j, last, k
      integer z(rowdim,coldim), zero
      character*8 asa(3), column, ctl, blank
      data column/'column'   /,asa/' ','00000000'  ,
     1     '--------'  /,blank/' '/
      data kcol/8/
      data zero/0/
      do 11 i=rowlow,rowhi
         do 10 j=collow,colhi
            if (z(i,j).ne.zero) go to 15
 10      continue
 11   continue
      write (6,3000)
 3000 format (/' zero matrix'/)
      go to 3
 15   continue
      ctl = blank
      if ((nctl.le.3).and.(nctl.gt.0)) ctl = asa(nctl)
      if (rowhi.lt.rowlow) go to 3
      if (colhi.lt.collow) go to 3
      last = min(colhi,collow+kcol-1)
      do 2 begin = collow,colhi,kcol
*         write (6,1000) (column,i,i = begin,last)
         write (6,1000) (i,i = begin,last)
         do 1 k = rowlow,rowhi
            do 4 i=begin,last
               if (z(k,i).ne.zero) go to 5
 4          continue
            go to 1
 5          write (6,2000) ctl,k,(z(k,i), i = begin,last)
 1       continue
         last = min(last+kcol,colhi)
 2    continue
 3    return
* kcol = 4
* 1000 format (/1h ,16x,3(a6,i3,2x),(a6,i3))
* 2000 format (a1,3hrow,i4,2x,4f17.11)
* kcol = 8
*
* if U like having rows and columns labelled with row and col
* use these
*
* 1000 format (/1h ,11x,7(a3,i3,3x),(a3,i3))
* 2000 format (a1,'row',i4,1x,8f9.4)
c
 1000 format (/1h ,8x,7('   ',i3,3x),('   ',i3))
 2000 format (a1,i4,1x,8i9)
      end
