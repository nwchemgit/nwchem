      subroutine int_so_sifs(ibas, oskel, aoints, nbft, max1e, mem1,
     &  l1rec, n1max, 
     &  sifbuf, sifval, ninfo, info, clab, fcore, ibvtyp, ibitv,
     &  ifmt1)
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "cscfps.fh"
#include "sym.fh"
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... 
c     
c     arguments
c     
      integer ibas
      logical oskel             ! [input] If true generate symmetry unique list
      double precision sifval(n1max)
      double precision sifbuf(l1rec)
      integer*4 ninfo
      integer*4 info(ninfo)  ! not to be confused with NWChem info
      integer*4 clab(2,max1e)
      real*8 fcore

c     
c     local variables
c     
      integer nshell
      integer ishell, jshell, mem1, max1e, lrecal, n1mx, maxint
      integer ilo, ihi, jlo, jhi, max1eso, mem1so
      integer l_buf, l_scr, 
     &        l_sox, l_soy, l_soz, l_iso, l_jso
      integer k_buf, k_scr, k_ilab, k_jlab, 
     &        k_sox, k_soy, k_soz, k_iso, k_jso
      integer type, noffsetx, noffsety, noffsetz
      logical odoit
      double precision q2
      external block_int1e      ! For T3D
      integer i, noffset,g_loc, j, ijmap, symmap, symints
c
      logical ocache_save
c
      logical odbug
      logical osome
      character*8 blabtmp
      integer nbft
c  
c     sifs parameters
c
      integer*4 aoints
      integer*4 ifmt1, l1rec, n1max, ierr 
      integer*4 ibvtyp
      integer*4 ibitv
      integer*4 mxbf
c     make sure this is compatible with columbus
      parameter (mxbf=1000)
c     header 1
      integer*4 ntitle, nsym, nbas, mxenrgy
      parameter(ntitle=1,nsym=1,mxenrgy=1)
c     header 2
      integer otype, nprim, ngen, sphcart, iat, shdim, igen, ibf
      integer numtot
      integer*4 msame,nmsame,nomore
      parameter(msame=0, nmsame=1, nomore= 2)
      integer*4 nrec, num
      integer itypea, itypeb, last
      integer*4 ibuf
      logical basok
*     character*3 sotype(3)
*     data sotype /'sox','soy','soz'/
      real*8 soval
      

cgk provisional
*     change thresh to the appropriate user supplied zero tolerance
      real*8 thresh
      parameter (thresh=1d-12)
cgk end
      integer numints, icart, kcart
      
c
cgk debug
*      write(*,*)'gk: **** entered int_so_sifs'
cgk end
      basok=.false.
      odbug=.true.
      osome=.true.
      osome=osome.or.odbug
      odbug=odbug.and.(ga_nodeid().eq.0)
      osome=osome.and.(ga_nodeid().eq.0)
c
      call ga_sync()
      if (oscfps) call pstat_on(ps_int_1e)
c     
c     
c
c     ----- save ocache logical variable -----
c
      if(type.eq.6 .or. type.eq.7 ) then
         ocache_save=ocache
         ocache     =.false.
      endif
c     
c     Get info about the basis sets
c     
      if (.not. bas_numcont(ibas, nshell)) call errquit
     $     ('int_so_sifs: bas_numcont failed for ibas', ibas,
     &       BASIS_ERR)
      if (nbft.gt.mxbf) call errquit
     $     ('int_so_sifs: nbft gt maximum aoints basis functions', nbft,
     &       BASIS_ERR)
c     
c     allocate necessary local temporary arrays on the stack
c     
c     l_buf ... buffer to hold shell block of matrix
c     l_s   ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c     
c     k_* are the offsets corrsponding to the l_* handles
c     
cgk debug
*      write(*,*)'gk: max1e=', max1e
*      write(*,*)'gk: mem1=', mem1
*      write(*,*)'gk: nbft=', nbft
cgk end
      max1eso=3*max1e
      mem1so=3*mem1
      maxint=nbft*(nbft+1)/2 
cgk debug
*      write(*,*)'gk: max1eso=', max1eso
*      write(*,*)'gk: mem1so=', mem1so
*      write(*,*)'gk: maxint=', maxint
*      WRITE(*,*)"LB: start MA_push_get"
cgk end

C LB
*      WRITE(*,*)"LB: integral buffer"
C LB
c     get memory for integral buffer
      if (.not. MA_push_get(MT_DBL,max1eso,'int_so_sifs:buf',l_buf,
     &  k_buf))
     $     call errquit('int_so_sifs: ma failed', max1eso, MA_ERR)
c     get memory for scratch space
      if (.not. MA_push_get(MT_DBL, mem1so,'int_so_sifs:scr',l_scr,
     &  k_scr))
     $     call errquit('int_so_sifs: ma failed', mem1so, MA_ERR)
c     get memory for x spin orbit
      if (.not. MA_push_get(MT_DBL,maxint,'int_so_sifs:sox',l_sox,
     &        k_sox))
     $     call errquit('int_so_sifs: ma failed for so x', maxint, 
     &        MA_ERR)
c     get memory for y spin orbit
      if (.not. MA_push_get(MT_DBL,maxint,'int_so_sifs:soy',l_soy,
     &        k_soy))
     $     call errquit('int_so_sifs: ma failed for so y', maxint, 
     &        MA_ERR)
c     get memory for z spin orbit
      if (.not. MA_push_get(MT_DBL,maxint,'int_so_sifs:soz',l_soz,
     &        k_soz))
     $     call errquit('int_so_sifs: ma failed for so z', maxint, 
     &        MA_ERR)
c     get memory for full bra label array
      if (.not. MA_push_get(MT_Int,maxint,'int_so_sifs:iso',l_iso,
     &        k_iso))
     $     call errquit('int_so_sifs: ma failed for iso', maxint, 
     &        MA_ERR)
c     get memory for full ket label array
      if (.not. MA_push_get(MT_Int,maxint,'int_so_sifs:jso',l_jso,
     &        k_jso))
     $     call errquit('int_so_sifs: ma failed for iso', maxint, 
     &        MA_ERR)

      if (0 .eq. ga_nodeid()) then

c        start loop over so integrals here.
         numints=0
         do jshell = 1, nshell
            do ishell = 1, nshell
cgk debug
*      write(*,*)'gk: ishell=',ishell,' jshell=',jshell
cgk end
               odoit = .true.
               if (oskel) 
     $            odoit = sym_shell_pair(ibas, ishell, jshell, q2)
 
               if (odoit) then
                  if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) 
     $               call errquit('int_1e_ga: bas_cn2bfr ishell', 
     &                   ishell, BASIS_ERR)
                  if (.not. bas_cn2bfr(ibas, jshell, jlo, jhi)) 
     $               call errquit('int_1e_ga: bas_cn2bfr jshell', 
     &                   jshell, BASIS_ERR)
                  call intso_1e (ibas, jshell, ibas, ishell, 
     $               mem1, dbl_mb(k_scr), max1eso, dbl_mb(k_buf))

                  noffsetz = 0
                  noffsety = (ihi-ilo+1)*(jhi-jlo+1)   
                  noffsetx = (ihi-ilo+1)*(jhi-jlo+1)*2   
                  do j=jlo,jhi
                    do i=ilo,ihi
cgk debug
*     write(*,'(a,i4,a,i4,a,i4,a,i4)')
*    &  'gk: ilo=',ilo,' ihi=',ihi,' jlo=',jlo,' jhi=',jhi
*      write(*,*)'gk: numints=',numints
cgk end
                       ijmap=(j-jlo)*(ihi-ilo+1)+(i-ilo)
                       if(i.ge.j) then
                         symmap=nbft*(j-1)-((j-1)*j)/2+i
*     write(*,*)'gk: symmap=',symmap
                         dbl_mb(k_sox+symmap-1)=
     &                     dbl_mb(k_buf+noffsetx+ijmap)
                         dbl_mb(k_soy+symmap-1)=
     &                     dbl_mb(k_buf+noffsety+ijmap)
                         dbl_mb(k_soz+symmap-1)=
     &                     dbl_mb(k_buf+noffsetz+ijmap)
                         int_mb(k_iso+symmap-1)=i
                         int_mb(k_jso+symmap-1)=j 
cgk debug
*                      write(*,'(a,3i3,f18.12)')
*    &                   'gk:  i,j,ij,dbl_mb(so_x):',
*    &                    i,j,numints+ijmap,dbl_mb(k_sox+symmap-1)
*                      write(*,'(a,3i3,f18.12)')
*    &                   'gk:  i,j,ij,dbl_mb(so_y):',
*    &                    i,j,numints+ijmap,dbl_mb(k_soy+symmap-1)
*                      write(*,'(a,3i3,f18.12)')
*    &                   'gk:  i,j,ij,dbl_mb(so_z):',
*    &                    i,j,numints+ijmap,dbl_mb(k_soz+symmap-1)
cgk end
                       endif
                    enddo
                  enddo
                  numints=numints+ijmap+1
               end if
            end do
         end do
c call sifs
      endif
      call ga_sync()            ! So that no nasty races can result
cgk debug
*      write(*,*)'gk: numints = ', numints
*      write(*,*)'gk: maxint = ', maxint
*     symmap=0
*     do i=1,nbft
*       do j = 1, i
*         write(*,'(a,2i8,3f18.12)')
*    &    'gk:  i,j,so_x,so_y,so_z:', int_mb(k_iso+symmap),
*    &    int_mb(k_jso+symmap),dbl_mb(k_sox+symmap), 
*    &    dbl_mb(k_soy+symmap), dbl_mb(k_soz+symmap)
*         symmap=symmap+1
*       enddo
*     enddo
      do icart = 1, 3
        itypea=2
        if(icart.eq.1) then
          kcart=k_soz
          itypeb=2
        elseif(icart.eq.2) then
          kcart=k_soy
          itypeb=1
        elseif(icart.eq.3) then
          kcart=k_sox
          itypeb=0
        endif
        symmap=0
        numtot=0
        ibuf=0
        do i=1,nbft
          do j = 1, i
            if(ibuf.eq.n1max) then
              numtot = numtot + ibuf
cgk debug
*      WRITE(*,*)"int_so_sifs calling sifew1 1"
*      write(*,*)'gk: nrec=', nrec
*      write(*,*)'gk:  ibuf=', ibuf
*      write(*,*)'gk:  numtot=', numtot
*      WRITE(*,*)"last=",last,"ibvtyp=",ibvtyp
cgk end
              call sifew1(aoints, info, 2, ibuf, msame, 
     &          itypea, itypeb, ifmt1, ibvtyp, 
     &          sifval, clab, fcore, ibitv, 
     &          sifbuf, nrec, ierr)
              numtot = numtot - ibuf
            endif

            soval=dbl_mb(kcart+symmap)
            
            if(abs(soval).gt.thresh) then
               ibuf=ibuf+1
               clab(1,ibuf)=int_mb(k_iso+symmap)
               clab(2,ibuf)=int_mb(k_jso+symmap)
               sifval(ibuf)=soval
            endif
            symmap=symmap+1
          enddo
        enddo
*       if (icart.eq.3) then
*         last=nomore 
*       else
          last=nmsame
*       endif
        numtot=numtot+ibuf
cgk debug
*      WRITE(*,*)"int_so_sifs calling sifew1 2"
*     write(*,*)'gk: nrec=', nrec
*      write(*,*)'gk:  ibuf=', ibuf
*     write(*,*)'gk:  numtot=', numtot
*      WRITE(*,*)"last=",last,"ibvtyp=",ibvtyp
cgk end
        call sifew1(aoints, info, 2, ibuf, last, 
     &    itypea, itypeb, ifmt1, ibvtyp, 
     &    sifval, clab, fcore, ibitv, 
     &    sifbuf, nrec, ierr)
*       write(6,'(a,i10,1x,a,a)') 'Wrote ',numtot, sotype(icart), 
*    &     ' integrals to aoints'
      enddo
cgk end
cgk debug
*     write(*,*)'gk:  l_jso = ', l_jso
*     write(*,*)'gk:  k_jso = ', k_jso
cgk end
c     
c     chop stack at first item allocated
c     
*     if (.not. MA_pop_stack(l_sifval)) call errquit
*    $     ('int_so_sifs: pop failed at l_sifval', 0, GA_ERR)
*     if (.not. MA_pop_stack(l_sifbuf)) call errquit
*    $     ('int_so_sifs: pop failed at l_sifbuf', 0, GA_ERR)
      if (.not. MA_pop_stack(l_jso)) call errquit
     $     ('int_so_sifs: pop failed at l_jso', 0, GA_ERR)
      if (.not. MA_pop_stack(l_iso)) call errquit
     $     ('int_so_sifs: pop failed at l_iso', 0, GA_ERR)
      if (.not. MA_pop_stack(l_soz)) call errquit
     $     ('int_so_sifs: pop failed at l_soz', 0, GA_ERR)
      if (.not. MA_pop_stack(l_soy)) call errquit
     $     ('int_so_sifs: pop failed at l_soy', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sox)) call errquit
     $     ('int_so_sifs: pop failed at l_sox', 0, GA_ERR)
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_so_sifs: pop failed at l_scr', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_so_sifs: pop failed at l_buf', 0, GA_ERR)
 
      call ga_sync()            ! So that no nasty races can result
      
      if (oscfps) call pstat_off(ps_int_1e)
c     
cgk debug
*     write(*,*)'gk: leaving int_so_sifs at end of routine'
cgk end
      end
