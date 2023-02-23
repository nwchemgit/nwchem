      subroutine int_1e_sifs(ibas, aoints, energy, nenrgy, nbft,
     &   nmap,map,imtype,ifmt1,ibvtyp,ibitv,l1rec,n1max,
     &   clab, ninfo, info)
      implicit none
#include "errquit.fh"
#include "cint1cache.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "cdft.fh"
c     
c     Compute the desired type of integrals (kinetic, potential, overlap)
c     
c     Oskel indicates that the skeleton (petite-list symmetry) matrix should be
c     built ... 
c     
c     arguments
c     
      integer ibas
      character*(10) integ_type  ! [input] Name of integrals to compute
*     logical oskel             ! [input] If true generate symmetry unique list
c     
c     local variables
c     
      integer nshell
      integer ishell, jshell, mem1, max1e, lrecal, n1mx
      integer ijshell, ilo, ihi, idim
      integer l_buf, l_scr, l_ilab, l_jlab, l_info, l_bfnlab, l_sifbuf,
     &        l_sifval
      integer k_buf, k_scr, k_ilab, k_jlab, k_info, k_bfnlab, k_sifbuf, 
     &        k_sifval
      integer type
      logical odoit
      double precision q2
      external block_int1e      ! For T3D
      integer i, noffset,g_loc, j, ijmap
c
      logical ocache_save
c
      logical odbug, Lresult
      logical osome
      character*8 blabtmp
      integer nbft
c  
c     sifs parameters
c
      integer*4 aoints
      integer*4 ifmt1, ibvtyp, n1max, ierr 
      INTEGER l1rec, ntitle, ibuf, nsym
      integer*4 ibitv
      integer*4 mxbf
      parameter (mxbf=1000)
c     header 1
      integer*4 nbas, ninfo, mxenrgy, nenrgy, nmap
      parameter(ntitle=1,nsym=1,mxenrgy=1)
c     header 2
      integer*4 info(ninfo)  ! not to be confused with NWChem info
      integer*4 ietype(mxenrgy)
      integer*4 imtype(nmap) ! map will be taken from dynamic memory 
      integer*4 map(nbft,nmap)
      integer*4 clab(2,*)
      real*8 energy(nenrgy)
      character*80 stitle(ntitle)
      character*4 slabel(nsym)
      integer*4 nbpsy(nsym)
      character*1 shtypes(-1:7)
      data shtypes/'l','s', 'p', 'd', 'f', 'g', 'h', 'i', 'k'/
      integer*4 shmap(-1:7)
      data shmap  /  0,  1,   2,   4,   6,   9,  12,  16,  0 /
      integer otype, nprim, ngen, sphcart, iat, shdim, igen, ibf
      integer numtot
      integer*4 msame,nmsame,nomore
      parameter(msame=0, nmsame=1, nomore= 2)
      integer*4 itypea, itypeb, nrec, num,  last
      real*8 fcore
      logical basok
      

cgk provisional
*     change thresh to the appropriate user supplied zero tolerance
      real*8 thresh
      parameter (thresh=1d-12)
cgk end
      integer numints
      
c
cgk debug
!      write(*,*)'gk: **** entered int_1e_sifs'
!      WRITE(*,*)'LB: info=',info
cgk end
      basok=.false.
      odbug=.true.
      osome=.true.
      osome=osome.or.odbug
      odbug=odbug.and.(ga_nodeid().eq.0)
      osome=osome.and.(ga_nodeid().eq.0)


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
     $     ('int_1e_sifs: bas_numcont failed for ibas', ibas,
     &       BASIS_ERR)
      if (nbft.gt.mxbf) call errquit
     $     ('int_1e_sifs: nbft gt maximum aoints basis functions', nbft,
     &       BASIS_ERR)
c     
c     allocate necessary local temporary arrays on the stack
c     
c     k_* are the offsets corrsponding to the l_* handles
c     
      call int_mem_1e(max1e, mem1)
cgk debug
!      write(*,*)'gk: in int_mem_1e, max1e=', max1e, ' mem1=',mem1
cgk end
c     get memory for integral buffer
      if (.not. MA_push_get(MT_DBL,max1e,'int_1e_sifs:buf',l_buf,k_buf))
     $     call errquit('int_1e_sifs: ma failed', max1e, MA_ERR)
c     get memory for scratch space
      if (.not. MA_push_get(MT_DBL, mem1,'int_1e_sifs:scr',l_scr,k_scr))
     $     call errquit('int_1e_sifs: ma failed', mem1, MA_ERR)
c     get memory for bra labels
      if (.not. MA_push_get(MT_Int,max1e,'int_1e_sifs:ilab',l_ilab,
     &        k_ilab))
     $     call errquit('int_1e_sifs: ma failed for ilab', max1e, 
     &        MA_ERR)
c     get memory for ket labels
      if (.not. MA_push_get(MT_Int,max1e,'int_1e_sifs:jlab',l_jlab,
     &        k_jlab))
     $     call errquit('int_1e_sifs: ma failed for jlab', max1e, 
     &        MA_ERR)
c     get memory for basis info (character*16)
      if (.not. ma_push_get(mt_byte,nbft*16,'int_1e_sifs:info',
     $        l_info,k_info))
     $     call errquit('int_1e_sifs:  error getting info mem', 
     &        nbft*16, MA_ERR)
c     get memory for sifs bfnlab (character*8)
      if (.not. ma_push_get(mt_byte,nbft*8,'int_1e_sifs:bfnlab',
     $        l_bfnlab,k_bfnlab))
     $     call errquit('int_1e_sifs:  error getting bfnlab mem', 
     &        nbft*8, MA_ERR)
c     get memory for sifs buffer space
      lrecal=l1rec
      if (.not. MA_push_get(MT_DBL, lrecal,'int_1e_sifs:sifbuf',
     &     l_sifbuf, k_sifbuf))
     $     call errquit('int_1e_sifs: ma failed getting sifbuf', lrecal,
     &        MA_ERR)
      n1mx=n1max
      if (.not. MA_push_get(MT_DBL, n1mx,'int_1e_sifs:sifval',
     &     l_sifval, k_sifval))
     $     call errquit('int_1e_sifs: ma failed getting sifval', n1mx, 
     &        MA_ERR)

c
c     Set up SIFS header iformation
c
      call bas_vec_info(ibas, byte_mb(k_info))
cgk debug
*     write(*,*)'gk: basis set info'
*     do i=0,nbft-1
*        write(6,'(i4,2x,16a1)') i+1, (byte_mb(k_info+i*16+j), j=0,15)
*     enddo
cgk end
c     gk: instead of doing this, use 16 byte labels in SIFS/COLUMBUS.
c     gk: consider this a temporary hack.  8 characters is not enough.
      do i=0,nbft-1
         do j=0,7
            write(byte_mb(k_bfnlab+i*8+j),'(a1)') 
     &         byte_mb(k_info+i*16+j) 
         enddo
      enddo
cgk debug
!      write(*,*)'gk: bfnlab:'
!      do i=0,nbft-1
!         write(6,'(i4,2x,8a1)') i+1, (byte_mb(k_bfnlab+i*8+j), j=0,7)
!      enddo
cgk end

c     put effective nuclear repulsion energy in energy(1)
      ietype(1) = -1 ! nuclear repulsion=1 ... for now, but core later
      LResult = geom_nuc_rep_energy(geom, energy(1))

      imtype(1)=3 ! bfn-to-center map vector
      imtype(2)=4 ! bfn-to-orbital_type map vector
      fcore=0d0 

cgk debug
!       write(*,*)'gk:  i, iat, shmap(otype)'
!       write(*,*)'gk:  nmap=',nmap
cgk end
      i=0
      do ishell=1,nshell
         basok=bas_continfo(ibas, ishell, otype, nprim, ngen, sphcart)
         basok=bas_cn2bfr(ibas, ishell, ilo, ihi)
         basok=bas_cn2ce(ibas, ishell, iat)
         if (otype.gt.6 .or. otype.lt.0) call errquit
     $        ('int_1e_sifs: unsupported sifs basis otype?',otype,
     &       BASIS_ERR)
c
         shdim = (ihi-ilo+1)
         if (otype.ge.0) then
            shdim = shdim / ngen
         else     ! Fudge for SP shells
            ngen = 1
         endif
         if (otype .le. 1) sphcart = 0
          
         do igen = 1, ngen 
            do ibf = 1, shdim
 
               i = i + 1 

               if (i.gt.nbft) call errquit
     $            ('int_1e_sifs: error counting basis functions',nbft,
     &             BASIS_ERR)

               map(i,1)=iat
               map(i,2)=shmap(otype)
cgk debug
*       write(*,'(a,i,i,i)')'gk: ', i, iat, shmap(otype)
*       write(*,'(a,i,i,i)')'gk: ', i, map(i,1), map(i,2)
cgk end
*              int_mb(k_map+i-1) = iat
*              int_mb(k_map+nbft+i-1) = shmap(otype)

            enddo
         enddo
      enddo

      slabel(1)='a1'
      stitle(1)='AO integrals from NWChem'
      nbpsy(1)=nbft

*     call sifwh( aoints,  ntitle,  nsym,    nbft,
*    &   ninfo,   nenrgy,  nmap,    stitle,   nbpsy,   slabel,  info,    
*    &   bfnlab,  ietype,  energy,  imtype,  map, ierr )
      
      ijshell = 0
      q2 = 1.0d0
      if (0 .eq. ga_nodeid()) then

cgk debug
!       write(*,*)'gk: calling sifwh with aoints=',aoints
!       write(*,*)'gk:    nenrgy=',nenrgy
!       write(*,*)'gk:    ntitle=',ntitle
!       write(*,*)'gk:    ninfo=',ninfo
!       WRITE(*,*)'LB: l2rec=',info(4)
cgk end
         call sifwh( aoints,  ntitle,  nsym,    nbft,
     &      ninfo,   nenrgy,  nmap,    stitle,   nbpsy,   slabel,  info,
     &      byte_mb(k_bfnlab),  ietype,  energy,  imtype,  map, 
     &      ierr )
cgk debug
!       write(*,*)'gk: back from sifwh with ierr=',ierr
cgk end

c        start loop over integral types here.
         do type=1,3
         nrec=0
         num=0
         numtot=0
         ibuf=0
         last=msame
         if (type.eq.1) then
            integ_type='kinetic'
            itypea=0
            itypeb=1
         else if (type.eq.2) then
            integ_type='potential'
            itypea=0
            itypeb=2
         else if (type.eq.3) then
            integ_type='overlap'
            itypea=0
            itypeb=0
         end if

cgk debug
!      write(*,*)'gk: processing integral type: ', integ_type
cgk end

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
     $               call errquit('int_1e_ga: bas_cn2bfr ?', ibas,
     &                   BASIS_ERR)
                  idim = ihi - ilo + 1

                  call ifill(max1e,0,int_mb(k_ilab),1)
                  call ifill(max1e,0,int_mb(k_jlab),1)
      
                  if (type .eq. 1) then
                     call int_l1eke (ibas, jshell, ibas, ishell, 
     &                  thresh, int_mb(k_ilab),int_mb(k_jlab), max1e,
     $                  dbl_mb(k_buf), mem1, dbl_mb(k_scr), numints)
                  else if (type .eq. 2) then
cgk debug
*       write(*,*)'gk: calling int_l1epe'
*       write(*,*)'gk: ibas = ', ibas
*       write(*,*)'gk: ishell = ', ishell
*       write(*,*)'gk: jshell = ', jshell
*       write(*,*)'gk: thresh = ', thresh
*       write(*,*)'gk: k_ilab = ', k_ilab
*       write(*,*)'gk: k_jlab = ', k_jlab
*       write(*,*)'gk: max1e = ', max1e
*       write(*,*)'gk: k_buf = ', k_buf
*       write(*,*)'gk: mem1 = ', mem1
*       write(*,*)'gk: k_scr = ', k_scr
*       write(*,*)'gk: numints = ', numints
cgk end
                     call int_l1epe (ibas, jshell, ibas, ishell, 
     &                  thresh, int_mb(k_ilab),int_mb(k_jlab), max1e,
     $                  dbl_mb(k_buf), mem1, dbl_mb(k_scr), numints)
cgk debug
*       write(*,*)'gk: back from int_l1epe'
cgk end
                  else if (type .eq. 3) then
c                    ECP is summed here, but COLUMBUS should not care.
                     call int_l1eov (ibas, ishell, ibas, jshell, thresh,
     &                  int_mb(k_ilab),int_mb(k_jlab), max1e,
     &                  dbl_mb(k_buf), mem1, dbl_mb(k_scr), numints)
                  end if
c     
cgk debug
*      write(*,*)'gk: this batch, numints = ', numints
cgk end

                  do i=1,numints

                     if(ibuf.eq.n1max) then
                       numtot = numtot + ibuf
C LB
!                       WRITE(*,*)"int_1e_sifs calling sifew1 1,
!     &last=",last,"ibvtyp=",ibvtyp
!                   WRITE(*,*)"sizeof(dbl_mb)=",sizeof(dbl_mb(k_sifval)),
!     &             "sizeof(buf)=",sizeof(dbl_mb(k_sifbuf))

C LB
                       call sifew1(aoints, info, 2, ibuf, last, 
     &                   itypea, itypeb, ifmt1, ibvtyp, 
     &                   dbl_mb(k_sifval), clab, fcore, ibitv, 
     &                   dbl_mb(k_sifbuf), nrec, ierr)
c                      ibuf on return has the number of unwritten
c                      integrals.  dbl_mb(k_sifval+0:(ibuf-1)) 
c                      contains the unwritten values, ditto for labels 
                       numtot = numtot - ibuf
                     endif
                      
                     ibuf=ibuf+1
                     clab(1,ibuf)=int_mb(k_ilab+i-1)
                     clab(2,ibuf)=int_mb(k_jlab+i-1)
                     dbl_mb(k_sifval+ibuf-1)=dbl_mb(k_buf+i-1)      
cgk debug
*                     write(*,*)'gk: ibuf=', ibuf
*                     write(*,'(a,2i4,f20.12)') 
*     &               'gk:  ilab,jlab,dbl_mb:',
*     &               int_mb(k_ilab+i-1),int_mb(k_jlab+i-1),
*     &               dbl_mb(k_buf+i-1)
cgk end

                  enddo

               end if
               ijshell = ijshell + 1
            end do
         end do
         last=nmsame
         numtot=numtot+ibuf
cgk debug
*       write(*,*)'gk: nrec=', nrec
*       write(*,*)'gk:  ibuf=', ibuf
*       write(*,*)'gk:  last=', last,"ibvtyp=",ibvtyp
*       write(*,*)'gk:  numtot=', numtot
*       WRITE(*,*)"int_mom_sifs calling sifew1 2"
cgk end
         call sifew1(aoints, info, 2, ibuf, last, 
     &     itypea, itypeb, ifmt1, ibvtyp, 
     &     dbl_mb(k_sifval), clab, fcore, ibitv, 
     &     dbl_mb(k_sifbuf), nrec, ierr)

         write(6,'(a,i10,1x,a,a)') 'Wrote ',numtot, integ_type, 
     &     ' integrals to aoints'

         enddo
      endif
      call ga_sync()            ! So that no nasty races can result

cgk debug
!      write(*,*)'gk: calling int_so_sifs from int_1e_sifs'
*      write(*,*)'gk:  mem1=', mem1
cgk end
      call int_so_sifs(ibas, oskel, aoints, nbft, max1e, mem1, l1rec,
     &  n1max, dbl_mb(k_sifbuf), dbl_mb(k_sifval), ninfo, info, clab,
     &  fcore, ibvtyp, ibitv, ifmt1)

cgk debug
!      write(*,*)'gk: calling int_mom_sifs from int_1e_sifs'
*      write(*,*)'gk:  mem1=', mem1
cgk end
      call int_mom_sifs(ibas, oskel, aoints, nbft, l1rec,
     &  n1max, dbl_mb(k_sifbuf), dbl_mb(k_sifval), ninfo, info, 
     &  fcore, ibvtyp, ibitv, ifmt1)
cgk debug
!      write(*,*)'gk: back from int_so_sifs in int_1e_sifs'
cgk end
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_sifval)) call errquit
     $     ('int_1e_sifs: pop failed at l_bfnlab', 0, GA_ERR)
      if (.not. MA_pop_stack(l_sifbuf)) call errquit
     $     ('int_1e_sifs: pop failed at l_bfnlab', 0, GA_ERR)
      if (.not. MA_pop_stack(l_bfnlab)) call errquit
     $     ('int_1e_sifs: pop failed at l_bfnlab', 0, GA_ERR)
      if (.not. MA_pop_stack(l_info)) call errquit
     $     ('int_1e_sifs: pop failed at l_info', 0, GA_ERR)
      if (.not. MA_pop_stack(l_jlab)) call errquit
     $     ('int_1e_sifs: pop failed at l_jlab', 0, GA_ERR)
      if (.not. MA_pop_stack(l_ilab)) call errquit
     $     ('int_1e_sifs: pop failed at l_ilab', 0, GA_ERR)
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('int_1e_sifs: pop failed at l_scr', 0, GA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('int_1e_sifs: pop failed at l_buf', 0, GA_ERR)
c
c     ----- restore ocache -----
c
      if(type.eq.6 .or. type.eq.7) then
         ocache=ocache_save
      endif
c     
*     if (.not.inp_compare(.false., integ_type, 'so')) then
*        call ga_dadd(1.0d0, g_loc, 1.0d0, g, g)
*        if (.not. ga_destroy(g_loc)) call errquit('i1ega: gad?',0,
*    &       GA_ERR)
*     endif
      call ga_sync()            ! So that no nasty races can result
c     
      if (oscfps) call pstat_off(ps_int_1e)
c     
cgk debug
!      write(*,*)'gk: leaving int_1e_sifs at end of routine'
cgk end
      end
