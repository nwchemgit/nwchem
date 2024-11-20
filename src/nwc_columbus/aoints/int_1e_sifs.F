      subroutine int_1e_sifs(ibas, aoints, energy, nenrgy, nbft,
     &   nmap,map,imtype,ibvtyp,ibitv,l1rec,n1max,
     &   clab, ninfo, info, csolab)
      USE nwc_sym, ONLY: mxsym,mxbf,zir,nir,nbpsy
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
     &        l_sifval,l_SOval,l_SOlab
      integer k_buf, k_scr, k_ilab, k_jlab, k_info, k_bfnlab, k_sifbuf, 
     &        k_sifval,k_SOval,k_SOlab
      integer type
      logical odoit
      double precision q2
      external block_int1e      ! For T3D
      integer i, noffset,g_loc, j, ijmap, k
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
      integer aoints
      integer n1max
      INTEGER l1rec, ntitle, ibuf, nsym, ibvtyp, ierr
      integer ibitv
!      integer mxbf
!      parameter (mxbf=1000)
c     header 1
      integer nbas, mxenrgy, nenrgy, nmap
      INTEGER ninfo
      !parameter(ntitle=1,nsym=1,mxenrgy=1)
      parameter(ntitle=1,mxenrgy=1)
c     header 2
      integer info(ninfo)  ! not to be confused with NWChem info
      integer ietype(mxenrgy)
      integer imtype(nmap) ! map will be taken from dynamic memory 
      integer map(nbft,nmap)
      integer clab(2,*)
      integer cSOlab(2,*)
      real*8 energy(nenrgy)
      character*80 stitle(ntitle)
      character*4 slabel(mxsym)
      character*1 shtypes(-1:7)
      data shtypes/'l','s', 'p', 'd', 'f', 'g', 'h', 'i', 'k'/
      integer shmap(-1:7)
      data shmap  /  0,  1,   2,   4,   6,   9,  12,  16,  0 /
      integer otype, nprim, ngeno, sphcart, iat, shdim, igen, ibf
      integer numtot
      integer msame,nmsame,nomore
      parameter(msame=0, nmsame=1, nomore= 2)
      integer nrec
      INTEGER itypea, itypeb, num, last
      real*8 fcore
      logical basok
      
C LB
       INTEGER :: nsoints ! number of symmetry orbital integrals
C LB

cgk provisional
*     change thresh to the appropriate user supplied zero tolerance
      real*8 thresh
      parameter (thresh=1d-12)
cgk end
      integer numints
      
c
cgk debug
      write(*,*)'gk: **** entered int_1e_sifs'
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
c     LB type is never initialized, ocache never used
      !if(type.eq.6 .or. type.eq.7 ) then
      !   ocache_save=ocache
      !   ocache     =.false.
      !endif
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
c     get memory for SO buffer space
      if (.not. MA_push_get(MT_DBL, n1mx,'int_1e_sifs:SOval',
     &     l_SOval, k_SOval))
     $     call errquit('int_1e_sifs: ma failed getting SOval', n1mx, 
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
!        WRITE(*,*)"C LB, nshell=",nshell
cgk end
      i=0
      do ishell=1,nshell
         basok=bas_continfo(ibas, ishell, otype, nprim, ngeno, sphcart)
!        WRITE(*,*)"C LB, ibas=",ibas
!        WRITE(*,*)"C LB, ishell=",ishell
!        WRITE(*,*)"C LB, otype=",otype
!        WRITE(*,*)"C LB, nprim=",nprim
!        WRITE(*,*)"C LB, ngeno=",ngeno
!        WRITE(*,*)"C LB, sphcart=",sphcart
         basok=bas_cn2bfr(ibas, ishell, ilo, ihi)
!        WRITE(*,*)"C LB, ilo=",ilo
!        WRITE(*,*)"C LB, ihi=",ihi
         basok=bas_cn2ce(ibas, ishell, iat)
!        WRITE(*,*)"C LB, iat=",iat
         if (otype.gt.6 .or. otype.lt.0) call errquit
     $        ('int_1e_sifs: unsupported sifs basis otype?',otype,
     &       BASIS_ERR)
c
         shdim = (ihi-ilo+1)
         if (otype.ge.0) then
            shdim = shdim / ngeno
         else     ! Fudge for SP shells
            ngeno = 1
         endif
         if (otype .le. 1) sphcart = 0
          
         do igen = 1, ngeno 
            do ibf = 1, shdim
 
               i = i + 1 

               if (i.gt.nbft) call errquit
     $            ('int_1e_sifs: error counting basis functions',nbft,
     &             BASIS_ERR)

               map(i,1)=iat
               map(i,2)=shmap(otype)
cgk debug
!       write(*,'(a,i,i,i)')'gk: ', i, iat, shmap(otype)
!       write(*,'(a,i,i,i)')'gk: ', i, map(i,1), map(i,2)
cgk end
*              int_mb(k_map+i-1) = iat
*              int_mb(k_map+nbft+i-1) = shmap(otype)

            enddo
         enddo
      enddo

c      slabel(1)='a1'
C LB: Add symmetry labels to aoints files
!      status = sym_char_table(grp_name,nop,nir,iclass,zir,zclass,chars)
!      WRITE(*,*)"C LB, nir=",nir
!      WRITE(*,*)"C LB, ibas=",ibas
!      WRITE(*,*)"C LB, zir=",zir(1:nir)
      nsym=nir
!      WRITE(*,*)"C LB, nsym=",nsym
      slabel(1:nsym)=zir(1:nsym)
!      CALL sym_bas_irreps(ibas,.true.,nbpsy(1:nsym))
!      WRITE(*,*)"----Geometry and Symmetry Information----"
!      write(*,'(a,i3)')"geom= ",geom
!      write(*,'(a,a)')"sym group= ",grp_name
!      write(*,*)"slabel= ",slabel(1:nsym)
!      WRITE(*,'(a,i3)')"ibas= ",ibas
!      WRITE(*,'(a,i3)')"nbft= ",nbft
      WRITE(*,'(a,8i3)')"nbpsy= ",nbpsy(1:nsym)
!      WRITE(*,'(a,i3)')"sym_num_ops= ",sym_num_ops(geom)
!      WRITE(*,'(a,i3)')"nop= ",nop
C LB
      stitle(1)='AO integrals from NWChem'
      !nbpsy(1)=nbft

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
!       WRITE(*,*)'LB: ietype',ietype
!        WRITE(*,*)"LB before sifwh, info=",info(1:5)

cgk end
         call sifwh( aoints,  ntitle,  nsym,    nbft,
     &               ninfo,   nenrgy,  nmap,    stitle,   
     &               nbpsy,   slabel(1:nsym),  info, 
     &      byte_mb(k_bfnlab),  ietype,  energy,  imtype,  map, 
     &      ierr )
      WRITE(*,*)"C LB, O"
cgk debug
!       write(*,*)'gk: back from sifwh with ierr=',ierr
!        WRITE(*,*)"LB after sifwh, info=",info(1:5)
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
      write(*,*)'gk: processing integral type: ', integ_type
cgk end

         do jshell = 1, nshell
            do ishell =1, nshell
cgk debug
      write(*,*)'gk: ishell=',ishell,' jshell=',jshell
cgk end
 
               odoit = .true.
!               if (oskel) 
!    &           odoit = sym_shell_pair(ibas, ishell, jshell, q2)
!               WRITE(*,*)"C LB, odoit= ", odoit
!               WRITE(*,*)"C LB, q2= ", q2
 
               if (odoit) then
                  if (.not. bas_cn2bfr(ibas, ishell, ilo, ihi)) 
     $               call errquit('int_1e_ga: bas_cn2bfr ?', ibas,
     &                   BASIS_ERR)
                  idim = ihi - ilo + 1

                  call ifill(max1e,0,int_mb(k_ilab),1)
                  call ifill(max1e,0,int_mb(k_jlab),1)
      
                  if      (type .eq. 1) then !Kinetic Energy
                     call int_l1eke (ibas, jshell, ibas, ishell, 
     &                  thresh, int_mb(k_ilab),int_mb(k_jlab), max1e,
     $                  dbl_mb(k_buf), mem1, dbl_mb(k_scr), numints)

                  else if (type .eq. 2) then !Potential Energy
                     call int_l1epe (ibas, jshell, ibas, ishell, 
     &                  thresh, int_mb(k_ilab),int_mb(k_jlab), max1e,
     $                  dbl_mb(k_buf), mem1, dbl_mb(k_scr), numints)

                  else if (type .eq. 3) then !Overlap
c                    ECP is summed here, but COLUMBUS should not care.
                     call int_l1eov (ibas, ishell, ibas, jshell,
     &                  thresh, int_mb(k_ilab),int_mb(k_jlab), max1e,
     &                  dbl_mb(k_buf), mem1, dbl_mb(k_scr), numints)
                  end if
       WRITE(*,*)"LB, numints=",numints
c     
                  do i=1,numints

                     if(ibuf.eq.n1max) then !dump full record to file
                       numtot = numtot + ibuf
C LB
!                       WRITE(*,*)"C LB, calling sym_1int 1"
!                       CALL sym_1int(ibuf,nsoints,
!     &                               dbl_mb(k_sifval), clab,
!     &                               dbl_mb(k_SOval), cSOlab)
!                       WRITE(*,'(a,1i4)')"C LB, nsoints =",nsoints
!                       WRITE(*,*)"C LB, cSOlab ="
!                       WRITE(*,'(2i4)') cSOlab(1:2,1:nsoints)
C LB
                       WRITE(*,*)"C LB, calling sifew1"
!                       call sifew1(aoints, info, 2, nsoints, last, 
                       call sifew1(aoints, info, 2, ibuf, last, 
     &                   itypea, itypeb, ibvtyp, 
!     &                   dbl_mb(k_SOval), cSOlab, fcore, ibitv,
     &                   dbl_mb(k_sifval), clab, fcore, ibitv,
     &                   dbl_mb(k_sifbuf), nrec, ierr)
c                      ibuf on return has the number of unwritten
c                      integrals.  dbl_mb(k_sifval+0:(ibuf-1)) 
c                      contains the unwritten values, ditto for labels 
                       numtot = numtot - ibuf
                     endif
                      
                     IF ( ABS(dbl_mb(k_buf+i-1)) .gt. thresh) THEN
                       ibuf=ibuf+1
                       clab(1,ibuf)=int_mb(k_ilab+i-1)
                       clab(2,ibuf)=int_mb(k_jlab+i-1)
                       dbl_mb(k_sifval+ibuf-1)=dbl_mb(k_buf+i-1)      
                     ENDIF

                  enddo ! i

               end if
               ijshell = ijshell + 1
            end do !i_shell
         end do !j_shell
         last=nmsame
         numtot=numtot+ibuf
!         WRITE(*,*)"C LB, calling sym_1int 2, ibuf=",ibuf
!         CALL sym_1int(ibuf,nsoints,
!     &     dbl_mb(k_sifval), clab,
!     &     dbl_mb(k_SOval), cSOlab)
         WRITE(*,*)"C LB, calling sifew1"
!         call sifew1(aoints, info, 2, nsoints, last, 
         call sifew1(aoints, info, 2, ibuf, last, 
     &     itypea, itypeb, ibvtyp, 
!     &     dbl_mb(k_SOval), cSOlab, fcore, ibitv, 
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
     &  fcore, ibvtyp, ibitv,dbl_mb(k_SOval),cSOlab)

cgk debug
*      write(*,*)'gk: calling int_mom_sifs from int_1e_sifs'
*      write(*,*)'gk:  mem1=', mem1
cgk end
      call int_mom_sifs(ibas, oskel, aoints, nbft, l1rec,
     &  n1max, dbl_mb(k_sifbuf), dbl_mb(k_sifval), ninfo, info, 
     &  fcore, ibvtyp, ibitv,dbl_mb(k_SOval),cSOlab)
cgk debug
*      write(*,*)'gk: back from int_so_sifs in int_1e_sifs'
cgk end
c     
c     chop stack at first item allocated
c     
      if (.not. MA_pop_stack(l_SOval)) call errquit
     $     ('int_1e_sifs: pop failed at l_SOval', 0, GA_ERR)
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
c     LB type is never initialized, ocache never used
      !if(type.eq.6 .or. type.eq.7) then
      !   ocache=ocache_save
      !endif
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
