c---------------- get_chi_centers_ga() ------------- START
      subroutine get_chi_centers(chi_cntr, ! out
     &                           basis,    ! in  : basis    handle
     &                           nbf,      ! in  : nr basis functions
     &                           geom,     ! in  : geometry handle
     &                           mcenters) ! in  : nr. atoms
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none

#include "rtdb.fh" 
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "msgids.fh"
      integer basis,geom,nbf
      double precision chi_cntr(3,nbf) ! OUTPUT
      double precision cnt(3),valZ
      integer ictr,ic1,ic2,icset
      integer l,nprim,ncontr,isphere,nshbf
      integer mcenters,i,n1
      integer iniz,ifin
      character*16 at_tag  
      integer lo1(3),hi1(3),ld(2)
      logical status     
      n1=0
        do ictr=1,mcenters
           iniz=0
           ifin=0
         if (.not.bas_ce2cnr(basis,ictr,ic1,ic2))
     &       call errquit('Exiting in get_chi_centers_ga.',
     &                    11, BASIS_ERR)
         do icset = ic1,ic2       
c ----- get info about current contraction set      
          if (.not. bas_continfo(basis,icset,l,nprim,
     &         ncontr,isphere))
     &         call errquit('Exiting in get_chi_centers_ga.',
     &                       5, BASIS_ERR)
          nshbf=ncontr*(((l+1)*(l+2))/2)
          if(isphere.eq.1) then
            nshbf=ncontr*(2*l+1)
          endif
          if (iniz.eq.0) iniz=n1+1
          n1=n1+nshbf
         enddo ! end loop icset
         ifin= n1
         status=geom_cent_get(geom,ictr,at_tag,
     &                        cnt,valZ)
c if iniz=0, must be a bq or bare ecp
         if(iniz.ne.0) then
         do i=iniz,ifin
           chi_cntr(1,i)=cnt(1)
           chi_cntr(2,i)=cnt(2)
           chi_cntr(3,i)=cnt(3)
         enddo ! end loop i 
         endif
        enddo ! end loop ictr
      return
      end
      subroutine get_chi_centers_ga(g_chi_cntr, ! out
     &                              basis,      ! in  : basis handle
     &                              nbf,geom,mcenters)
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none

#include "rtdb.fh" 
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "msgids.fh"
      integer basis,geom,nbf
      integer g_chi_cntr(3) ! OUTPUT
      double precision cnt(3),valZ
      integer ictr,ic1,ic2,icset
      integer l,nprim,ncontr,isphere,nshbf
      integer mcenters,i,n1
      integer iniz,ifin
      character*16 at_tag  
      integer l_buf,k_buf
      integer lo1(3),hi1(3),ld(2)
      logical status     
c ----- allocate array to store centers ---- START  
       if(.not.MA_push_get(MT_DBL,3*nbf,'get_chi_centers_ga:buf',
     &                    l_buf,k_buf))
     $     call errquit('get_chi_centers_ga: ma failed',
     &                  3*nbf, MA_ERR) 
c ----- allocate array to store centers ---- END
        n1=0
        do ictr=1,mcenters
           iniz=0
           ifin=0
         if (.not.bas_ce2cnr(basis,ictr,ic1,ic2))
     &       call errquit('Exiting in get_chi_centers_ga.',
     &                    11, BASIS_ERR)
         do icset = ic1,ic2       
c ----- get info about current contraction set      
          if (.not. bas_continfo(basis,icset,l,nprim,
     &         ncontr,isphere))
     &         call errquit('Exiting in get_chi_centers_ga.',
     &                       5, BASIS_ERR)
          nshbf=ncontr*(((l+1)*(l+2))/2)
          if(isphere.eq.1) then
            nshbf=ncontr*(2*l+1)
          endif
          if (iniz.eq.0) iniz=n1+1
          n1=n1+nshbf
         enddo ! end loop icset
         ifin= n1
         status=geom_cent_get(geom,ictr,at_tag,
     &                        cnt,valZ)
c if iniz=0, must be a bq or bare ecp
         if(iniz.ne.0) then
         do i=iniz,ifin
           dbl_mb(k_buf      +i-1)=cnt(1)
           dbl_mb(k_buf+nbf  +i-1)=cnt(2)
           dbl_mb(k_buf+2*nbf+i-1)=cnt(3)
         enddo ! end loop i 
         endif
        enddo ! end loop ictr
c ----- store in g_chi_cntr() --- START
c       dbl_mb() ---> g_chi_cntr()
        do i=1,3
         ld(1)=nbf
         lo1(1)=1
         hi1(1)=nbf
         lo1(2)=i
         hi1(3)=i
         call nga_put(g_chi_cntr(i),
     &                lo1,hi1,dbl_mb(k_buf+(i-1)*nbf),ld)
        enddo ! end-loop-i
c ----- store in g_chi_cntr() --- END
c --- Free memory 
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('get_chi_centers_ga: pop failed', 0, GA_ERR)
      return
      end
c---------------- get_chi_centers_ga() ------------- END
      subroutine get_3rdterm_R(g_N,     ! to be scaled
     &                         g_R,     ! scaling
     &                         ind_a,   ! from kab=123,231,312
     &                         ind_b,   ! from kab=123,231,312
     &                         g_tmp2,  ! scratch
     &                         g_N_scld)! output
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
      integer g_N,g_R(3),g_M
      integer g_tmp2,g_N_scld
      integer ind_a,ind_b
      call ga_copy(g_N,g_tmp2)
      call ga_copy(g_N,g_N_scld)
      call ga_scale_cols(g_tmp2,g_R(ind_b))   ! R_{nu,b} g_N
      call ga_scale_rows(g_tmp2,g_R(ind_a))   ! R_{mu,a} [R_{nu,b} g_N] -> g_tmp2
      call ga_scale_cols(g_N_scld,g_R(ind_a)) ! R_{nu,a} g_N
      call ga_scale_rows(g_N_scld,g_R(ind_b)) ! R_{mu,b} [R_{nu,a} g_N] -> g_N_scld
      call ga_add(1.0d0,g_tmp2,-1.0d0,g_N_scld,g_N_scld)
      return
      end

      subroutine get_scld_A(g_A,  ! ga-arr to scale - OUT
     &                      g_R,  ! scaling arr
     &                      g_tmp)! scratch arr 
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
      integer g_A,g_R
      integer g_tmp
      integer nbf
c     Purpose: Compute R_{nu} U_{munu} - R_{mu} U_{munu}
c              g_R ->  R_{mu}
c              g_A -> U_{munu}
      call ga_copy(g_A,g_tmp)
      call ga_scale_cols(g_A  ,g_R)
      call ga_scale_rows(g_tmp,g_R)
      call ga_add(1.0d0,g_A,-1.0d0,g_tmp,g_A)
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++ READ/WRITE NMR-ZORA data +++++++++++++ START
c Note.- Using modified versions of
c        dft_zora_read() and dft_zoraNMR_write()
c        --> located in dft_zora_utils.F
czora...Write out the zora NMR shieldings to disk

      logical function dft_zoraNMR_write(
     &           filename, ! in: filename
     &       type_nmrdata, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                nbf, ! in: number of basis functions
     &              nlist, ! in: number of selected atoms
     &            g_AtNr1, ! in: list of atoms to calc. shieldings
     &              g_dia, ! in: dia   tensor
     &            g_para1, ! in: para1 tensor
     &              g_h01, ! in: h01 AO matrix
     &              g_Fji) ! in: Perturbed Fock matrix
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer nbf               ! [input] No. of functions in basis
      integer nlist ! nr. slc atoms
      integer g_AtNr1, ! in: list of atoms to calc. shieldings 
     &        g_dia,   ! in: dia   tensor
     &        g_para1, ! in: para1 tensor
     &        g_h01,   ! in: h01 AO matrix
     &        g_Fji
      integer unitno
      parameter (unitno = 77)
      integer l_AtNr,k_AtNr,
     &        l_tens,k_tens,
     &        l_h01 ,k_h01,
     &        l_Fji ,k_Fji
      integer ok, iset, i, j
      integer inntsize
      integer n9,nxyz,nh01,nFji,type_nmrdata
      integer alo(3),ahi(3),ld(2)
      nxyz=3 ! =x,y,z
      l_tens = -1   ! An invalid MA handle
      l_h01  = -1   ! An invalid MA handle
      l_Fji  = -1   ! An invalid MA handle
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Open the file
       open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c     Write out the number of sets and basis functions
       write(unitno, err=1001) nbf
       write(unitno, err=1001) nxyz
       write(unitno, err=1001) nlist
c     Allocate the temporary buffer
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       if (.not. ma_alloc_get(mt_dbl,nlist,'dft_zoraNMR_write',
     &                     l_AtNr,k_AtNr))
     $  call errquit('dft_zoraNMR_write: ma failed', 
     &               nlist, MA_ERR)
      n9=nxyz*nxyz*nlist
       if (.not. ma_alloc_get(mt_dbl,n9,'dft_zoraNMR_write',
     &                     l_tens,k_tens))
     $  call errquit('dft_zoraNMR_write: ma failed', 
     &               n9, MA_ERR)
      nh01=nbf*nbf*nxyz*nlist
       if (.not. ma_alloc_get(mt_dbl,nh01,'dft_zoraNMR_write',
     &                     l_h01,k_h01))
     $  call errquit('dft_zoraNMR_write: ma failed', 
     &               nh01, MA_ERR)
      nFji=nbf*nbf*nxyz
       if (.not. ma_alloc_get(mt_dbl,nFji,'dft_zoraNMR_write',
     &                     l_Fji,k_Fji))
     $  call errquit('dft_zoraNMR_write: ma failed', 
     &               nFji, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.2) then
       call ga_get(g_AtNr1,1,1,1,nlist,dbl_mb(k_AtNr),1)
       call swrite(unitno,dbl_mb(k_AtNr),nlist)
      endif
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=nlist
      ld(1)=3
      ld(2)=3
      call nga_get(g_dia,alo,ahi,dbl_mb(k_tens),ld)
      call swrite(unitno,dbl_mb(k_tens),n9)
      if (type_nmrdata .ne. 1 .and. type_nmrdata .ne. 2 .and.
     &    type_nmrdata .ne. 3) then
       write(*,*) 'Error in dft_zoraNMR_read:: ',
     &            'type_nmrdata not correct.',
     &            'It should be equal 1 or 2 or 3'
      endif
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.3) then
       call nga_get(g_para1,alo,ahi,dbl_mb(k_tens),ld)
       call swrite(unitno,dbl_mb(k_tens),n9)
      endif
      alo(1)=1
      ahi(1)=nbf
      alo(2)=1
      ahi(2)=nbf
      alo(3)=1
      ahi(3)=nxyz*nlist
      ld(1)=nbf
      ld(2)=nbf
      call nga_get(g_h01,alo,ahi,dbl_mb(k_h01),ld)
      call swrite(unitno,dbl_mb(k_h01),nh01)
      alo(1)=1
      ahi(1)=nbf
      alo(2)=1
      ahi(2)=nbf
      alo(3)=1
      ahi(3)=nxyz
      ld(1)=nbf
      ld(2)=nbf
      call nga_get(g_Fji,alo,ahi,dbl_mb(k_Fji),ld)
      call swrite(unitno,dbl_mb(k_Fji),nFji)
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_AtNr))
     $  call errquit('dft_zoraNMR_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_tens))
     $  call errquit('dft_zoraNMR_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_h01))
     $  call errquit('dft_zoraNMR_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_Fji))
     $  call errquit('dft_zoraNMR_write: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c     Close the file
      close(unitno,err=1002)
      ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraNMR_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote ZORA NMR data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_zoraNMR_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNMR_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNMR_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_zoraNMR_read(
     &           filename, ! in: filename
     &       type_nmrdata, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                nbf, ! in: number of basis functions
     &              nlist, ! in: number of selected atoms
     &            g_AtNr1, ! out: list of atoms to calc. shieldings
     &              g_dia, ! out: dia   tensor
     &            g_para1, ! out: para1 tensor
     &              g_h01, ! out: h01 AO matrix
     &              g_Fji) ! out: Perturbed Fock matrix
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer nbf               ! [input] No. of functions in basis
      integer nlist ! nr. slc atoms
      integer g_AtNr1, ! in: list of atoms to calc. shieldings 
     &        g_dia,   ! in: dia   tensor
     &        g_para1, ! in: para1 tensor
     &        g_h01,   ! in: h01 AO matrix
     &        g_Fji
      integer unitno
      parameter (unitno = 77)
      integer l_AtNr,k_AtNr,
     &        l_tens,k_tens,
     &        l_h01 ,k_h01,
     &        l_Fji ,k_Fji
      integer n9,nxyz,nh01,nFji,type_nmrdata
      integer alo(3),ahi(3),ld(2)
      integer ok,inntsize
      integer nxyz_read,nlist_read,
     &        nbf_read
      if (type_nmrdata .ne. 1 .and. type_nmrdata .ne. 2 .and.
     &    type_nmrdata .ne. 3) then
       write(*,*) 'Error in dft_zoraNMR_read::',
     &            ' type_nmrdata not correct.',
     &            'It should be equal 1 or 2 or 3'
       stop
      endif
      nxyz=3 ! =x,y,z   
c     Initialise to invalid MA handle
      l_tens = -1   ! An invalid MA handle
      l_h01  = -1   ! An invalid MA handle
      l_Fji  = -1   ! An invalid MA handle
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c---- Create GA arrays: g_AtNr1,g_dia,g_para1,g_h01,g_Fji --- START
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.2) then
         if (.not. ga_create(mt_dbl,1,nlist,
     &   'dft_zoraNMR_read: g_AtNr1',0,0,g_AtNr1)) 
     $   call errquit('gCSSR: g_AtNr1',0,GA_ERR)
        call ga_zero(g_AtNr1)   
      endif
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) = nlist
      if (.not.nga_create(MT_DBL,3,ahi,'g_DIA matrix',alo,g_dia)) call 
     &    errquit('dft_zoraNMR_read: nga_create failed g_dia',
     &            0,GA_ERR)
      call ga_zero(g_dia)
      if (type_nmrdata .eq. 1 .or. type_nmrdata .eq. 3) then
       if (.not.nga_create(MT_DBL,3,ahi,'gPAR1 matrix',
     &                    alo,g_para1)) 
     &    call errquit('dft_zoraNMR_read: nga_create failed gpar1',
     &            0,GA_ERR)
       call ga_zero(g_para1)
      endif
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = nxyz
      if (.not.nga_create(MT_DBL,3,ahi,'Fji matrix',alo,g_Fji)) 
     &    call 
     &    errquit('dft_zoraNMR_read: nga_create failed g_Fji',
     &            0,GA_ERR)
      call ga_zero(g_Fji)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*nlist
      if (.not.nga_create(MT_DBL,3,ahi,'h01 matrix',alo,g_h01)) call 
     &    errquit('dft_zoraNMR_read: nga_create failed g_h01_num',
     &            0,GA_ERR)
      call ga_zero(g_h01)
c---- Create GA arrays: g_AtNr1,g_dia,g_para1,g_h01,g_Fji --- END
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read ZORA NMR data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) nxyz_read
       read(unitno, err=1001, end=1001) nlist_read 
c      Error checks
       if ((nxyz_read  .ne. nxyz) .or.
     &     (nbf_read   .ne. nbf)  .or.
     &     (nlist_read .ne. nlist) ) goto 1003
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       if (.not. ma_alloc_get(mt_dbl,nlist,'dft_zoraNMR_read',
     &                     l_AtNr,k_AtNr))
     $  call errquit('dft_zoraNMR_read: ma failed', 
     &               nlist, MA_ERR)
      n9=nxyz*nxyz*nlist
       if (.not. ma_alloc_get(mt_dbl,n9,'dft_zoraNMR_read',
     &                     l_tens,k_tens))
     $  call errquit('dft_zoraNMR_read: ma failed', 
     &               n9, MA_ERR)
      nh01=nbf*nbf*nxyz*nlist
       if (.not. ma_alloc_get(mt_dbl,nh01,'dft_zoraNMR_read',
     &                     l_h01,k_h01))
     $  call errquit('dft_zoraNMR_read: ma failed', 
     &               nh01, MA_ERR)
      nFji=nbf*nbf*nxyz
       if (.not. ma_alloc_get(mt_dbl,nFji,'dft_zoraNMR_read',
     &                     l_Fji,k_Fji))
     $  call errquit('dft_zoraNMR_read: ma failed', 
     &               nFji, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.2) then
       call sread(unitno,dbl_mb(k_AtNr),nlist)
       call ga_put(g_AtNr1,1,1,1,nlist,dbl_mb(k_AtNr),1)
      endif
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=nlist
      ld(1)=3
      ld(2)=3
      call sread(unitno,dbl_mb(k_tens),n9)
      call nga_put(g_dia,alo,ahi,dbl_mb(k_tens),ld)
      if (type_nmrdata .eq. 1 .or. type_nmrdata .eq. 3) then
       call sread(unitno,dbl_mb(k_tens),n9)
       call nga_put(g_para1,alo,ahi,dbl_mb(k_tens),ld)
      endif
      alo(1)=1
      ahi(1)=nbf
      alo(2)=1
      ahi(2)=nbf
      alo(3)=1
      ahi(3)=nxyz*nlist
      ld(1)=nbf
      ld(2)=nbf
      call sread(unitno,dbl_mb(k_h01),nh01)
      call nga_put(g_h01,alo,ahi,dbl_mb(k_h01),ld)
      alo(1)=1
      ahi(1)=nbf
      alo(2)=1
      ahi(2)=nbf
      alo(3)=1
      ahi(3)=nxyz
      ld(1)=nbf
      ld(2)=nbf
      call sread(unitno,dbl_mb(k_Fji),nFji)
      call nga_put(g_Fji,alo,ahi,dbl_mb(k_Fji),ld)
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_AtNr))
     $  call errquit('dft_zoraNMR_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_tens))
     $  call errquit('dft_zoraNMR_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_h01))
     $  call errquit('dft_zoraNMR_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_Fji))
     $  call errquit('dft_zoraNMR_read: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraNMR_read = ok .eq. 1
      return
 1000 write(6,*) 'dft_zoraNMR_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNMR_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 'dft_zshield_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNMR_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
c --- 05-02-11 ------- writing/reading A,B contributions ----- START
      logical function dft_zoraNMR_write_AB(
     &           filename, ! in: filename
     &       type_nmrdata, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                nbf, ! in: number of basis functions
     &              nlist, ! in: number of selected atoms
     &            g_AtNr1, ! in: list of atoms to calc. shieldings
     &              g_dia, ! in: dia A,B tensor
     &            g_para1) ! in: par A,B tensor
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer nbf               ! [input] No. of functions in basis
      integer nlist ! nr. slc atoms
      integer g_AtNr1, ! in: list of atoms to calc. shieldings 
     &        g_dia,   ! in: dia   tensor
     &        g_para1  ! in: para1 tensor
      integer unitno
      parameter (unitno = 77)
      integer l_AtNr,k_AtNr,
     &        l_tens,k_tens
      integer ok, iset, i, j
      integer inntsize
      integer n9,nxyz,nh01,nFji,type_nmrdata
      integer alo(3),ahi(3),ld(2)
      nxyz=3 ! =x,y,z
      l_tens = -1   ! An invalid MA handle
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Open the file
       open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c     Write out the number of sets and basis functions
       write(unitno, err=1001) nbf
       write(unitno, err=1001) nxyz
       write(unitno, err=1001) nlist
c     Allocate the temporary buffer
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       if (.not. ma_alloc_get(mt_dbl,nlist,'dft_zoraNMR_write',
     &                     l_AtNr,k_AtNr))
     $  call errquit('dft_zoraNMR_write: ma failed', 
     &               nlist, MA_ERR)
      n9=2*nxyz*nxyz*nlist
       if (.not. ma_alloc_get(mt_dbl,n9,'dft_zoraNMR_write',
     &                     l_tens,k_tens))
     $  call errquit('dft_zoraNMR_write: ma failed', 
     &               n9, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.2) then
       call ga_get(g_AtNr1,1,1,1,nlist,dbl_mb(k_AtNr),1)
       call swrite(unitno,dbl_mb(k_AtNr),nlist)
      endif
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=2*nlist
      ld(1)=3
      ld(2)=3
      call nga_get(g_dia,alo,ahi,dbl_mb(k_tens),ld)
      call swrite(unitno,dbl_mb(k_tens),n9)
      if (type_nmrdata .ne. 1 .and. type_nmrdata .ne. 2 .and.
     &    type_nmrdata .ne. 3) then
       write(*,*) 'Error in dft_zoraNMR_read:: ',
     &            'type_nmrdata not correct.',
     &            'It should be equal 1 or 2 or 3'
      endif
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.3) then
       call nga_get(g_para1,alo,ahi,dbl_mb(k_tens),ld)
       call swrite(unitno,dbl_mb(k_tens),n9)
      endif
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_AtNr))
     $  call errquit('dft_zoraNMR_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_tens))
     $  call errquit('dft_zoraNMR_write: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c     Close the file
      close(unitno,err=1002)
      ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraNMR_write_AB = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote ZORA NMR data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_zoraNMR_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNMR_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNMR_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_zoraNMR_read_AB(
     &           filename, ! in: filename
     &       type_nmrdata, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                nbf, ! in: number of basis functions
     &              nlist, ! in: number of selected atoms
     &            g_AtNr1, ! out: list of atoms to calc. shieldings
     &              g_dia, ! out: dia-A,B   tensor
     &            g_para1) ! out: par-A,B   tensor    
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer nbf               ! [input] No. of functions in basis
      integer nlist ! nr. slc atoms
      integer g_AtNr1, ! in: list of atoms to calc. shieldings 
     &        g_dia,   ! in: dia   tensor
     &        g_para1  ! in: para1 tensor
      integer unitno
      parameter (unitno = 77)
      integer l_AtNr,k_AtNr,
     &        l_tens,k_tens
      integer n9,nxyz,nh01,nFji,type_nmrdata
      integer alo(3),ahi(3),ld(2)
      integer ok,inntsize
      integer nxyz_read,nlist_read,
     &        nbf_read
      if (type_nmrdata .ne. 1 .and. type_nmrdata .ne. 2 .and.
     &    type_nmrdata .ne. 3) then
       write(*,*) 'Error in dft_zoraNMR_read::',
     &            ' type_nmrdata not correct.',
     &            'It should be equal 1 or 2 or 3'
       stop
      endif
      nxyz=3 ! =x,y,z   
c     Initialise to invalid MA handle
      l_tens = -1   ! An invalid MA handle
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c---- Create GA arrays: g_AtNr1,g_dia,g_para1,g_h01,g_Fji --- START
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.2) then
         if (.not. ga_create(mt_dbl,1,nlist,
     &   'dft_zoraNMR_read: g_AtNr1',0,0,g_AtNr1)) 
     $   call errquit('gCSSR: g_AtNr1',0,GA_ERR)
        call ga_zero(g_AtNr1)   
      endif
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) = 2*nlist
      if (.not.nga_create(MT_DBL,3,ahi,'g_DIA matrix',alo,g_dia)) call 
     &    errquit('dft_zoraNMR_read: nga_create failed g_dia',
     &            0,GA_ERR)
      call ga_zero(g_dia)
      if (type_nmrdata .eq. 1 .or. type_nmrdata .eq. 3) then
       if (.not.nga_create(MT_DBL,3,ahi,'gPAR1 matrix',
     &                    alo,g_para1)) 
     &    call errquit('dft_zoraNMR_read: nga_create failed gpar1',
     &            0,GA_ERR)
       call ga_zero(g_para1)
      endif
c---- Create GA arrays: g_AtNr1,g_dia,g_para1,g_h01,g_Fji --- END
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read ZORA NMR data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) nxyz_read
       read(unitno, err=1001, end=1001) nlist_read 
c      Error checks
       if ((nxyz_read  .ne. nxyz) .or.
     &     (nbf_read   .ne. nbf)  .or.
     &     (nlist_read .ne. nlist) ) goto 1003
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       if (.not. ma_alloc_get(mt_dbl,nlist,'dft_zoraNMR_read',
     &                     l_AtNr,k_AtNr))
     $  call errquit('dft_zoraNMR_read: ma failed', 
     &               nlist, MA_ERR)
      n9=2*nxyz*nxyz*nlist
       if (.not. ma_alloc_get(mt_dbl,n9,'dft_zoraNMR_read',
     &                     l_tens,k_tens))
     $  call errquit('dft_zoraNMR_read: ma failed', 
     &               n9, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
      if (type_nmrdata.eq.1 .or. type_nmrdata.eq.2) then
       call sread(unitno,dbl_mb(k_AtNr),nlist)
       call ga_put(g_AtNr1,1,1,1,nlist,dbl_mb(k_AtNr),1)
      endif
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=2*nlist
      ld(1)=3
      ld(2)=3
      call sread(unitno,dbl_mb(k_tens),n9)
      call nga_put(g_dia,alo,ahi,dbl_mb(k_tens),ld)
      if (type_nmrdata .eq. 1 .or. type_nmrdata .eq. 3) then
       call sread(unitno,dbl_mb(k_tens),n9)
       call nga_put(g_para1,alo,ahi,dbl_mb(k_tens),ld)
      endif
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_AtNr))
     $  call errquit('dft_zoraNMR_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_tens))
     $  call errquit('dft_zoraNMR_read: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraNMR_read_AB = ok .eq. 1
      return
 1000 write(6,*) 'dft_zoraNMR_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNMR_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 'dft_zshield_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNMR_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
c --- 05-02-11 ------- writing/reading A,B contributions ----- END
c +++++++++++ READ/WRITE NMR-ZORA data +++++++++++++ END
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ========================================================
c =========== READ/WRITE CPHF (g_rhs) data ==========START
      logical function dft_zoraCPHF_write(
     &           filename, ! in: filename
     &           npol,     ! in: nr polarization
     &           nocc,     ! in: nr occupied MOs
     &           nvirt,    ! in: nr virtual  MOs
     &           nbf,      ! in: nr basis functions
     &           vectors,  ! in: MOs
     &           g_rhs0,   ! in: (ntot,3)       GA matrix
     &           g_rhs)    ! in: (nocc*nvirt,3) GA matrix
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer npol,nbf,
     &        nocc(npol),nvirt(npol),
     &        vectors(npol),
     &        ispin,ntot,
     &        g_rhs0,g_rhs
      integer unitno
      parameter (unitno = 77)
      integer l_rhs0,k_rhs0,
     &        l_rhs,k_rhs,
     &        l_mo,k_mo
      integer ok,i,j
      integer inntsize
      integer nrhs,nxyz

      nxyz=3 ! =x,y,z
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Open the file
       open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000)
c     Write out the number of sets and basis functions
       write(unitno, err=1001) npol
       do i=1,npol
        write(unitno, err=1001) nocc(i)
       enddo
       do i=1,npol
        write(unitno, err=1001) nvirt(i)
       enddo
       write(unitno, err=1001) nxyz
       write(unitno, err=1001) nbf
c     Allocate the temporary buffer
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
c ----- Add MOs in file ----- START
       if (.not. ma_alloc_get(
     &        mt_dbl,nbf,'dft_zoraNLMO_writehyp',
     &        l_mo,k_mo))
     $  call errquit('dft_zoraCPHF_write: k_mo failed', 
     &               nbf,MA_ERR)
        do i=1,npol
         do j=1,nbf
         call dcopy(nbf,0.0d0,0,dbl_mb(k_mo),1) ! reset
         call ga_get(vectors(i),1,nbf,j,j,dbl_mb(k_mo),nbf)
         call swrite(unitno,dbl_mb(k_mo),nbf)         
         enddo ! end-loop-j
        enddo ! end-loop-i
c ----- Add MOs in file ----- END
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
       ntot=0
       do ispin=1,npol
         ntot=ntot+nocc(ispin)*nocc(ispin)
       enddo
       write(unitno, err=1001) ntot
       if (.not. ma_alloc_get(mt_dbl,ntot,
     &           'dft_zoraCPHF_write',l_rhs0,k_rhs0))
     $  call errquit('dft_zoraCPHF_write: ma failed', 
     &               ntot, MA_ERR)
       do i=1,nxyz
        call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs0),1)
        call ga_get(g_rhs0,1,ntot,i,i,dbl_mb(k_rhs0),ntot)
        call swrite(unitno,dbl_mb(k_rhs0),ntot)
       enddo
       ntot=0
       do ispin=1,npol
        ntot=ntot+nocc(ispin)*nvirt(ispin)
       enddo
       write(unitno, err=1001) ntot
       if (.not. ma_alloc_get(mt_dbl,ntot,
     &           'dft_zoraCPHF_write',l_rhs,k_rhs))
     $  call errquit('dft_zoraCPHF_write: ma failed', 
     &               ntot, MA_ERR)
       do i=1,nxyz
        call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
        call ga_get(g_rhs,1,ntot,i,i,dbl_mb(k_rhs),ntot)
        call swrite(unitno,dbl_mb(k_rhs),ntot)
       enddo
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_rhs0))
     $  call errquit('dft_zoraCPHF_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_rhs))
     $  call errquit('dft_zoraCPHF_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_mo))
     $  call errquit('dft_zoraCPHF_write: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c     Close the file
      close(unitno,err=1002)
      ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraCPHF_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote CPHF data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_zoraCPHF_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraCPHF_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraCPHF_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_zoraCPHF_read(
     &           filename, !  in: filename
     &           npol,     !  in: nr polarization
     &           nocc,     !  in: nr occupied MOs
     &           nvirt,    !  in: nr virtual  MOs
     &           nbf,      !  in: nr basis functions
     &           vectors,  ! out: MOs
     &           g_rhs0,   ! out: (ntot,3)       GA matrix
     &           g_rhs)    ! out: (nocc*nvirt,3) GA matrix
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer npol,nbf,
     &        nocc(npol),nvirt(npol),
     &        vectors(npol),
     &        ispin,ntot,
     &        g_rhs0,g_rhs
      integer unitno
      parameter (unitno = 77)
      integer l_rhs0,k_rhs0,
     &        l_rhs,k_rhs,
     &        l_mo,k_mo
      integer ok,i,j
      integer inntsize
      integer nrhs,nxyz
      integer npol_read,nxyz_read,ntot_read,
     &        nbf_read,
     &        nocc_read(2),nvirt_read(2)

      nxyz=3 ! =x,y,z   
c     Initialise to invalid MA handle
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read ZORA NMR data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) npol_read
       do i=1,npol_read
        read(unitno, err=1001, end=1001) nocc_read(i)
       enddo
       do i=1,npol_read
        read(unitno, err=1001, end=1001) nvirt_read(i)
       enddo
       read(unitno, err=1001, end=1001) nxyz_read
       read(unitno, err=1001, end=1001) nbf_read
c      Error checks
       if ((nxyz_read  .ne. nxyz) .or.
     &     (npol_read  .ne. npol) .or.
     &     (nbf_read   .ne. nbf) ) goto 1003
c ----- Read MOs ----- START
       if (.not. ma_alloc_get(
     &        mt_dbl,nbf,'dft_zoraNLMO_readhyp',
     &        l_mo,k_mo))
     $  call errquit('dft_zoraNLMO_readhyp: ma failed', 
     &               nbf,MA_ERR)
        do i=1,npol
         do j=1,nbf
          call dcopy(nbf,0.0d0,0,dbl_mb(k_mo),1) ! reset
          call sread(unitno,dbl_mb(k_mo),nbf)   
          call ga_put(vectors(i),1,nbf,j,j,dbl_mb(k_mo),nbf)
         enddo ! end-loop-j
        enddo ! end-loop-i
c ----- Read MOs ----- END
       ntot=0
       do ispin=1,npol
         ntot=ntot+nocc(ispin)*nocc(ispin)
       enddo
       read(unitno, err=1001, end=1001) ntot_read
       if (.not. ma_alloc_get(mt_dbl,ntot,
     &           'dft_zoraCPHF_write',l_rhs0,k_rhs0))
     $  call errquit('dft_zoraCPHF_write: ma failed', 
     &               ntot, MA_ERR)
       do i=1,nxyz
        call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs0),1)
        call sread(unitno,dbl_mb(k_rhs0),ntot)
        call ga_put(g_rhs0,1,ntot,i,i,dbl_mb(k_rhs0),ntot)
       enddo
       ntot=0
       do ispin=1,npol
        ntot=ntot+nocc(ispin)*nvirt(ispin)
       enddo
       read(unitno, err=1001, end=1001) ntot_read
       if (.not. ma_alloc_get(mt_dbl,ntot,
     &           'dft_zoraCPHF_write',l_rhs,k_rhs))
     $  call errquit('dft_zoraCPHF_write: ma failed', 
     &               ntot, MA_ERR)
       do i=1,nxyz
        call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
        call sread(unitno,dbl_mb(k_rhs),ntot)
        call ga_put(g_rhs,1,ntot,i,i,dbl_mb(k_rhs),ntot)
       enddo
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_mo))       ! deallocate memory
     $  call errquit('dft_zoraCPHF_read: ma free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_rhs0))
     $  call errquit('dft_zoraCPHF_read: ma free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_rhs))
     $  call errquit('dft_zoraCPHF_read: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraCPHF_read = ok .eq. 1
      return
 1000 write(6,*) 'dft_zoraCPHF_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraCPHF_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 
     & 'dft_zoraCPHF_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraCPHF_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
c =========== READ/WRITE CPHF (g_rhs) data ==========END
c $Id$
c =========== READ/WRITE CPHF-1 (g_rhs) data ==========START
c To be used in aoresponse module: fiao_f1_movecs
      logical function dft_CPHF1_write(
     &           filename, ! in: filename
     &           npol,     ! in: nr polarization
     &           nocc,     ! in: nr occupied MOs
     &           nvirt,    ! in: nr virtual  MOs
     &           ncomp,    ! in: nr. components
     &           g_rhs_re, ! in: (nocc*nvirt,3)(ipm) GA matrix
     &           g_rhs_im, ! in: (nocc*nvirt,3)(ipm) GA matrix
     &           lifetime) ! in: =T if (RE,IM) =F if RE
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c
c Note.- nmo ne nbf if linear dependencies appear.

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer npol,
     &        nocc(npol),nvirt(npol),
     &        ispin,ntot,
     &        ipm,ncomp,
     &        g_rhs_re(ncomp),
     &        g_rhs_im(ncomp)
      integer unitno
      parameter (unitno = 77)
      integer l_rhs0,k_rhs0,
     &        l_rhs,k_rhs
      integer ok,i,j
      integer inntsize
      integer nrhs,nxyz
      logical lifetime

      nxyz=3 ! =x,y,z
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Open the file
       open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000)
c     Write out the number of sets and basis functions
       write(unitno, err=1001) npol
       do i=1,npol
        write(unitno, err=1001) nocc(i)
       enddo
       do i=1,npol
        write(unitno, err=1001) nvirt(i)
       enddo
       write(unitno, err=1001) nxyz
c     Allocate the temporary buffer
       ntot=0
       do ispin=1,npol
        ntot=ntot+nocc(ispin)*nvirt(ispin)
       enddo
       write(unitno, err=1001) ntot
       if (.not. ma_alloc_get(mt_dbl,ntot,
     &           'dft_CPHF_write',l_rhs,k_rhs))
     $  call errquit('dft_CPHF_write: ma failed', 
     &               ntot, MA_ERR)
       do ipm=1,ncomp
        do i=1,2*nxyz  ! write (g_b,g_rhs_sol)
         call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
         call ga_get(g_rhs_re(ipm),1,ntot,i,i,dbl_mb(k_rhs),ntot)
         call swrite(unitno,dbl_mb(k_rhs),ntot)
        enddo ! end-loop-i
        if (lifetime) then
        do i=1,2*nxyz ! write (g_b,g_rhs_sol)
         call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
         call ga_get(g_rhs_im(ipm),1,ntot,i,i,dbl_mb(k_rhs),ntot)
         call swrite(unitno,dbl_mb(k_rhs),ntot)
        enddo ! end-loop-i
        endif ! end-if-lifetime
       enddo ! end-loop-ipm
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_rhs))
     $  call errquit('dft_CPHF_write: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c     Close the file
      close(unitno,err=1002)
      ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_CPHF1_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' dft_CPHF_write: Wrote aoresponse g_rhs data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_CPHF_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_CPHF_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_CPHF_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_CPHF1_read(
     &           filename, !  in: filename
     &           npol,     !  in: nr polarization
     &           nocc,     !  in: nr occupied MOs
     &           nvirt,    !  in: nr virtual  MOs
     &           ncomp,    ! out: nr. components
     &           g_rhs_re, ! out: (nocc*nvirt,3)(ipm) GA matrix
     &           g_rhs_im, ! out: (nocc*nvirt,3)(ipm) GA matrix
     &           lifetime) ! out: =T if (RE,IM) =F if RE
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer npol,
     &        nocc(npol),nvirt(npol),
     &        ispin,ntot,
     &        ipm,ncomp,
     &        g_rhs_re(ncomp),
     &        g_rhs_im(ncomp)
      integer unitno
      parameter (unitno = 77)
      integer l_rhs,k_rhs
      integer ok,i,j
      integer inntsize
      integer nrhs,nxyz
      integer npol_read,nxyz_read,ntot_read,
     &        nocc_read(2),nvirt_read(2)
      logical lifetime
      nxyz=3 ! =x,y,z   
c     Initialise to invalid MA handle
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' dft_CPHF_read:Read aoresponse g_rhs data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) npol_read
       do i=1,npol_read
        read(unitno, err=1001, end=1001) nocc_read(i)
       enddo
       do i=1,npol_read
        read(unitno, err=1001, end=1001) nvirt_read(i)
       enddo
       read(unitno, err=1001, end=1001) nxyz_read
c      Error checks
       if ((nxyz_read  .ne. nxyz) .or.
     &     (npol_read  .ne. npol)) goto 1003
       ntot=0
       do ispin=1,npol
        ntot=ntot+nocc(ispin)*nvirt(ispin)
       enddo
       read(unitno, err=1001, end=1001) ntot_read
       if (.not. ma_alloc_get(mt_dbl,ntot,
     &           'dft_CPHF_read',l_rhs,k_rhs))
     $  call errquit('dft_CPHF_read: ma failed', 
     &               ntot, MA_ERR)

       do ipm=1,ncomp
        do i=1,nxyz ! skip 1st subspace
         call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
         call sread(unitno,dbl_mb(k_rhs),ntot)
        enddo 
        do i=nxyz+1,2*nxyz ! read 2nd subspace and copy to g_rhs_re
         call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
         call sread(unitno,dbl_mb(k_rhs),ntot)
         call ga_put(g_rhs_re(ipm),1,ntot,i,i,dbl_mb(k_rhs),ntot)
        enddo ! end-loop-i
        if (lifetime) then
        do i=1,nxyz ! skip 1st subspace
         call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
         call sread(unitno,dbl_mb(k_rhs),ntot)
        enddo ! end-loop-i
        do i=nxyz+1,2*nxyz ! read 2nd subspace and copy to g_rhs_im
         call dcopy(ntot,0.0d0,0,dbl_mb(k_rhs),1)
         call sread(unitno,dbl_mb(k_rhs),ntot)
         call ga_put(g_rhs_im(ipm),1,ntot,i,i,dbl_mb(k_rhs),ntot)
        enddo ! end-loop-i
        endif ! end-if-lifetime
       enddo ! end-loop-ipm
c
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
       if (.not. ma_free_heap(l_rhs))
     $  call errquit('dft_CPHF_read: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_CPHF1_read = ok .eq. 1
      return
 1000 write(6,*) 'dft_CPHF_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_CPHF_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 
     & 'dft_CPHF_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_CPHF_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
c 000000000000000000000000000000000000000000000000000000
c 000000000000000000000000000000000000000000000000000000
      logical function dft_CPHF2_write(
     &           filename, ! in: filename
     &           n,        ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,    ! in: nr. components
     &           nvec,     ! in: nr. of directions = 3
     &           n1,       ! in: =n*ncomp
     &           nsub,     ! in: last subspace index (nsub+1)= nr of subspaces stored
     &           nsub_file,! ub: subspace counter
     &           g_z1,     ! in: history matrix z
     &           g_Az1)    ! in: history matrix Az
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c
c Note.- nsub = cc * nvec (nvec=3=x,y,z) multiple of nvec
c        dim(g_z1) =(n1,maxsub)   maxsub=nvec*11 default
c        dim(g_Az1)=(n1,maxsub)   maxsub=nvec*11 default
c        nsub <= maxsub
c -->It will write only subscape nsub for (z1,Az1)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      character*255 filename_mini ! only to store nblocks
      integer g_z1,g_Az1
      integer unitno
      parameter (unitno = 77)
      integer l_dat,k_dat,
     &        l_z,k_z
      integer ok,i,j,nblock,nblock_file,idat
      integer inntsize,g_xre,g_xim,m1,iter
      integer n,ncomp,nvec,n1,nsub,nsub_file
      double precision val_re,val_im
      external conv2reim_z1 ! defined in ga_lkain_2cpl3.F

      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c +++++++++ store nsub +++++++++ START
      write(filename_mini,30) trim(filename),'_nblock'
 30   format(a,4a)
c      write(*,*) 'Writing mini:',filename_mini
      nblock=nsub/3+1
      nblock_file=nsub_file/3+1
c      write(*,2) nblock,nblock_file
c 2    format('Writing (nblock,nblock_file)=(',i5,',',i5,')')
       open(unitno, status='unknown', form='unformatted',
     $      file=filename_mini, err=1000)
       write(unitno, err=1001) nblock_file
       close(unitno,err=1002)  
c +++++++++ store nsub +++++++++ END
c     Open the file
       open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000,position='append')
c     Write out the number of sets and basis functions
       write(unitno, err=1001) n
       write(unitno, err=1001) ncomp
       write(unitno, err=1001) nvec
       write(unitno, err=1001) n1
       write(unitno, err=1001) nsub_file
c     Allocate the temporary buffer
       if (.not. ma_alloc_get(mt_dcpl,n1,
     &           'dft_CPHF2_write',l_z,k_z))
     $  call errquit('dft_CPHF2_write: ma failed', 
     &               n1, MA_ERR)
       if (.not. ma_alloc_get(mt_dbl,n1,
     &           'dft_CPHF2_write',l_dat,k_dat))
     $  call errquit('dft_CPHF2_write: ma failed', 
     &               n1, MA_ERR)
c ------- write g_z1 ---------------------- START
        do i=1,nvec 
         m1=(nblock-1)*nvec+i
         call ga_get(g_z1,1,n1,m1,m1,dcpl_mb(k_z),n1)
         do idat=1,n1
          val_re=dreal(dcpl_mb(k_z+idat-1))
          dbl_mb(k_dat+idat-1)=val_re
         enddo ! end-loop-idat
         call swrite(unitno,dbl_mb(k_dat),n1)
         do idat=1,n1
          val_im=dimag(dcpl_mb(k_z+idat-1))
          dbl_mb(k_dat+idat-1)=val_im
         enddo ! end-loop-idat
         call swrite(unitno,dbl_mb(k_dat),n1)
        enddo ! end-loop-i
c ------- write g_z1 ---------------------- END 
c ------- write g_Az1 ---------------------- START
        do i=1,nvec 
         m1=(nblock-1)*nvec+i
         call ga_get(g_Az1,1,n1,m1,m1,dcpl_mb(k_z),n1)
         do idat=1,n1
          val_re=dreal(dcpl_mb(k_z+idat-1))
          dbl_mb(k_dat+idat-1)=val_re
         enddo ! end-loop-idat
         call swrite(unitno,dbl_mb(k_dat),n1)
         do idat=1,n1
          val_im=dimag(dcpl_mb(k_z+idat-1))
          dbl_mb(k_dat+idat-1)=val_im
         enddo ! end-loop-idat
         call swrite(unitno,dbl_mb(k_dat),n1)
        enddo ! end-loop-i
c ------- write g_Az1 ---------------------- END
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not. ma_free_heap(l_dat))
     $  call errquit('dft_CPHF2_write: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_z))
     $  call errquit('dft_CPHF2_write: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
c     Close the file
      close(unitno,err=1002)
      ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_CPHF2_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
        iter=nsub/3+1 ! estimate iter from nsub
         write(6,22) filename(1:inp_strlen(filename)),
     &               iter,nsub,nsub_file             
 22      format(' dft_CPHF2_write: Wrote ',a,
     &           ' (iter,nsub,nsub_file)=(',
     &          i4,',',i5,',',i5,')')
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_CPHF2_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_CPHF2_write: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_CPHF2_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_CPHF2_read(
     &           filename, ! in: filename
     &           n,        ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,    ! in: nr. components
     &           nvec,     ! in: nr. of directions = 3
     &           n1,       ! in: =n*ncomp
     &           nsub,     ! ou: last subspace index (nsub+1)= nr of subspaces stored
     &           nsub_read,! ou: last subspace read from file
     &           maxsub,   ! in: max subspace of (g_z1,g_Az1)
     &           g_z1,     ! ou: history matrix z
     &           g_Az1)    ! ou: history matrix Az
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      character*255 filename_mini ! only to store nsub
      integer ivec,idat,
     &        g_z1,g_Az1
      integer unitno
      parameter (unitno = 77)
      integer l_zre,k_zre,
     &        l_zim,k_zim
      integer ok,i,j,m1,m2,
     &        imin,imax,iskip,nskip
      integer inntsize
      integer n,ncomp,nvec,n1,nsub,maxsub,iblock,
     &        n_read,n1_read,nvec_red,ncomp_read,
     &        nsub_read,nvec_read,
     &        nblocks_ga,nblocks_file
      double precision val_re,val_im
      double complex val_cmplx
      external conv2reim_z1 ! defined in ga_lkain_2cpl3.F

      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      imin=0
      imax=0
      if (ga_nodeid() .eq. 0) then ! ----- if-ga_nodeid-eq-0 ----- START
c +++++++++ read nsub +++++++++ START
      write(filename_mini,30) trim(filename),'_nblock'
 30   format(a,4a)
c      write(*,*) 'Reading mini:',filename_mini
       open(unitno, status='old', form='unformatted',
     $      file=filename_mini, err=1000)
       read(unitno, err=1001, end=1001) nblocks_file
       close(unitno,err=1002)  
       nblocks_ga=maxsub/nvec
       write(*,2) nblocks_file,nblocks_ga
 2     format('(nblocks_file,nblocks_ga)=(',i4,',',i4,')')
c --- Compare (nblocks_ga,nblocks_file)
       if (nblocks_file .le. nblocks_ga-1) then
        imin=1
        imax=nblocks_file
        nskip=0
       else
        imin=1
        imax=nblocks_ga-1
        nskip=nblocks_file-(nblocks_ga-1)
       endif
       write(*,3) imin,imax,nskip
 3     format('(imin,imax,nskip)=(',i4,',',i4,',',i4,')')
c +++++++++ read nsub +++++++++ END
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' dft_CPHF_read:Read aoresponse g_rhs data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
      if (.not.MA_Push_Get(mt_dbl,n1,'hessv jfacs',l_zre,k_zre))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  n1, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,n1,'hessv kfacs',l_zim,k_zim))
     &     call errquit('conv2complex: cannot allocate zim',
     &                  n1, MA_ERR)
c ======= skip blocks ======================= START
       do iskip=1,nskip
        read(unitno, err=1001, end=1001) n_read
        read(unitno, err=1001, end=1001) ncomp_read
        read(unitno, err=1001, end=1001) nvec_read
        read(unitno, err=1001, end=1001) n1_read
        read(unitno, err=1001, end=1001) nsub_read
c ------- skip g_z1 ---------------------- START
         do ivec=1,nvec 
          call sread(unitno,dbl_mb(k_zre),n1)
          call sread(unitno,dbl_mb(k_zim),n1)
         enddo ! end-loop-i
c ------- skip g_z1 ---------------------- END
c ------- skip g_Az1 ---------------------- START
         do ivec=1,nvec
          call sread(unitno,dbl_mb(k_zre),n1)
          call sread(unitno,dbl_mb(k_zim),n1)
         enddo ! end-loop-i
c ------- skip g_Az1 ---------------------- END
       enddo ! end-loop-iskip
c ======= skip blocks ======================= END
c ======= Loop in subspaces ===== START
       do iblock=imin,imax
        read(unitno, err=1001, end=1001) n_read
        read(unitno, err=1001, end=1001) ncomp_read
        read(unitno, err=1001, end=1001) nvec_read
        read(unitno, err=1001, end=1001) n1_read
        read(unitno, err=1001, end=1001) nsub_read
        write(*,14) iblock,nsub_read
 14     format('(iblock,nsub_read)=(',i4,',',i4,')')
c ------- read g_z1 ---------------------- START
        m1=(iblock-1)*nvec+1
        m2=m1+nvec-1
c        write(*,4) m1,m2,nvec
c 4      format('(m1,m2,nvec)=(',i4,',',i4,',',i4,')')
         do ivec=m1,m2 
         call dcopy(n1,0.0d0,0,dbl_mb(k_zre),1)
         call sread(unitno,dbl_mb(k_zre),n1)
         call dcopy(n1,0.0d0,0,dbl_mb(k_zim),1)
         call sread(unitno,dbl_mb(k_zim),n1)
          do idat=1,n1
           val_cmplx=dcmplx(dbl_mb(k_zre+idat-1),
     &                      dbl_mb(k_zim+idat-1))
           call ga_put(g_z1,idat,idat,ivec,ivec,val_cmplx,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
c ------- read g_z1 ---------------------- END
c ------- read g_Az1 ---------------------- START
         do ivec=m1,m2
          call dcopy(n1,0.0d0,0,dbl_mb(k_zre),1)
          call sread(unitno,dbl_mb(k_zre),n1)
          call dcopy(n1,0.0d0,0,dbl_mb(k_zim),1)
          call sread(unitno,dbl_mb(k_zim),n1)
          do idat=1,n1
           val_cmplx=dcmplx(dbl_mb(k_zre+idat-1),
     &                      dbl_mb(k_zim+idat-1))
           call ga_put(g_Az1,idat,idat,ivec,ivec,val_cmplx,1)
          enddo ! end-loop-idat
         enddo ! end-loop-i
c ------- read g_Az1 ---------------------- END
        enddo ! end-loop-iblock
        nsub=(imax-1)*3
        write(*,5) nsub,nblocks_file,nblocks_ga
 5      format('dft_CPHF2_read: (nsub,nblocks_file,nblocks_ga)=(',
     &         i12,',',i12,',',i12,')')
c ======= Loop in subspaces ===== END
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not.ma_pop_stack(l_zim))
     $  call errquit('dft_CPHF2_read: pop problem with l_zim',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_zre))
     $  call errquit('dft_CPHF2_read: pop problem with l_zre',
     &               555,MA_ERR)
c ----- Using ma_free_heap ------------ END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if ! ----- if-ga_nodeid-eq-0 ----- END
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_CPHF2_read = ok .eq. 1
      return
 1000 write(6,*) 'dft_CPHF2_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_CPHF2_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 
     & 'dft_CPHF2_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_CPHF2_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end

      logical function dft_CPHF2_read2fix(
     &           filename, ! in: filename
     &           n,        ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,    ! in: nr. components
     &           nvec,     ! in: nr. of directions = 3
     &           n1,       ! in: =n*ncomp
     &           nsub,     ! ou: last subspace index (nsub+1)= nr of subspaces stored
     &           maxsub,   ! in: max subspace of (g_z1,g_Az1)
     &           g_z1,     ! ou: history matrix z
     &           g_Az1)    ! ou: history matrix Az
c
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      character*255 filename_mini   ! only to store nsub
      character*255 filename_mini_1 ! only to store nsub
      character*255 filename_1      ! only to store nsub
      integer ivec,idat,
     &        g_z1,g_Az1
      integer unitno
      parameter (unitno = 77)

      integer unitno1
      parameter (unitno1 = 78)

      integer l_zre,k_zre,
     &        l_zim,k_zim
      integer ok,i,j,m1,m2,
     &        imin,imax,iskip,nskip,nsub1
      integer inntsize
      integer n,ncomp,nvec,n1,nsub,maxsub,iblock,
     &        n_read,n1_read,nvec_red,ncomp_read,
     &        nsub_read,nvec_read,
     &        nblocks_ga,nblocks_file,
     &        ntotblock_true
      double complex val_cmplx
      external conv2reim_z1 ! defined in ga_lkain_2cpl3.F
c 00000000000000000000000000
      ntotblock_true=75 ! FePc
c 00000000000000000000000000
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
      write(filename_1,31) trim(filename),'_1'
 31   format(a,2a)
      write(filename_mini_1,32) trim(filename),'_nblock_1'
 32   format(a,9a)
      write(*,33) filename_1(1:inp_strlen(filename_1)),
     &            filename_mini_1(1:inp_strlen(filename_mini_1))
 33   format('Creating files: filename_1=',a,
     &       ' filename_mini_1=',a)
c +++++++++ read nsub +++++++++ START
      write(filename_mini,30) trim(filename),'_nblock'
 30   format(a,4a)
c      write(*,*) 'Reading mini:',filename_mini
       open(unitno, status='unknown', form='unformatted',
     $      file=filename_mini, err=1000)
       read(unitno, err=1001, end=1001) nblocks_file
       close(unitno,err=1002)  

      write(*,*) 'Writing nblock=',ntotblock_true
       open(unitno1, status='unknown', form='unformatted',
     $      file=filename_mini_1, err=1000)
       write(unitno1, err=1001) ntotblock_true
       close(unitno1,err=1002)  

       nblocks_ga=maxsub/nvec
       write(*,2) nblocks_file,nblocks_ga
 2     format('(nblocks_file,nblocks_ga)=(',i4,',',i4,')')
c --- Compare (nblocks_ga,nblocks_file)
       if (nblocks_file .le. nblocks_ga) then
        imin=1
        imax=nblocks_file
        nskip=0
       else
        imin=1
        imax=nblocks_ga
        nskip=nblocks_file-nblocks_ga
       endif
       write(*,3) imin,imax,nskip
 3     format('(imin,imax,nskip)=(',i4,',',i4,',',i4,')')
c +++++++++ read nsub +++++++++ END
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' dft_CPHF_read:Read aoresponse g_rhs data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)

       open(unitno1, status='unknown', form='unformatted',
     $      file=filename_1, err=1000,position='append')

c      Read in some basics to check if they are consistent with the calculation
      if (.not.MA_Push_Get(mt_dbl,n1,'hessv jfacs',l_zre,k_zre))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  n1, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,n1,'hessv kfacs',l_zim,k_zim))
     &     call errquit('conv2complex: cannot allocate zim',
     &                  n1, MA_ERR)
c ======= Loop in subspaces ===== START
       nsub1=0
       do iblock=1,ntotblock_true
        read(unitno, err=1001, end=1001) n_read
        read(unitno, err=1001, end=1001) ncomp_read
        read(unitno, err=1001, end=1001) nvec_read
        read(unitno, err=1001, end=1001) n1_read
        read(unitno, err=1001, end=1001) nsub_read

        write(unitno1, err=1001) n_read
        write(unitno1, err=1001) ncomp_read
        write(unitno1, err=1001) nvec_read
        write(unitno1, err=1001) n1_read
        write(unitno1, err=1001) nsub1
        nsub1=nsub1+nvec
        write(*,14) iblock,nsub1,nsub_read
 14     format('(iblock,nsub1,nsub_read)=(',
     &         i4,',',i4,',',i4,')')
c ------- read g_z1 ---------------------- START
        m1=(iblock-1)*nvec+1
        m2=m1+nvec-1
c        write(*,4) m1,m2,nvec
c 4      format('(m1,m2,nvec)=(',i4,',',i4,',',i4,')')
         do ivec=m1,m2 
          call dcopy(n1,0.0d0,0,dbl_mb(k_zre),1)
          call sread(unitno,dbl_mb(k_zre),n1)
          call dcopy(n1,0.0d0,0,dbl_mb(k_zim),1)
          call sread(unitno,dbl_mb(k_zim),n1)
          call swrite(unitno1,dbl_mb(k_zre),n1)
          call swrite(unitno1,dbl_mb(k_zim),n1)
         enddo ! end-loop-i
c ------- read g_z1 ---------------------- END
c ------- read g_Az1 ---------------------- START
         do ivec=m1,m2
          call dcopy(n1,0.0d0,0,dbl_mb(k_zre),1)
          call sread(unitno,dbl_mb(k_zre),n1)
          call dcopy(n1,0.0d0,0,dbl_mb(k_zim),1)
          call sread(unitno,dbl_mb(k_zim),n1)
          call swrite(unitno1,dbl_mb(k_zre),n1)
          call swrite(unitno1,dbl_mb(k_zim),n1)
         enddo ! end-loop-i
c ------- read g_Az1 ---------------------- END
        enddo ! end-loop-iblock
c ======= Loop in subspaces ===== END
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
      if (.not.ma_pop_stack(l_zim))
     $  call errquit('dft_CPHF2_read: pop problem with l_zim',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_zre))
     $  call errquit('dft_CPHF2_read: pop problem with l_zre',
     &               555,MA_ERR)
c ----- Using ma_free_heap ------------ END
c      Close the file
       close(unitno,err=1002)
       close(unitno1,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_CPHF2_read2fix = ok .eq. 1
      return
 1000 write(6,*) 'dft_CPHF2_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_CPHF2_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 
     & 'dft_CPHF2_read: file inconsistent with calculation',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_CPHF2_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
c =========== READ/WRITE CPHF-1 g_rhs) data ==========END
