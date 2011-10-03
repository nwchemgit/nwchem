c---------------- get_chi_centers_ga() ------------- START
      subroutine get_chi_centers(chi_cntr, ! out
     &                           basis,    ! in  : basis    handle
     &                           nbf,      ! in  : nr basis functions
     &                           geom,     ! in  : geometry handle
     &                           mcenters) ! in  : nr. atoms
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
      integer iniz(mcenters),ifin(mcenters)
      character*16 at_tag  
      integer lo1(3),hi1(3),ld(2)
      logical status     
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
      n1=0
        do ictr=1,mcenters
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
          if (iniz(ictr).eq.0) iniz(ictr)=n1+1
          n1=n1+nshbf
         enddo ! end loop icset
         ifin(ictr)= n1
         status=geom_cent_get(geom,ictr,at_tag,
     &                        cnt,valZ)
         do i=iniz(ictr),ifin(ictr)
           chi_cntr(1,i)=cnt(1)
           chi_cntr(2,i)=cnt(2)
           chi_cntr(3,i)=cnt(3)
         enddo ! end loop i 
        enddo ! end loop ictr
      return
      end
      subroutine get_chi_centers_ga(g_chi_cntr, ! out
     &                              basis,      ! in  : basis handle
     &                              nbf,geom,mcenters)
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
      integer iniz(mcenters),ifin(mcenters)
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
        call ifill(mcenters,0,iniz,1)
        call ifill(mcenters,0,ifin,1)
        n1=0
        do ictr=1,mcenters
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
          if (iniz(ictr).eq.0) iniz(ictr)=n1+1
          n1=n1+nshbf
         enddo ! end loop icset
         ifin(ictr)= n1
         status=geom_cent_get(geom,ictr,at_tag,
     &                        cnt,valZ)
         do i=iniz(ictr),ifin(ictr)
           dbl_mb(k_buf      +i-1)=cnt(1)
           dbl_mb(k_buf+nbf  +i-1)=cnt(2)
           dbl_mb(k_buf+2*nbf+i-1)=cnt(3)
         enddo ! end loop i 
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
         call ga_get(vectors(i),1,nbf,j,j,dbl_mb(k_mo),1)
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
        call ga_get(g_rhs0,1,ntot,i,i,dbl_mb(k_rhs0),1)
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
        call ga_get(g_rhs,1,ntot,i,i,dbl_mb(k_rhs),1)
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
 22      format(/' Wrote ZORA CPHF data to ',a/)
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
          call ga_put(vectors(i),1,nbf,j,j,dbl_mb(k_mo),1)
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
        call ga_put(g_rhs0,1,ntot,i,i,dbl_mb(k_rhs0),1)
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
        call ga_put(g_rhs,1,ntot,i,i,dbl_mb(k_rhs),1)
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
