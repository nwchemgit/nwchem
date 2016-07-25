c
c     Calculates dipole spectrum between occupied orbitals
c     specified by the vspec keyword and all virtuals
c
c     Calculates dipole integrals once
c
      subroutine dft_vspec_so(rtdb, ao_bas_han,filename,nbf_ao,nbf_mo,
     &  nocc,ntotel,g_moso,k_eval,owstart,owend)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "msgtypesf.h"
#include "tcgmsg.fh"
#include "util.fh"
#include "inp.fh"
#include "util_params.fh"
c
      integer rtdb               ! [input]
      integer ao_bas_han         ! [input]
      character*(*) filename     ! [input]
      integer nbf_ao             ! [input]
      integer nbf_mo             ! [input]
      integer nocc(2)            ! [input]
      integer g_moso(2)          ! [input]
      integer k_eval(2)          ! [input]
      integer owstart            ! [input]
      integer owend              ! [input]
      integer g_dpx_ao        ! dipole x-component ao basis
      integer g_dpy_ao        ! dipole y-component ao basis
      integer g_dpz_ao        ! dipole z-component ao basis
      integer g_dpx_mo(2)     ! dipole x-component mo basis
      integer g_dpy_mo(2)     ! dipole y-component mo basis
      integer g_dpz_mo(2)     ! dipole z-component mo basis
      integer g_dp_wk         ! work array
      integer nvirt           ! number of virtuals
      integer nexcite         ! spectrum size
      integer twonexcite      ! nexcite*2
      integer mpone           ! num mo's + 1
      double precision ea,ei  ! orbital energies
      logical nodezero
      integer i,j,k,l,nn,jsav,kend,nindex
      integer nelec
      double precision esav, sort
      double precision two, three
      double precision tmpx1, tmpx2, tmpx12
      double precision tmpy1, tmpy2, tmpy12
      double precision tmpz1, tmpz2, tmpz12
      double precision fosc, exen
      integer ntotel
c
      double precision scal 
      parameter (two=2.d0,three=3.d0)
c
      integer ok
c
      integer unitno
      parameter (unitno = 77)
      integer inntsize,ddblsize
c
c     Preliminaries
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      nodezero=(ga_nodeid().eq.0)
c
c     Defaults
c
      if (owstart.le.0) owstart = 1 ! avoids failing
      if (owend.le.0) owend = 1  ! avoids failing
c
c     Get oscillator strength scale factor
c
      scal = 1.d0
      if (.not. rtdb_get(rtdb, 'vspec:scal', mt_dbl, 1, scal))
     &     scal = 1.d0
c
c     Determine size of spectrum
c
      nelec = ntotel  ! total number of electrons without fractional occs
      nvirt = nbf_mo-nelec
      nexcite = (owend-owstart+1)*nvirt
c
c     Check if there are no excitations present
c
      if (nexcite.eq.0) then
        if (nodezero) then
          write(luout, *) "No excitations present"
          write(luout, *) "nexcite :",nexcite
          return
        endif
      end if
c
c     Allocate ga memory for the AO and MO dipole components
c
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'ao x dipoles',nbf_ao,0,
     &  g_dpx_ao)) call errquit
     &  ('dft_vspec_so: error creating ao x dipole mat',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'ao y dipoles',nbf_ao,0,
     &  g_dpy_ao)) call errquit
     &  ('dft_vspec_so: error creating ao y dipole mat',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'ao z dipoles',nbf_ao,0,
     &  g_dpz_ao)) call errquit
     &  ('dft_vspec_so: error creating ao z dipole mat',0,GA_ERR)
      do i=1,2
        if (.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'mo dpx',nbf_mo,0,
     &      g_dpx_mo(i))) call errquit
     &      ('dft_vspec_so: error creating mo dipole mat',0,GA_ERR)
        if (.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'mo dpy',nbf_mo,0,
     &      g_dpy_mo(i))) call errquit
     &      ('dft_vspec_so: error creating mo dipole mat',0,GA_ERR)
        if (.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'mo dpz',nbf_mo,0,
     &      g_dpz_mo(i))) call errquit
     &      ('dft_vspec_so: error creating mo dipole mat',0,GA_ERR)
      enddo
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_mo,'mo dpx',nbf_mo,0,
     &  g_dp_wk)) call errquit
     &  ('dft_vspec_so: error creating mo work array',0,GA_ERR)
c
c     Compute dipole integrals in the AO basis
c
      call ga_sync()
      call ga_zero(g_dpx_ao)
      call ga_zero(g_dpy_ao)
      call ga_zero(g_dpz_ao)
      call int_dip_ga(ao_bas_han,ao_bas_han,g_dpx_ao,
     &  g_dpy_ao,g_dpz_ao)
c
c     Transform to MO basis - Right now this code assumes that
c     nbf_mo = 2*nbf_ao
c
      mpone=nbf_mo/2+1
c
c     X-component: Term 1
c
c     alpha real * M * alpha real
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,0.0d0,g_moso(1),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c     alpha imag * M * alpha imag
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(2),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c     beta real * M * beta real
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(1),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c     beta imag * M * beta imag
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(2),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
c     X-component: Term 2
c     
c     + alpha real^t * M * alpha imag
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,0.0d0,g_moso(1),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c     - alpha imag * M * alpha real
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,-1.0d0,g_moso(2),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c     + beta real * M * beta imag
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(1),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c     - beta imag * M * beta real
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpx_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,-1.0d0,g_moso(2),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpx_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
c     Y-component: Term 1
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,0.0d0,g_moso(1),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(2),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(1),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(2),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
c     Y-component: Term 2
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,0.0d0,g_moso(1),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,-1.0d0,g_moso(2),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(1),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpy_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,-1.0d0,g_moso(2),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpy_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
c     Z-component: Term 1
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,0.0d0,g_moso(1),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(2),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(1),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(2),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(1),1,nbf_mo,
     &  1,nbf_mo)
c
c     Z-component: Term 2
c     
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,0.0d0,g_moso(1),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c 
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),1,nbf_ao,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,-1.0d0,g_moso(2),1,nbf_mo,
     &  1,nbf_ao,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(2),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,1.0d0,g_moso(1),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
      call ga_matmul_patch('n','n',1.0d0,0.0d0,g_dpz_ao,1,nbf_ao,
     &  1,nbf_ao,g_moso(1),mpone,nbf_mo,1,nbf_mo,g_dp_wk,1,nbf_ao,
     &  1,nbf_mo)
      call ga_matmul_patch('t','n',1.0d0,-1.0d0,g_moso(2),1,nbf_mo,
     &  mpone,nbf_mo,g_dp_wk,1,nbf_ao,1,nbf_mo,g_dpz_mo(2),1,nbf_mo,
     &  1,nbf_mo)
c
c     Calculate Excitation energies and get arrays with dipole moments
c
      if (nodezero) then
        open(unitno, status='unknown', form='formatted',file=filename)
        write(unitno, 1000)
        write(unitno, 1400) nbf_ao,nbf_mo,owstart,owend,nelec,nexcite
      endif
c
      if (nodezero) then
         write (unitno, "(1x,a)")
     $        "      N   Occ  Virt       Energy [eV]        "//
     $        "<x> [au]          <y> [au]          <z> [au]    "//
     $        "     fosc [au] "
         write (unitno, "(a,i6)") "<START>",nexcite
      end if
c
      nindex=1
      do i=owstart,owend
        ei = dbl_mb(k_eval(1)+i-1) ! orbital energy of the initial state
        do j=nelec+1,nbf_mo
          ea = dbl_mb(k_eval(1)+j-1) ! orbital energy of the final state
          exen = (ea-ei)*cau2ev
          call ga_get(g_dpx_mo(1),i,i,j,j,tmpx1,1)
          call ga_get(g_dpx_mo(2),i,i,j,j,tmpx2,1)
          call ga_get(g_dpy_mo(1),i,i,j,j,tmpy1,1)
          call ga_get(g_dpy_mo(2),i,i,j,j,tmpy2,1)
          call ga_get(g_dpz_mo(1),i,i,j,j,tmpz1,1)
          call ga_get(g_dpz_mo(2),i,i,j,j,tmpz2,1)
c         tmpx1 = (tmpx1*tmpx1+tmpx2*tmpx2)*scal
c         tmpy1 = (tmpy1*tmpy1+tmpy2*tmpy2)*scal
c         tmpz1 = (tmpz1*tmpz1+tmpz2*tmpz2)*scal
          tmpx12 = (tmpx1*tmpx1)*scal   ! imag part doesn't contribute
          tmpy12 = (tmpy1*tmpy1)*scal   ! imag part doesn't contribute
          tmpz12 = (tmpz1*tmpz1)*scal   ! imag part doesn't contribute
          fosc = two/three*exen*(tmpx12+tmpy12+tmpz12)
          if (nodezero) then
            write (unitno, "(i8,2i6,5e18.6)")
     $         nindex, i, j, exen, tmpx12, tmpy12, tmpz12, fosc
          endif
          nindex=nindex+1
        enddo
      enddo
c
      if (nodezero) then
        write (unitno, "(a)") "<END>"
        close(unitno)
      end if
c
c     Broadcast status to other nodes
c
      ok = 1
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
 1000 format('== Spin-Orbit DFT Virtual Spectrum ==')
c
 1400 format(  'Number of AO basis functions: ',i5,
     1       /,'Total Number of MOs:          ',i5,
     3       /,'Occupied Win Start:           ',i5,
     4       /,'Occupied Win End:             ',i5,
     5       /,'Number of Occupied MOs:       ',i5,
     6       /,'Number of Excitations:        ',i5)
c
 1600 format(/,'Dipole spectrum for Koopman Excitation from orbital '
     &        ,i5,/)
 1700 format(3x,'N',7x,'Energy',11x,'<x>   ',11x,'<y>   ',11x,'<z>   ',
     &       11x,'fosc',/,'-------------------------------------------',
     &       '--------------------------------------------')
 1800 format(i5,2x,f12.7,5x,f12.7,5x,f12.7,5x,f12.7,5x,f12.7)
 2200 format(2x,i4,2x,f12.7)
c
c     Clean Up
c
      if (.not.ga_destroy(g_dpx_ao)) call errquit
     &  ('dft_vspec_so: failed to destroy g_dpx_ao',0,GA_ERR)
      if (.not.ga_destroy(g_dpy_ao)) call errquit
     &  ('dft_vspec_so: failed to destroy g_dpy_ao',0,GA_ERR)
      if (.not.ga_destroy(g_dpz_ao)) call errquit
     &  ('dft_vspec_so: failed to destroy g_dpz_ao',0,GA_ERR)
c
      do i=1,2
          if (.not.ga_destroy(g_dpx_mo(i))) call errquit
     &      ('dft_vspec_so: failed to destroy g_dpx_mo',0,GA_ERR)
          if (.not.ga_destroy(g_dpy_mo(i))) call errquit
     &      ('dft_vspec_so: failed to destroy g_dpy_mo',0,GA_ERR)
          if (.not.ga_destroy(g_dpz_mo(i))) call errquit
     &      ('dft_vspec_so: failed to destroy g_dpz_mo',0,GA_ERR)
      enddo
      if (.not.ga_destroy(g_dp_wk)) call errquit
     &  ('dft_vspec_so: failed to destroy g_dp_wk',0,GA_ERR)
c
      end 
c $Id$
