c
c     Calculates dipole spectrum between occupied orbitals
c     specified by the vspec keyword and all virtuals
c
      subroutine dft_vspec(rtdb, ao_bas_han,filename,ipol,nmo,nbf_ao,
     &  nocc, g_movecs,k_eval,owstart,owend)
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
c
      integer rtdb               ! [input]
      integer ao_bas_han         ! [input]
      character*(*) filename     ! [input] 
      integer ipol               ! [input]
      integer nmo(2)             ! [input]
      integer nbf_ao             ! [input]
      integer nocc(2)            ! [input]
      integer owstart(2)         ! [input]
      integer owend(2)           ! [input]
      integer g_movecs(2)        ! [input]
      integer k_eval(2)          ! [input]
      integer g_dpx_ao           ! dipole x-component ao basis
      integer g_dpy_ao           ! dipole y-component ao basis
      integer g_dpz_ao           ! dipole z-component ao basis
      integer g_dpx_mo(2)        ! dipole x-component mo basis
      integer g_dpy_mo(2)        ! dipole y-component mo basis
      integer g_dpz_mo(2)        ! dipole z-component mo basis
      integer g_dp_wk(2)         ! work array
      integer l_exen,k_exen      ! MA handles for excitation energies 
      integer l_amo,k_amo        ! MA handle of virtual orbs
      integer l_imo,k_imo        ! MA handle of occupied orb
      integer l_ex_occ, k_ex_occ
      integer l_ex_virt, k_ex_virt
      integer l_dpx,k_dpx
      integer l_dpy,k_dpy
      integer l_dpz,k_dpz
      integer l_fosc,k_fosc
      integer nvirt(2)           ! number of virtuals
      integer nexcite            ! spectrum size
      double precision ea,ei     ! orbital energies
      logical nodezero
      integer i,j,k,l,nn,jsav,kend,nindex
      double precision esav, sort
      integer sort_int
      integer ok
c
      integer unitno
      parameter (unitno = 77)
      integer inntsize,ddblsize
c
      double precision scal
      double precision two, three
      parameter (two=2.d0,three=3.d0)
      double precision au2ev
      parameter (au2ev=27.211396)
c
c     Preliminaries
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      nodezero=(ga_nodeid().eq.0)
c
c     Defaults
c
      do i=1,ipol
        if (owstart(i).le.0) owstart(i) = 1 ! avoids failing
        if (owend(i).le.0) owend(i) = 1  ! avoids failing
      enddo
C
C     Check that user didnt ask for excitations from virtuals, i.e.,
C     that owstart and owend both correspond to occupied orbitals.
C
      do i = 1,ipol
         if ( (owstart(i).gt.nocc(i)) .or.
     $        (owend(i).gt.nocc(i)) )
     $        call errquit("dft_vspec: cannot calculate excitations "//
     $        "from virtual orbitals, please adjust vspec input",0,0)
      enddo
C
C     Check that starting orbital less than ending one.
C
      do i = 1, ipol
         if (owstart(i).gt.owend(i))
     $        call errquit ("dft_vspec: invalid orbital range: "//
     $        "start > end?",0,0)
      enddo
c
c     Get oscillator strength scale factor
c
      scal = 1.d0
      if (.not. rtdb_get(rtdb, 'vspec:scal', mt_dbl, 1, scal))
     &     scal = 1.d0
c
c     Determine size of spectrum
c
      nexcite=0
      do i=1,ipol
        if (owstart(i).gt.0) then
          nexcite = nexcite+(owend(i)-owstart(i)+1)*(nmo(i)-nocc(i))
        endif
        nvirt(i) = nmo(i)-nocc(i)
      enddo
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
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'AO X dipole',nbf_ao,0,
     &  g_dpx_ao)) call errquit
     &  ('dft_vspec: error creating AO X dipole matrix',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'AO Y dipole',nbf_ao,0,
     &  g_dpy_ao)) call errquit
     &  ('dft_vspec: error creating AO Y dipole matrix',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'AO Z dipole',nbf_ao,0,
     &  g_dpz_ao)) call errquit
     &  ('dft_vspec: error creating AO Z dipole matrix',0,GA_ERR)
      do i=1,ipol
        if (owstart(i).gt.0) then
         if (.not.ga_create(mt_dbl,nmo(i),nmo(i),'MO X dipole',nmo(i),0,
     &      g_dpx_mo(i))) call errquit
     &      ('dft_vspec: error creating MO X dipole matrix',0,GA_ERR)
         if (.not.ga_create(mt_dbl,nmo(i),nmo(i),'MO Y dipole',nmo(i),0,
     &      g_dpy_mo(i))) call errquit
     &      ('dft_vspec: error creating MO Y dipole matrix',0,GA_ERR)
         if (.not.ga_create(mt_dbl,nmo(i),nmo(i),'MO Z dipole',nmo(i),0,
     &      g_dpz_mo(i))) call errquit
     &      ('dft_vspec: error creating MO Z dipole matrix',0,GA_ERR)
         if (.not.ga_create(mt_dbl,nbf_ao,nmo(i),'MO work',nbf_ao,0,
     &      g_dp_wk(i))) call errquit
     &      ('dft_vspec: error creating work array',0,GA_ERR)
        endif
      enddo
c
c     Create local arrays based on spectrum size
c
      if (.not.ma_push_get
     &  (mt_dbl,nexcite,'excitation energies',l_exen,k_exen))
     &  call errquit('dft_vspec: failed exen alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_int,nexcite,'occupied',l_ex_occ,k_ex_occ))
     &  call errquit('dft_vspec: failed ex_occ',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_int,nexcite,'virtuals',l_ex_virt,k_ex_virt))
     &  call errquit('dft_vspec: failed ex_virt',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite,'x dipoles',l_dpx,k_dpx))
     &  call errquit('dft_vspec: failed dpx alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite,'y dipoles',l_dpy,k_dpy))
     &  call errquit('dft_vspec: failed dpy alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite,'z dipoles',l_dpz,k_dpz))
     &  call errquit('dft_vspec: failed dpz alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite,'oscillator strength',l_fosc,k_fosc))
     &  call errquit('dft_vspec: failed fosc alloc',0,MA_ERR)
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
c     Transform to MO basis
c
      do i=1,ipol
        if (owstart(i).gt.0) then
c         x-comp
          call ga_dgemm('n','n',nbf_ao,nmo(i),nbf_ao,1.0d0,g_dpx_ao,
     &      g_movecs(i),0.d0,g_dp_wk(i))
          call ga_dgemm('t','n',nmo(i),nmo(i),nbf_ao,1.0d0,g_movecs(i),
     &      g_dp_wk(i),0.d0,g_dpx_mo(i))
c         y-comp
          call ga_dgemm('n','n',nbf_ao,nmo(i),nbf_ao,1.0d0,g_dpy_ao,
     &      g_movecs(i),0.d0,g_dp_wk(i))
          call ga_dgemm('t','n',nmo(i),nmo(i),nbf_ao,1.0d0,g_movecs(i),
     &      g_dp_wk(i),0.d0,g_dpy_mo(i))
c         z-comp
          call ga_dgemm('n','n',nbf_ao,nmo(i),nbf_ao,1.0d0,g_dpz_ao,
     &      g_movecs(i),0.d0,g_dp_wk(i))
          call ga_dgemm('t','n',nmo(i),nmo(i),nbf_ao,1.0d0,g_movecs(i),
     &      g_dp_wk(i),0.d0,g_dpz_mo(i))
        endif
      enddo
c
c     Calculate excitation energies
c
      nn = 0
      do i=1,ipol
c       account for case when beta but no alpha excitations 
        do j=owstart(i),owend(i)
          if (j.gt.0) then
            ei = dbl_mb(k_eval(i)+j-1)
            do k=nocc(i)+1,nmo(i)
              nn=nn+1
              ea = dbl_mb(k_eval(i)+k-1)
              dbl_mb(k_exen+nn-1) = ea-ei
              int_mb(k_ex_occ+nn-1) = j  !occupied MO for excitation i
              int_mb(k_ex_virt+nn-1) = k !virtual MO for excitation i
            enddo
          endif
        enddo
      enddo
c
c     Start print
c
      if (nodezero) then
        open(unitno, status='unknown', form='formatted',file=filename)
        write(unitno, 1000)
        if (ipol.eq.1) then
        write(unitno,1401) nbf_ao,nmo(1),owstart(1),
     &    owend(1),nocc(1),nexcite
        else
        write(unitno,1400) nbf_ao,nmo(1),nmo(2),owstart(1),owstart(2),
     &    owend(1),owend(2),nocc(1),nocc(2),nexcite
        end if
      endif
c
c     Get arrays with transition dipoles moments
c
      nindex=0
      do i=1,ipol
        do j=owstart(i),owend(i)
c         account for openshell with no alpha excitations
          if (j.gt.0) then
            call ga_get(g_dpx_mo(i),j,j,nocc(i)+1,nmo(i),
     &        dbl_mb(k_dpx+nindex),1)
            call ga_get(g_dpy_mo(i),j,j,nocc(i)+1,nmo(i),
     &        dbl_mb(k_dpy+nindex),1)
            call ga_get(g_dpz_mo(i),j,j,nocc(i)+1,nmo(i),
     &        dbl_mb(k_dpz+nindex),1)
            nindex=nindex+nvirt(i)
          endif
        enddo
      enddo
c
c     Calculate Oscillator Strengths
c
      do i=1,nexcite
        dbl_mb(k_fosc+i-1)=two/three*au2ev*dbl_mb(k_exen+i-1)*
     &    (dbl_mb(k_dpx+i-1)**2+dbl_mb(k_dpy+i-1)**2+
     &    dbl_mb(k_dpz+i-1)**2)
        dbl_mb(k_fosc+i-1)=scal*dbl_mb(k_fosc+i-1)  ! scale oscillator strength
      enddo
c
c     Sort by energy
c
      if (((owend(1)-owstart(1)).gt.0).or.(ipol.gt.1)) then
        do i=1,nexcite-1 
          esav = dbl_mb(k_exen+i-1)
          jsav = i
          do j=i+1,nexcite
            if (dbl_mb(k_exen+j-1) .lt. esav) then
              esav = dbl_mb(k_exen+j-1)
              jsav = j
            endif
          enddo
          sort=dbl_mb(k_exen+i-1)
          dbl_mb(k_exen+i-1)=dbl_mb(k_exen+jsav-1)
          dbl_mb(k_exen+jsav-1)=sort
c         sort transition moments and fosc
c         x
          sort=dbl_mb(k_dpx+i-1)
          dbl_mb(k_dpx+i-1)=dbl_mb(k_dpx+jsav-1)
          dbl_mb(k_dpx+jsav-1)=sort
c         y
          sort=dbl_mb(k_dpy+i-1)
          dbl_mb(k_dpy+i-1)=dbl_mb(k_dpy+jsav-1)
          dbl_mb(k_dpy+jsav-1)=sort
c         z
          sort=dbl_mb(k_dpz+i-1)
          dbl_mb(k_dpz+i-1)=dbl_mb(k_dpz+jsav-1)
          dbl_mb(k_dpz+jsav-1)=sort
c         fosc
          sort=dbl_mb(k_fosc+i-1)
          dbl_mb(k_fosc+i-1)=dbl_mb(k_fosc+jsav-1)
          dbl_mb(k_fosc+jsav-1)=sort
c         sort occupied orbitals
          sort=int_mb(k_ex_occ+i-1)
          int_mb(k_ex_occ+i-1)=int_mb(k_ex_occ+jsav-1)
          int_mb(k_ex_occ+jsav-1)=sort
c         sort virtual orbitals
          sort=int_mb(k_ex_virt+i-1)
          int_mb(k_ex_virt+i-1)=int_mb(k_ex_virt+jsav-1)
          int_mb(k_ex_virt+jsav-1)=sort
        enddo ! sorting
      endif
c
c     Output results
c
      if (nodezero) then
         write (unitno, "(1x,a)")
     $        "      N   Occ  Virt       Energy [eV]        "//
     $        "<x> [au]          <y> [au]          <z> [au]    "//
     $        "     fosc [au] "
         write (unitno, "(a,i6)") "<START>",nexcite
       do i=1,nexcite
          write (unitno, "(i8,2i6,5e18.6)")
     $         i, int_mb(k_ex_occ+i-1), int_mb(k_ex_virt+i-1),
     $         dbl_mb(k_exen+i-1)*au2ev,
     $         dbl_mb(k_dpx+i-1), dbl_mb(k_dpy+i-1), dbl_mb(k_dpz+i-1),
     $         dbl_mb(k_fosc+i-1)
       enddo
       write (unitno, "(a)") "<END>"
       close(unitno)
      endif
c
c     Broadcast status to other nodes
c
      ok = 1
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
c
 1000 format('== Scalar DFT Virtual Spectrum ==')
c
 1400 format(  'Number of basis functions in AO basis: ',i5,
     1       /,'Number of Alpha MOs:                   ',i5,
     2       /,'Number of Beta Mos:                    ',i5,
     3       /,'Occupied Win Start Alpha:              ',i5,
     4       /,'Occupied Win Start Beta:               ',i5,
     5       /,'Occupied Win End   Alpha:              ',i5,
     5       /,'Occupied Win End   Beta:               ',i5,
     6       /,'Number of Occupied Alpha MOs:          ',i5,
     6       /,'Number of Occupied Beta MOs :          ',i5,
     7       /,'Number of Excitations:                 ',i5)
c
 1401 format(  'Number of basis functions in AO basis: ',i5,
     1       /,'Number of Alpha MOs:                   ',i5,
     3       /,'Occupied Win Start Alpha:              ',i5,
     5       /,'Occupied Win End   Alpha:              ',i5,
     6       /,'Number of Occupied Alpha MOs:          ',i5,
     7       /,'Number of Excitations:                 ',i5)
c
 1600 format(/,x,'N',7x,'Energy',13x,'<x>',14x,'<y>',14x,'<z>',14x,
     &       'fosc')
c
 1700 format(i3,2x,f12.7,5x,f12.7,5x,f12.7,5x,f12.7,5x,f12.7)
c
c     Clean Up
c
      if (.not.ma_pop_stack(l_fosc)) call errquit
     &  ('dft_vspec: failed to deallocate fosc',0,MA_ERR)
      if (.not.ma_pop_stack(l_dpz)) call errquit
     &  ('dft_vspec: failed to deallocate dpz',0,MA_ERR)
      if (.not.ma_pop_stack(l_dpy)) call errquit
     &  ('dft_vspec: failed to deallocate dpy',0,MA_ERR)
      if (.not.ma_pop_stack(l_dpx)) call errquit
     &  ('dft_vspec: failed to deallocate dpx',0,MA_ERR)
      if (.not.ma_pop_stack(l_ex_virt)) call errquit
     &  ('dft_vspec: failed to deallocate ex_virt',0,MA_ERR)
      if (.not.ma_pop_stack(l_ex_occ)) call errquit
     &  ('dft_vspec: failed to deallocate ex_occ',0,MA_ERR)
      if (.not.ma_pop_stack(l_exen)) call errquit
     &  ('dft_vspec: failed to deallocate exen',0,MA_ERR)
c
      if (.not.ga_destroy(g_dpx_ao)) call errquit
     &  ('dft_vspec: failed to destroy g_dpx_ao',0,GA_ERR)
      if (.not.ga_destroy(g_dpy_ao)) call errquit
     &  ('dft_vspec: failed to destroy g_dpy_ao',0,GA_ERR)
      if (.not.ga_destroy(g_dpz_ao)) call errquit
     &  ('dft_vspec: failed to destroy g_dpz_ao',0,GA_ERR)
c
      do i=1,ipol
        if (owstart(i).gt.0) then
          if (.not.ga_destroy(g_dpx_mo(i))) call errquit
     &      ('dft_vspec: failed to destroy g_dpx_mo',0,GA_ERR)
          if (.not.ga_destroy(g_dpy_mo(i))) call errquit
     &      ('dft_vspec: failed to destroy g_dpy_mo',0,GA_ERR)
          if (.not.ga_destroy(g_dpz_mo(i))) call errquit
     &      ('dft_vspec: failed to destroy g_dpz_mo',0,GA_ERR)
          if (.not.ga_destroy(g_dp_wk(i))) call errquit
     &      ('dft_vspec: failed to destroy g_dp_wk',0,GA_ERR)
        endif
      enddo
c
      end 
c $Id: dft_vspec.F 22700 2012-08-12 05:20:52Z niri $
