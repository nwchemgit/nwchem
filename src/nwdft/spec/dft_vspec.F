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
#include "util_params.fh"
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
      integer nexcite_tot        ! total spectrum size
      integer nexcite_win        ! user specified window
      double precision ea,ei     ! orbital energies
      logical nodezero,polfromsos
      integer i,j,k,l,nn,jsav,kend,nindex
      double precision esav, sort
      integer sort_int
      integer l_scr,k_scr,l_mape,k_mape
c
      integer unitno
      parameter (unitno = 77)
      integer inntsize,ddblsize
c
      double precision scal
      double precision two, three
      parameter (two=2.d0,three=3.d0)
      double precision alpha(3,3) ! polariz
      double precision t0
      logical write_vspec
      t0=util_wallsec()

c
c     Preliminaries
c
      if (.not. rtdb_get(rtdb,'dft:write_vspec',mt_log,1,write_vspec))
     &     write_vspec=.true.
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      nodezero=(ga_nodeid().eq.0)
      if (.not. rtdb_get(rtdb,'dft:polfromsos',mt_log,1,polfromsos))
     &     polfromsos=.false.
      if(polfromsos) then
         write_vspec=.false.
         do i=1,ipol
            owstart(i)=1
            owend(i)=nocc(i)
         enddo
      endif
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
c
c     Check if the number of excitations have been specified
      nexcite_win = 0
      if (.not. rtdb_get(rtdb,'dft:nexcite_win',mt_int,1,nexcite_win))
     &         nexcite_win = 0
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
      nexcite_tot=0
      do i=1,ipol
        if (owstart(i).gt.0) then
          nexcite_tot = nexcite_tot +
     &     (owend(i)-owstart(i)+1)*(nmo(i)-nocc(i))
        endif
        nvirt(i) = nmo(i)-nocc(i)
      enddo
      if (nexcite_win.eq.0) nexcite_win = nexcite_tot
      if (nexcite_win.gt.nexcite_tot) then
        if (nodezero) then
          write(luout, *) "nexcite_win:", nexcite_win
          write(luout, *) "nexcite_tot:", nexcite_tot
        end if
        call errquit ("dft_vspec: invalid excitation window: "//
     $        "nexcite_win > nexcite_tot?",0,0)
      end if
c
c     Check if there are no excitations present
c
      if (nexcite_tot.eq.0) then
        if (nodezero) then
          write(luout, *) "No excitations present"
          write(luout, *) "nexcite_tot :",nexcite_tot
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
     &  (mt_dbl,nexcite_tot,'excitation energies',l_exen,k_exen))
     &  call errquit('dft_vspec: failed exen alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_int,nexcite_tot,'occupied',l_ex_occ,k_ex_occ))
     &  call errquit('dft_vspec: failed ex_occ',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_int,nexcite_tot,'virtuals',l_ex_virt,k_ex_virt))
     &  call errquit('dft_vspec: failed ex_virt',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite_tot,'x dipoles',l_dpx,k_dpx))
     &  call errquit('dft_vspec: failed dpx alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite_tot,'y dipoles',l_dpy,k_dpy))
     &  call errquit('dft_vspec: failed dpy alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite_tot,'z dipoles',l_dpz,k_dpz))
     &  call errquit('dft_vspec: failed dpz alloc',0,MA_ERR)
      if (.not.ma_push_get
     &  (mt_dbl,nexcite_tot,'oscillator strength',l_fosc,k_fosc))
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
         call dft_spece(ipol,owstart,owend,nocc,nbf_ao,nmo,
     W     dbl_mb(k_eval(1)),dbl_mb(k_exen),
     I     int_mb(k_ex_occ),int_mb(k_ex_virt))

c
c     Start print
c
      if (nodezero.and.write_vspec) then
        open(unitno, status='unknown', form='formatted',file=filename)
        write(unitno, 1000)
        if (ipol.eq.1) then
        write(unitno,1401) nbf_ao,nmo(1),owstart(1),
     &    owend(1),nocc(1),nexcite_tot
        else
        write(unitno,1400) nbf_ao,nmo(1),nmo(2),owstart(1),owstart(2),
     &    owend(1),owend(2),nocc(1),nocc(2),nexcite_tot
        end if
      endif
c
c      only nodezero from here on
      if(nodezero) then
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
      do i=1,nexcite_tot
        dbl_mb(k_fosc+i-1)=two/three*cau2ev*dbl_mb(k_exen+i-1)*
     &    (dbl_mb(k_dpx+i-1)**2+dbl_mb(k_dpy+i-1)**2+
     &    dbl_mb(k_dpz+i-1)**2)
        dbl_mb(k_fosc+i-1)=scal*dbl_mb(k_fosc+i-1)  ! scale oscillator strength
      enddo
c
c     Sort by energy
c
      if (((owend(1)-owstart(1)).gt.0).or.(ipol.gt.1)) then
         if (.not.ma_push_get
     &        (mt_int,nexcite_tot,'mape',l_mape,k_mape))
     &        call errquit('dft_vspec: failed mape alloc',0,MA_ERR)
         
         call dft_gsorte(nexcite_tot,
     E        dbl_mb(k_exen), int_mb(k_mape))

      if (.not.ma_push_get
     &     (mt_dbl,nexcite_tot,'scr',l_scr,k_scr))
     &  call errquit('dft_vspec: failed scr alloc',0,MA_ERR)
      call dft_sortmap(nexcite_tot,
     D     dbl_mb(k_scr),dbl_mb(k_fosc),int_mb(k_mape))
      call dft_sortmap(nexcite_tot,
     D     dbl_mb(k_scr),dbl_mb(k_dpx),int_mb(k_mape))
      call dft_sortmap(nexcite_tot,
     D     dbl_mb(k_scr),dbl_mb(k_dpy),int_mb(k_mape))
      call dft_sortmap(nexcite_tot,
     D     dbl_mb(k_scr),dbl_mb(k_dpz),int_mb(k_mape))
      if (.not.ma_chop_stack(l_scr)) call errquit
     &  ('dft_vspec: failed to deallocate scr',0,MA_ERR)
      if (.not.ma_push_get
     &     (mt_int,nexcite_tot,'scr',l_scr,k_scr))
     &  call errquit('dft_vspec: failed to alloc  scre',0,MA_ERR)
      call dft_sortmapi(nexcite_tot,
     I     int_mb(k_scr),int_mb(k_ex_occ),int_mb(k_mape))
      call dft_sortmapi(nexcite_tot,
     I     int_mb(k_scr),int_mb(k_ex_virt),int_mb(k_mape))
      if (.not.ma_chop_stack(l_mape)) call errquit
     &  ('dft_vspec: failed to deallocate mape',0,MA_ERR)

      endif
      endif ! nodezero
c
c     Output results
c
      if (nodezero.and.write_vspec) then
         write (unitno, "(1x,a)")
     $        "      N   Occ  Virt       Energy [eV]        "//
     $        "<x> [au]          <y> [au]          <z> [au]    "//
     $        "     fosc [au] "
         write (unitno, "(a,i10)") "<START>",nexcite_win
       do i=1,nexcite_tot
          if (i.le.nexcite_win) then
           write (unitno, "(i8,2i6,5e18.6)")
     $         i, int_mb(k_ex_occ+i-1), int_mb(k_ex_virt+i-1),
     $         dbl_mb(k_exen+i-1)*cau2ev,
     $         dbl_mb(k_dpx+i-1), dbl_mb(k_dpy+i-1), dbl_mb(k_dpz+i-1),
     $         dbl_mb(k_fosc+i-1)
          end if ! nexcite_win
       enddo
       write (unitno, "(a)") "<END>"
       close(unitno)
      endif
c
c     Broadcast status to other nodes
c
c
      call ga_sync()
 1000 format('== Scalar DFT Virtual Spectrum ==')
c
 1400 format(  'Number of basis functions in AO basis: ',i10,
     1       /,'Number of Alpha MOs:                   ',i10,
     2       /,'Number of Beta Mos:                    ',i10,
     3       /,'Occupied Win Start Alpha:              ',i10,
     4       /,'Occupied Win Start Beta:               ',i10,
     5       /,'Occupied Win End   Alpha:              ',i10,
     5       /,'Occupied Win End   Beta:               ',i10,
     6       /,'Number of Occupied Alpha MOs:          ',i10,
     6       /,'Number of Occupied Beta MOs :          ',i10,
     7       /,'Total Number of Excitations:           ',i10)
c
 1401 format(  'Number of basis functions in AO basis: ',i10,
     1       /,'Number of Alpha MOs:                   ',i10,
     3       /,'Occupied Win Start Alpha:              ',i10,
     5       /,'Occupied Win End   Alpha:              ',i10,
     6       /,'Number of Occupied Alpha MOs:          ',i10,
     7       /,'Total Number of Excitations:           ',i10)
c
 1600 format(/,x,'N',7x,'Energy',13x,'<x>',14x,'<y>',14x,'<z>',14x,
     &       'fosc')
c
 1700 format(i3,2x,f12.7,5x,f12.7,5x,f12.7,5x,f12.7,5x,f12.7)
c
c     call to compute polariz. with SOS
c
      if(polfromsos.and.nodezero)
     A     call dft_polfromsos(alpha,nexcite_tot,
     D     dbl_mb(k_dpx),dbl_mb(k_dpy),dbl_mb(k_dpz),
     E     dbl_mb(k_exen),ipol)

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
      if(nodezero) write(luout,'(a,f10.1,a)') ' VSPEC done in ',
     c     util_wallsec()-t0,'secs'
c
      end 
c $Id$
      subroutine dft_spece(ipol,owstart,owend,nocc,nbf,nmo,
     W     eval,exen,ex_occ,ex_virt)
      implicit none
      integer ipol,owstart(*),owend(*),nocc(*),nmo(*),nbf ! [in]
      double precision eval(nbf,*)  ! [in]
      integer ex_occ(*)
      integer ex_virt(*)
      double precision exen(*)
c
      integer nn,i,j,k
      double precision ea,ei
c
c
c     Calculate excitation energies
c
      nn = 0
      do i=1,ipol
c       account for case when beta but no alpha excitations 
        do j=owstart(i),owend(i)
          if (j.gt.0) then
            ei = eval(j,i)
            do k=nocc(i)+1,nmo(i)
              nn=nn+1
              ea = eval(k,i)
              exen(nn) = ea-ei
              ex_occ(nn) = j  !occupied MO for excitation i
              ex_virt(nn) = k !virtual MO for excitation i
            enddo
          endif
        enddo
      enddo
      return
      end
      subroutine dft_sorte(nexcite_tot,
     E     exen,mape)
      implicit none
c optimized bubble sort
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
      integer nexcite_tot ! [in]
      double precision   exen(*)
      integer mape(*)
      
      
c
      double precision esav,sort,t0
      integer i,j,jsav,isort,new_lgt,lgt
      t0=util_wallsec()
c      
      do i=1,nexcite_tot
         mape(i)=i
      enddo

c
c     Sort by energy
c
      lgt=nexcite_tot     
      do i=1,nexcite_tot-1 
         esav = exen(i)
         jsav = i
         do j=1,lgt-1
            if (exen(j) .gt. esav) then
               esav = exen(j)
               new_lgt=j
               jsav = j
            endif
         enddo
         lgt=new_lgt
         sort=exen(i)
         exen(i)=exen(jsav)
         exen(jsav)=sort
c         
         isort=mape(i)
         mape(i)=mape(jsav)
         mape(jsav)=isort
      enddo                     ! sorting
      if(ga_nodeid().eq.0) then
         write(luout,'(a,f10.3,a)') '  vspec sorting done in ',
     U        util_wallsec()-t0,'secs'
         call util_flush(luout)
      endif
      return
      end
      subroutine dft_gsorte(nexcite_tot,
     E     exen,mape)
      implicit none
c gnome sort
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer nexcite_tot ! [in]
      double precision   exen(*)
      integer mape(*)
      
c
      double precision esav,sort,t0
      integer i,jsav,isort,lgt
      t0=util_wallsec()

      do i=1,nexcite_tot
         mape(i)=i
      enddo
c
c     Sort by energy
c
      lgt=nexcite_tot
      i=1
      do while (i.le.lgt)
         if (i.eq.1.or.exen(i).ge.exen(i-1)) then
            i = i +1
         else
            jsav=i-1
c     swap i and jsav
         sort=exen(i)
         exen(i)=exen(jsav)
         exen(jsav)=sort
c         
         isort=mape(i)
         mape(i)=mape(jsav)
         mape(jsav)=isort
         endif
      enddo                     ! sorting
      if(ga_nodeid().eq.0) then
         write(luout,'(a,f10.3,a)') '  vspec sorting done in ',
     U        util_wallsec()-t0,'secs'
         call util_flush(luout)
       endif
      return
      end
      subroutine dft_sortmap(n,scre,newe,mape)
      implicit none
#include "mafdecls.fh"
      integer n
      integer mape(*)
      double precision newe(*),scre(*)
c
      integer i
c
      do i=1,n
         scre(i)=newe(i)
      enddo
      do i=1,n
         newe(i)=scre(mape(i))
      enddo
      return
      end
      subroutine dft_sortmapi(n,scre,newe,mape)
      implicit none
#include "mafdecls.fh"
      integer n
      integer mape(*)
      integer newe(*),scre(*)
c
      integer i
c
      do i=1,n
         scre(i)=newe(i)
      enddo
      do i=1,n
         newe(i)=scre(mape(i))
      enddo
      return
      end
