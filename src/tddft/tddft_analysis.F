      subroutine tddft_analysis(rtdb,geom,ao_bas_han,nroots,nbf_ao,
     1  ipol,nocc,nmo,ntrials,tda,g_trials,g_prod_p,g_apbvec,apbval,
     2  g_ambvec,g_movecs,evl,singlet,triplet,target,mult)
c
c $Id: tddft_analysis.F,v 1.2 2002-03-22 23:10:48 sohirata Exp $
c
c Calculate various properties of excitations such as
c (1) transition moments and oscillator strengths,
c (2) expectation values of r^2,
c (3) spatial symmetry and spin expectation values.
c
c !! CAUTION !! au2ev conversion factor is hardwired here.
c
c Written by So Hirata, Mar 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      double precision thresh  ! Threshold for printing
      parameter (thresh=5.0d-2)
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)
c
      integer rtdb             ! RTDB handle
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      integer nroots           ! Number of roots sought
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! = 1 (RDFT); =2 (UDFT)
      integer nocc(2)          ! Number of occupied orbitals
      integer nmo(2)           ! Number of orbitals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_prod_p(2)      ! GA handle for (A+B) product vectors
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      integer g_apbvec         ! GA handle for (X+Y)
      integer g_ambvec         ! GA handle for (X-Y)
      integer g_movecs(2)      ! GA handle for MO coefficients
      double precision apbval(ntrials)  ! (A-B)^(1/2)(A+B)(A-B)^(1/2) eigenvalues
      double precision evl(nbf_ao,2)    ! Orbital energies
      logical singlet          ! True if singlet excited state calculation
      logical triplet          ! True if triplet excited state calculation
      integer g_x(2)           ! GA handle for X vector
      integer g_y(2)           ! GA handle for Y vector
      character*4 oname        ! Irrep name
      character*4 vname        ! Irrep name
      integer l_irs(2),k_irs(2)! Irreps of MO's
      integer nov(2)           ! Number of occupied virtual pairs
      integer nmix(2)          ! Number of symmetry contaminations
      double precision tmom(4) ! Transition moments
      double precision cntr(3) ! Center of dipole (arbitrary)
      integer target           ! Target root
      double precision energy  ! DFT/TDDFT energy
      integer mult             ! Ground state spin multiplicity
c
      integer i,j,k,l,n
      integer g_work
      integer g_tmomd(2)       ! Transition density matrix
      double precision r,maxr
      integer occ,vir,pol
      integer isym,jsym
      character*5 spin(2)
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c --------------------------------
c Number of occupied virtual pairs
c --------------------------------
c
      spin(1)='alpha'
      spin(2)='beta '
      do i=1,ipol
        nov(i)=(nmo(i)-nocc(i))*nocc(i)
      enddo
      if (ipol.eq.1) nov(2)=0
c
c ---------------------------------------
c Determine symmetry (irreps) of orbitals
c ---------------------------------------
c
      do i=1,ipol
c
c We re-symmetry-adapt the MO's to get irreps information
c Make sure not to alter g_movecs
c
        if (.not.ma_push_get(mt_int,nbf_ao,'irreps',
     1    l_irs(i),k_irs(i))) call errquit
     2    ('tddft_analysis: failed to allocate irs',0)
        if (.not.ga_duplicate(g_movecs(i),g_work,'Work')) call errquit
     1    ('tddft_analysis: failed to duplicate g_movecs',0)
        call ga_copy(g_movecs(i),g_work)
        call sym_movecs_adapt(ao_bas_han,1.0d-8,g_work,
     1    int_mb(k_irs(i)),nmix(i))
        if (.not.ga_destroy(g_work)) call errquit
     1    ('tddft_analysis: failed to destroy g_work',0)
      enddo
c
c -------------------------------------
c Recover X and Y vectors from Z vector
c -------------------------------------
c
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nov(i),nroots,'X vector',
     1    -1,-1,g_x(i))) call errquit
     2    ('tddft_analysis: failed to create g_x',0) 
        if (.not.tda) then
          if (.not.ga_create(mt_dbl,nov(i),nroots,'Y vector',
     1      -1,-1,g_y(i))) call errquit
     2      ('tddft_analysis: failed to create g_y',0) 
        endif
      enddo
c
c ------------
c Tamm-Dancoff
c ------------
c
      if (tda) then
        do i=1,ipol
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_trials(i),g_apbvec,0.0d0,g_x(i))
        enddo
c
c --------------------
c Full linear response
c --------------------
c
      else
        do i=1,ipol
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_trials(i),g_apbvec,0.0d0,g_x(i))
          call ga_dgemm('N','N',nov(i),nroots,ntrials,
     1      1.0d0,g_trials(i),g_ambvec,0.0d0,g_y(i))
          call ga_add(0.5d0,g_x(i), 0.5d0,g_y(i),g_x(i))
          call ga_add(1.0d0,g_x(i),-1.0d0,g_y(i),g_y(i))
        enddo
      endif
c
c ---------
c All roots
c ---------
c
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'Transition density',
     1    -1,-1,g_tmomd(i))) call errquit
     2    ('tddft_analysis: failed to create g_tmomd',0) 
      enddo
c
c Ground state
c
      if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,energy)) 
     1  call errquit('tddft_analysis: failed to get dft energy',0)
      isym=1
      if (ipol.eq.2) then
        do i=1,ipol
          do j=1,nocc(i)
            call sym_product(geom,isym,int_mb(k_irs(i)+j-1),isym)
          enddo
        enddo
      endif
      call sym_irrepname(geom,isym,oname)
      if (nodezero.and.util_print('ground state',print_default)) then
        write(LuOut,*)
        write(LuOut,9100) oname,energy
        call util_flush(LuOut)
      endif
 9100 format(2x,'Ground state ',a4,f16.9,' a.u.')
c
c Excited states
c
      if (util_print('excited state',print_default)) then
        if (nodezero) write(LuOut,*)
        do n=1,nroots
c
c Identify symmetry of the excited state
c
          maxr=0.0d0
          occ=0
          vir=0
          do i=1,ipol
            l=0
            do j=1,nocc(i)
              do k=nocc(i)+1,nmo(i)
                l=l+1
                call ga_get(g_x(i),l,l,n,n,r,1)
                if (dabs(r).gt.maxr) then
                  maxr=dabs(r)
                  pol=ipol
                  occ=j
                  vir=k
                endif
              enddo
            enddo
          enddo
          call sym_product(geom,isym,
     1      int_mb(k_irs(pol)+occ-1),jsym)
          call sym_product(geom,jsym,
     1      int_mb(k_irs(pol)+vir-1),jsym)
          call sym_irrepname(geom,jsym,oname)
          if (nodezero) then
            if (ipol.eq.2) then
              write(LuOut,9110) n,'       ',oname,apbval(n),
     1          apbval(n)*au2ev
            else if (singlet) then
              write(LuOut,9110) n,'singlet',oname,apbval(n),
     1          apbval(n)*au2ev
            else if (triplet) then
              write(LuOut,9110) n,'triplet',oname,apbval(n),
     1          apbval(n)*au2ev
            endif
          endif
c
c Compute transition moments and oscillator strength
c
          do i=1,3
            cntr(i)=0.0d0
          enddo
          do i=1,4
            tmom(i)=0.0d0
          enddo
          call tddft_transfm(n,g_x,g_movecs,nbf_ao,
     1      nocc,nmo,ipol,g_tmomd,1)
          do i=1,ipol
            call multipole_density(ao_bas_han,cntr,1,
     1        g_tmomd(i),tmom,4)
          enddo
          if (.not.tda) then
            call tddft_transfm(n,g_y,g_movecs,nbf_ao,
     1        nocc,nmo,ipol,g_tmomd,1)
            do i=1,ipol
              call multipole_density(ao_bas_han,cntr,1,
     1          g_tmomd(i),tmom,4)
            enddo
          endif
          if (ipol.eq.1) then
            do i=2,4
              tmom(i)=tmom(i)*dsqrt(2.0d0)
            enddo
          endif
          if (nodezero) then
            if ((singlet.and.(mult.ne.1)).or.
     1          (triplet.and.(mult.ne.3))) then
              write(LuOut,9210)
            else
              write(LuOut,9200) tmom(2),tmom(3),tmom(4),
     1        2.0d0/3.0d0*(tmom(2)**2+tmom(3)**2+tmom(4)**2)*apbval(n)
          endif
          endif
c
c Print dominant components of X & Y
c
          if (nodezero) write(LuOut,*)
          if (ipol.eq.1) then
            l=0
            do j=1,nocc(1)
              do k=nocc(1)+1,nmo(1)
                l=l+1
                if (tda) then
                  call ga_get(g_x(1),l,l,n,n,r,1)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(1)+j-1),oname)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(1)+k-1),vname)
                  if (nodezero.and.(dabs(r).gt.thresh))
     1              write(LuOut,9120) j,oname,k,vname,r
                else
                  call ga_get(g_x(1),l,l,n,n,r,1)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(1)+j-1),oname)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(1)+k-1),vname)
                  if (nodezero.and.(dabs(r).gt.thresh))
     1              write(LuOut,9130) j,oname,k,vname,r
                  call ga_get(g_y(1),l,l,n,n,r,1)
                  if (nodezero.and.(dabs(r).gt.thresh))
     1              write(LuOut,9140) j,oname,k,vname,r
                endif
              enddo
            enddo
          else
            do i=1,ipol
              l=0
              do j=1,nocc(i)
                do k=nocc(i)+1,nmo(i)
                  l=l+1
                  if (tda) then
                    call ga_get(g_x(i),l,l,n,n,r,1)
                    call sym_irrepname
     1                (geom,int_mb(k_irs(i)+j-1),oname)
                    call sym_irrepname
     1                (geom,int_mb(k_irs(i)+k-1),vname)
                    if (nodezero.and.(dabs(r).gt.thresh))
     1                write(LuOut,9150) j,spin(i),oname,
     2                k,spin(i),vname,r
                  else
                    call ga_get(g_x(i),l,l,n,n,r,1)
                    call sym_irrepname
     1                (geom,int_mb(k_irs(i)+j-1),oname)
                    call sym_irrepname
     1                (geom,int_mb(k_irs(i)+k-1),vname)
                    if (nodezero.and.(dabs(r).gt.thresh))
     1                write(LuOut,9160) j,spin(i),oname,
     2                k,spin(i),vname,r
                    call ga_get(g_y(i),l,l,n,n,r,1)
                    if (nodezero.and.(dabs(r).gt.thresh))
     1                write(LuOut,9170) j,spin(i),oname,
     2                k,spin(i),vname,r
                  endif
                enddo
              enddo
            enddo      
          endif
          call util_flush(LuOut)
        enddo
        if (nodezero) write(LuOut,*)
      endif
 9110 format(2x,
     1  '-------------------------------------------------------',
     2  /,2x,'Root',i3,1x,a7,1x,a4,f13.9,' a.u. (',f11.7,' eV)',/,
     3  2x,'-------------------------------------------------------')
 9200 format(5x,'Transition Moments  X',f9.5,' Y',f9.5,' Z',f9.5
     1    ,/,5x,'Oscillator Strength',21x,f12.5)
 9210 format(5x,'Transition Moments                    Spin forbidden'
     1    ,/,5x,'Oscillator Strength                   Spin forbidden')
 9120 format(5x,'Occ.',i3,2x,a4,'---  Virt.',i3,2x,a4,f10.5)
 9130 format(5x,'Occ.',i3,2x,a4,'---  Virt.',i3,2x,a4,f8.5,' X')
 9140 format(5x,'Occ.',i3,2x,a4,'---  Virt.',i3,2x,a4,f8.5,' Y')
 9150 format(5x,'Occ.',i3,1x,a5,1x,a4,'---  Virt.',
     1  i3,1x,a5,1x,a4,f10.5)
 9160 format(5x,'Occ.',i3,1x,a5,1x,a4,'---  Virt.',
     1  i3,1x,a5,1x,a4,f8.5,' X')
 9170 format(5x,'Occ.',i3,1x,a5,1x,a4,'---  Virt.',
     1  i3,1x,a5,1x,a4,f8.5,' Y')
c
c -----------
c Target root
c -----------
c
      if (nodezero.and.util_print('excited state',print_low)) then
        write(LuOut,9300) target
        write(LuOut,9310) energy
        write(LuOut,9320) apbval(target)
      endif
      energy=energy+apbval(target)
      if (nodezero.and.util_print('excited state',print_low)) then
        write(LuOut,9330) energy
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      if (.not.rtdb_put(rtdb,'tddft:energy',mt_dbl,1,energy))
     1  call errquit('tddft_analysis: failed to put tddft energy',0)
 9300 format('              Target root =',i7)
 9310 format('      Ground state energy =',f20.12)
 9320 format('        Excitation energy =',f20.12)
 9330 format('     Excited state energy =',f20.12)
c
c ------
c Return
c ------
c
      do i=1,ipol
        if (.not.ga_destroy(g_x(i))) call errquit
     1    ('tddft_analysis: failed to destroy g_x',0) 
        if (.not.tda) then
          if (.not.ga_destroy(g_y(i))) call errquit
     1      ('tddft_analysis: failed to destroy g_y',0) 
        endif
      enddo
      do i=1,ipol
        if (.not.ga_destroy(g_tmomd(i))) call errquit
     2    ('tddft_analysis: failed to destroy g_tmomd',0) 
      enddo
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_irs(i))) call errquit
     1    ('tddft_analysis: failed to deallocate irs',0)
      enddo
      return
      end
