      subroutine tddft_analysis(rtdb,geom,ao_bas_han,nroots,nbf_ao,
     1  ipol,nocc,nmo,nfc,nfv,ntrials,tda,oskel,g_trials,g_buffer1,
     2  d_trials,apbvec,apbval,ambvec,g_movecs,evl,singlet,triplet,
     4  target,targetsym,mult,algorithm,l_irs,k_irs)
c
c $Id: tddft_analysis.F,v 1.18 2009-03-18 01:32:06 jhammond Exp $
c
c Calculate various properties of excitations such as
c (1) transition moments and oscillator strengths,
c (2) spatial symmetry and spin expectation values.
c
c !! CAUTION !! au2ev conversion factor is hardwired here.
c
c Written by So Hirata, Mar 2002. (c) Battelle, PNNL, 2002.
c
c Circular dichroism and velocity-representation dipole moments:
c J. Autschbach, SUNY Buffalo (2009).
c reference: J.Chem.Phys. 116 (2002), pages 891 and 6930
c
      implicit none
#include "errquit.fh"
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
#include "dra.fh"
c
      double precision thresh  ! Threshold for printing
      parameter (thresh=5.0d-2)
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)

      integer rtdb             ! RTDB handle
      integer geom             ! Geometry handle
      integer ao_bas_han       ! AO basis set handle
      integer nroots           ! Number of roots sought
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! = 1 (RDFT); =2 (UDFT)
      integer nocc(2)          ! Number of occupied orbitals
      integer nmo(2)           ! Number of orbitals
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer g_trials(2)      ! GA handle for trial vectors
      integer g_buffer1(2)     ! GA handle for buffer vectors
      integer d_trials(2)      ! DRA handle for trial vectors
      integer ntrials          ! Current number of trial vectors
      logical tda              ! True if Tamm-Dancoff approximation
      logical oskel            ! True if symmetry is used in integral evaluation
      double precision apbvec(ntrials,ntrials) ! (X+Y)
      double precision ambvec(ntrials,ntrials) ! (X-Y)
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
      double precision tmom(20) ! Transition moments
      double precision cntr(3) ! Center of dipole (arbitrary)
      integer target           ! Target root
      character*4 targetsym    ! Target symmetry
      double precision energy  ! DFT/TDDFT energy
      double precision tenergy ! Excitation energy of target root
      integer mult             ! Ground state spin multiplicity
      integer algorithm        ! Algorithm
c
      integer g_ovlp           ! GA handle for AO overlap matrix
      integer g_corr           ! GA handle for alpha-beta MO correlation matrix
      integer g_work
      integer g_tmomd(2)       ! Transition density matrix
      double precision r,maxr
      integer occ,vir,pol
      integer isym,jsym
      character*5 spin(2)
      logical nodezero
      integer i,j,k,l,m,n
      double precision na,nb,s2
      integer ja1,ja2,ka1,ka2,jb1,jb2,kb1,kb2,la1,la2,lb1,lb2,ja3,jb3
      double precision xa1,xa2,ya1,ya2,xb1,xb2,yb1,yb2,rc1,rc2
      integer l_x(2),k_x(2)
      integer l_y(2),k_y(2)
      integer l_corr,k_corr
      integer dummy,request

c     jochen 11/2009: add calculation of rotatory strengths (CD spectra)    
      logical cdspectrum ! true if CDSPECTRUM in TDDFT input
      logical velocity   ! true if VELOCITY in TDDFT input
      logical lvelok     ! check for division by zero
c     We will declare the variables needed for the g_dip{x,y,z}
c     This was modified from dft_main0d.F
      integer g_dipole(3)  ! ga for dipole integrals
      integer g_dipmag ! ga for magnetic dipole matrices (all 3)
      integer g_dipvel ! ga for dipole-velocity integrals
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer icomp ! field components
c     caution: speed of light and rotatory strength conversion factor
c     are hardwired here:
      double precision cl, rau2cgs
      parameter (cl = 137.0359895d0, rau2cgs = 235.726327d0)
      integer g_temp, g_tmpdip
      double precision el_transdip(3)  ! length-gauge electric dipoles
      double precision mag_transdip(3) ! magnetic transition dipoles (rxp)
      double precision vel_transdip(3) ! velocity electric dipoles
      double precision ga_trace_diag
      external ga_trace_diag
      character*32 pname
      double precision origin(3)
      data origin/0d0,0d0,0d0/
c     /jochen

C     character*80 vector
C     logical vector_exists,fileerror
c
      nodezero=(ga_nodeid().eq.0)

c     jochen:

c     check for CDSPECTRUM input key in RTDB
      if (.not.rtdb_get(rtdb,'tddft:cdspectrum',mt_log,1,cdspectrum))
     $  cdspectrum=.false.

c     check for VELOCITY input key in RTDB
      if (.not.rtdb_get(rtdb,'tddft:velocity',mt_log,1,velocity))
     $  velocity=.false.

      pname="tddft_analysis: "
      if (cdspectrum) then         
        
c     We will create arrays for some dipole matrix elements here.
c     adapted from dft_main0d.F
c     == Create and get the dipole matrix element ga: g_dipole ==

         do icomp = 1,3
            if (.not. ga_create(MT_DBL, nbf_ao, nbf_ao, 'g_dipole',
     $       nbf_ao, 0, g_dipole(icomp)))
     $        call errquit(pname//'ga_create failed g_dipole',0, GA_ERR)
            call ga_zero(g_dipole(icomp))
         end do ! icomp
c     
c     == Compute the dipole integrals for the 3 components ==
         call int_dip_ga(ao_bas_han, ao_bas_han,
     $        g_dipole(1), g_dipole(2), g_dipole(3))

c     magnetic dipole moment integrals (r x nabla), times 0.5:
         alo(1) = nbf_ao
         alo(2) = -1
         alo(3) = -1
         ahi(1) = nbf_ao
         ahi(2) = nbf_ao
         ahi(3) = 3         
         if (.not.nga_create(MT_DBL,3,ahi,'mag-dipole',alo,g_dipmag))
     $       call errquit(pname//'nga_create failed g_dipmag',0,GA_ERR)
         call ga_zero(g_dipmag)

         call int_giao_1ega(ao_bas_han,ao_bas_han,g_dipmag,'angmom',
     $        origin,1, .false.)
         call ga_scale (g_dipmag, 0.5d0)
        endif ! cdspectrum

        if (cdspectrum .or. velocity) then
           if (.not. ga_create(MT_DBL, nbf_ao, nbf_ao, 'g_temp',
     $       nbf_ao, 0, g_temp)) call errquit(pname//'g_temp',0, GA_ERR)
             call ga_zero(g_temp)

           if (.not. ga_create(MT_DBL, nbf_ao, nbf_ao, 'g_tmpdip',
     $       nbf_ao, 0, g_tmpdip))
     $          call errquit(pname//'g_tmpdip',0, GA_ERR)
           call ga_zero(g_tmpdip)
        end if ! cdspectrum .or. velocity

        if (velocity) then
c     allocate workspace for dipole-velocity integrals
         alo(1) = nbf_ao
         alo(2) = -1
         alo(3) = -1
         ahi(1) = nbf_ao
         ahi(2) = nbf_ao
         ahi(3) = 3   
           if (.not.nga_create(MT_DBL,3,ahi,'e-dipvel',alo,g_dipvel))
     $        call 
     &        errquit(pname//' nga_create failed g_dipvel',0,GA_ERR)
           call ga_zero(g_dipvel)
c     compute  dipole-velocity integrals
           call int_giao_1ega(ao_bas_han,ao_bas_han,g_dipvel,'velocity',
     $          origin,1,.false.)
c     invert sign of velocity integrals for compatibility 
c     w/ dipole-length formalism
           call ga_scale (g_dipvel, -1d0) 
        end if                  ! velocity
c
c     /jochen

c
c --------------------------------
c Number of occupied virtual pairs
c --------------------------------
c
      spin(1)='alpha'
      spin(2)='beta '
      do i=1,ipol
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
      enddo
      if (ipol.eq.1) nov(2)=0
c
c -------------------------------------
c Recover X and Y vectors from Z vector
c -------------------------------------
c
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nov(i),nroots,'X vector',
     1    -1,-1,g_x(i))) call errquit
     2    ('tddft_analysis: failed to create g_x',0, GA_ERR)
        if (.not.tda) then
          if (.not.ga_create(mt_dbl,nov(i),nroots,'Y vector',
     1      -1,-1,g_y(i))) call errquit
     2      ('tddft_analysis: failed to create g_y',0, GA_ERR)
        endif
      enddo
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nov(i),'X vector',
     1    l_x(i),k_x(i))) call errquit
     2    ('tddft_analysis: failed to allocate x',0, MA_ERR)
        if (.not.ma_push_get(mt_dbl,nov(i),'Y vector',
     1    l_y(i),k_y(i))) call errquit
     2    ('tddft_analysis: failed to allocate y',0, MA_ERR)
      enddo
c
c ----------------------------------
c Get alpha-beta orbital correlation
c ----------------------------------
c
      if (ipol.eq.2) then
        if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,
     1    'AO overlap',-1,-1,g_ovlp))
     2    call errquit('tddft_analysis: failed to create g_ovlp',0,
     &       GA_ERR)
        if (.not.ga_create(mt_dbl,nmo(1),nmo(2),
     1    'MO correlation',-1,-1,g_corr))
     2    call errquit('tddft_analysis: failed to create g_corr',0,
     &       GA_ERR)
        if (.not.ga_create(mt_dbl,nbf_ao,nmo(2),'work',-1,-1,g_work))
     1    call errquit('tddft_analysis: failed to create g_work',0,
     &       GA_ERR)
        if (.not.ma_push_get(mt_dbl,nmo(1)*nmo(2),'corr',
     1    l_corr,k_corr)) call errquit
     2    ('tddft_analysis: failed to allocate corr',0, MA_ERR)
        call ga_zero(g_ovlp)
        call int_1e_ga(ao_bas_han,ao_bas_han,g_ovlp,'overlap',oskel)
        if (oskel) call sym_symmetrize
     1    (geom,ao_bas_han,.false.,g_ovlp)
        call ga_dgemm('N','N',nbf_ao,nmo(2),nbf_ao,1.0d0,
     1    g_ovlp,g_movecs(2),0.0d0,g_work)
        call ga_dgemm('T','N',nmo(1),nmo(2),nbf_ao,1.0d0,
     1    g_movecs(1),g_work,0.0d0,g_corr)
        if (.not.ga_destroy(g_work)) call errquit
     1    ('tddft_analysis: failed to destroy g_work',0, GA_ERR)
        if (.not.ga_destroy(g_ovlp)) call errquit
     1    ('tddft_analysis: failed to destroy g_ovlp',0, GA_ERR)
        if (util_print('excited state',print_debug))
     1    call ga_print(g_corr)
        call ga_get(g_corr,1,nmo(1),1,nmo(2),dbl_mb(k_corr),nmo(1))
      endif
c
c --------------------
c Solution vector file
c --------------------
c
C     if (.not.rtdb_cget(rtdb,'tddft:vector',1,vector))
C    1  call errquit('tddft_analysis: failed to read vector',0)
C     inquire(file=vector,exist=vector_exists)
C     fileerror=.false.
C     if (vector_exists) then
C       if (nodezero)
C    1    open(unit=69,file=vector,form='formatted',status='old',
C    2    err=90000,access='sequential')
C         goto 90010
C     else
C       if (nodezero)
C    1    open(unit=69,file=vector,form='formatted',status='new',
C    2    err=90000,access='sequential')
C         goto 90010
C     endif
C0000 fileerror=.true.
C0010 continue
C     if (fileerror) then
C       if (nodezero.and.util_print('CI vector',print_default))
C    1    write(LuOut,2000)
C     else
C       if (nodezero.and.util_print('CI vector',print_default)) then
C         write(LuOut,2010) vector
C         rewind(69)
C         write(69,*) tda
C         write(69,*) ipol
C         write(69,*) nroots
C         write(69,*) nov(1),nov(2)
C       endif
C     endif
C2000 format(/,2x,'No CI vector file is created')
C2010 format(/,2x,'CI vectors are stored in ',a20)
c
c ------------
c Tamm-Dancoff
c ------------
c
      if (tda) then
        do i=1,ipol
          call ga_zero(g_x(i))
          do n=1,nroots
            do m=1,ntrials
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                call ga_add_patch(apbvec(m,n),g_trials(i),1,nov(i),m,m,
     1                            1.0d0,g_x(i),1,nov(i),n,n,
     2                                  g_x(i),1,nov(i),n,n)
              else
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,d_trials(i),1,nov(i),m,m,
     2            request)
                dummy=dra_wait(request)
                call ga_add_patch(apbvec(m,n),g_buffer1(i),1,nov(i),1,1,
     1                            1.0d0,g_x(i),1,nov(i),n,n,
     2                                  g_x(i),1,nov(i),n,n)
              endif 
            enddo
C           if ((.not.fileerror).and.nodezero) then
C             do m=1,nov(i)
C               call ga_get(g_x(i),m,m,n,n,r,1)
C               write(69,*) r
C             enddo
C           endif
          enddo
        enddo
c
c --------------------
c Full linear response
c --------------------
c
      else
        do i=1,ipol
          call ga_zero(g_x(i))
          call ga_zero(g_y(i))
          do n=1,nroots
            do m=1,ntrials
              if ((algorithm.eq.1).or.(algorithm.eq.2)) then
                call ga_add_patch(apbvec(m,n),g_trials(i),1,nov(i),m,m,
     1                            1.0d0,g_x(i),1,nov(i),n,n,
     2                                  g_x(i),1,nov(i),n,n)
                call ga_add_patch(ambvec(m,n),g_trials(i),1,nov(i),m,m,
     1                            1.0d0,g_y(i),1,nov(i),n,n,
     2                                  g_y(i),1,nov(i),n,n)
              else
                dummy=dra_read_section(.false.,
     1            g_buffer1(i),1,nov(i),1,1,d_trials(i),1,nov(i),m,m,
     2            request)
                dummy=dra_wait(request)
                call ga_add_patch(apbvec(m,n),g_buffer1(i),1,nov(i),1,1,
     1                            1.0d0,g_x(i),1,nov(i),n,n,
     2                                  g_x(i),1,nov(i),n,n)
                call ga_add_patch(ambvec(m,n),g_buffer1(i),1,nov(i),1,1,
     1                            1.0d0,g_y(i),1,nov(i),n,n,
     2                                  g_y(i),1,nov(i),n,n)
              endif
            enddo
          enddo
          call ga_add(0.5d0,g_x(i), 0.5d0,g_y(i),g_x(i))
          call ga_add(1.0d0,g_x(i),-1.0d0,g_y(i),g_y(i))
C         if ((.not.fileerror).and.nodezero) then
C           do n=1,nroots
C             do m=1,nov(i)
C               call ga_get(g_x(i),m,m,n,n,r,1)
C               write(69,*) r
C             enddo
C           enddo
C           do n=1,nroots
C             do m=1,nov(i)
C               call ga_get(g_y(i),m,m,n,n,r,1)
C               write(69,*) r
C             enddo
C           enddo
C         endif
        enddo
      endif
C     if (.not.fileerror) close(unit=69,status='keep')
c
c ---------
c All roots
c ---------
c
      do i=1,ipol
        if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'Transition density',
     1    -1,-1,g_tmomd(i))) call errquit
     2    ('tddft_analysis: failed to create g_tmomd',0, GA_ERR)
      enddo
c
c ------------
c Ground state
c ------------
c
      if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,energy)) 
     1  call errquit('tddft_analysis: failed to get dft energy',0,
     &       RTDB_ERR)
      isym=1
      if (ipol.eq.2) then
        do i=1,ipol
c here we are just getting ground state symmetry ...
c ... frozen cores are irrelevant
          do j=1,nocc(i)
            isym=ieor(isym-1,int_mb(k_irs(i)+j-1)-1)+1
          enddo
        enddo
      endif
      call sym_irrepname(geom,isym,oname)
      if (nodezero.and.util_print('ground state',print_default)) then
        write(LuOut,*)
        write(LuOut,9100) oname,energy
        call util_flush(LuOut)
      endif
      if (ipol.eq.2) then
        na=dfloat(nocc(1))
        nb=dfloat(nocc(2))
        s2=na-(na-nb)/2.0d0+((na-nb)/2.0d0)**2
        do j=1,nocc(1)
          do k=1,nocc(2)
            r=dbl_mb(k_corr+(k-1)*nmo(1)+j-1)
            s2=s2-r**2
          enddo
        enddo
        if (nodezero.and.util_print('ground state',print_default)) then
          write(LuOut,9180) s2
          call util_flush(LuOut)
        endif
        na=dfloat(nocc(1)-nfc(1))
        nb=dfloat(nocc(2)-nfc(2))
      endif
c
c --------------
c Excited states
c --------------
c
      if (nodezero.and.util_print('excited state',print_default))
     1  write(LuOut,*)
c m is a counter of excited states in a certain irrep 
c for the purpose of finding the target root
      m=0
      do n=1,nroots
        do i=1,ipol
          call ga_get(g_x(i),1,nov(i),n,n,dbl_mb(k_x(i)),1)
          if (.not.tda)
     1      call ga_get(g_y(i),1,nov(i),n,n,dbl_mb(k_y(i)),1)
        enddo
c
c --------------------------------------
c Identify symmetry of the excited state
c --------------------------------------
c
        maxr=0.0d0
        occ=0
        vir=0
        do i=1,ipol
          l=0
          do j=nfc(i)+1,nocc(i)
            do k=nocc(i)+1,nmo(i)-nfv(i)
              l=l+1
              r=dbl_mb(k_x(i)+l-1)
              if (dabs(r).gt.maxr) then
                maxr=dabs(r)
                pol=i
                occ=j
                vir=k
              endif
            enddo
          enddo
        enddo
        jsym=ieor(isym-1,int_mb(k_irs(pol)+occ-1)-1)+1
        jsym=ieor(jsym-1,int_mb(k_irs(pol)+vir-1)-1)+1
        call sym_irrepname(geom,jsym,oname)
        if ((targetsym.eq.'none').or.(targetsym.eq.oname)) then
          m=m+1
          if (m.eq.target) tenergy=apbval(n)
        endif
        if (nodezero.and.
     1    util_print('excited state',print_default)) then
          if (ipol.eq.2) then
            write(LuOut,9110) n,'       ',oname,apbval(n),
     1        apbval(n)*au2ev
          else if (singlet) then
            write(LuOut,9110) n,'singlet',oname,apbval(n),
     1        apbval(n)*au2ev
          else if (triplet) then
            write(LuOut,9110) n,'triplet',oname,apbval(n),
     1        apbval(n)*au2ev
          endif
        endif
c
c --------------------------------
c Compute <S**2> expectation value
c --------------------------------
c
        if (ipol.eq.2) then
          s2=na-(na-nb)/2.0d0+((na-nb)/2.0d0)**2
c
c alpha excitation ket
c
          la1=0
          do ja1=nfc(1)+1,nocc(1)
            do ka1=nocc(1)+1,nmo(1)-nfv(1)
              la1=la1+1
              xa1=dbl_mb(k_x(1)+la1-1)
              ya1=0.0d0
              if (.not.tda) ya1=dbl_mb(k_y(1)+la1-1)
              if (dabs(xa1+ya1).gt.1.0d-5) then
c a alpha and j beta exchange
                do jb1=nfc(2)+1,nocc(2)
c alpha excitation bra
                  la2=0
                  do ja2=nfc(1)+1,nocc(1)
                    do ka2=nocc(1)+1,nmo(1)-nfv(1)
                      la2=la2+1
                      if (ja1.eq.ja2) then
                        xa2=dbl_mb(k_x(1)+la2-1)
                        ya2=0.0d0
                        if (.not.tda) ya2=dbl_mb(k_y(1)+la2-1)
                        rc1=dbl_mb(k_corr+(jb1-1)*nmo(1)+ka2-1)
                        rc2=dbl_mb(k_corr+(jb1-1)*nmo(1)+ka1-1)
                        s2=s2-(xa2-ya2)*(xa1+ya1)*rc1*rc2
                      endif
                    enddo
                  enddo
c beta excitation bra
                  lb2=0
                  do jb2=nfc(2)+1,nocc(2)
                    do kb2=nocc(2)+1,nmo(2)-nfv(2)
                      lb2=lb2+1
                      if (jb1.eq.jb2) then
                        xb2=dbl_mb(k_x(2)+lb2-1)
                        yb2=0.0d0
                        if (.not.tda) yb2=dbl_mb(k_y(2)+lb2-1)
                        rc1=dbl_mb(k_corr+(jb2-1)*nmo(1)+ja1-1)
                        rc2=dbl_mb(k_corr+(kb2-1)*nmo(1)+ka1-1)
                        s2=s2-(xb2-yb2)*(xa1+ya1)*rc1*rc2
                      endif
                    enddo
                  enddo
                enddo
c i alpha and j beta exchange
                do ja3=nfc(1)+1,nocc(1)
                  if (ja3.ne.ja1) then
                    do jb1=nfc(2)+1,nocc(2)
c same alpha excitation bra
                      rc1=dbl_mb(k_corr+(jb1-1)*nmo(1)+ja3-1)
                      s2=s2-(xa1-ya1)*(xa1+ya1)*rc1*rc1
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
c
c beta excitation ket
c
          lb1=0
          do jb1=nfc(2)+1,nocc(2)
            do kb1=nocc(2)+1,nmo(2)-nfv(2)
              lb1=lb1+1
              xb1=dbl_mb(k_x(2)+lb1-1)
              yb1=0.0d0
              if (.not.tda) yb1=dbl_mb(k_y(2)+lb1-1)
              if (dabs(xb1+yb1).gt.1.0d-5) then
c i alpha and b beta exchange
                do ja1=nfc(1)+1,nocc(1)
c beta excitation bra
                  lb2=0
                  do jb2=nfc(2)+1,nocc(2)
                    do kb2=nocc(2)+1,nmo(2)-nfv(2)
                      lb2=lb2+1
                      if (jb1.eq.jb2) then
                        xb2=dbl_mb(k_x(2)+lb2-1)
                        yb2=0.0d0
                        if (.not.tda) yb2=dbl_mb(k_y(2)+lb2-1)
                        rc1=dbl_mb(k_corr+(kb1-1)*nmo(1)+ja1-1)
                        rc2=dbl_mb(k_corr+(kb2-1)*nmo(1)+ja1-1)
                        s2=s2-(xb2-yb2)*(xb1+yb1)*rc1*rc2
                      endif
                    enddo
                  enddo
c alpha excitation bra
                  la2=0
                  do ja2=nfc(1)+1,nocc(1)
                    do ka2=nocc(1)+1,nmo(1)-nfv(1)
                      la2=la2+1
                      if (ja1.eq.ja2) then
                        xa2=dbl_mb(k_x(1)+la2-1)
                        ya2=0.0d0
                        if (.not.tda) ya2=dbl_mb(k_y(1)+la2-1)
                        rc1=dbl_mb(k_corr+(jb1-1)*nmo(1)+ja2-1)
                        rc2=dbl_mb(k_corr+(kb1-1)*nmo(1)+ka2-1)
                        s2=s2-(xa2-ya2)*(xb1+yb1)*rc1*rc2
                      endif
                    enddo
                  enddo
                enddo
c i alpha and j beta exchange
                do jb3=nfc(2)+1,nocc(2)
                  if (jb3.ne.jb1) then
                    do ja1=nfc(1)+1,nocc(1)
c same beta excitation bra
                      rc1=dbl_mb(k_corr+(jb3-1)*nmo(1)+ja1-1)
                      s2=s2-(xb1-yb1)*(xb1+yb1)*rc1*rc1
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
          if (nodezero.and.util_print('excited state',print_default)
     1      .and.(nfc(1).eq.nfc(2)))
     2      write(LuOut,9180) s2
        endif
c
c --------------------------------------------------
c Compute transition moments and oscillator strength
c --------------------------------------------------
c
        do i=1,3
           cntr(i)=0.0d0
        enddo
        do i=1,20
          tmom(i)=0.0d0
        enddo

c     jochen:
        if (cdspectrum) then           
c     initialize transition dipoles for this transition (mainly debug)
           do icomp = 1,3
              el_transdip (icomp) = 0d0
              mag_transdip(icomp) = 0d0 ! 3 components
           end do
        end if ! cdspectrum
        if (velocity) then
           do icomp = 1,3
              vel_transdip(icomp) = 0d0 ! 3 components
           end do
           lvelok = (apbval(n) .gt. 1d-6) ! make sure excit. E is >0
        end if ! velocity

        if (cdspectrum .or. velocity) then
           alo(1) = 1
           ahi(1) = nbf_ao
           alo(2) = 1
           ahi(2) = nbf_ao
           alo(3) = 1
           ahi(3) = 1
           blo(1) = 1
           bhi(1) = nbf_ao
           blo(2) = 1
           bhi(2) = nbf_ao
        end if ! cdpectrum .or. velocity
c     /jochen

c ... jochen - comment: the following code 
c     calculates X(omega) in the AO basis, then contracts with multipole
c     integrals:

        call tddft_transfm(n,g_x,g_movecs,nbf_ao,
     1    nocc,nmo,nfc,nfv,ipol,g_tmomd,1)

        do i=1,ipol

c     jochen:

          if (cdspectrum) then

c     The following is adapted from rttdf_utils.F:
c     recompute the electric transition dipoles (mainly to
c     debug the code) and also compute the magnetic transition dip.

             do icomp = 1,3
                call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &               g_tmomd(i), g_dipole(icomp), 0.d0, g_temp)
                el_transdip(icomp) = el_transdip(icomp)
     $               - ga_trace_diag(g_temp) 

                blo(3) = icomp 
                bhi(3) = icomp       
                
              call nga_copy_patch('n',g_dipmag,blo,bhi,g_tmpdip,alo,ahi)
                call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &               g_tmomd(i), g_tmpdip, 0.d0, g_temp)
                mag_transdip(icomp) = mag_transdip(icomp)  +
     $               ga_trace_diag(g_temp)
                
             end do ! icomp

          endif                 ! cdspectrum

          if (velocity .and. lvelok) then
             do icomp = 1,3
                blo(3) = icomp 
                bhi(3) = icomp        
     
              call nga_copy_patch('n',g_dipvel,blo,bhi,g_tmpdip,alo,ahi)
                call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &               g_tmomd(i), g_tmpdip, 0.d0, g_temp)
                vel_transdip(icomp) = vel_transdip(icomp)  +
     $               ga_trace_diag(g_temp) / apbval(n)
             end do             ! icomp
          end if                ! velocity
c     /jochen

          call multipole_density(ao_bas_han,cntr,3,
     1      g_tmomd(i),tmom,20)

        enddo ! ipol

        if (.not.tda) then ! we have a also matrix Y to consider

c ... jochen - comment: repeat the preceeding step with Y:
c     calculate Y(omega) in AO basis, then contract with multipole
c     integrals, and add results to the transition moments

          call tddft_transfm(n,g_y,g_movecs,nbf_ao,
     1      nocc,nmo,nfc,nfv,ipol,g_tmomd,1)

          do i=1,ipol

c     jochen: 
            if (cdspectrum) then

             do icomp = 1,3

                call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &               g_tmomd(i), g_dipole(icomp), 0.d0, g_temp)
                el_transdip(icomp) = el_transdip(icomp)
     $               - ga_trace_diag(g_temp) 

                blo(3) = icomp  ! there are 3 components of the
                bhi(3) = icomp  ! responding magnetic moment
     
              call nga_copy_patch('n',g_dipmag,blo,bhi,g_tmpdip,alo,ahi)
                call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &               g_tmomd(i), g_tmpdip, 0.d0, g_temp)
                mag_transdip(icomp) = mag_transdip(icomp)  -
     $               ga_trace_diag(g_temp)
                
             end do ! icomp

            endif ! cdspectrum

          if (velocity .and. lvelok) then
             do icomp = 1,3
                blo(3) = icomp 
                bhi(3) = icomp 
     
              call nga_copy_patch('n',g_dipvel,blo,bhi,g_tmpdip,alo,ahi)
                call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     &               g_tmomd(i), g_tmpdip, 0.d0, g_temp)
                vel_transdip(icomp) = vel_transdip(icomp)  -
     $               ga_trace_diag(g_temp) / apbval(n)
             end do             ! icomp
          end if                ! velocity
c     /jochen

            call multipole_density(ao_bas_han,cntr,3,
     1        g_tmomd(i),tmom,20)

         enddo ! ipol

        endif ! .not.tda 

        if (ipol.eq.1) then
          do i=1,20
            tmom(i)=tmom(i)*dsqrt(2.0d0)
          enddo

c     jochen:
          if (cdspectrum) then
            do icomp = 1,3
              el_transdip(icomp) = el_transdip(icomp) *dsqrt(2.0d0)
              mag_transdip(icomp) = mag_transdip(icomp) *dsqrt(2.0d0)
            end do
          end if
          if (velocity) then
            do icomp = 1,3
               if (.not.cdspectrum) ! keep the old stuff for now
     $          el_transdip(icomp) = tmom(icomp+1)   

              vel_transdip(icomp) = vel_transdip(icomp) *dsqrt(2.0d0)
            end do
         end if ! velocity
c     /jochen

        endif ! ipol.eq.1 ?


        if (nodezero.and.util_print('excited state',print_default))
     1    write(LuOut,9190)
        if (nodezero.and.
     1    util_print('excited state',print_default)) then
          if ((singlet.and.(mult.ne.1)).or.
     1        (triplet.and.(mult.ne.3))) then
            write(LuOut,9210)
          else
            write(LuOut,9200) tmom(2),tmom(3),tmom(4),
     1                        tmom(5),tmom(6),tmom(7),
     2                        tmom(8),tmom(9),tmom(10),
     2                        tmom(11),tmom(12),tmom(13),
     2                        tmom(14),tmom(15),tmom(16),
     2                        tmom(17),tmom(18),tmom(19),
     2                        tmom(20),
     1      2.0d0/3.0d0*(tmom(2)**2+tmom(3)**2+tmom(4)**2)*apbval(n)

c     jochen:
            if (cdspectrum) then
               write (LuOut,*)
            write (LuOut,9250) 'Electric Transition Dipole:',
     $           el_transdip(1), el_transdip(2), el_transdip(3)
            write (LuOut,9250) 'Magnetic Transition Dipole:',
     $           mag_transdip(1)/cl,
     $           mag_transdip(2)/cl,
     $           mag_transdip(3)/cl
            write (LuOut,9251)
     $          'Rotatory Strength (1E-40 esu**2cm**2):',
     $          (el_transdip(1)*mag_transdip(1) +
     $           el_transdip(2)*mag_transdip(2) +
     $           el_transdip(3)*mag_transdip(3))
     $           * 2d0 * rau2cgs 
            end if ! cdspectrum
            if (velocity .and. lvelok) then
               write (LuOut,*)
            write (LuOut,9250)
     $      'Electric Transition Dipole (velocity representation):',
     $           vel_transdip(1),
     $           vel_transdip(2),
     $           vel_transdip(3)
            write (LuOut,9251)
     $           'Oscillator Strength (velocity repr.) :',
     $           (vel_transdip(1)*vel_transdip(1) +
     $            vel_transdip(2)*vel_transdip(2) +
     $            vel_transdip(3)*vel_transdip(3))
     $           * apbval(n) * 2.0d0/3.0d0
            write (LuOut,9251)
     $           'Oscillator Strength (mixed repr.   ) :',
     $           (vel_transdip(1)*el_transdip(1) +
     $            vel_transdip(2)*el_transdip(2) +
     $            vel_transdip(3)*el_transdip(3))
     $           * apbval(n) * 2.0d0/3.0d0 ! it can be negative ...
            if (cdspectrum) then
               write (LuOut,9251)
     $           'Rotatory Strength   (velocity repr.) :',
     $          (vel_transdip(1)*mag_transdip(1) +
     $           vel_transdip(2)*mag_transdip(2) +
     $           vel_transdip(3)*mag_transdip(3)) * 2d0 * rau2cgs 
            end if
            else if(velocity .and. .not.lvelok) then
              write (LuOut,*) 
     $        'Velocity: E <= 1e-6. skipping this excitation'
            end if ! velocity
c     /jochen
          endif
        endif

c
c ----------------------------------
c Print dominant components of X & Y
c ----------------------------------
c
        if (nodezero.and.
     1    util_print('excited state',print_default)) write(LuOut,*)
        if (ipol.eq.1) then
          l=0
          do j=nfc(1)+1,nocc(1)
            do k=nocc(1)+1,nmo(1)-nfv(1)
              l=l+1
              if (tda) then
                r=dbl_mb(k_x(1)+l-1)
                call sym_irrepname
     1            (geom,int_mb(k_irs(1)+j-1),oname)
                call sym_irrepname
     1            (geom,int_mb(k_irs(1)+k-1),vname)
                if (nodezero.and.util_print('excited state',
     1            print_default).and.(dabs(r).gt.thresh))
     2            write(LuOut,9120) j,oname,k,vname,r
              else
                r=dbl_mb(k_x(1)+l-1)
                call sym_irrepname
     1            (geom,int_mb(k_irs(1)+j-1),oname)
                call sym_irrepname
     1            (geom,int_mb(k_irs(1)+k-1),vname)
                if (nodezero.and.util_print('excited state',
     1            print_default).and.(dabs(r).gt.thresh))
     2            write(LuOut,9130) j,oname,k,vname,r
                r=dbl_mb(k_y(1)+l-1)
                if (nodezero.and.util_print('excited state',
     1            print_default).and.(dabs(r).gt.thresh))
     2            write(LuOut,9140) j,oname,k,vname,r
              endif
            enddo
          enddo
        else
          do i=1,ipol
            l=0
            do j=nfc(i)+1,nocc(i)
              do k=nocc(i)+1,nmo(i)-nfv(i)
                l=l+1
                if (tda) then
                  r=dbl_mb(k_x(i)+l-1)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(i)+j-1),oname)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(i)+k-1),vname)
                  if (nodezero.and.util_print('excited state',
     1              print_default).and.(dabs(r).gt.thresh))
     2              write(LuOut,9150) j,spin(i),oname,
     3              k,spin(i),vname,r
                else
                  r=dbl_mb(k_x(i)+l-1)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(i)+j-1),oname)
                  call sym_irrepname
     1              (geom,int_mb(k_irs(i)+k-1),vname)
                  if (nodezero.and.util_print('excited state',
     1              print_default).and.(dabs(r).gt.thresh))
     2              write(LuOut,9160) j,spin(i),oname,
     3              k,spin(i),vname,r
                  r=dbl_mb(k_y(i)+l-1)
                  if (nodezero.and.util_print('excited state',
     1              print_default).and.(dabs(r).gt.thresh))
     2              write(LuOut,9170) j,spin(i),oname,
     3              k,spin(i),vname,r
                endif
              enddo
            enddo
          enddo      
        endif
        if (nodezero) call util_flush(LuOut)
      enddo
      if (m.lt.target) 
     1  call errquit('tddft_analysis: target root not found',m,
     &       UNKNOWN_ERR)
      if (nodezero.and.
     1  util_print('excited state',print_default)) write(LuOut,*)

c       jochen:
        if (cdspectrum) then
c     Clean up arrays used for CD spectra
           do icomp = 1,3
              if (.not. ga_destroy(g_dipole(icomp))) call errquit
     &             (pname//'ga_destroy failed g_dipole', 0, GA_ERR)
           end do
           if (.not.ga_destroy(g_dipmag)) call
     &          errquit(pname//'ga_destroy failed g_dipmag',0,GA_ERR)

        end if

        if (cdspectrum .or. velocity) then
           if (.not. ga_destroy(g_temp)) call errquit
     &          (pname//'Could not destroy g_temp', 0, GA_ERR)
           if (.not. ga_destroy(g_tmpdip)) call errquit
     &          (pname//'Could not destroy g_tmpdip', 0, GA_ERR)
        endif ! cdspectrum .or. velocity
        
        if (velocity) then
c     clean up arrays used for dipole-velocity integrals
           if (.not.ga_destroy(g_dipvel)) call
     &       errquit(pname//'ga_destroy failed g_dipvel',0,GA_ERR)
        endif                   ! velocity
c      /jochen


 9100 format(2x,'Ground state ',a4,f16.9,' a.u.')
 9110 format(2x,
     1  '-------------------------------------------------------',
     2  /,2x,'Root',i3,1x,a7,1x,a4,f13.9,' a.u. (',f11.7,' eV)')
 9190 format(2x,
     1  '-------------------------------------------------------')
 9200 format(5x,'Transition Moments    X',f9.5,'   Y',f9.5,'   Z',f9.5
     1    ,/,5x,'Transition Moments   XX',f9.5,'  XY',f9.5,'  XZ',f9.5
     2    ,/,5x,'Transition Moments   YY',f9.5,'  YZ',f9.5,'  ZZ',f9.5
     3    ,/,5x,'Transition Moments  XXX',f9.5,' XXY',f9.5,' XXZ',f9.5
     4    ,/,5x,'Transition Moments  XYY',f9.5,' XYZ',f9.5,' XZZ',f9.5
     5    ,/,5x,'Transition Moments  YYY',f9.5,' YYZ',f9.5,' YZZ',f9.5
     6    ,/,5x,'Transition Moments  ZZZ',f9.5
     7    ,/,5x,'Dipole Oscillator Strength',21x,f12.5)
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
 9180 format(2x,'<S2> = ',f8.4)
 9250 format(5x,a
     $    ,/,5x,'       X ',f12.7,'   Y',f12.7,'   Z',f12.7)
 9251 format( 5x,a,1x,f20.7)
c
c -----------
c Target root
c -----------
c
      if (nodezero.and.util_print('excited state',print_low)) then
        write(LuOut,9300) target
        write(LuOut,9340) targetsym
        write(LuOut,9310) energy
        write(LuOut,9320) tenergy
      endif
      energy=energy+tenergy
      if (nodezero.and.util_print('excited state',print_low)) then
        write(LuOut,9330) energy
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      if (.not.rtdb_put(rtdb,'tddft:energy',mt_dbl,1,energy))
     1  call errquit('tddft_analysis: failed to put tddft energy',0,
     &       RTDB_ERR)
 9300 format('              Target root =',i7)
 9340 format('          Target symmetry =',1x,a4)
 9310 format('      Ground state energy =',f20.12)
 9320 format('        Excitation energy =',f20.12)
 9330 format('     Excited state energy =',f20.12)
c
c ------
c Return
c ------
c
      if (ipol.eq.2) then
        if (.not.ma_pop_stack(l_corr)) call errquit
     1    ('tddft_analysis: failed to deallocate corr',0, MA_ERR)
      endif
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_y(i))) call errquit
     1    ('tddft_analysis: failed to deallocate y',0, MA_ERR)
        if (.not.ma_pop_stack(l_x(i))) call errquit
     1    ('tddft_analysis: failed to deallocate x',0, MA_ERR)
      enddo
      if (ipol.eq.2) then
        if (.not.ga_destroy(g_corr)) call errquit
     1    ('tddft_analysis: failed to destroy g_corr',0, GA_ERR)
      endif
      do i=1,ipol
        if (.not.ga_destroy(g_x(i))) call errquit
     1    ('tddft_analysis: failed to destroy g_x',0, GA_ERR)
        if (.not.tda) then
          if (.not.ga_destroy(g_y(i))) call errquit
     1      ('tddft_analysis: failed to destroy g_y',0, GA_ERR)
        endif
      enddo
      do i=1,ipol
        if (.not.ga_destroy(g_tmomd(i))) call errquit
     2    ('tddft_analysis: failed to destroy g_tmomd',0, GA_ERR)
      enddo
      return
      end
