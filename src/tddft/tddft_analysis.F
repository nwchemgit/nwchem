      subroutine tddft_analysis(rtdb,geom,ao_bas_han,nroots,nbf_ao,
     1  ipol,nocc,nmo,nfc,nfv,ntrials,tda,oskel,g_trials,g_buffer1,
     2  d_trials,apbvec,apbval,ambvec,g_movecs,evl,singlet,triplet,
     4  target,targetsym,mult,algorithm)
c
c $Id: tddft_analysis.F,v 1.9 2002-04-28 19:09:13 sohirata Exp $
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
#include "dra.fh"
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
      integer nmix(2)          ! Number of symmetry contaminations
      double precision tmom(4) ! Transition moments
      double precision cntr(3) ! Center of dipole (arbitrary)
      integer target           ! Target root
      integer targetsym        ! Target symmetry
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
      double precision corr(nbf_ao,nbf_ao)
      integer dummy,request
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
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
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
      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nov(i),'X vector',
     1    l_x(i),k_x(i))) call errquit
     2    ('tddft_analysis: failed to allocate x',0)
        if (.not.ma_push_get(mt_dbl,nov(i),'Y vector',
     1    l_y(i),k_y(i))) call errquit
     2    ('tddft_analysis: failed to allocate y',0)
      enddo
c
c ----------------------------------
c Get alpha-beta orbital correlation
c ----------------------------------
c
      if (ipol.eq.2) then
        if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,
     1    'AO overlap',-1,-1,g_ovlp))
     2    call errquit('tddft_analysis: failed to create g_ovlp',0)
        if (.not.ga_create(mt_dbl,nmo(1),nmo(2),
     1    'MO correlation',-1,-1,g_corr))
     2    call errquit('tddft_analysis: failed to create g_corr',0)
        if (.not.ga_create(mt_dbl,nbf_ao,nmo(2),'work',-1,-1,g_work))
     1    call errquit('tddft_analysis: failed to create g_work',0)
        call ga_zero(g_ovlp)
        call int_1e_ga(ao_bas_han,ao_bas_han,g_ovlp,'overlap',oskel)
        if (oskel) call sym_symmetrize
     1    (geom,ao_bas_han,.false.,g_ovlp)
        call ga_dgemm('N','N',nbf_ao,nmo(2),nbf_ao,1.0d0,
     1    g_ovlp,g_movecs(2),0.0d0,g_work)
        call ga_dgemm('T','N',nmo(1),nmo(2),nbf_ao,1.0d0,
     1    g_movecs(1),g_work,0.0d0,g_corr)
        if (.not.ga_destroy(g_work)) call errquit
     1    ('tddft_analysis: failed to destroy g_work',0)
        if (.not.ga_destroy(g_ovlp)) call errquit
     1    ('tddft_analysis: failed to destroy g_ovlp',0)
        if (util_print('excited state',print_debug))
     1    call ga_print(g_corr)
        call ga_get(g_corr,1,nmo(1),1,nmo(2),corr,nbf_ao)
      endif
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
c here we are just getting ground state symmetry ...
c ... frozen cores are irrelevant
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
      if (ipol.eq.2) then
        na=dfloat(nocc(1))
        nb=dfloat(nocc(2))
        s2=na-(na-nb)/2.0d0+((na-nb)/2.0d0)**2
        do j=1,nocc(1)
          do k=1,nocc(2)
            r=corr(j,k)
            s2=s2-r**2
          enddo
        enddo
        if (nodezero.and.util_print('ground state',print_default)) then
          write(LuOut,9180) s2
          call util_flush(LuOut)
        endif
      endif
 9180 format(2x,'<S2> = ',f8.4)
c
c Excited states
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
c Identify symmetry of the excited state
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
                pol=ipol
                occ=j
                vir=k
              endif
            enddo
          enddo
        enddo
        call sym_product(geom,isym,
     1    int_mb(k_irs(pol)+occ-1),jsym)
        call sym_product(geom,jsym,
     1    int_mb(k_irs(pol)+vir-1),jsym)
        call sym_irrepname(geom,jsym,oname)
        if ((targetsym.eq.0).or.(targetsym.eq.jsym)) then
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
c Compute <S**2> expectation value
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
              if (dabs(xa1+ya1).lt.1.0d-5) continue
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
                      rc1=corr(ka2,jb1)
                      rc2=corr(ka1,jb1)
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
                      rc1=corr(ja1,jb2)
                      rc2=corr(ka1,kb2)
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
                    rc1=corr(ja3,jb1)
                    s2=s2-(xa1-ya1)*(xa1+ya1)*rc1*rc1
                  enddo
                endif
              enddo
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
              if (dabs(xb1+yb1).lt.1.0d-5) continue
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
                      rc1=corr(ja1,kb1)
                      rc2=corr(ja1,kb2)
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
                      rc1=corr(ja2,jb1)
                      rc2=corr(ka2,kb1)
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
                    rc1=corr(ja1,jb3)
                    s2=s2-(xb1-yb1)*(xb1+yb1)*rc1*rc1
                  enddo
                endif
              enddo
            enddo
          enddo
          if (nodezero.and.util_print('excited state',print_default))
     1      write(LuOut,9180) s2
        endif
        if (nodezero.and.util_print('excited state',print_default))
     1    write(LuOut,9190)
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
     1    nocc,nmo,nfc,nfv,ipol,g_tmomd,1)
        do i=1,ipol
          call multipole_density(ao_bas_han,cntr,1,
     1      g_tmomd(i),tmom,4)
        enddo
        if (.not.tda) then
          call tddft_transfm(n,g_y,g_movecs,nbf_ao,
     1      nocc,nmo,nfc,nfv,ipol,g_tmomd,1)
          do i=1,ipol
            call multipole_density(ao_bas_han,cntr,1,
     1        g_tmomd(i),tmom,4)
          enddo
        endif
        if (ipol.eq.1) then
          do i=2,4
            tmom(i)=tmom(i)*dsqrt(2.0d0)
          enddo
        endif
        if (nodezero.and.
     1    util_print('excited state',print_default)) then
          if ((singlet.and.(mult.ne.1)).or.
     1        (triplet.and.(mult.ne.3))) then
            write(LuOut,9210)
          else
            write(LuOut,9200) tmom(2),tmom(3),tmom(4),
     1      2.0d0/3.0d0*(tmom(2)**2+tmom(3)**2+tmom(4)**2)*apbval(n)
          endif
        endif
c
c Print dominant components of X & Y
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
        call util_flush(LuOut)
      enddo
      if (m.lt.target) 
     1  call errquit('tddft_analysis: target root not found',m)
      if (nodezero.and.
     1  util_print('excited state',print_default)) write(LuOut,*)
 9110 format(2x,
     1  '-------------------------------------------------------',
     2  /,2x,'Root',i3,1x,a7,1x,a4,f13.9,' a.u. (',f11.7,' eV)')
 9190 format(2x,
     1  '-------------------------------------------------------')
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
        if (targetsym.eq.0) then
          write(LuOut,9350)
        else
          write(LuOut,9340) targetsym
        endif
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
     1  call errquit('tddft_analysis: failed to put tddft energy',0)
 9300 format('              Target root =',i7)
 9340 format('          Target symmetry =',i7)
 9350 format('          Target symmetry = not specified')
 9310 format('      Ground state energy =',f20.12)
 9320 format('        Excitation energy =',f20.12)
 9330 format('     Excited state energy =',f20.12)
c
c ------
c Return
c ------
c
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_y(i))) call errquit
     1    ('tddft_analysis: failed to deallocate y',0)
        if (.not.ma_pop_stack(l_x(i))) call errquit
     1    ('tddft_analysis: failed to deallocate x',0)
      enddo
      if (ipol.eq.2) then
        if (.not.ga_destroy(g_corr)) call errquit
     1    ('tddft_analysis: failed to destroy g_corr',0)
      endif
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
