      subroutine dft_mpole( geom, i_basis, j_basis, lmax_pole, 
     &     iga_dens, ndens, rmoments, len_rmoments)
C     $Id: dft_mpole.F,v 1.5 1995-05-08 19:46:44 og845 Exp $
      implicit none

#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
      double precision au2Debye
      parameter(au2Debye=2.541d0)
      double precision ZERO
      parameter(ZERO=0.d0)
C     
      Integer NxTask ! From TCGMSG
      External NxTask
c     
c     Compute the multipole integrals and contract against the input
c     densities.  Return the n-pole moments, n=1,...,lmax_pole, in the
c     array rmoments
c     
c     This version computes the full square of integrals and should work
c     OK even if i_basis != j_basis.
c-----------------------------------------------------------------------
c     
c     arguments
c     
      Integer geom ! geometry handle
      integer i_basis, j_basis ! [input] bra and ket basis sets
      integer ndens, len_rmoments
      integer iga_dens(ndens) ! [input] array of GA handles to densities
      integer lmax_pole ! [input] up to 2^{lmax_pole} multipoles
      Double Precision centerl(3) ! [input] expansion center for multipoles
      Double Precision rmoments(len_rmoments)
c     ! [output] cartesian 2^{lmax_pole}-pole moments
c-----------------------------------------------------------------------
C     
C     Parameters
C     
      Integer Task_Chunk, MAX_LMAX_POLE,
     &     MAX_MPOLE_CARTESIAN_COMPONENTS
      Double Precision ONE, TWO

      Parameter (Task_Chunk = 10, ! Large enough to minimize overhead?
     &     MAX_LMAX_POLE=1, ! dipoles only
     &     MAX_MPOLE_CARTESIAN_COMPONENTS=
     &     (MAX_LMAX_POLE+1)*(MAX_LMAX_POLE+2)/2 )

      Parameter ( ONE=1.d0, TWO=2.d0 )
c-----------------------------------------------------------------------
c     
c     local variables
c     
      Integer nmvals
      Integer nshell_i, nshell_j
      Integer ishell, jshell, iproc, nproc, Next
      Integer ijshell, ilo, ihi, jlo, jhi
      Integer  l_scr
      Integer  k_scr
      Character*8 mpole_ints_name 
      Character*1 coord_label(MAX_MPOLE_CARTESIAN_COMPONENTS)

      Integer isphere, jsphere, nshbf_j, nprimo_i, l_i, k_mp
      Integer nbf_shell_i,  nbf_shell_j
      Integer l_mp, nprimo_j, l_j, ioff
      Integer idens, n_ij_pairs, nshbf_i
      Integer lmpint, k_mp_reord
      Integer mscratch_3ov
      Integer maxbuf 
      Integer max_ij_pairs, lval, l_mp_reord, lmpmax,  lmax_i
      Integer i_nbf_max, j_nbf_max, nbf_ao
      Integer n_operators, lmax_j
      
c**** Integer iga_mpole(MAX_MPOLE_CARTESIAN_COMPONENTS) ! GA for m-pole integrals
      Integer iga_mpole(3)
      Integer ncenter ! number of atomic centers in geometry object

      Logical LResult

      Integer ixyz

C      Integer ii  ! debug
c-----------------------------------------------------------------------
c**** 
c**** Check input:  only dipoles for now
c**** 
      if ( lmax_pole .gt. MAX_LMAX_POLE) then
        call errquit(
     &       'dft_mpole: only dipoles (lmax_pole.eq.1) implemented',
     &       3)        
      endif 
c****
c**** Nuclear contribution to moments
c****
      if ( .not. geom_ncent(geom,ncenter)) call errquit(
     &     'dft_mpole:  geom_ncent failed',90)

      call dzero(rmoments,len_rmoments)
      call dft_nuc_mpole(geom, ncenter, lmax_pole, centerl, 
     &                   rmoments, len_rmoments)
C      write (6,'(a,3e15.5)')
C     &     'dft_mpole: nuclear part of rmoments: ',
C     &     (rmoments(ii),ii=1,3)
c**** 
c**** number of contractions in left and right bases
c**** 
      if (.not. bas_numcont(i_basis, nshell_i)) call errquit
     $     ('dft_mpole: bas_numcont failed for i_basis', i_basis)
      if (.not. bas_numcont(j_basis, nshell_j)) call errquit
     $     ('dft_mpole: bas_numcont failed for j_basis', j_basis)
c**** 
c**** allocate GA arrays for multipole integrals (dipole for now)
c**** 
      coord_label(1) = 'x'
      coord_label(2) = 'y'
      coord_label(3) = 'z'

      do ixyz = 1,3
        mpole_ints_name = 'dipole_'//coord_label(ixyz)
        if (.not. bas_numbf(i_basis, nbf_ao)) call errquit
     &       ('dft_moments: bas_nbf failed', 13)

        call gacrea( geom, i_basis,
     &               nbf_ao, nbf_ao, mpole_ints_name,
     &               iga_mpole(ixyz),'atom')
        call ga_zero(iga_mpole(ixyz))

      enddo 
c**** 
c**** estimate memory requirements for int_mpole call; allocate
c**** 2 buffers (need to sort)
c**** 
      call int_mem_3ov(maxbuf,mscratch_3ov)

      if (.not. bas_high_angular(i_basis,lmax_i)) call errquit(
     &     'dft_mpole: bas_high_angular call failed',15)
      if (.not. bas_high_angular(j_basis,lmax_j)) call errquit(
     &     'dft_mpole: bas_high_angular call failed',16)
c**** 
c**** length of int_mpole integral output for full square list
c**** includes l_pole = 0,...,lmax_pole, where l_pole = 0 is simply
c**** the 2-c overlap matrix.  (Cartesian components).
c**** 
      n_operators = 0 ! number of cartesian 2^{l-pole} operators
      do lval = 0,lmax_pole
        n_operators = n_operators + (lval+1)*(lval+2)/2
      enddo 
      
      if ( .not. bas_nbf_cn_max(i_basis,i_nbf_max)) call errquit(
     &     'dft_mpole: bas_nbf_cn_max failed',20)
      if ( .not. bas_nbf_cn_max(j_basis,j_nbf_max)) call errquit(
     &     'dft_mpole: bas_nbf_cn_max failed',22)

      max_ij_pairs = i_nbf_max * j_nbf_max

      lmpmax = n_operators * max_ij_pairs

c     
c     allocate necessary local temporary arrays on the stack
c     
      if (.not. MA_push_get( MT_DBL, lmpmax,'dft_mpole:MP_reord',
     &     l_MP_reord, k_MP_reord))
     &     call errquit('dft_mpole: ma failed', lmpmax)

      if (.not. MA_push_get( MT_DBL, lmpmax,'dft_mpole:MP',
     &     l_MP, k_MP))
     &     call errquit('dft_mpole: ma failed', lmpmax)

      if (.not. MA_push_get( MT_DBL, mscratch_3ov,'dft_mpole:scr30v',
     &     l_scr, k_scr))
     &     call errquit('dft_mpole: ma failed', mscratch_3ov)
c     
c     Loop thru shells with static parallel work decomposition
c     
      LResult = .true.
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      Next = NxTask( nproc, Task_Chunk)
C     
      do ishell = 1, nshell_i
c****   
c****   get i_basis info
c****   
        if (.not. bas_cn2bfr(i_basis, ishell, ilo, ihi)) call errquit
     &       ('dft_mpole: bas_cn2bfr failed for i_basis', i_basis)
        if (.not. bas_continfo(i_basis, ishell,
     &       l_i, nprimo_i, nshbf_i, isphere))call errquit
     &       ('dft_mpole: bas_continfo failed for i_basis', i_basis)
        nbf_shell_i = ihi - ilo + 1

        if ( isphere .ne. 0 ) call errquit(
     &       'dft_mpole: only coded for cartesian bases',42)

        do jshell = 1, nshell_j
          if ( ijshell .eq. Next) then
c****       
c****       get j_basis info
c****       
            if (.not. bas_cn2bfr(j_basis, jshell, jlo, jhi))
     &           call errquit('dft_mpole: bas_cn2bfr', j_basis)
            if (.not. bas_continfo(j_basis, jshell,
     &           l_j, nprimo_j, nshbf_j, jsphere))
     &           call errquit('dft_mpole: bas_cn2bfr', j_basis)
            if ( jsphere .ne. 0 ) call errquit(
     &           'dft_mpole: only coded for cartesian bases',43)

            nbf_shell_j = jhi - jlo + 1
c****       
c****       No sp-orbitals yet!
c****       
            if (l_i .lt. 0 .or. l_j .lt. 0) then
              call errquit(
     &             'int_eval_gbsets: sp-type orbital not coded', 1)
            endif 
c****       
c****       Generate the integrals
c****       subroutine int_mpole(i_basis, ish, j_basis, jsh, lmax, centerl,
c****       &              lscr, scr, lmpint, MP)
c****       
            n_ij_pairs = nbf_shell_i * nbf_shell_j
            lmpint = n_ij_pairs * n_operators

            call int_mpole(i_basis, ishell, j_basis, jshell,
     &           lmax_pole, centerl, mscratch_3ov, dbl_mb(k_scr),
     &           lmpint, dbl_mb(k_MP))

c****       
c****       output from int_mpole is: overlap, dipole, q-pole, ...
c****       within a multipole block, the order is
c****       MP order <i|j|k>  k fastest, then j, then i
c****       sort integrals in dipole, q-pole, etc. blocks from i
c****       nt_mpole order to operator-index-slowest order
c****       
            ioff = n_ij_pairs                       ! skip 2-c overlap elements
            call dzero(dbl_mb(k_MP_reord), lmpmax)

c****       copy overlap block
            call dcopy( n_ij_pairs, dbl_mb(k_MP), 1,
     &                              dbl_mb(k_MP_reord), 1)

            do lval = 1,lmax_pole
              nmvals = (lval+1)*(lval+2)/2
              call dft_sort23_3( dbl_mb(k_MP + ioff),
     &             dbl_mb(k_MP_reord + ioff), nbf_shell_j,
     &             nmvals, nbf_shell_i)
              ioff = ioff + n_ij_pairs*nmvals
            enddo 
c****       
c****       Sum the integrals into the global array (one GA for each operator)
c****       
            ioff = n_ij_pairs ! skip l_pole=0 part of buffer

            do ixyz = 1,3
              call ga_acc( iga_mpole(ixyz), jlo, jhi, ilo, ihi,
     &                     dbl_mb(k_MP_reord+ioff), nbf_shell_j, ONE)
              ioff = ioff + n_ij_pairs
            enddo 

            Next = NxTask( NProc, Task_Chunk)
          endif
          ijshell = ijshell + 1
        enddo
      enddo
      Next = NxTask( -NProc, Task_Chunk)
c**** 
c**** clean up stack
c**** 
      if (.not. MA_pop_stack(l_scr)) call errquit
     &     ('dft_mpole: pop failed', 0)
c     
      if (.not. MA_pop_stack(l_MP)) call errquit
     &     ('dft_mpole: pop failed', 0)
c     
      if (.not. MA_pop_stack(l_MP_reord)) call errquit
     &     ('dft_mpole: pop failed', 0)
c     
      call ga_sync() ! So that no nasty races can result
c     
c**** 
c**** contract multipole integrals against densities
c**** double precision function ga_ddot(g_a, g_b)
c**** 
c**** debug
c****      write (6,*) 
c****      write (6,*) 'iga_mpole(1):'
c****      call ga_print(iga_mpole(1))
c****      write (6,*) 
c****      write (6,*) 'iga_mpole(2):'
c****      call ga_print(iga_mpole(2))
c****      write (6,*) 
c****      write (6,*) 'iga_mpole(3):'
c****      call ga_print(iga_mpole(3))
c****      write (6,*) 
c****      write (6,*) 'iga_dens(1):'
c****      call ga_print(iga_dens(1))
c**** debug
      do ixyz = 1,3
        rmoments(ixyz) = ZERO
        do idens = 1,ndens
          rmoments(ixyz) = rmoments(ixyz) +
     &         ga_ddot(iga_mpole(ixyz),iga_dens(idens))
C          WRITE (6,'(A,I5,I5,3E15.5)') 'IXYZ, IDENS, RMOMENTS: ',
C     &         ixyz, idens, (rmoments(ii),ii=1,3)
        enddo 
      enddo 
c**** 
c**** destroy multipole matrices
c**** 
      do ixyz = 1,3
        if (.not. ga_destroy(iga_mpole(ixyz)) ) call errquit(
     &       'dft_mpole:  ga_destroy failed ',ixyz)
      enddo
 
      if(ga_nodeid().eq.0) then
        call dft_header('Dipole Moments (debye):       ')
        write(LuOut,'(3A12)') 'x','y','z'
        write (LuOut,'(3f12.5)') (rmoments(ixyz)*au2Debye,ixyz=1,3)
        write(LuOut,*)
      endif
        
      return
      end
c**** 
c**** 
c**** 
      subroutine dft_sort23_3( MP, MP_reord, imax, jmax, kmax)
c**** 
c**** reorder a three-indexed rectangular array by exchanging the
c**** two slowest indices
c**** 
      Implicit none
c-----------------------------------------------------------------------
c**** arguments

      Integer imax, jmax, kmax
      Double Precision MP(imax,jmax,kmax), MP_reord(imax,kmax,jmax)
c-----------------------------------------------------------------------
c**** locals

      Integer i, j, k
c-----------------------------------------------------------------------

      do k = 1,kmax
        do j = 1,jmax
          do i = 1,imax
            MP_reord(i,k,j) = MP(i,j,k)
          enddo 
        enddo 
      enddo 
      return
      end
c**** 
c**** 
c**** 
      subroutine  dft_nuc_mpole(geom, ncenter, lmax_pole, centerl, 
     &                          rmoments, len_rmoments)
c**** 
c**** nuclear contribution to multipole moments (dipole for now)
c**** 
      Implicit none

#include "bas.fh"
#include "geom.fh"
c-----------------------------------------------------------------------
c**** 
c**** parameters
c**** 
      Double Precision ZERO, ONE, THREE
      Parameter ( ZERO=0.d0, ONE=1.d0, THREE=3.d0 )
c-----------------------------------------------------------------------
c**** 
c**** arguments
c**** 
      Integer geom ! [input] geometry handle
      Integer ncenter ! number of centers
      Integer lmax_pole
      Double Precision centerl(3)
      Integer len_rmoments
      Double Precision rmoments(len_rmoments)
c-----------------------------------------------------------------------
c**** 
c**** locals
c**** 
      Integer iquad, ixyz, jxyz, kxyz, icent
      Character*16 tag
      Double Precision coord(3), q, rsquare, delta, total_nuc_charge

c-----------------------------------------------------------------------
c**** 
c**** compute center of charge
c**** 
      do ixyz = 1,3
        centerl(ixyz) = ZERO
      enddo 
      do icent = 1,ncenter
        if (.not. geom_cent_get(geom,icent,tag,coord,q))
     &       call errquit(
     &       'dft_mpole:dft_nuc_mpole:  geom_cent_get failed',10)
c****   
c****   contribution to center of charge vector
c****   
        do ixyz = 1,3
          centerl(ixyz) = centerl(ixyz) + q*coord(ixyz)
        enddo 
      enddo 
      if (.not. geom_nuc_charge(geom,total_nuc_charge)) call errquit(
     &     'dft_mpole:dft_nuc_mpole: geom_nuc_charge failed',20)
      do ixyz = 1,3
        centerl(ixyz) = centerl(ixyz)/total_nuc_charge
      enddo 
C      write (6,*) 'dft_mpole: center of charge: ',centerl
c**** 
c**** loop over centers, summing multipole moments
c**** 
      do icent = 1,ncenter
        if (.not. geom_cent_get(geom,icent,tag,coord,q))
     &       call errquit(
     &       'dft_mpole:dft_nuc_mpole:  geom_cent_get failed',10)
c****   
c****   dipole moment contributions
c****   
        do ixyz = 1,3
          rmoments(ixyz) = rmoments(ixyz) +
     &         q * (coord(ixyz)-centerl(ixyz))
        enddo 
c****   
c****   traceless quadrupole term (if required)
c****   ordering is xx xy xz yy yz zz
c****   
        if ( lmax_pole .gt. 1 ) then
          iquad = 0
          do ixyz = 1,3
            do jxyz = ixyz+1,3
              iquad = iquad+1

              if (ixyz.eq.jxyz ) then
                delta = ONE
                rsquare = ZERO
                do kxyz = 1,3
                  rsquare = rsquare +
     &                 (coord(kxyz)-centerl(kxyz))*
     &                 (coord(kxyz)-centerl(kxyz))
                enddo 
              else 
                delta = ZERO
              endif 

              rmoments(iquad) = rmoments(iquad) +
     &             q* (THREE*(coord(ixyz)-centerl(ixyz))
     &             * (coord(jxyz)-centerl(jxyz))
     &             - rsquare*delta)
            enddo 
          enddo 
        endif 
      enddo 
      return
      end 




