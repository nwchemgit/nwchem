      subroutine aoresponse_giao_rhs (rtdb, basis, geom, 
     &   ncomp, g_rhs,
     &   g_vectors, g_smat1, g_sket1, g_dens, g_amat,
     &   nocc, nvir, nmo, nbf, froct, epst, omega, lstatic)

c $Id$ 

c     ==================================================================
c     
c     purpose: compute giao terms for the right-hand-side of
c              the Ax = b linear response equation system
c
c     called from: aoresponse, aoresponse_beta
c     
c     note: it is assumed that the integral code is initialized
c           prior to calling this routine.
c
c     note: most of this was taken from hnd_giaox without significant
c           modifications

c
c     ==================================================================
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
c
c     subroutine arguments:

      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer ncomp  ! [input] number of Fourier components
      integer g_rhs(ncomp)   ! [output] part of CPKS rhs vectors
      integer g_vectors ! [input] MO vectors for given spin
      integer g_dens ! [input] density matrix
      integer g_amat(ncomp) ! [in/out] A-matrices
      integer g_smat1 ! [input] perturbed S-matrix
      integer g_sket1  ! [input] "half perturbed" S-matrix

      integer nbf, nocc, nvir, nmo
      double precision froct(nbf), epst(nbf), omega

      logical lstatic ! we can skip all S-ket related terms if true

c     local variables:

      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)

      integer g_u, g_fock, g_d1, g_temp

      integer iocc, ifld, ipm

      double precision tol2e
      double precision jfac(3),kfac(3), xfac
      character*3 scftyp

      logical debug, oskel, dbgmat

c     ==================================================================

      debug = .false. .and. ga_nodeid().eq.0 ! .true. during development
      dbgmat = .false. .and. ga_nodeid().eq.0 ! debug large matrices

c     THERE WAS A BUG HERE BECAUSE ONLY IF DEBUG WAS SET TO TRUE DID I
c     GET THE CORRECT RESULTS FOR ORBETA. the code got stuck where the
c     fock matrix is printed if I ran on more than one processor;
c     another sign of a problem.
c     The problem seems to have been with the Fock matrix from
c     the occ-occ block, not the other pieces. This Fock matrix ended
c     up being zero

c     For the time being (01/2011) compile with FOPTIMIZE=-O1 or -O0
c     GNUmakefile has bee changed so that by default this routine
c     is not optimized

c      if (debug) write (luout,*) 'hello from aoresponse_giao_rhs'
      
c     Current CPHF does not handle symmetry 
      
      oskel = .false.
      
c     If DFT get part of the exact exchange defined
      
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()

      if (debug) write (luout,*) 'use_theory, xfac =', use_theory, xfac
      
      call scf_get_fock_param(rtdb,tol2e)

      if (debug) write (luout,*) 'nocc, nvir, nmo, nbf',
     &   nocc, nvir, nmo, nbf
      
c     Create workspace matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c     
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nocc
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'aor U(1)',alo,g_u)) call 
     &    errquit('aorgiao: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u)

c     transform Smat1 and S-ket to MO basis

      call giao_aotomo(g_smat1,g_vectors,nocc,nvir,1,3,nbf)
      if (.not.lstatic)
     &   call giao_aotomo(g_sket1,g_vectors,nocc,nvir,1,3,nbf)

c     ga_rhs(a,i) = ga_rhs(a,i) - e(i) * S10(a,i)
c     Scale (occ,virt) block of g_smat1 with - (minus) MO energies 

      alo(1) = nocc+1
      ahi(1) = nmo
      alo(3) = 1
      ahi(3) = 3
      do iocc = 1, nocc
        alo(2) = iocc
        ahi(2) = iocc
        call nga_scale_patch(g_smat1,alo,ahi,-epst(iocc)) 
      enddo

c     Copy to ga_rhs 
c     alo(1) and ahi(1) the same as before
c
      alo(2) = 1
      ahi(2) = nocc

c     these dimensions will be used repeatedly:
      blo(1) = 1
      bhi(1) = nocc * nvir
      blo(2) = 1
      bhi(2) = 3

      do ipm = 1,ncomp
        call nga_copy_patch('n',g_smat1,alo,ahi,g_rhs(ipm),blo,bhi)
      enddo

      if (debug) write (luout,*) 'added Smat1 term to rhs'


c     Construct occ-occ part of the three A matrices
c     Occ-occ blocks for each field direction are defined as -1/2 S10
c     Scale (occ,occ) block of g_smat1 with -1/2 and copy to g_u
c
c     alo(2) and ahi(2) will stay as 1 and nocc for a while
c
      alo(1) = 1
      ahi(1) = nocc
      call nga_scale_patch(g_smat1,alo,ahi,-0.5d0)
      call nga_copy_patch('n',g_smat1,alo,ahi,g_u,alo,ahi)

c     note: the calling routine needs the same information in g_amat
c     occ-occ blocks

      do ipm = 1,ncomp
        call nga_copy_patch('n',g_smat1,alo,ahi,g_amat(ipm),alo,ahi)
      enddo

      if (debug) write(luout,*) 'initial setup of g_u, g_amat  complete'

c     We also need the occupied-occupied contribution of
c     g_u contributing
c     to the first order density matrix. As this block does not change 
c     during the CPHF we can calculate it once and subtract it from the
c     RHS. We will reuse g_s10 as scratch space.

c     Note: For consistency, in this case the density matrix build
c     in the CPHF must also not include the occ-occ block
c
      call ga_zero(g_smat1)

      clo(1) = 3
      clo(2) = nbf
      clo(3) = nbf
      chi(1) = 1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'aor_giao F1',chi,g_fock)) call 
     &    errquit('aorgiao: nga_create failed g_fock',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,clo,'D10 matrix',chi,g_d1)) call 
     &    errquit('aorgiao: nga_create failed g_d1',0,GA_ERR)

      call ga_zero(g_fock)
      call ga_zero(g_d1)


      alo(1) = 1
      alo(2) = 1
      blo(1) = 1
      blo(2) = 1
      bhi(1) = nbf
      clo(2) = 1
      clo(3) = 1
      chi(2) = nbf
      chi(3) = nbf
      dlo(1) = 1
      dlo(2) = 1
      dhi(1) = nbf
      dhi(2) = nocc

c
c     Create perturbed density matrix for closed-closed g_u block.
c     We do this in a loop over the directions of the perturbing field
c
      do ifld = 1, 3
        
        jfac(ifld) = 0.0d0
        kfac(ifld) = -1.0d0*xfac
        alo(3) = ifld
        ahi(3) = ifld
        clo(1) = ifld
        chi(1) = ifld
        dlo(3) = ifld
        dhi(3) = ifld
        ahi(1) = nmo
        ahi(2) = nocc
        bhi(2) = nmo 
        call nga_matmul_patch('n','n',1.0d0,0.0d0,g_vectors,blo,bhi,  
     &     g_u,alo,ahi,g_smat1,dlo,dhi)  
        if (debug) write (6,*) 'occ-occ perturbed density matrix step1'
        ahi(2) = nbf
        ahi(1) = nocc
        bhi(2) = nocc
        
c       
c       Minus sign as we will subtract it from the RHS and not consider
c       it in the LHS
c      
c        if (debug) then
c          write (6,*) blo(1:3)
c          write (6,*) bhi(1:3)
c          write (6,*) alo(1:3)
c          write (6,*) ahi(1:3)
c          write (6,*) clo(1:3)
c          write (6,*) chi(1:3)
c        endif 
        call nga_matmul_patch('n','t',-1.0d0,0.0d0,g_vectors,blo,bhi,
     &     g_smat1,alo,ahi,g_d1,clo,chi)  

      enddo                     ! ifld = 1,3

      if (debug) write(6,*) 'occ-occ perturbed density matrix built'
      
c     
c     Build first-order fock matrix in AO basis due to occ-occ terms
c     



c from old (2006) hnd_giaox:

c$$$      if (use_theory.eq.'dft') then
c$$$        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
c$$$     $     call errquit('aorgiao: rtdb_put of xc_active failed',0,
c$$$     &     RTDB_ERR)
c$$$        if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
c$$$     $     call errquit('aorgiao: rtdb_put of calc_type failed',0,
c$$$     &     RTDB_ERR)
c$$$        if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
c$$$     $     call errquit('aorgiao: rtdb_put of j_derfit failed',0,
c$$$     &     RTDB_ERR)
c$$$      endif
c$$$      call shell_fock_build(geom, basis, 0, 3,
c$$$     $   jfac, kfac,tol2e, g_d1, g_fock,.false.)
c$$$      if(use_theory.eq.'dft') then
c$$$        ifld = 0
c$$$        if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
c$$$     $     call errquit('aorgiao: rtdb_put failed',0,RTDB_ERR)
c$$$        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
c$$$     $     call errquit('aorgiao: rtdb_put of xc_active failed',0,
c$$$     &     RTDB_ERR)
c$$$      endif

c from new (2011) hnd_giaox:
      if (use_theory.eq.'dft') then
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('aor_giao: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('aor_giao: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('aor_giao: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
      endif
      call shell_fock_build(geom, basis, 0, 3,
     $     jfac, kfac,tol2e, g_d1, g_fock,.false.)
      if(use_theory.eq.'dft') then
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 0))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif

      if (debug) then
        write(6,*) 'call to shell_fock_build complete'
        write (luout,*) 'g_fock =',g_fock
        call ga_check_handle(g_fock, 'check g_fock')
        call ga_print(g_fock)
      end if
c     
      
c     
c     Transform to the occ-virt MO basis and add to RHS
c     
      call ga_zero(g_smat1)
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nocc
      do ifld = 1, 3
        alo(3) = ifld
        ahi(3) = ifld
        clo(1) = ifld
        chi(1) = ifld
        call nga_matmul_patch('n','n',2.0d0,0.0d0,
     $     g_fock,clo,chi,
     $     g_vectors,alo,ahi,
     $     g_smat1,alo,ahi)
      enddo
      
      call ga_zero(g_fock)
      clo(2) = nocc+1
      clo(3) = 1
      chi(2) = nmo
      chi(3) = nocc
      do ifld = 1, 3
        blo(1) = nocc+1
        blo(2) = 1
        bhi(1) = nmo
        bhi(2) = nbf
        alo(3) = ifld
        ahi(3) = ifld
        clo(1) = ifld
        chi(1) = ifld
        call nga_matmul_patch('t','n',1.0d0,0.0d0,
     $     g_vectors, blo,bhi,
     $     g_smat1, alo,ahi,
     $     g_fock, clo,chi )

      if (debug) write(6,*) 'occ-vir F(1) build complete'

        blo(1) = 1
        bhi(1) = nocc*nvir
        blo(2) = ifld
        bhi(2) = ifld

        do ipm = 1,ncomp
        call nga_add_patch(1.0d0,g_rhs(ipm),blo,bhi,
     &       1.0d0,g_fock,clo,chi,
     &     g_rhs,blo,bhi)
        enddo

        if (debug) write(6,*) 'added F(1) occocc to rhs'
        if (debug) then
          write (luout,*) 'printing Fock matrix ifld=,',ifld
          call nga_print_patch(g_fock,clo,chi,1)
        end if

      enddo                     ! ifld
c     
c     Cleanup of g_d1 and g_fock, not needed for now
c     
      if (.not.ga_destroy(g_d1)) call 
     &   errquit('aorgiao: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(g_fock)) call 
     &   errquit('aorgiao: ga_destroy failed g_fock',0,GA_ERR)
      
c     The Remaining term is the Perturbed (GIAO) two-electron term
c     times Unperturbed density matrix
c     Calculate Sum(r,s) D0(r,s) * G10(m,n,r,s) in AO basis


c old code from 2006 giaox: 
     
c$$$      alo(1) = nbf
c$$$      alo(2) = -1
c$$$      alo(3) = -1
c$$$      ahi(1) = nbf
c$$$      ahi(2) = nbf
c$$$      ahi(3) = 3
c$$$      blo(1) = 1
c$$$      bhi(1) = nocc*nvir
c$$$      blo(2) = 1
c$$$      bhi(2) = 3
c$$$
c$$$
c$$$
c$$$      if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',alo,g_fock)) call 
c$$$     &   errquit('aorgiao: nga_create failed g_fock',0,GA_ERR)
c$$$      call ga_zero(g_fock)
c$$$      if(use_theory.eq.'dft') then
c$$$        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
c$$$     $     call errquit('aorgiao: rtdb_put of xc_active failed',0,
c$$$     &     RTDB_ERR)
c$$$        ifld = 4
c$$$        if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
c$$$     $     call errquit('aorgiao: rtdb_put failed',0,RTDB_ERR)
c$$$      endif
c$$$      call new_giao_2e(geom, basis, nbf, tol2e, g_dens, g_fock, xfac)
c$$$      if(use_theory.eq.'dft') then
c$$$        ifld = 0
c$$$        if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
c$$$     $     call errquit('aorgiao: rtdb_put failed',0,RTDB_ERR)
c$$$        if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
c$$$     $     call errquit('aorgiao: rtdb_put of xc_active failed',0,
c$$$     &     RTDB_ERR)
c$$$      endif


c new code from 2011 giaox:

      alo(1) = -1 
      alo(2) = -1
      alo(3) = 1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'aor_giao F2',alo,g_fock)) call 
     &    errquit('aor_giao: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
      if(use_theory.eq.'dft') then
         ifld = 4
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('aor_giao: rtdb_put failed',0,RTDB_ERR)
      endif
c     call new_giao_2e(geom, basis, nbf, tol2e, g_dens, g_fock, xfac)
      call new_giao_2e(geom, basis, nbf, tol2e, g_dens, g_fock, xfac, 1) ! FA restrict calc
      if(use_theory.eq.'dft') then
         ifld = 0
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('aor_giao: rtdb_put failed',0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $       call errquit('aor_giao: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      endif

      if (debug) write(6,*) 'Fmat derivatives calculated'
      
c     
c     Transform to MO basis and add to right-hand-side
c     
      call giao_aotomo(g_fock,g_vectors,nocc,nvir,1,3,nbf)
      if (debug) write(6,*) 'Fmat derivatives transformed to MO basis'
      alo(1) = nocc+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nocc
      alo(3) = 1
      ahi(3) = 3

      blo(1) = 1
      bhi(1) = nocc * nvir
      blo(2) = 1
      bhi(2) = 3

      do ipm = 1,ncomp
      call nga_add_patch(1.0d0,g_rhs(ipm),blo,bhi,1.0d0,g_fock,alo,ahi,
     &   g_rhs(ipm),blo,bhi)
      enddo

      if (.not.ga_destroy(g_fock)) call 
     &   errquit('aorgiao: ga_destroy failed g_fock',0,GA_ERR)
      if (.not.ga_destroy(g_u)) call 
     &   errquit('aorgiao: ga_destroy failed g_u',0,GA_ERR)

c     Add S-ket matrix contribution (now in MO basis) to rhs
c     Note that in JCP 123, 114103 S-Ket contains a factor of (+/-)w
c     which we account for here
c     THE PREFACTORS NEED TO BE DEBUGGED YET, THIS IS NOT WORKING

c     we will use g_smat1 as scratch space since it is not needed
c     anymore


      if (.not.lstatic) then
        call ga_zero(g_smat1)
        
        do ifld = 1,3
          
          alo(1) = 1
          ahi(1) = nmo
          alo(2) = 1
          ahi(2) = nmo
          alo(3) = ifld
          ahi(3) = ifld
          
          call nga_copy_patch('t',g_sket1,alo,ahi,g_smat1,alo,ahi)
        end do                  ! ifld
        
        
        alo(1) = nocc+1
        ahi(1) = nmo
        alo(2) = 1
        ahi(2) = nocc
        alo(3) = 1
        ahi(3) = 3
        
        blo(1) = 1
        bhi(1) = nocc * nvir
        blo(2) = 1
        bhi(2) = 3
        
        if (ncomp.ne.2) then
          call errquit('aor_giao_rhs: adding S1ket but ncomp.ne.2',
     &       0, INPUT_ERR)
        end if
        
c       component 1: use sket1 vir-occ block unmodified
        
        call nga_add_patch(-1.0d0*omega,g_sket1,
     &     alo,ahi,1d0,g_rhs(1),blo,bhi,
     &     g_rhs(1),blo,bhi)
        
c       component 2: use the transpose of sket 
        
c       jochen: below used to be smat1 , not sket1.
c       try what happens with sket1
        
        if (ncomp.eq.2) then
          call nga_add_patch(-1.0d0*omega,g_sket1,
     &       alo,ahi,1d0,g_rhs(2),blo,bhi,
     &       g_rhs(2),blo,bhi)             
        endif
        
        if (debug) write (luout,*) 'added Sket1 term to rhs'
        
      end if ! lstatic

      if (debug) write(6,*) 'rhs build complete'

      if (debug) then
        do ipm = 1,ncomp
          call ga_print(g_rhs(ipm)) 
          call ga_print(g_amat(ipm))
        end do
        call ga_print( g_vectors)
        call ga_print( g_dens) 
        call ga_print( g_smat1)
        call ga_print(g_sket1)
      end if


c     ==================================================================
      return

      end
