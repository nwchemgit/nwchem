c $Id: et_calc.F,v 1.4 2001-12-31 17:50:28 windus Exp $
      logical function et_calc(rtdb,Vba)
      implicit none
c
c     This function actually calculates Vba
*
* $Id: et_calc.F,v 1.4 2001-12-31 17:50:28 windus Exp $
c
*
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
c
c     Generic NWChem interface to compute Vba. SCF of DFT determinants
c     can be used.
c
c     RTDB output parameters
c     ----------------------
c     task:status (logical)- T/F for success/failure
c     if (status) then
c     .  task:et (real)   - Vba (electron transfer matrix element)
c     .  task:cputime (real)  - cpu time to execute the task
c     .  task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
      integer rtdb
      integer i, j, k
      integer l_vals,k_vals
      integer g_tmp,g_s,g_ahandle(2),g_bhandle(2),g_d,g_temp
      integer nocc(2)
      integer g_u,g_vt,g_v,g_phandle(2)
      integer g_av,g_bu,g_t,g_h1
      integer g_avecsa,g_avecsb,g_bvecsa,g_bvecsb
c
c     logical et_movecs_read
c     external et_movecs_read
      logical ga_copy
      external ga_copy
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      character*2 pname(2)
      double precision Vba
      double precision prod(2)
      double precision Sba,pa_val,pb_val,omega1,h1_val
      double precision detu(2),detv(2),detprod,omega2
      double precision haa,hbb,hba,val, enrep, dum
c
c========================
c     for printing:
c========================
c
      call util_print_push
      call util_print_rtdb_load(rtdb,"et")
c
c
c     load geometry and symmetry info
c     -------------------------------
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0)
c
c     load the basis set and get info about it
c     ----------------------------------------
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $        call errquit('scf_init: no ao basis set', 0)
c
c     get nbf
c     -------
c
      if (.not.bas_numbf(basis,nbf))
     $     call errquit('bas_numbf failed',0)
c
c     initialize the API 
c     ------------------
c
      call int_init(rtdb, 1, basis)
      call int_acc_std()
c
c     normalize the basis
c     -------------------
c
      call int_normalize(rtdb,basis)
c
c     calculate the overlap integrals
c     -------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,g_s))
     $        call errquit('task_et ga_create failed', 0)
       call int_1e_ga(basis, basis, g_s, 'overlap', .false.)  
c
      if (util_print("et",print_debug)) then
       call ga_print(g_s)
      endif
c
c      terminate API
c      -------------
c       call int_terminate
c
c======================================================================= 
c     Read MO vectors from files.
c     Read mo vectors for state A (initial state) from file 'movecsa.mo' 
c     read mo vectors for state B (final state)   from file 'movecsb.mo' 
c=======================================================================
c
c     Allocate occupation vector. Used to calculate nocc(j),j=1,2 for alpha,
c     beta occupied orbitals.
c     ------------------------------------------------------------------------------
      if (.not.ma_push_get(mt_dbl,nbf,'occ',l_occ,k_occ))
     $         call errquit('ma_push_get failed for occ',k_occ)
c
c     Allocate evals, vector containing eigenvalues. Only used for
c     calling et_movecs_read. Eigenvalues themselves are not needed.
c     ----------------------------------------------------------
      if (.not.ma_push_get(mt_dbl,nbf,'evals',l_eval,k_eval))
     $        call errquit('ma_push_get failed for evals',k_eval)
c
c     allocate arrays for MO coefficients, for state A
c     -------------------------------------------------------
      if (.not.ga_create(mt_dbl,nbf,nbf,'A alpha',0,0,g_avecsa))
     $         call errquit('ga_create failed for A alpha',g_avecsa)

      if (.not.ga_create(mt_dbl,nbf,nbf,'A beta',0,0,g_avecsb))
     $         call errquit('ga_create failed for A beta',g_avecsb)

c     read occ, evals, mo coefficients for alpha orbitals from movecsa.mo
c     ------------------------------------------------------------------
c     if(.not.et_movecs_read('movecsa.mo',1,dbl_mb(k_occ),dbl_mb(k_eval)
c    $   ,g_avecsa,dum))
c    $ call errquit('et_movecs_read failed for alpha mos',0)
c
c     get nocc(1): # of occupied alpha mo's
c     -------------------------------------
      nocc(1)=0
      do i = 1,nbf
       nocc(1) = nocc(1) + dbl_mb(k_occ+i-1)
      enddo
c
      if (util_print("et",print_high)) then
      write(6,*) nocc(1), " occupied alpha orbitals"
      write(6,*)
      call ga_print(g_avecsa)
      write(6,*)
      endif
c
c     read occ, evals, mo coefficients for beta orbitals from movecsa.mo
c     ------------------------------------------------------------------
c     if(.not.et_movecs_read('movecsa.mo',2,dbl_mb(k_occ),dbl_mb(k_eval)
c    $     ,g_avecsb, dum))
c    $ call errquit('et_movecs_read failed for beta mos',0)
c
c     get nocc(2): # of occupied alpha mo's
c     -------------------------------------
      nocc(2)=0
      do i = 1,nbf
       nocc(2) = nocc(2) + dbl_mb(k_occ+i-1)
      enddo
c
      if (util_print("et",print_high)) then
      write(6,*) nocc(2), " occupied  beta orbitals"
      write(6,*)
      call ga_print(g_avecsb)
      write(6,*)
      endif
c
c     allocate arrays for MO coefficients, for state B
c     -------------------------------------------------------
      if (.not.ga_create(mt_dbl,nbf,nbf,'B alpha',0,0,g_bvecsa))
     $         call errquit('ga_create failed for B alpha',g_bvecsa)

      if (.not.ga_create(mt_dbl,nbf,nbf,'B beta',0,0,g_bvecsb))
     $         call errquit('ga_create failed for B beta',g_bvecsb)

c
c     read occ, evals, mo coefficients from movecsb.mo
c     ------------------------------------------------
c     if(.not.et_movecs_read('movecsb.mo',1,dbl_mb(k_occ),dbl_mb(k_eval)
c    $   ,g_bvecsa, dum))
c    $   call errquit('et_movecs_read failed for alpha mos',0)
c
c     if(.not.et_movecs_read('movecsb.mo',2,dbl_mb(k_occ),dbl_mb(k_eval)
c    $  ,g_bvecsb, dum))
c    $   call errquit('et_movecs_read failed for beta mos',0)
c
      if (util_print("et",print_high)) then
       call ga_print(g_bvecsa)
       call ga_print(g_bvecsb)
      endif
c
c==================================================
c     calculate D=B'*S*A for alpha and beta electrons
c==================================================
c
      g_ahandle(1)=g_avecsa
      g_ahandle(2)=g_avecsb
      g_bhandle(1)=g_bvecsa
      g_bhandle(2)=g_bvecsb
c
c
c     detprod is the product of det(U)*det(V) for alpha and beta electrons,
c     and is used later.
c
      detprod=1.0d0
c
c  for generalized density matrices:
c  Pa is for alpha electrons, Pb is for beta electrons.
c
      pname(1)='pa'
      pname(2)='pb'
c
c j is the alpha/beta loop:
c
      do j = 1,2 
c
      if (util_print("et",print_high)) then
       write(6,*) 
       write(6,*) '---------------'
       if(j.eq.1) write(6,*) 'ALPHA ELECTRONS'
       if(j.eq.2) write(6,*) ' BETA ELECTRONS'
       write(6,*) '---------------'
      else if (util_print("et",print_debug)) then 
       write(6,*) 
       write(6,*) '---------------'
       if(j.eq.1) write(6,*) 'ALPHA ELECTRONS'
       if(j.eq.2) write(6,*) ' BETA ELECTRONS'
       write(6,*) '---------------'
      endif
c
c  Perform matrix multiplication to get D = B'*S*A
c  set size of D to be only nocc*nocc
c  -----------------------------------------------
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'d',0,0,g_d))
     $        call errquit('task_et ga_create failed', g_d)
c
c D only contains elements for occupied orbitals
c
c g_tmp: handle for temporary array B'*S
      g_tmp = ga_create_atom_blocked(geom, basis, 'g_tmp')    
c
      call ga_dgemm('T','N',nocc(j),nbf,nbf,1.d0,g_bhandle(j)
     &               ,g_s,0.d0,g_tmp)
c
      call ga_dgemm('N','N',nocc(j),nocc(j),nbf,1.d0,g_tmp
     &               ,g_ahandle(j),0.d0,g_d)
c
c  tmp array is no longer needed.
      if (.not.ga_destroy(g_tmp))
     &  call errquit('could not destroy tmp array',g_tmp)
c
c
      if (util_print("et",print_debug)) then
       call ga_print(g_d)
      endif
c
c Perform SVD to get U, V', and singluar values
c The largest possible size of U and V' is nocc*nocc
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'u',0,0,g_u))
     $        call errquit('task_et ga_create failed', g_u)
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'vt',0,0,g_vt))
     $        call errquit('task_et ga_create failed', g_vt)
c
      if (.not.ma_push_get(mt_dbl,nocc(j),'values',l_vals,k_vals))
     $        call errquit('ma_push_get failed for values ',k_vals)
c
      call ga_svd_seq(g_d,g_u,g_vt,dbl_mb(k_vals))
c
c  D array is no longer needed.
      if (.not.ga_destroy(g_d))
     &  call errquit('could not destroy D array',g_d)
c
      if (util_print("et",print_debug)) then
      call ga_print(g_u)
      call ga_print(g_vt)
      endif
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'v',0,0,g_v))
     $        call errquit('task_et ga_create failed', g_v)
      call ga_transpose(g_vt,g_v)
c
      if (util_print("et",print_debug)) then
      call ga_print(g_v) 
      endif
c
c
c      ===========================
c      calculate det(U) and det(V)
c      ===========================
c
       call et_det('N',g_u,nocc(j),detu(j))
       call et_det('T',g_vt,nocc(j),detv(j))
       detprod = detprod*detu(j)*detv(j)
c
c calculate UU',U'U,VV',V'V to check unitarity.
c
c      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'temp',0,0,g_temp))
c     $        call errquit('task_et ga_create failed', g_temp)
cc
c      call ga_dgemm('N','T',nocc(j),nocc(j),nocc(j),1.d0,g_u,g_u
c     &         ,0.d0,g_temp)
c      call ga_print (g_temp)
cc
c      call ga_dgemm('T','N',nocc(j),nocc(j),nocc(j),1.d0,g_u,g_u
c     &         ,0.d0,g_temp)
c      call ga_print (g_temp)
cc
c      call ga_dgemm('T','N',nocc(j),nocc(j),nocc(j),1.d0,g_vt,g_vt
c     &         ,0.d0,g_temp)
c      call ga_print (g_temp)
cc
c      call ga_dgemm('N','T',nocc(j),nocc(j),nocc(j),1.d0,g_vt,g_vt
c     &         ,0.d0,g_temp)
c      call ga_print (g_temp)
c
c=============================================================
c                                               ^       ^
c perform matrix mult to get transformed A & B: A = AV, B = BU
c=============================================================
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'av',0,0,g_av))
     $        call errquit('task_et ga_create failed', g_av)
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'bu',0,0,g_bu))
     $        call errquit('task_et ga_create failed', g_bu)
c
      call ga_zero(g_av)
      call ga_dgemm('N','T',nbf,nocc(j),nocc(j),1.d0,g_ahandle(j),g_vt
     &         ,0.d0,g_av)
c
      call ga_zero(g_bu)
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_bhandle(j),g_u
     &         ,0.d0,g_bu)
c
c  vectors are no longer needed
      if (.not.ga_destroy(g_ahandle(j))) 
     &  call errquit('could not destroy mo vectors array',g_ahandle(j))
      if (.not.ga_destroy(g_bhandle(j))) 
     &  call errquit('could not destroy mo vectors array',g_bhandle(j))
c
c  U, V, VT are no longer needed
      if (.not.ga_destroy(g_u))
     &  call errquit('could not destroy U array',g_u)
      if (.not.ga_destroy(g_v))
     &  call errquit('could not destroy V array',g_v)
      if (.not.ga_destroy(g_vt))
     &  call errquit('could not destroy VT array',g_vt)
c
      if (util_print("et",print_high)) then
      write(6,*) 
      write(6,*) 'SVD transformed MO vectors for state A'
      call ga_print (g_av)
      write(6,*) 
      write(6,*) 'SVD transformed MO vectors for state B'
      call ga_print (g_bu)
      endif
c
c
C                                   N
C                                -------
C     Calculates   T(i) =     PRODUCT  d(j) where N = nalpha or nbeta
C                                -------
C                                j.NE.i  
c
c    also calculates prod, the product of diagonal elements in T(i), 
c    including both alpha and beta contributions.
c    --------------------------------------------------------------
c    
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'t',0,0,g_t))
     $        call errquit('task_et ga_create failed', g_t)
      call ga_zero(g_t)
c
      prod(j)=1.0d0
c
      if(nocc(j).gt.1) then
c
        do i = 1,nocc(j)
          val = 1.0d0
            do k = 1,nocc(j)
              if(k.ne.i) val = val * dbl_mb(k_vals+i-1)
            enddo
          call ga_put(g_t,i,i,i,i,val,1)
          prod(j) = prod(j)*val
        enddo
      else
         val = 1.0d0
         call ga_put(g_t,1,1,1,1,val,1)
         prod(j) = dbl_mb(k_vals)
c
      endif
c
c      write(6,*)
c      write(6,*)'prod',prod(j)
c      write(6,*)
c
      if (util_print("et",print_debug)) then
      call ga_print(g_t) 
      endif
c
c =====================================================
c                                            ^ ^
c  calculate generalized density matrix, P = ATB'
c  Pa is for alpha electrons, Pb is for beta electrons.
c =====================================================
c
      if(.not.ga_create(mt_dbl,nbf,nbf,pname(j),0,0,g_phandle(j)))
     $        call errquit('task_et ga_create failed', g_phandle(j))
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'temp',0,0,g_temp))
     $        call errquit('task_et ga_create failed', g_temp)
c
c
      call ga_zero(g_phandle(j))
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_av,g_t
     &         ,0.d0,g_temp)
c
      call ga_dgemm('N','T',nbf,nbf,nocc(j),1.0d0,g_temp,g_bu
     &         ,0.0d0,g_phandle(j))
c
c
c  transformed MO's are no longer needed
      if (.not.ga_destroy(g_av)) 
     &  call errquit('could not destroy transformed MO vectors',g_av)
      if (.not.ga_destroy(g_bu)) 
     &  call errquit('could not destroy transformed MO vectors',g_bu)
c
c  T and temp arrays are no longer needed
      if (.not.ga_destroy(g_t)) 
     &  call errquit('could not destroy T array',g_t)
      if (.not.ga_destroy(g_temp)) 
     &  call errquit('could not destroy temp array',g_temp)
c
c
c
c
      if (util_print("et",print_debug)) then
      call ga_print(g_phandle(j))
      endif
c
c
c end alpha/beta loop:
      enddo    
c
c
c  overlap integrals are not longer needed
      if (.not.ga_destroy(g_s)) 
     &  call errquit('could not destroy overlap array',g_s)
c
c
c
c
c
c
c why is Pa scaled by prod(b), Pb scaled by prod(a)?
c -------------------------------------------------
c   
       call ga_dscal(g_phandle(1),prod(2))
       call ga_dscal(g_phandle(2),prod(1))
c
      if (util_print("et",print_debug)) then
       write(6,*) 
       write(6,*) 'scaled generalized density matrices:' 
       write(6,*) 
       call ga_print(g_phandle(1))
       call ga_print(g_phandle(2))
      endif
c
c      =============
c      calculate SBA
c      =============
c
       Sba = prod(1)*prod(2)*detprod
c
c      ============================================
c      calculate the 1e contribution to Vba, omega1
c      ============================================
c
c     get hamiltonian integrals 
c     -------------------------------
      if(.not.ga_create(mt_dbl,nbf,nbf,'h1',0,0,g_h1))
     $        call errquit('task_et ga_create failed', g_h1)
c
      call ga_zero(g_h1)
      call int_1e_ga(basis, basis, g_h1, 'kinetic', .false.)  
      call int_1e_ga(basis, basis, g_h1, 'potential', .false.)  
c
      if (util_print("et",print_debug)) then
      call ga_print(g_h1)
      endif
c
c
       omega1 = 0.0d0
       do i = 1,nbf
       do j = 1,nbf
         call ga_get(g_phandle(1), i, i, j, j, pa_val, 1)
         call ga_get(g_phandle(2), i, i, j, j, pb_val, 1)
         call ga_get(g_h1, i, i, j, j, h1_val, 1)
         omega1 = omega1 + (pa_val+pb_val)*h1_val
c         write(6,*) i,j,pa_val+pb_val,h1_val
       enddo
       enddo
       omega1 = omega1*detprod
c
c  H1 array is no longer needed.
      if (.not.ga_destroy(g_h1))
     &  call errquit('could not destroy H1 array',g_h1)
c
c
c
c      ============================================
c      calculate the 2e contribution to Vba, omega2
c      ============================================
c
       call et_2e(g_phandle(1),g_phandle(2),omega2)
c
c  pa and pb arrays are no longer needed.
      if (.not.ga_destroy(g_phandle(1)))
     &  call errquit('could not destroy P array',g_phandle(1))
      if (.not.ga_destroy(g_phandle(2)))
     &  call errquit('could not destroy P array',g_phandle(2))
c
       omega2 = omega2 * (detprod/ ( prod(1)*prod(2) ) ) 
       omega2 = omega2 / 2.0d0
c
c      =============
c      calculate VBA
c      =============
c
c     Read converged scf or dft energy from movecs file
c     vectors, occupation vector, eigenvalues will be re-written,
c     but they are no longer needed at this point.
c     Filename and Haa or Hbb need to be correct.
c     -------------------------------------------------
c     if(.not.et_movecs_read('movecsa.mo',3,dbl_mb(k_occ),dbl_mb(k_eval)
c    $   ,g_bvecsa, Haa))
c    $   call errquit('et_movecs_read failed for Haa',0)
c
c     if(.not.et_movecs_read('movecsb.mo',3,dbl_mb(k_occ),dbl_mb(k_eval)
c    $   ,g_bvecsa, Hbb))
c    $   call errquit('et_movecs_read failed for Hbb',0)
c
c
c     get the nuclear repulsion energy. It needs to be subtracted
c     from the total energy to get the electronic energy.
c     -----------------------------------------------------------
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $  call errquit('unable to get nuclear rep energy',geom)
c
      Haa = Haa - enrep
      Hbb = Hbb - enrep
c
       Hba = omega1+omega2
       Vba =  Hba - Sba*(Haa+Hbb)*0.5d0
       Vba = Vba/(1.0d0-Sba*Sba)
       Vba = dabs (Vba)
c
c      =============
c      Print results
c      =============
c
 1111 FORMAT(1x,a34,1x,F15.10)
       write(6,*) 
       write(6,*) '====================================='
       write(6,*) 'Summary of Electron Transport Results'
       write(6,*) '====================================='
       write(6,*) 
       write(6,1111)  'Overlap, SBA                     ',Sba
       write(6,1111)  'Electronic energy of state A, HAA',Haa
       write(6,1111)  'Electronic energy of state B, HBB',Haa
       write(6,*) 
       write(6,*) 'Interaction Energy, HBA'
       write(6,*) '-----------------------'
       write (6,1111) '1e contribution, HBA(1)          ',omega1
       write (6,1111) '2e contribution, HBA(2)          ',omega2
       write (6,1111) 'Total HBA (Hartrees)             ',Hba
       write (6,*)
       write(6,*) 'Electron Transfer Matrix Element, VBA'
       write(6,*) '-------------------------------------'
       write (6,1111) '|VBA| (Hartrees)                 ',Vba
       write (6,1111) '      (cm-1)                     ',Vba*2.194746d5
       write (6,1111) '      (eV)                       ',Vba*27.2114d0
       write (6,1111) '      (kcal/mol)                 ',Vba*627.510d0
       write (6,1111) 
       write (6,1111) 'ln(|VBA|)                        ',log(Vba)
       write (6,*) 
       write (6,*) 
       write (6,*) 
       write (6,*) 
c
c
c
      et_calc=.true.
c
      call util_print_pop
      return
      end

