      logical function task_et(rtdb)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "cscf.fh" 
      integer rtdb
c
c  11/26/01 DMAS: this was a copy of task_energy.F 
c  changed 'energy' tasks to 'et', changed 'energy' variable to 'Vab',
c  removed ONIOM stuff, and theories other than SCF and DFT.
c  Also removed 'rtdb_clean' subroutine, and added task_et.o to GNUMakefile.
c
c      logical status, ignore
      logical status
c      character*32 theory
      double precision cpu, wall, Vab
      logical task_et_doit
      external task_et_doit
c
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c
c      ignore = rtdb_delete(rtdb, 'task:et')
c      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
c     $     call errquit('task_et: failed to invalidate status',0)
cc
c      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
c     $     call errquit('task:et: theory not specified',0)
c
      call util_flush(6)
      status = task_et_doit(rtdb,Vab)
c      status = task_et_doit(rtdb,theory,Vab)
      call util_flush(6)
      if (.not. rtdb_put(rtdb, 'task:et', mt_dbl, 1, Vab))
     $     call errquit('task_et: failed storing Vab',0)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_et: failed to store status',0)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_et: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_et: failed storing walltime',0)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call geom_print_rtdb_ecce(rtdb)
         call ecce_print1('Vab', mt_dbl, Vab, 1)
         call ecce_print_module_exit('task et', 'ok')
      else
         call ecce_print_module_exit('task et', 'failed')
      endif
c         
      task_et = status
c
      end
c===========================================================================
      logical function task_et_doit(rtdb,Vab)
c      logical function task_et_doit(rtdb,theory,Vab)
      implicit none
c
c     This function actually calculates Vab
*
* $Id: task_et.F,v 1.5 2001-12-14 23:03:03 d3k958 Exp $
c
*
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
c
c     Generic NWChem interface to compute Vab. SCF of DFT determinants
c     can be used.
c
c     RTDB output parameters
c     ----------------------
c     task:status (logical)- T/F for success/failure
c     if (status) then
c     .  task:et (real)   - Vab (electron transfer matrix element)
c     .  task:cputime (real)  - cpu time to execute the task
c     .  task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
      integer rtdb
      integer i, j, k
      integer l_evals,k_evals
      integer g_vecsa,g_vecsb
      integer g_tmp,g_s,g_ahandle(2),g_bhandle(2),g_d,g_temp
      integer nocc(2)
      integer g_u,g_vt,g_phandle(2)
      integer g_av,g_bu,g_t,g_pa,g_pb,g_h1
      integer g_avecsa,g_avecsb,g_bvecsa,g_bvecsb
c      integer k_val,l_val,g_av,g_bu
c
      logical movecs_read
      external movecs_read
      logical ga_copy
      external ga_copy
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical status
c
      character*2 pname(2)
      double precision scr(mem_1e),h1(nbf,nbf)
c      double precision c_alpha(nbf,nbf),c_beta(nbf,nbf)
      double precision occ_alpha(nbf),occ_beta(nbf)
      double precision evals_alpha(nbf),evals_beta(nbf)
      double precision Vab,values(nbf*nbf)
      double precision t(nbf,nbf),prod(2)
      double precision sba
c
      task_et_doit = .false.
c
c=======================================================
c     use the integral API to get the 1e integrals
c     S and H1:
c=======================================================
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
c     calculate the overlap integrals and store in s(nbf,nbf)
c     ------------------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,g_s))
     $        call errquit('task_et ga_create failed', 0)
       call int_1e_ga(basis, basis, g_s, 'overlap', .false.)  
       call ga_print(g_s)
c
c     calculate hamiltonian integrals and store in h1(nbf,nbf)
c     ------------------------------------------------------
      if(.not.ga_create(mt_dbl,nbf,nbf,'h1',0,0,g_h1))
     $        call errquit('task_et ga_create failed', g_h1)
c
       do i = 1,nbf
       do j = 1,nbf
         call int_1eh1(basis,i,basis,j,mem_1e,scr,isz_1e,h1(i,j))
         call ga_put(g_h1, i, i, j, j, h1(i,j), 1)
       enddo
c         write(6,1112) i,(h1(i,j),j=1,nbf)
       enddo
       call ga_print(g_h1)
c
c      terminate API
c      -------------
c       call int_terminate
c
c 1112 FORMAT (i2,6x,4(F8.5,5x))
c
c======================================================================= 
c     Read MO vectors from files.
c     Read mo vectors for state A (initial state) from file 'movecsa.mo' 
c     read mo vectors for state B (final state)   from file 'movecsb.mo' 
c=======================================================================
c
c     allocate occ_alpha & beta, containing 0 or 1 for virtual and occupied orbitals
c     ------------------------------------------------------------------------------
c
      if (.not.ma_push_get(mt_dbl,nbf,'occ_alpha',l_occ,k_occ))
     $         call errquit('ma_push_get failed for occ_alpha',k_occ)
c
      if (.not.ma_push_get(mt_dbl,nbf,'occ_beta',l_occ,k_occ))
     $         call errquit('ma_push_get failed for occ_beta',k_occ)
c
c
c     allocate evals_alpha & beta, vector containing eigenvalues
c     ----------------------------------------------------------
      if (.not.ma_push_get(mt_dbl,nbf,'evals_alpha',l_evals,k_evals))
     $        call errquit('ma_push_get failed for evals_alpha',k_evals)
c
      if (.not.ma_push_get(mt_dbl,nbf,'evals_beta',l_evals,k_evals))
     $        call errquit('ma_push_get failed for evals_beta',k_evals)
c
c     allocate arrays for MO coefficients, c_alpha and c_beta
c     -------------------------------------------------------
      if (.not.ga_create(mt_dbl,nbf,nbf,'c_alpha',0,0,g_vecsa))
     $         call errquit('ga_create failed for c_alpha',g_vecsa)

      if (.not.ga_create(mt_dbl,nbf,nbf,'c_beta',0,0,g_vecsb))
     $         call errquit('ga_create failed for c_beta',g_vecsb)

c     read occ, evals, mo coefficients from movecsa.mo
c     ------------------------------------------------
      if(.not.movecs_read('movecsa.mo',1,occ_alpha,evals_alpha,g_vecsa))
     $   call errquit('movecs_read failed for alpha mos',0)
c
      if(.not.movecs_read('movecsa.mo',2,occ_beta,evals_beta,g_vecsb))
     $         call errquit('movecs_read failed for beta mos',0)
c
c     copy mo coefficients to arrays for state A
c     ------------------------------------------
c
      if (.not.ga_duplicate(g_vecsa,g_avecsa,'A alpha'))
     $         call errquit('ga_duplicate failed',g_avecsa)
      if (.not.ga_copy(g_vecsa,g_avecsa))
     $         call errquit('ga_copy failed',g_avecsa)
c
      if (.not.ga_duplicate(g_vecsb,g_avecsb,'A beta'))
     $         call errquit('ga_duplicate failed',g_avecsb)
      if (.not.ga_copy(g_vecsb,g_avecsb))
     $         call errquit('ga_copy failed',g_avecsb)
c
      call ga_print(g_avecsa)
      call ga_print(g_avecsb)
c
c      write(6,*) 
c      write(6,*) 'occupation numbers and eigenvectors for state A'
c      write(6,*) 'occ_alpha',(occ_alpha(i),i=1,nbf)
c      write(6,*) 'occ_beta ',(occ_beta(i),i=1,nbf)
c      write(6,*) 'evals alpha',(evals_alpha(i),i=1,nbf)
c      write(6,*) 'evals beta ',(evals_beta(i),i=1,nbf)
c      write(6,*) 
c
c     read occ, evals, mo coefficients from movecsb.mo
c     ------------------------------------------------
      if(.not.movecs_read('movecsb.mo',1,occ_alpha,evals_alpha,g_vecsa))
     $   call errquit('movecs_read failed for alpha mos',0)
c
      if(.not.movecs_read('movecsb.mo',2,occ_beta,evals_beta,g_vecsb))
     $         call errquit('movecs_read failed for beta mos',0)
c
c     copy mo coefficients to arrays for state B
c     ------------------------------------------
c
      if (.not.ga_duplicate(g_vecsa,g_bvecsa,'B alpha'))
     $          call errquit('ga_duplicate failed',g_bvecsa) 
      if (.not.ga_copy(g_vecsa,g_bvecsa))
     $         call errquit('ga_copy failed',g_bvecsa)

c
      if (.not.ga_duplicate(g_vecsb,g_bvecsb,'B beta'))
     $          call errquit('ga_duplicate failed',g_bvecsb)
      if (.not.ga_copy(g_vecsb,g_bvecsb))
     $         call errquit('ga_copy failed',g_bvecsb)
c
      call ga_print(g_bvecsa)
      call ga_print(g_bvecsb)
c
c      write(6,*) 
c      write(6,*) 'occupation numbers and eigenvectors for state B'
c      write(6,*) 'occ_alpha',(occ_alpha(i),i=1,nbf)
c      write(6,*) 'occ_beta ',(occ_beta(i),i=1,nbf)
c      write(6,*) 'evals alpha',(evals_alpha(i),i=1,nbf)
c      write(6,*) 'evals beta ',(evals_beta(i),i=1,nbf)
c      write(6,*) 
c
c
c==================================================
c     calculate D=B'*S*A for alpha and beta electrons
c==================================================
c
c g_tmp: handle for temporary array B'*S
c
      g_ahandle(1)=g_avecsa
      g_ahandle(2)=g_avecsb
      g_bhandle(1)=g_bvecsa
      g_bhandle(2)=g_bvecsb
c
      g_tmp = ga_create_atom_blocked(geom, basis, 'g_temp')    
c
c get nocc(1 and 2): # of occupied alpha and beta mo's
c
      nocc(1)=0
      do i = 1,nbf
       nocc(1) = nocc(1) + occ_alpha(i)
      enddo
c
      nocc(2)=0
      do i = 1,nbf
       nocc(2) = nocc(2) + occ_beta(i)
      enddo
c
      write(6,*)
      write(6,*)'# occupied alpha and beta orbitals ',nocc(1),nocc(2)
      write(6,*)
c
c  for generalized density matrices:
c  Pa is for alpha electrons, Pb is for beta electrons.
c
      pname(1)='pa'
      pname(2)='pb'
      g_phandle(1)=g_pa
      g_phandle(2)=g_pb
c
c j is the alpha/beta loop:
c
      do j = 1,2 
c
      write(6,*) 
      write(6,*) '---------------'
      if(j.eq.1) write(6,*) 'ALPHA ELECTRONS'
      if(j.eq.2) write(6,*) ' BETA ELECTRONS'
      write(6,*) '---------------'
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
      call ga_dgemm('T','N',nocc(j),nbf,nbf,1.d0,g_bhandle(j)
     &               ,g_s,0.d0,g_tmp)
      call ga_dgemm('N','N',nocc(j),nocc(j),nbf,1.d0,g_tmp
     &               ,g_ahandle(j),0.d0,g_d)
      call ga_print(g_d)
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
      call ga_svd_seq(g_d,g_u,g_vt,values)
c
      call ga_print(g_u)
      call ga_print(g_vt)
c
      write(6,*)
      write(6,*)'diagonal values',(values(i),i=1,nocc(j))
      write(6,*)
c
c calculate UU',U'U,VV',V'V to check unitarity.
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'temp',0,0,g_temp))
     $        call errquit('task_et ga_create failed', g_temp)
c
      call ga_dgemm('N','T',nocc(j),nocc(j),nocc(j),1.d0,g_u,g_u
     &         ,0.d0,g_temp)
      call ga_print (g_temp)
c
      call ga_dgemm('T','N',nocc(j),nocc(j),nocc(j),1.d0,g_u,g_u
     &         ,0.d0,g_temp)
      call ga_print (g_temp)
c
      call ga_dgemm('T','N',nocc(j),nocc(j),nocc(j),1.d0,g_vt,g_vt
     &         ,0.d0,g_temp)
      call ga_print (g_temp)
c
      call ga_dgemm('N','T',nocc(j),nocc(j),nocc(j),1.d0,g_vt,g_vt
     &         ,0.d0,g_temp)
      call ga_print (g_temp)
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
      call ga_print (g_av)
c
      call ga_zero(g_bu)
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_bhandle(j),g_u
     &         ,0.d0,g_bu)
      call ga_print (g_bu)
c
c
C                                   N
C                                -------
C     Calculates   T(i) =     PRODUCT  d(j) where N = nalpha or nbeta
C                                -------
C                                j.NE.i  
c
c
c    
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'t',0,0,g_t))
     $        call errquit('task_et ga_create failed', g_t)
c
      prod(j)=1.0d0
c
      if(nocc(j).gt.1) then
      do i = 1,nocc(j)
       t(i,i) = 1.0d0
        do k = 1,nocc(j)
         if(k.ne.i) t(i,i) = t(i,i) * values(k)
        enddo
        prod(j) = prod(j)*t(i,i)
      enddo
      else
       t(1,1) = 1.0d0
       prod(j) = values(1)
      endif
c
      write(6,*)
      write(6,*)'prod',prod(j)
      write(6,*)
c
      call ga_put(g_t, 1, nocc(j), 1, nocc(j),t, nocc(j)*nocc(j))
      call ga_print(g_t) 
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
      call ga_print(g_phandle(j))
c
c
c
c end alpha/beta loop:
      enddo    
c
c
c
c why is Pa scaled by prod(b), Pb scaled by prod(a)?
c -------------------------------------------------
c   
       write(6,*) 
       write(6,*) 'scaled generalized density matrices:' 
       write(6,*) 
       call ga_dscal(g_phandle(1),prod(2))
       call ga_print(g_phandle(1))
       call ga_dscal(g_phandle(2),prod(1))
       call ga_print(g_phandle(2))
c
c =============
c calcualte SBA
c =============
c
       sba = prod(1)*prod(2)
       write(6,*) 
       write(6,*) 'SBA',sba
       write(6,*) 
c
c      calculate Vab
c      -------------
c
       Vab=1.0d0
       write(6,*)
       write(6,*)'Vab = ',Vab
c
      status=.true.
      task_et_doit = status
c
      end


