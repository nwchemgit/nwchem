      logical function et_calc(rtdb,Vba)
      implicit none
c
c     This function actually calculates Vba
*
* $Id: et_calc.F,v 1.1 2001-12-20 18:47:06 d3k958 Exp $
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
      logical movecs_read
      external movecs_read
      logical ga_copy
      external ga_copy
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical status
c
      character*2 pname(2)
      double precision Vba
      double precision prod(2)
      double precision Sba,pa_val,pb_val,omega1,h1_val
      double precision detu(2),detv(2),detprod,omega2
      double precision haa,hbb,hba,val
c
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
c     calculate the overlap integrals
c     -------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,g_s))
     $        call errquit('task_et ga_create failed', 0)
       call int_1e_ga(basis, basis, g_s, 'overlap', .false.)  
       call ga_print(g_s)
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
c     calling movecs_read. Eigenvalues themselves are not needed.
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
       if(.not.movecs_read('movecsa.mo',1,dbl_mb(k_occ),dbl_mb(k_eval)
     $      ,g_avecsa))
     $ call errquit('movecs_read failed for alpha mos',0)
c
c     get nocc(1): # of occupied alpha mo's
c     -------------------------------------
      nocc(1)=0
      do i = 1,nbf
       nocc(1) = nocc(1) + dbl_mb(k_occ+i-1)
      enddo
      write(6,*) 'new occ',nocc(1)
c
c     read occ, evals, mo coefficients for beta orbitals from movecsa.mo
c     ------------------------------------------------------------------
      if(.not.movecs_read('movecsa.mo',2,dbl_mb(k_occ),dbl_mb(k_eval)
     $     ,g_avecsb))
     $ call errquit('movecs_read failed for beta mos',0)
c
c     get nocc(2): # of occupied alpha mo's
c     -------------------------------------
      nocc(2)=0
      do i = 1,nbf
       nocc(2) = nocc(2) + dbl_mb(k_occ+i-1)
      enddo
      write(6,*) 'new occ',nocc(2)
c
      call ga_print(g_avecsa)
      call ga_print(g_avecsb)
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
      if(.not.movecs_read('movecsb.mo',1,dbl_mb(k_occ),dbl_mb(k_eval)
     $   ,g_bvecsa))
     $   call errquit('movecs_read failed for alpha mos',0)
c
      if(.not.movecs_read('movecsb.mo',2,dbl_mb(k_occ),dbl_mb(k_eval)
     $  ,g_bvecsb))
     $   call errquit('movecs_read failed for beta mos',0)
c
      call ga_print(g_bvecsa)
      call ga_print(g_bvecsb)
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
      g_tmp = ga_create_atom_blocked(geom, basis, 'g_tmp')    
c
c      write(6,*)
c      write(6,*)'# occupied alpha and beta orbitals ',nocc(1),nocc(2)
c      write(6,*)
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
      write(6,*)'nocc',nocc(j)
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
      if (.not.ma_push_get(mt_dbl,nocc(j),'values',l_vals,k_vals))
     $        call errquit('ma_push_get failed for values ',k_vals)
c
      call ga_svd_seq(g_d,g_u,g_vt,dbl_mb(k_vals))
c
      call ga_print(g_u)
      call ga_print(g_vt)
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'v',0,0,g_v))
     $        call errquit('task_et ga_create failed', g_v)
      call ga_transpose(g_vt,g_v)
      call ga_print(g_v) 
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
      call ga_print (g_av)
c
      call ga_zero(g_bu)
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_bhandle(j),g_u
     &         ,0.d0,g_bu)
      call ga_print (g_bu)
c
c  check that D is diagonal after SVD
c
c      if(.not.ga_create(mt_dbl,nocc(j),nbf,'temp',0,0,g_temp))
c     $        call errquit('task_et ga_create failed', g_temp)
c      call ga_dgemm('T','N',nocc(j),nbf,nbf,1.d0,g_bu
c     &               ,g_s,0.d0,g_temp)
c      call ga_dgemm('N','N',nocc(j),nocc(j),nbf,1.d0,g_temp
c     &               ,g_av,0.d0,g_d)
c      call ga_print(g_d)
c      DO i = 1,nocc(j)
c      DO k = 1,nocc(j)
c         call ga_get(g_d, i, i, k, k, val, 1)
c         write (6,*) i,k,val
c      ENDDO
c      ENDDO
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
      call ga_print(g_h1)
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
c      ============================================
c      calculate the 2e contribution to Vba, omega2
c      ============================================
c
       call vab2e(g_phandle(1),g_phandle(2),omega2)
c
       omega2 = omega2 * (detprod/ ( prod(1)*prod(2) ) ) 
       omega2 = omega2 / 2.0d0
c
c      =============
c      calculate VBA
c      =============
c      Haa and Hbb will have to be read in somehow. Temporarily hard-wired.
c
       Haa = -4.811062033425d0 -  0.529177249000d0
       Hbb = -4.811062033425d0 -  0.529177249000d0
       Hba = omega1+omega2
       Vba =  Hba - Sba*(Haa+Hbb)*0.5d0
       Vba = Vba/(1.0d0-Sba*Sba)
       Vba = dabs (Vba)
c
c      =============
c      Print results
c      =============
c
 1111 FORMAT(1x,a15,1x,F15.10,1x,a10)
       write(6,*) 
       write(6,*) '-------------------------------------'
       write(6,*) 'Summary of Electron Transport Results'
       write(6,*) '-------------------------------------'
       write(6,*) 
       write(6,1111)  'overlap, Sba   ',Sba
       write(6,*) 
       write (6,1111) '1e contribution',omega1
       write (6,1111) '2e contribution',omega2
       write (6,1111) 'Hba            ',Hba
       write (6,*)
       write (6,1111) 'Vba            ',Vba*2.194746d5,'CM-1'
       write (6,1111) '               ',Vba*27.2114d0,'EV'
       write (6,1111) '               ',Vba*627.510d0,'KCAL/MOL'
       write (6,*) 
c
c
c
      et_calc=.true.
c
      return
      end

c==================================================================
      subroutine et_det(op,handle,n,detval)
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
c returns the determinant of the square matrix with file handle 'handle'.
c n is the dimension of the array, detval is the returned value of the
c determinant. If op='N', a=a. If op='T', a=a'.
c
      character*1 op
      integer handle
      integer n,i,j,ipvt(n),info
      double precision detval,a(n,n),work(n)
c
c
c
      do i = 1,n
      do j = 1,n
       if(op.eq.'N') call ga_get(handle, i, i, j, j, a(i,j), 1)  
       if(op.eq.'T') call ga_get(handle, j, j, i, i, a(i,j), 1)  
      enddo
      enddo
c
      call dgefa(a,n,n,ipvt,info)
      call dgedi (a,n,n,ipvt,detval,work,10)
c
      return
      end

c=======================================================================
       subroutine vab2e(g_pa,g_pb,omega2)
c
c      calculates the two-electron contribution to Vba
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
       integer g_pa,g_pb,g_pt,i,j,k,l
       double precision int2e1,int2e2,omega2
       double precision pa1,pa2,pb1,pb2,pt1,pt2
       double precision scr(mem_2e4c)
c
c      create array pt = pa + pb
c 
      if(.not.ga_create(mt_dbl,nbf,nbf,'pt',0,0,g_pt))
     $        call errquit('task_et ga_create failed', g_pt)
      call ga_dadd (1.0d0,g_pa,1.0d0,g_pb,g_pt)
c
      omega2 = 0.0d0
c
       DO i = 1,nbf
       DO j = 1,nbf
       DO k = 1,nbf
       DO l = 1,nbf
         call ga_get(g_pt, i, i, j, j, pt1, 1)
         call ga_get(g_pt, k, k, l, l, pt2, 1)
         call ga_get(g_pa, i, i, j, j, pa1, 1)
         call ga_get(g_pa, k, k, l, l, pa2, 1)
         call ga_get(g_pb, i, i, j, j, pb1, 1)
         call ga_get(g_pb, k, k, l, l, pb2, 1)
         call int_2e4c(basis,i,j,basis,k,l,mem_2e4c,scr,1,int2e1)
         call int_2e4c(basis,i,l,basis,k,j,mem_2e4c,scr,1,int2e2)
        omega2 = omega2 + ( pt1*pt2 )*int2e1
        omega2 = omega2 + ( -pa1*pa2-pb1*pb2 )*int2e2
       ENDDO
       ENDDO
       ENDDO
       ENDDO
c
       return
       end
c-----------------------------------------------------
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end




 
