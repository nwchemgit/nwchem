c $Id: et_calc.F,v 1.11 2002-03-27 22:16:51 d3k958 Exp $
      logical function et_calc(rtdb,Vba)
      implicit none
c
c     This function actually calculates Vba
*
* $Id: et_calc.F,v 1.11 2002-03-27 22:16:51 d3k958 Exp $
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
c
c
      integer rtdb, me
      integer i, j, k
      integer l_vals,k_vals
      integer g_tmp,g_s,g_ahandle(2),g_bhandle(2),g_d,g_temp
      integer nocc(2)
      integer g_u,g_vt,g_v,g_phandle(2)
      integer g_av,g_bu,g_t,g_h1
      integer g_avecsa,g_avecsb,g_bvecsa,g_bvecsb
      integer mem2, max2e, lena, lenb
      integer method_2e
c
      logical et_movecs_read
      external et_movecs_read
      logical ga_copy
      external ga_copy
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      character*255 a_file, b_file
      character*2 pname(2)
      double precision Vba
      double precision prod(2)
      double precision Sba,pa_val,pb_val,omega1,h1_val
      double precision detu(2),detv(2),detprod,omega2
      double precision Haa,Hbb,Hba,val, dum
c
      me=ga_nodeid()
c
      if (me.eq.0) then
      Call Util_Print_Centered(6,'Electron Transfer Calculation'
     &     ,40,.true.)
      write(6,*)
      endif
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
     $     call errquit('et_calc: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $        call errquit('et_calc: no ao basis set', 0)
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
      call int_normalize(rtdb,basis)
c
c     calculate the overlap integrals
c     -------------------------------
      call ga_sync()
      if(.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,g_s))
     $        call errquit('task_et ga_create failed', 0)
      call ga_zero(g_s)
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)  
c
      if (util_print("et",print_debug)) then
       call ga_print(g_s)
      endif
c
c======================================================================= 
c     Read MO vectors from files.
c     First, get the A and B filenames from the rtdb
c=======================================================================
c
      if (.not. rtdb_cget(rtdb, 'a_vecs', 1, a_file))
     $     call errquit('et_calc: couldnt get a_vecs filename',0)
c
      if (.not. rtdb_cget(rtdb, 'b_vecs', 1, b_file))
     $     call errquit('et_calc: couldnt get a_vecs filename',0)
c
      lena = inp_strlen(a_file)
      lenb = inp_strlen(b_file)
c
      if (me.eq.0) then
      write(6,*)'MO vectors for reactants: ',a_file(1:lena)
      write(6,*)'MO vectors for products : ',b_file(1:lenb)
      write(6,*)
      endif
c
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

c     read occ, evals, mo coefficients for alpha orbitals from a_file
c     ------------------------------------------------------------------
      if(.not.et_movecs_read(a_file(1:lena),1,dbl_mb(k_occ)
     $   ,dbl_mb(k_eval),g_avecsa,dum))
     $ call errquit('et_movecs_read failed for alpha mos',0)
c
c     get nocc(1): # of occupied alpha mo's
c     -------------------------------------
      nocc(1)=0
      do i = 1,nbf
       nocc(1) = nocc(1) + dbl_mb(k_occ+i-1)
      enddo
c
      if (util_print("et",print_high)) then
      if(me.eq.0) then
        write(6,*) nocc(1), " occupied alpha orbitals"
        write(6,*)
      endif
        call ga_print(g_avecsa)
      endif
c
c     read occ, evals, mo coefficients for beta orbitals from a_file
c     ------------------------------------------------------------------
      if(.not.et_movecs_read(a_file(1:lena),2,dbl_mb(k_occ)
     $     ,dbl_mb(k_eval),g_avecsb, dum))
     $ call errquit('et_movecs_read failed for beta mos',0)
c
c     get nocc(2): # of occupied alpha mo's
c     -------------------------------------
      nocc(2)=0
      do i = 1,nbf
       nocc(2) = nocc(2) + dbl_mb(k_occ+i-1)
      enddo
c
      if (util_print("et",print_high)) then
      if(me.eq.0) then
        write(6,*) nocc(2), " occupied  beta orbitals"
        write(6,*)
      endif
        call ga_print(g_avecsb)
      endif
c
c     allocate arrays for MO coefficients, for state B
c     -------------------------------------------------------
      if (.not.ga_create(mt_dbl,nbf,nbf,'B alpha',0,0,g_bvecsa))
     $         call errquit('ga_create failed for B alpha',g_bvecsa)

      if (.not.ga_create(mt_dbl,nbf,nbf,'B beta',0,0,g_bvecsb))
     $         call errquit('ga_create failed for B beta',g_bvecsb)

c
c     read occ, evals, mo coefficients from b_file
c     ------------------------------------------------
      if(.not.et_movecs_read(b_file(1:lenb),1,dbl_mb(k_occ)
     $   ,dbl_mb(k_eval) ,g_bvecsa, dum))
     $   call errquit('et_movecs_read failed for alpha mos',0)
c
      if(.not.et_movecs_read(b_file(1:lenb),2,dbl_mb(k_occ)
     $  ,dbl_mb(k_eval) ,g_bvecsb, dum))
     $   call errquit('et_movecs_read failed for beta mos',0)
c
      if (util_print("et",print_high)) then
       call ga_print(g_bvecsa)
       call ga_print(g_bvecsb)
      endif
c
c     Read converged scf or dft energy from movecs file
c     vectors, occupation vector, eigenvalues will be re-written,
c     but they are no longer needed at this point.
c     ==========================================================
c
      if(.not.et_movecs_read(a_file(1:lena),3
     $   ,dbl_mb(k_occ),dbl_mb(k_eval),g_bvecsa, Haa)) 
     $   call errquit('et_movecs_read failed for Haa',0)
c
      if(.not.et_movecs_read(b_file(1:lenb),3
     $   ,dbl_mb(k_occ),dbl_mb(k_eval),g_bvecsa, Hbb)) 
     $   call errquit('et_movecs_read failed for Hbb',0)
c
      if (me.eq.0) then
      write(6,9995) 'Electronic energy of reactants   ','H(rr)',Haa
      write(6,9995) 'Electronic energy of products    ','H(pp)',Hbb
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
      if(me.eq.0) then
       write(6,*) 
       write(6,*) '---------------'
       if(j.eq.1) write(6,*) 'ALPHA ELECTRONS'
       if(j.eq.2) write(6,*) ' BETA ELECTRONS'
       write(6,*) '---------------'
      endif
      endif
c
c  Perform matrix multiplication to get D = B'*S*A
c  set size of D to be only nocc*nocc
c  -----------------------------------------------
c
      call ga_sync()
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'d',0,0,g_d))
     $        call errquit('task_et ga_create failed', g_d)
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'temp',0,0,g_tmp))
     $        call errquit('task_et ga_create failed', g_d)
c
      call ga_dgemm('T','N',nocc(j),nbf,nbf,1.d0,g_bhandle(j)
     &               ,g_s,0.d0,g_tmp)
c
      call ga_dgemm('N','N',nocc(j),nocc(j),nbf,1.d0,g_tmp
     &               ,g_ahandle(j),0.d0,g_d)
c
      if (.not.ga_destroy(g_tmp))
     &  call errquit('could not destroy tmp array',g_tmp)
c
      call ga_sync()
      if (util_print("et",print_debug)) then
       call ga_print(g_d)
      endif
c
c Perform SVD to get U, V', and singluar values
c Size of U and V' is nocc*nocc
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'u',0,0,g_u))
     $        call errquit('task_et ga_create failed', g_u)
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'vt',0,0,g_vt))
     $        call errquit('task_et ga_create failed', g_vt)
c
      call ga_zero(g_u)
      call ga_zero(g_vt)
      if (.not.ma_push_get(mt_dbl,nocc(j),'values',l_vals,k_vals))
     $        call errquit('ma_push_get failed for values ',k_vals)
c
      call ga_svd_seq(g_d,g_u,g_vt,dbl_mb(k_vals))
c       call et_svd(nocc(j),g_d,g_u,g_vt,dbl_mb(k_vals))
c
c prod is the product of singular values.
c
      prod(j) = 1.0d0
      do k = 1,nocc(j)
         prod(j) = prod(j)*dbl_mb(k_vals+k-1)
      enddo
c      write(6,*)'prod',j,prod(j)
c      write(6,*)'singular values',(dbl_mb(k_vals+i-1),i=1,nocc(j))
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
c       write(6,*) j,'detu,detv',detu(j),detv(j)
c
c=============================================================
c                                               ^       ^
c perform matrix mult to get transformed A & B: A = AV, B = BU
c=============================================================
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'av',0,0,g_av))
     $        call errquit('task_et ga_create failed', g_av)
c
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
      if(me.eq.0) then
        write(6,*) 
        write(6,*) 'SVD transformed MO vectors for state A'
      endif
        call ga_print (g_av)
      if(me.eq.0) then
        write(6,*) 
        write(6,*) 'SVD transformed MO vectors for state B'
      endif
        call ga_print (g_bu)
      endif
c
c
C                                   N
C                                -------
C     Calculates   T(i) =     PRODUCT  d(j) where N = nalpha or nbeta
C                                -------
C                                j.NE.i  
c    --------------------------------------------------------------
c    
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'t',0,0,g_t))
     $        call errquit('task_et ga_create failed', g_t)
      call ga_zero(g_t)
c
      if(nocc(j).gt.1) then
c
        do i = 1,nocc(j)
          val = 1.0d0
            do k = 1,nocc(j)
              if(k.ne.i) val = val * dbl_mb(k_vals+k-1)
            enddo
          call ga_put(g_t,i,i,i,i,val,1)
        enddo
      else
         val = 1.0d0
         call ga_put(g_t,1,1,1,1,val,1)
c
      endif
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
      call ga_zero(g_phandle(j))
c
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_av,g_t
     &         ,0.d0,g_temp)
c
      call ga_dgemm('N','T',nbf,nbf,nocc(j),1.0d0,g_temp,g_bu
     &         ,0.0d0,g_phandle(j))
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
      if (util_print("et",print_debug)) then
      call ga_print(g_phandle(j))
      endif
c
c
c end alpha/beta loop:
      enddo    
c
c  overlap integrals are not longer needed
      if (.not.ga_destroy(g_s)) 
     &  call errquit('could not destroy overlap array',g_s)
c
c
c why is Pa scaled by prod(b), Pb scaled by prod(a)?
c -------------------------------------------------
c   
       call ga_dscal(g_phandle(1),prod(2))
       call ga_dscal(g_phandle(2),prod(1))
c
      if (util_print("et",print_debug)) then
      if(me.eq.0) then
       write(6,*) 
       write(6,*) 'scaled generalized density matrices:' 
       write(6,*) 
      endif
       call ga_print(g_phandle(1))
       call ga_print(g_phandle(2))
      endif
c
c      =============
c      calculate SBA
c      =============
c
       Sba = prod(1)*prod(2)*detprod
c       write(6,*) 'Srp debug:',prod(1),prod(2),detprod
c
c
      if(me.eq.0) then
       write(6,*)
       write(6,9995) 'Reactants/Products overlap       ','S(rp)',Sba
       write(6,*)
      endif
c
      call util_flush(6) 
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
       omega1 = 0.0d0
       do i = 1,nbf
       do j = 1,nbf
         call ga_get(g_phandle(1), i, i, j, j, pa_val, 1)
         call ga_get(g_phandle(2), i, i, j, j, pb_val, 1)
         call ga_get(g_h1, i, i, j, j, h1_val, 1)
         omega1 = omega1 + (pa_val+pb_val)*h1_val
       enddo
       enddo
c
       omega1 = omega1*detprod
c
      if(me.eq.0) then
       write(6,*) 'Reactants/Products interaction energy:    '
       write(6,*) '-------------------------------------           '
       write(6,9995) 'One-electron contribution        ','H1(rp)',omega1
      endif
c
      call util_flush(6) 
c
c  H1 array is no longer needed.
      if (.not.ga_destroy(g_h1))
     &  call errquit('could not destroy H1 array',g_h1)
c
c      ============================================
c      calculate the 2e contribution to Vba, omega2
c      ============================================
c
       call int_mem_2e4c(mem2,max2e)
       mem2 = max(mem2,max2e)
c
      if (.not. rtdb_get(rtdb, 'et:tol2e', MT_DBL, 1, tol2e)) then
       tol2e    =  min(1.0d-7,dabs(Sba)*1.0d-7)
       tol2e    =  max(1.0d-12,tol2e)
      endif
c
      if (.not.rtdb_get(rtdb,'et:method_2e',MT_INT,1,method_2e)) then
       method_2e = 1  !default uses fock builder for omega2
      endif
c
      if (method_2e.EQ.1) then
      call et_fock(rtdb,geom,basis,
     &  g_phandle(1),g_phandle(2),nbf,tol2e,omega2)
      endif
c
      if (method_2e.EQ.2) then
       call et_2e(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &           g_phandle(1),g_phandle(2),omega2)
      endif
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
       Hba = omega1+omega2
c
      if(me.eq.0) then
       write(6,9995) 'Two-electron contribution        ','H2(rp)',omega2
       write(6,9995) 'Total interaction energy         ','H(rp)',Hba
       write(6,*)
      endif
c
c
c      =============
c      calculate VBA
c      =============
c
c
       Vba =  Hba - Sba*(Haa+Hbb)*0.5d0
       Vba = Vba/(1.0d0-Sba*Sba)
       Vba = dabs (Vba)
c
       if (.not. rtdb_put(rtdb,'et:Vrp', MT_DBL, 1, Vba))
     $        call errquit('et_calc: rtdb_put of Vrp failed',0)
c
c
c      =============
c      Print results
c      =============
c
      if(me.eq.0) then
      write(6,9996) Vba,Vba*2.194746d5,
     &              Vba*27.2114d0,Vba*627.510d0
      write (6,*) 
      endif
c
 9995 format (a34,1x,a6,1x,F18.10)
 9996 format( 1x,'Electron Transfer Coupling Energy, |Vrp|',1x,F18.10,/, 
     &  49x,F11.3,' cm-1',/,
     &  49x,F11.6,' eV',/,
     &  49x,F11.3,' kcal/mol')
c
      if (.not.ma_chop_stack(l_vals)) 
     &    call errquit('unable to destroy MA',l_vals)
      if (.not.ma_chop_stack(l_eval)) 
     &    call errquit('unable to destroy MA',l_eval)
      if (.not.ma_chop_stack(l_occ)) 
     &    call errquit('unable to destroy MA',l_occ)
c
c
      et_calc=.true.
c
      call util_print_pop
c
         if (.not. geom_destroy(geom)) call errquit
     $        ('et_calc: geom destroy failed', 0)
         if (.not. bas_destroy(basis)) call errquit
     $        ('et_calc: basis destroy failed',0)
c
      call int_terminate()
c
c
      return
      end


