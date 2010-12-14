      logical function et_calc(rtdb,Vba)
      implicit none
#include "errquit.fh"
c
c     This function calculates Vba, the electron transfer matrix element.
*
* $Id$
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
#include "stdio.fh"
c
      integer rtdb, me
      integer i, j, k
      integer l_vals,k_vals
      integer g_tmp,g_s,g_ahandle(2),g_bhandle(2),g_d,g_temp
      integer nocc(2)
      integer g_u,g_vt,g_phandle(2)
      integer g_av,g_bu,g_t,g_h1
      integer g_avecsa,g_avecsb,g_bvecsa,g_bvecsb
      integer mem2, max2e, lena, lenb
      integer method_2e, anmo(2), bnmo(2), nbf_file, nset

      logical et_movecs_read
      external et_movecs_read
      logical ga_copy
      external ga_copy
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read_header
      external movecs_read_header

      character*255 a_file, b_file, basis_name,
     E     a_filedir,b_filedir
      character*2 pname(2)

      double precision Vba
      double precision prod(2)
      double precision Sba,omega1
      double precision detu(2),detv(2),detprod,omega2
      double precision Haa,Hbb,Hba,val, dum
      data Haa,Hbb,Hba/3*0.0d0/

      me=ga_nodeid()

      if (me.eq.0) then
      Call Util_Print_Centered(LuOut,'Electron Transfer Calculation'
     &     ,40,.true.)
      write(LuOut,*)
      endif
c
c -------------
c for printing:
c -------------
c
      call util_print_push
      call util_print_rtdb_load(rtdb,"et")
c
c -------------------------------
c load geometry and symmetry info
c -------------------------------
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('et_calc: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('et_calc: no geometry ', 0, RTDB_ERR)
c
c ----------------------------------------
c load the basis set and get info about it
c ----------------------------------------
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('et_calc: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $        call errquit('et_calc: no ao basis set', 0, RTDB_ERR)
c
c  -------
c  get nbf
c  -------
c
      if (.not.bas_numbf(basis,nbf))
     $     call errquit('bas_numbf failed',0, BASIS_ERR)
c
c ------------------
c initialize the API 
c ------------------
c
      call int_init(rtdb, 1, basis)
      call int_acc_std()
      call int_normalize(rtdb,basis)
c
c  -------------------------------
c  calculate the overlap integrals
c  -------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'s',0,0,g_s))
     $        call errquit('task_et ga_create failed', 0, GA_ERR)
      call ga_zero(g_s)
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)  

      if (util_print("et",print_debug)) then
       call ga_print(g_s)
      endif
c
c ---------------------------
c Read MO vectors from files
c ---------------------------
c
      if (.not. rtdb_cget(rtdb, 'a_vecs', 1, a_file))
     $     call errquit('et_calc: couldnt get a_vecs filename',0,
     &       RTDB_ERR)
c
      if (.not. rtdb_cget(rtdb, 'b_vecs', 1, b_file))
     $     call errquit('et_calc: couldnt get a_vecs filename',0,
     &       RTDB_ERR)
c
      lena = inp_strlen(a_file)
      lenb = inp_strlen(b_file)
c
      if (me.eq.0) then
      write(LuOut,*)'MO vectors for reactants: ',a_file(1:lena)
      write(LuOut,*)'MO vectors for products : ',b_file(1:lenb)
      write(LuOut,*)
      endif
c
c -----------------------------------------------------
c Allocate occupation and eigenvalue vectors. 
c 'occ' is used to calculate nocc, 'eval' is never used.
c ------------------------------------------------------
c
      if (.not.ma_push_get(mt_dbl,nbf,'occ',l_occ,k_occ))
     $         call errquit('ma_push_get failed for occ',k_occ, MA_ERR)

      if (.not.ma_push_get(mt_dbl,nbf,'evals',l_eval,k_eval))
     $        call errquit('ma_push_get failed for evals',k_eval,
     &       MA_ERR)
c
c ------------------------------------------------
c Read header to check that nset=2
c ------------------------------------------------
c
c a_file:
c -------
c
      call et_filename(a_file,a_filedir)
      if(.not.movecs_read_header(a_filedir(1:inp_strlen(a_filedir)),
     E     title,basis_name,
     &   scftype,nbf_file,nset,anmo,2)) 
     &   call errquit('et_calc: movecs_read_header failed',me, DISK_ERR)

      if (nset.ne.2) then
        call errquit('et_calc: only UHF/ODFT orbitals allowed',nset,
     &       INPUT_ERR)
      endif

      if (nbf_file.ne.nbf) then
        call errquit('et_calc: basis function mismatch',nbf_file,
     &       BASIS_ERR)
      endif

c b_file:
c -------
c
      call et_filename(b_file,b_filedir)
      if(.not.movecs_read_header(b_filedir(1:inp_strlen(b_filedir)),
     E     title,basis_name,
     &   scftype,nbf_file,nset,bnmo,2)) 
     &   call errquit('et_calc: movecs_read_header failed',me, DISK_ERR)

      if (nset.ne.2) then
        call errquit('et_calc: only UHF/ODFT orbitals allowed',nset,
     &       INPUT_ERR)
      endif

      if (nbf_file.ne.nbf) then
        call errquit('et_calc: basis function mismatch',nbf_file,
     &       BASIS_ERR)
      endif
c
c -------------------------------------------------------
c allocate arrays for MO coefficients, for state A
c -------------------------------------------------------
c
      if (.not.ga_create(mt_dbl,nbf,nbf,'A alpha',0,0,g_avecsa))
     $         call errquit('ga_create failed for A alpha',g_avecsa,
     &       GA_ERR)

      if (.not.ga_create(mt_dbl,nbf,nbf,'A beta',0,0,g_avecsb))
     $         call errquit('ga_create failed for A beta',g_avecsb,
     &       GA_ERR)

      call ga_zero (g_avecsa)
      call ga_zero (g_avecsb)
c 
c ------------------------------------------------------------------
c read occ, evals, mo coefficients for alpha orbitals from a_file
c ------------------------------------------------------------------
c
      if(.not.et_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     A     1,dbl_mb(k_occ)
     $   ,dbl_mb(k_eval),g_avecsa,dum))
     $ call errquit('et_movecs_read failed for alpha mos',0, DISK_ERR)
c
c -------------------------------------
c get nocc(1): # of occupied alpha mo's
c -------------------------------------
c
      nocc(1)=0
      do i = 1,nbf
       nocc(1) = nocc(1) + dbl_mb(k_occ+i-1)
      enddo

      if (util_print("et",print_high)) then
      if(me.eq.0) then
        write(LuOut,*) nocc(1), " occupied alpha orbitals"
        write(LuOut,*)
      endif
        call ga_print(g_avecsa)
      endif
c
c ------------------------------------------------------------------
c read occ, evals, mo coefficients for beta orbitals from a_file
c ------------------------------------------------------------------
c
      if(.not.et_movecs_read(a_filedir(1:inp_strlen(a_filedir)),
     E     2,dbl_mb(k_occ)
     $     ,dbl_mb(k_eval),g_avecsb, dum))
     $ call errquit('et_movecs_read failed for beta mos',0, DISK_ERR)
c
c -------------------------------------
c get nocc(2): # of occupied alpha mo's
c -------------------------------------
c
      nocc(2)=0
      do i = 1,nbf
       nocc(2) = nocc(2) + dbl_mb(k_occ+i-1)
      enddo

      if (util_print("et",print_high)) then
      if(me.eq.0) then
        write(LuOut,*) nocc(2), " occupied  beta orbitals"
        write(LuOut,*)
      endif
        call ga_print(g_avecsb)
      endif
c
c -------------------------------------------------------
c allocate arrays for MO coefficients, for state B
c -------------------------------------------------------
c
      if (.not.ga_create(mt_dbl,nbf,nbf,'B alpha',0,0,g_bvecsa))
     $         call errquit('ga_create failed for B alpha',g_bvecsa,
     &       GA_ERR)

      if (.not.ga_create(mt_dbl,nbf,nbf,'B beta',0,0,g_bvecsb))
     $         call errquit('ga_create failed for B beta',g_bvecsb,
     &       GA_ERR)

      call ga_zero (g_bvecsa)
      call ga_zero (g_bvecsb)
c
c -----------------------------------------------------------
c read occ, evals, mo coefficients from b_file.  It is 
c assumed that the number of occupied alpha and beta orbitals
c is the same as in the a_file.
c -----------------------------------------------------------
c
      if(.not.et_movecs_read(b_filedir(1:inp_strlen(b_filedir)),
     E     1,dbl_mb(k_occ)
     $   ,dbl_mb(k_eval) ,g_bvecsa, dum))
     $   call errquit('et_movecs_read failed for alpha mos',0, DISK_ERR)
c
      if(.not.et_movecs_read(b_filedir(1:inp_strlen(b_filedir)),
     E     2,dbl_mb(k_occ)
     $  ,dbl_mb(k_eval) ,g_bvecsb, dum))
     $   call errquit('et_movecs_read failed for beta mos',0, DISK_ERR)
c
      if (util_print("et",print_high)) then
       call ga_print(g_bvecsa)
       call ga_print(g_bvecsb)
      endif
c
c --------------------------------------------------
c Read converged scf or dft energy from movecs file.
c --------------------------------------------------
c
      if(.not.et_movecs_read(a_filedir(1:inp_strlen(a_filedir)),3
     $   ,dbl_mb(k_occ),dbl_mb(k_eval),1, Haa)) 
     $   call errquit('et_movecs_read failed for Haa',0, DISK_ERR)
c
      if(.not.et_movecs_read(b_filedir(1:inp_strlen(b_filedir)),3
     $   ,dbl_mb(k_occ),dbl_mb(k_eval),1, Hbb)) 
     $   call errquit('et_movecs_read failed for Hbb',0, DISK_ERR)
c
      if (me.eq.0) then
      write(LuOut,9995) 'Electronic energy of reactants   ','H(RR)',Haa
      write(LuOut,9995) 'Electronic energy of products    ','H(PP)',Hbb
      endif
c
c -----------------------------------------------
c calculate D=B'*S*A for alpha and beta electrons
c -----------------------------------------------
c
      g_ahandle(1)=g_avecsa
      g_ahandle(2)=g_avecsb
      g_bhandle(1)=g_bvecsa
      g_bhandle(2)=g_bvecsb
c
c
c detprod is the product of det(U)*det(V) for alpha and beta electrons,
c and is used later.
c ----------------------------------------------------------------------
c
      detprod=1.0d0
c
c -----------------------------------------------------
c  pname is needed for generalized density matrices:
c  Pa is for alpha electrons, Pb is for beta electrons.
c -----------------------------------------------------
c
      pname(1)='pa'
      pname(2)='pb'
c
c ------------------------
c begin alpha/beta loop j
c ------------------------
c
      do j = 1,2 
c
      if (util_print("et",print_high)) then
      if(me.eq.0) then
       write(LuOut,*) 
       write(LuOut,*) '---------------'
       if(j.eq.1) write(LuOut,*) 'ALPHA ELECTRONS'
       if(j.eq.2) write(LuOut,*) ' BETA ELECTRONS'
       write(LuOut,*) '---------------'
      endif
      endif
c
c  -----------------------------------------------
c  Perform matrix multiplication to get D = B'*S*A
c  -----------------------------------------------
c
      call ga_sync()
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'d',0,0,g_d))
     $        call errquit('task_et ga_create failed', g_d, GA_ERR)
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'temp',0,0,g_tmp))
     $        call errquit('task_et ga_create failed', g_d, GA_ERR)
c
      call ga_dgemm('T','N',nocc(j),nbf,nbf,1.d0,g_bhandle(j)
     &               ,g_s,0.d0,g_tmp)
c
      call ga_dgemm('N','N',nocc(j),nocc(j),nbf,1.d0,g_tmp
     &               ,g_ahandle(j),0.d0,g_d)
c
      if (.not.ga_destroy(g_tmp))
     &  call errquit('could not destroy tmp array',g_tmp, GA_ERR)
c
      if (util_print("et",print_debug)) then
       call ga_print(g_d)
      endif
c
c ----------------------------------------------
c Perform SVD to get U, V', and singluar values
c ----------------------------------------------
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'u',0,0,g_u))
     $        call errquit('task_et ga_create failed', g_u, GA_ERR)
c
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'vt',0,0,g_vt))
     $        call errquit('task_et ga_create failed', g_vt, GA_ERR)
c
      call ga_zero(g_u)
      call ga_zero(g_vt)
      if (.not.ma_push_get(mt_dbl,nocc(j),'values',l_vals,k_vals))
     $        call errquit('ma_push_get failed for values ',k_vals,
     &       MA_ERR)
c
      call ga_svd_seq(g_d,g_u,g_vt,dbl_mb(k_vals))
c
c prod is the product of singular values.
c ---------------------------------------
c
      prod(j) = 1.0d0
      do k = 1,nocc(j)
         prod(j) = prod(j)*dbl_mb(k_vals+k-1)
      enddo
c      write(LuOut,*)'prod',j,prod(j)
c      write(LuOut,*)'singular values',(dbl_mb(k_vals+i-1),i=1,nocc(j))
c
c
      if (.not.ga_destroy(g_d))
     &  call errquit('could not destroy D array',g_d, GA_ERR)
c
      if (util_print("et",print_debug)) then
      call ga_print(g_u)
      call ga_print(g_vt)
      endif
c
c
c --------------------------
c calculate det(U) and det(V)
c --------------------------
c
       call et_det('N',g_u,nocc(j),detu(j))
       call et_det('T',g_vt,nocc(j),detv(j))
       detprod = detprod*detu(j)*detv(j)
c       write(LuOut,*) j,'detu,detv',detu(j),detv(j)
c
c ------------------------------------------------------------
c                                               ^       ^
c perform matrix mult to get transformed A & B: A = AV, B = BU
c ------------------------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'av',0,0,g_av))
     $        call errquit('task_et ga_create failed', g_av, GA_ERR)
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'bu',0,0,g_bu))
     $        call errquit('task_et ga_create failed', g_bu, GA_ERR)
c
      call ga_zero(g_av)
      call ga_dgemm('N','T',nbf,nocc(j),nocc(j),1.d0,g_ahandle(j),g_vt
     &         ,0.d0,g_av)
c
      call ga_zero(g_bu)
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_bhandle(j),g_u
     &         ,0.d0,g_bu)
c
c
      if (.not.ga_destroy(g_ahandle(j))) 
     &  call errquit('could not destroy mo vectors array',g_ahandle(j),
     &       GA_ERR)
      if (.not.ga_destroy(g_bhandle(j))) 
     &  call errquit('could not destroy mo vectors array',g_bhandle(j),
     &       GA_ERR)
c
      if (.not.ga_destroy(g_u))
     &  call errquit('could not destroy U array',g_u,
     &       GA_ERR)
      if (.not.ga_destroy(g_vt))
     &  call errquit('could not destroy VT array',g_vt,
     &       GA_ERR)
c
      if (util_print("et",print_high)) then
      if(me.eq.0) then
        write(LuOut,*) 
        write(LuOut,*) 'SVD transformed MO vectors for state A'
      endif
        call ga_print (g_av)
      if(me.eq.0) then
        write(LuOut,*) 
        write(LuOut,*) 'SVD transformed MO vectors for state B'
      endif
        call ga_print (g_bu)
      endif
c
c
c    --------------------------------------------------------------
c                                   N
c                                -------
c     Calculates   T(i) =     PRODUCT  d(j) where N = nalpha or nbeta
c                                -------
c                                j.NE.i  
c    --------------------------------------------------------------
c    
      if(.not.ga_create(mt_dbl,nocc(j),nocc(j),'t',0,0,g_t))
     $        call errquit('task_et ga_create failed', g_t, GA_ERR)
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
      call ga_sync()
c
      if (util_print("et",print_debug)) then
      call ga_print(g_t) 
      endif
c
c -----------------------------------------------------
c                                            ^ ^
c  calculate generalized density matrix, P = ATB'
c  Pa is for alpha electrons, Pb is for beta electrons.
c -----------------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,pname(j),0,0,g_phandle(j)))
     $        call errquit('task_et ga_create failed', g_phandle(j),
     &       GA_ERR)
c
      if(.not.ga_create(mt_dbl,nbf,nocc(j),'temp',0,0,g_temp))
     $        call errquit('task_et ga_create failed', g_temp, GA_ERR)
c
      call ga_zero(g_phandle(j))
c
      call ga_dgemm('N','N',nbf,nocc(j),nocc(j),1.d0,g_av,g_t
     &         ,0.d0,g_temp)
c
      call ga_dgemm('N','T',nbf,nbf,nocc(j),1.0d0,g_temp,g_bu
     &         ,0.0d0,g_phandle(j))
c
c
      if (.not.ga_destroy(g_av)) 
     &  call errquit('could not destroy transformed MO vectors',g_av,
     &       GA_ERR)
      if (.not.ga_destroy(g_bu)) 
     &  call errquit('could not destroy transformed MO vectors',g_bu,
     &       GA_ERR)
c
      if (.not.ga_destroy(g_t)) 
     &  call errquit('could not destroy T array',g_t, GA_ERR)
      if (.not.ga_destroy(g_temp)) 
     &  call errquit('could not destroy temp array',g_temp, GA_ERR)
c
c
      if (util_print("et",print_debug)) then
      call ga_print(g_phandle(j))
      endif
c
c
c --------------------
c end alpha/beta loop:
c --------------------
      enddo    
c
      if (.not.ga_destroy(g_s)) 
     &  call errquit('could not destroy overlap array',g_s, GA_ERR)

c
c Pa is scaled by prod(b), Pb is scaled by prod(a)
c ------------------------------------------------
c
       call ga_dscal(g_phandle(1),prod(2))
       call ga_dscal(g_phandle(2),prod(1))

      if (util_print("et",print_debug)) then
      if(me.eq.0) then
       write(LuOut,*) 
       write(LuOut,*) 'scaled generalized density matrices:' 
       write(LuOut,*) 
      endif
       call ga_print(g_phandle(1))
       call ga_print(g_phandle(2))
      endif
c
c
c -----------------------------------------
c calculate overlap of state A with state B
c -----------------------------------------
c
       Sba = prod(1)*prod(2)*detprod
c
c
      if(me.eq.0) then
       write(LuOut,*)
       write(LuOut,9997) Sba
       write(LuOut,*)
      endif
c
 9997 format(1x,'Reactants/Products overlap S(RP) : ',1Pd9.2)
c
c --------------------------------------------
c calculate the 1e contribution to Vba, omega1
c --------------------------------------------
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'h1',0,0,g_h1))
     $        call errquit('task_et ga_create failed', g_h1, GA_ERR)
c
c get H1 integrals
c ----------------
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
       omega1 = omega1 + ga_ddot (g_h1,g_phandle(1))
       omega1 = omega1 + ga_ddot (g_h1,g_phandle(2))
c
       call ga_sync()
       omega1 = omega1*detprod
c
      if(me.eq.0) then
       write(LuOut,*) 'Reactants/Products interaction energy:    '
       write(LuOut,*) '-------------------------------------           '
       write(LuOut,9995) 'One-electron contribution        ','H1(RP)',
     &    omega1
      endif
c
      if (.not.ga_destroy(g_h1))
     &  call errquit('could not destroy H1 array',g_h1, GA_ERR)
c
c --------------------------------------------
c calculate the 2e contribution to Vba, omega2
c --------------------------------------------
c
      if (.not. rtdb_get(rtdb, 'et:tol2e', MT_DBL, 1, tol2e)) then
       tol2e    =  min(1.0d-7,dabs(Sba)*1.0d-7)
       tol2e    =  max(1.0d-12,tol2e)
      endif

c ------------------------------------------------------------------------
c Two methods are available. One uses the scf fock builder (method_2e=1),
c the other uses et_2e (method_2e=2). et_2e is more accurate, but much more
c time consuming and should only be used for very small systems, or very 
c weakly interacting systems (with small overlap), or for debugging purposes.
c ------------------------------------------------------------------------
c
      if (.not.rtdb_get(rtdb,'et:method_2e',MT_INT,1,method_2e)) then
          method_2e = 1
      endif
c
      if (method_2e.EQ.1) then
c for some unknown reason I have to terminate and re-init the API for
c et_fock to work with ECP's
      call int_terminate()
      call et_fock(rtdb,geom,basis,
     &  g_phandle(1),g_phandle(2),nbf,tol2e,omega2)
      endif
c
      if (method_2e.EQ.2) then
       call int_mem_2e4c(mem2,max2e)
       mem2 = max(mem2,max2e)
       call et_2e(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &           g_phandle(1),g_phandle(2),omega2)
      endif
c
c ----------------------------------------------------
c destroy density arrays and calculate omega2 and Hba
c ----------------------------------------------------
c
      if (.not.ga_destroy(g_phandle(1)))
     &  call errquit('could not destroy P array',g_phandle(1), GA_ERR)
      if (.not.ga_destroy(g_phandle(2)))
     &  call errquit('could not destroy P array',g_phandle(2), GA_ERR)
c
       omega2 = omega2 * (detprod/ ( prod(1)*prod(2) ) ) 
       omega2 = omega2 / 2.0d0
c
       Hba = omega1+omega2
c
      if(me.eq.0) then
       write(LuOut,9995) 'Two-electron contribution        ','H2(RP)',
     &    omega2
       write(LuOut,9995) 'Total interaction energy         ','H(RP)',
     &    Hba
       write(LuOut,*)
      endif
c
c
c ---------------------------------------------
c Finally, calculate the ET matrix element, Vba
c ---------------------------------------------
c
c
       Vba =  Hba - Sba*(Haa+Hbb)*0.5d0
       Vba = Vba/(1.0d0-Sba*Sba)
       Vba = dabs (Vba)
c
       if (.not. rtdb_put(rtdb,'et:V(RP)', MT_DBL, 1, Vba))
     $        call errquit('et_calc: rtdb_put of V(RP) failed',0,
     &       RTDB_ERR)
c
c
c -------------
c Print results
c -------------
c
      if(me.eq.0) then
      write(LuOut,9996) Vba,Vba*2.194746d5,
     &              Vba*27.2114d0,Vba*627.510d0
      write(LuOut,*) 
      endif
c
 9995 format (a34,1x,a6,1x,F18.10)
 9996 format( 1x,'Electron Transfer Coupling Energy |V(RP)|',F18.10,/, 
     &  49x,F11.3,' cm-1',/,
     &  49x,F11.6,' eV',/,
     &  49x,F11.3,' kcal/mol')
c


      if (.not.ma_chop_stack(l_vals)) 
     &    call errquit('unable to destroy MA',l_vals, MA_ERR)
      if (.not.ma_chop_stack(l_eval)) 
     &    call errquit('unable to destroy MA',l_eval, MA_ERR)
      if (.not.ma_chop_stack(l_occ)) 
     &    call errquit('unable to destroy MA',l_occ, MA_ERR)
c
c
      et_calc=.true.
c
      call util_print_pop
c
         if (.not. geom_destroy(geom)) call errquit
     $        ('et_calc: geom destroy failed', 0, GEOM_ERR)
         if (.not. bas_destroy(basis)) call errquit
     $        ('et_calc: basis destroy failed',0, BASIS_ERR)
c
      call int_terminate()
c
c
      return
      end
      subroutine et_filename(file,dirfile)
      implicit none
#include "inp.fh"
#include "global.fh"
      character*(*) file,dirfile
c
      character*255 dirss
      integer lenn0,lenn1,i
c

      call util_directory_name(dirss, .false., ga_nodeid())
      dirfile(1:inp_strlen(dirss))=dirss(1:inp_strlen(dirss))
      dirfile(inp_strlen(dirss)+1:inp_strlen(dirss)+1)='/'
      lenn0=inp_strlen(dirss)+2
      lenn1=inp_strlen(dirss)+inp_strlen(file)+1
      dirfile(lenn0:lenn1)=file(1:inp_strlen(file))
      do i=lenn1+1,255
      dirfile(i:i)=' '
      enddo
      
      return
      end
