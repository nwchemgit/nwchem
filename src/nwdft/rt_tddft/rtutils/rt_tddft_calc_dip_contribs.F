C
C
C     (i) Movecs in ao basis and the inverse of this matrix.
C     C = X C' 
C
C     need C^-1  and (C+)^-1
C     
C     (ii) Convert P and D into the *true* MO basis (not the orthogonal basis):
C      
C     Pmo = C^-1 P (C+)^-1
C     Dmo = C+ D C
C
C     (iii) Compute all the MO dipole contributions:
C
C     mu_{ia} = Dmo_{ia}*Pmo_{ai} + Dmo_{ai}*Pmo_{ia}
C
C     print Re[mu_{ia}]
C
CC      subroutine rt_tddft_calc_dip_contribs (params,
CC     $     g_movecs_mo_gs, g_zdens_ao, tt)
      subroutine rt_tddft_calc_dip_contribs (params,
     $     g_movecs_mo_gs, g_zdens_ao, tt, abc)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"
#include "matutils.fh"
      
C     inputs
CCC      type(rt_params_t), intent(in) :: params
      type(rt_params_t), intent(inout) :: params !INTENT
      integer, intent(in)           :: g_movecs_mo_gs   ! C', the eigenvectors of the Fock matrix in the orthogonal basis
      integer, intent(in)           :: g_zdens_ao    ! complex density matrix in the AO basis
      double precision, intent(in)  :: tt      

      character(*), intent(in)      :: abc    !Alpha, Beta, Total tag

      character(len=15) :: str

      character(len=*), dimension(3), parameter ::
     $     diptags = [
     $     "#Mu-X i-a Contrib ",
     $     "#Mu-Y i-a Contrib ",
     $     "#Mu-Z i-a Contrib " ]

C     local variables.  Note "gsmo" denotes quantity in the "true" ground state MO basis.
      integer :: n, m, i, j, id
      integer :: me
      integer :: g_zXloc
      integer :: g_zCprime
      integer :: g_zCao, g_zCao_inv, g_zCao_dag, g_zCao_dag_inv
      integer :: g_zDao(3)
      integer :: g_zDgsint(3)   !intermediate ground state dipole tensor
      integer :: g_zDgsmo(3), g_zPgsmo(3)
      integer :: g_zPgsint(3)   !intermediate ground state density matrix
      integer :: g_zCao_inv_int, g_zCao_inv_int2   !intermediate parts of the L and R operators
      integer :: g_zMuia_ia, g_zMuia_ai
      integer :: g_zDgsmo_trans(3), g_zPgsmo_trans(3), g_zMu_ia(3)
      integer :: g_densre_ao, g_densim_ao !real and imaginary parts of the density matix in the AO basis
      double complex :: zMu(3)

      integer :: istart, iend, jstart, jend
      integer :: ipol_in, idstart, idend
      
      me = ga_nodeid()


C
C     allocate GAs
C
      n = params%nbf_ao
      m = params%nbf_mo


      if (.not. ga_create(mt_dcpl, m, m, "zCprime", 0, 0, g_zCprime))
     $     call errquit("Failed to allocate g_zCprime", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, n, m, "zCao", 0, 0, g_zCao))
     $     call errquit("Failed to allocate g_zCao", 0, GA_ERR)


      if (.not. ga_create(mt_dcpl, m, n, "zCao_inv", 0, 0, g_zCao_inv))
     $     call errquit("Failed to allocate g_zCao_inv", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, m, m, "zCao_inv_int", 
     $     0, 0, g_zCao_inv_int))
     $     call errquit("Failed to allocate g_zCao_inv_int", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, m, m, "zCao_inv_int2", 
     $     0, 0, g_zCao_inv_int2))
     $     call errquit("Failed to allocate g_zCao_inv_int2", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, m, n, "zCao_dag", 0, 0, g_zCao_dag))
     $     call errquit("Failed to allocate g_zCao_dag", 0, GA_ERR)
      
      if (.not. ga_create(mt_dcpl, n, m, "zCao_dag_inv",  
     $     0, 0, g_zCao_dag_inv))
     $     call errquit("Failed to allocate g_zCao_dag_inv", 0, GA_ERR)

      do id = 1, 3
         if (.not. ga_create(mt_dcpl, n, n, "zDao", 0, 0, g_zDao(id)))
     $        call errquit("Failed to allocate g_zDao", 0, GA_ERR)

         if (.not. ga_create(mt_dcpl, m, m, "zDgsmo",
     $        0, 0, g_zDgsmo(id)))
     $        call errquit("Failed to allocate g_zDgsmo", 0, GA_ERR)

         if (.not. ga_create(mt_dcpl, m, m, "zPgsmo",
     $        0, 0, g_zPgsmo(id)))
     $        call errquit("Failed to allocate g_zPgsmo", 0, GA_ERR)

         if (.not. ga_create(mt_dcpl, m, m, "zDgsmo_trans",
     $        0, 0, g_zDgsmo_trans(id))) 
     $        call errquit("Failed to allocate g_zDgsmo_trans", 0, 
     $        GA_ERR)

         if (.not. ga_create(mt_dcpl, m, m, "zPgsmo_trans",
     $        0, 0, g_zPgsmo_trans(id))) 
     $        call errquit("Failed to allocate g_zPgsmo_trans", 0, 
     $        GA_ERR)

         if (.not. ga_create(mt_dcpl, n, m, "zDgsint", 0, 0, 
     $        g_zDgsint(id))) 
     $        call errquit("Failed to allocate g_zDgsint", 0, GA_ERR) 

         if (.not. ga_create(mt_dcpl, n, m, "zPgsint", 0, 0, 
     $        g_zPgsint(id))) 
     $        call errquit("Failed to allocate g_zPgsint", 0, GA_ERR)

         if (.not. ga_create(mt_dcpl, m, m, "zMu_ia",
     $        0, 0, g_zMu_ia(id))) 
     $        call errquit("Failed to allocate g_zMu_ia", 0, GA_ERR)
      enddo


      if (.not. ga_create(mt_dbl, n, n, "g_densre_ao", 0, 0, 
     $     g_densre_ao)) 
     $     call errquit("Failed to allocate g_densre_ao", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, n, n, "g_densim_ao", 0, 0, 
     $     g_densim_ao)) 
     $     call errquit("Failed to allocate g_densim_ao", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, m, m, "g_zMuia_ia", 0, 0, 
     $     g_zMuia_ia)) 
     $     call errquit("Failed to allocate g_zMuia_ia", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, m, m, "g_zMuia_ai", 0, 0, 
     $     g_zMuia_ai)) 
     $     call errquit("Failed to allocate g_zMuia_ai", 0, GA_ERR)


C
C     (i) Convert movecs to MO basis and invert
C
C     cast orthogonal movecs to complex data type
      call convert_d2z (1d0, g_movecs_mo_gs,
     $     0d0, g_movecs_mo_gs, g_zCprime)

C     cast dipole integrals to cmplx dbl
      call convert_d2z(1d0, params%g_dip%x,
     $     0d0, params%g_dip%x, g_zDao(1)) !assuming x is 1
C     dipole y to cmplx dbl
      call convert_d2z(1d0, params%g_dip%y,
     $     0d0, params%g_dip%y, g_zDao(2)) !assuming y is 2
C     dipole z to cmplx dbl
      call convert_d2z(1d0, params%g_dip%z,
     $     0d0, params%g_dip%z, g_zDao(3)) !assuming z is 3
      

C
C     (ii)
C

C     C = X C'
      call ga_zero(g_zCao)
      call ga_zgemm("N","N", n, m, m,
     $        (1d0,0d0), params%g_zX, g_zCprime, (0d0,0d0), g_zCao)
      
C     C^+ = conj(transpose(C)) (Adjoint / Hermitian Conjugate)
      call ga_zero(g_zCao_dag)
C      call ga_print(g_zCao)
      call mat_hermconj(g_zCao,g_zCao_dag)
CC      call ga_print(g_zCao_dag)

C     C^-1 and (C^+)^-1 (inverses)
      call ga_zgemm("N","N", m, m, n, (1d0,0d0),
     Z     g_zCao_dag, g_zCao, (0d0,0d0), g_zCao_inv_int)
      call zmat_inv(g_zCao_inv_int, g_zCao_inv_int2)
      call ga_zgemm("N","N", m, n, m, (1d0,0d0),
     Z     g_zCao_inv_int2, g_zCao_dag, (0d0,0d0), g_zCao_inv) !Left Operator
      call ga_zgemm("N","N", n, m, m, (1d0,0d0),
     Z     g_zCao, g_zCao_inv_int2, (0d0,0d0), g_zCao_dag_inv) !Right Operator

C     Pmo = (C^-1) P ((C^+)^-1)
      call ga_zero(g_zPgsint(1)) 
      call ga_zero(g_zPgsmo(1))
      call ga_Zgemm("N","N", n, m, n,(1d0,0d0), g_zdens_ao, 
     $     g_zCao_dag_inv, (0d0,0d0), g_zPgsint(1))
      
      call ga_Zgemm("N","N", m, m, n, (1d0,0d0), g_zCao_inv, 
     $     g_zPgsint(1), 0d0, g_zPgsmo(1))

      call ga_zero(g_zPgsint(2)) 
      call ga_zero(g_zPgsmo(2))
      call ga_Zgemm("N","N", n, m, n, (1d0,0d0), g_zdens_ao, 
     $     g_zCao_dag_inv, (0d0,0d0), g_zPgsint(2))
      
      call ga_Zgemm("N","N", m, m, n, (1d0,0d0), g_zCao_inv, 
     $     g_zPgsint(2), (0d0,0d0), g_zPgsmo(2))

      call ga_zero(g_zPgsint(3))
      call ga_zero(g_zPgsmo(3))
      call ga_Zgemm("N","N", n, m, n, (1d0,0d0), g_zdens_ao, 
     $     g_zCao_dag_inv, (0d0,0d0), g_zPgsint(3))
      
      call ga_Zgemm("N","N", m, m, n, (1d0,0d0), g_zCao_inv, 
     $     g_zPgsint(3), (0d0,0d0), g_zPgsmo(3))

C      call halt ()

C     Dmo = (C^+) D C (only in the ground state (gs) / for weak-field)
      call ga_zero(g_zDgsint(1)) !clear the variables
      call ga_zero(g_zDgsmo(1))
      call ga_Zgemm("N","N", n, m, n, (1d0,0d0), g_zDao(1), 
     $     g_zCao, (1d0,0d0), g_zDgsint(1))
      call ga_Zgemm("N","N", m, m, n, (1d0,0d0), g_zCao_dag, 
     $     g_zDgsint(1), (1d0,0d0), g_zDgsmo(1))

      call ga_zero(g_zDgsint(2)) 
      call ga_zero(g_zDgsmo(2))
      call ga_Zgemm("N","N", n, m, n, (1d0,0d0), g_zDao(2), 
     $     g_zCao, (1d0,0d0), g_zDgsint(2))
      call ga_Zgemm("N","N", m, m, n, (1d0,0d0), g_zCao_dag, 
     $     g_zDgsint(2), (1d0,0d0), g_zDgsmo(2))

      call ga_zero(g_zDgsint(3)) 
      call ga_zero(g_zDgsmo(3))
      call ga_Zgemm("N","N", n, m, n, (1d0,0d0), g_zDao(3), 
     $     g_zCao, (1d0,0d0), g_zDgsint(3))
      call ga_Zgemm("N","N", m, m, n, (1d0,0d0), g_zCao_dag, 
     $     g_zDgsint(3), (1d0,0d0), g_zDgsmo(3))


C     Calculate mu_ia values
C     (Dmo)^T and (Pmo)^T

      call ga_transpose(g_zDgsmo(1), g_zDgsmo_trans(1))
      call ga_transpose(g_zPgsmo(1), g_zPgsmo_trans(1))

      call ga_transpose(g_zDgsmo(2), g_zDgsmo_trans(2))
      call ga_transpose(g_zPgsmo(2), g_zPgsmo_trans(2))

      call ga_transpose(g_zDgsmo(3), g_zDgsmo_trans(3))
      call ga_transpose(g_zPgsmo(3), g_zPgsmo_trans(3))


C     mu_ia = Dia * Pai + Dai * Pia
      call ga_zero(g_zMuia_ia) !clear the variables
      call ga_zero(g_zMuia_ai)
      call ga_elem_multiply(g_zDgsmo(1), 
     $     g_zPgsmo_trans(1), g_zmuia_ia)
      call ga_elem_multiply(g_zDgsmo_trans(1), 
     $     g_zPgsmo(1), g_zMuia_ai)
      call ga_add((1d0,0d0), g_zMuia_ia, (1d0,0d0), 
     $     g_zMuia_ai, g_zMu_ia(1))

      call ga_zero(g_zMuia_ia) !clear the variables
      call ga_zero(g_zMuia_ai)
      call ga_elem_multiply(g_zDgsmo(2), 
     $     g_zPgsmo_trans(2), g_zmuia_ia)
      call ga_elem_multiply(g_zDgsmo_trans(2), 
     $     g_zPgsmo(2), g_zMuia_ai)
      call ga_add((1d0,0d0), g_zMuia_ia, (1d0,0d0), 
     $     g_zMuia_ai, g_zMu_ia(2))

      call ga_zero(g_zMuia_ia) !clear the variables
      call ga_zero(g_zMuia_ai)
      call ga_elem_multiply(g_zDgsmo(3), 
     $     g_zPgsmo_trans(3), g_zmuia_ia)
      call ga_elem_multiply(g_zDgsmo_trans(3), 
     $     g_zPgsmo(3), g_zMuia_ai)
      call ga_add((1d0,0d0), g_zMuia_ia, (1d0,0d0), 
     $     g_zMuia_ai, g_zMu_ia(3))

#if 0
C     Build Output
      if (me.eq.0) then
         write (*, FMT=1,advance="yes"), trim(params%tag), tt,
     Z         "### MO Dipole Contributions ###"
    1 FORMAT (A, 2X, 1F11.5, 8X, A)
    2 FORMAT (A, 2X, 1F11.5, 8X, A, 1X, I5, 1X, I5, 3X, 1ES14.7,
     Z       2X, 1ES14.7)
C     Print tag and finish line
         write(luout, *) "    "

         do i=1,n
            do j=i,m
               call ga_get(g_zMu_ia(1),i,i,j,j,zMu(1),0)
               write(*,FMT=3,advance="yes"), trim(params%tag), tt, 
C     Z         "#Mu-X i-a Contrib", i, j, zMu(1)
     Z         zMu(1), "#Mu-X i-a Contrib", i, j
               call ga_get(g_zMu_ia(2),i,i,j,j,zMu(2),0)
               write(*,FMT=3,advance="yes"), trim(params%tag), tt, 
C     Z         "#Mu-Y i-a Contrib", i, j, zMu(2)
     Z         zMu(2), "#Mu-Y i-a Contrib", i, j
               call ga_get(g_zMu_ia(3),i,i,j,j,zMu(3),0)
               write(*,FMT=3,advance="yes"), trim(params%tag), tt, 
     Z         zMu(3), "#Mu-Z i-a Contrib", i, j
            enddo
        enddo
    3 FORMAT (A, 2X, 1F11.5, 8X, 1ES14.7, 2x, 1ES14.7, 3X, A, 1X, 
     Z       I5, 1X, I5)
      endif
#endif

C
C     New compact output
C

C
C     Determine MO windows (which MOs to calculate contribs for).  If
C     not found in rtdb, then fall back to full MO range.


C	NEW WINDOWING; Turning on the ranges. [ASB] 
      if (.not. rtdb_get (params%rtdb, "rt_tddft:dipcontrib:mo1min",
     $     mt_int, 1, istart)) then
         istart = 1
      else
         istart = istart
      endif
      
      if (.not. rtdb_get (params%rtdb, "rt_tddft:dipcontrib:mo1max",
     $     mt_int, 1, iend)) then
         iend = m
      else
         iend = iend
      endif

      if (.not. rtdb_get (params%rtdb, "rt_tddft:dipcontrib:mo2min",
     $     mt_int, 1, jstart)) then
         jstart = 1
      else
         jstart = jstart
      endif
      
      if (.not. rtdb_get (params%rtdb, "rt_tddft:dipcontrib:mo2max",
     $     mt_int, 1, jend)) then
         jend = m
      else
         jend = jend
      endif

C     (debugging info)
      if (me.eq.0) then
         write (6,*) "istart = ", istart, "; iend = ", iend, 
     $     "; jstart = ", jstart, "; jend = ", jend
      endif

C    MO2 windows disabled for now; also note that loops over j (MO2) go
C    from i to m.  This will need to be fixed if doing windows for MO2.

C     (check the windows)
      if (istart > iend)
     $     call errquit("istart cannot be > iend",0,0)
      if (jstart > jend)
     $     call errquit("jstart cannot be > jend",0,0)

      if (istart > m)
     $     call errquit("istart cannot be > number of MOs",0,0)
      if (iend > m)
     $     call errquit("iend cannot be > number of MOs",0,0)
      if (jstart > m)
     $     call errquit("jstart cannot be > number of MOs",0,0)
      if (jend > m)
     $     call errquit("jend cannot be > number of MOs",0,0)

      if (istart < 1)
     $     call errquit("istart cannot be < 1",0,0)
      if (iend < 1)
     $     call errquit("iend cannot be < 1",0,0)
      if (jstart < 1)
     $     call errquit("jstart cannot be < 1",0,0)
      if (jend < 1)
     $     call errquit("jend cannot be < 1",0,0)


C
C     Which polarizations to print.  If user does not specify, it will
C     print all three.
C
      if (.not. rtdb_get (params%rtdb, "rt_tddft:dipcontrib:pol",
     $     mt_int, 1, ipol_in)) then
         idstart = 1  
         idend = 3
      else

         if ((ipol_in.ne.1).and.(ipol_in.ne.2).and.(ipol_in.ne.3))
     $        call errquit("pol must be 1, 2, or 3", 0, 0)
         
         idstart = ipol_in   
         idend = ipol_in   
      endif
      
      
      if (me.eq.0) then

         do id = idstart, idend           !x,y,z polarizations
            
C     (header)
            if (tt < params%dt) then ! print header at t=0
               write(luout, fmt="(a15)", advance="no") ""
               do i = istart, iend      
!                  do j = i, jend     !old way
                  do j = jstart, jend   
                     write(str,"(i0)") i
                     write(luout, fmt="(a15)", advance="no") str
                  enddo
               enddo
               write(luout,*) diptags(id), abc
               
               write(luout, fmt="(a15)", advance="no") ""
               do i = istart, iend      
!                  do j = i, jend     !old way
                  do j = jstart, jend   
                     write(str,"(i0)") j
                     write(luout, fmt="(a15)", advance="no") str
                  enddo
               enddo
               write(luout,*) diptags(id), abc
            endif
               
C     (dip contribs)
            write(luout, "(e15.7)", advance="no") dble(tt)
            do i = istart, iend
!               do j = i, jend     !old way
               do j = jstart, jend
                  call ga_get(g_zMu_ia(id),i,i,j,j,zMu(id),0)
                  write(luout, "(e15.7)", advance="no") dble(zMu(id))   ! print real part only
               enddo
            enddo
            write(luout,*) diptags(id), abc 
            
         enddo                  ! loop over dimensions
      endif                     ! proc 0

      call util_flush (luout)

C
C     Clean up
C
      if (.not.ga_destroy(g_zCprime))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zCao))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zCao_inv))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zCao_inv_int))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zCao_inv_int2))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zCao_dag))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zCao_dag_inv))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zMuia_ia))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_zMuia_ai))
     $     call errquit("Destroy failed", 0, GA_ERR)

      if (.not.ga_destroy(g_densre_ao)) 
     $     call errquit("Destroy failed", 0, GA_ERR)
      if (.not.ga_destroy(g_densim_ao)) 
     $     call errquit("Destroy failed", 0, GA_ERR)

      do id = 1, 3
         if (.not.ga_destroy(g_zDao(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)
         if (.not.ga_destroy(g_zDgsmo(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)
         if (.not.ga_destroy(g_zPgsmo(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)

         if (.not.ga_destroy(g_zDgsmo_trans(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)
         if (.not.ga_destroy(g_zPgsmo_trans(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)

         if (.not.ga_destroy(g_zDgsint(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)
         if (.not.ga_destroy(g_zPgsint(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)

         if (.not.ga_destroy(g_zMu_ia(id)))
     $        call errquit("Destroy failed", 0, GA_ERR)
      enddo
      
      end subroutine

      
