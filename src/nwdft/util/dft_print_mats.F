C
C     Prints most DFT matrices to stdout.  This likely
C     replicates other "print" functionality but I've
C     included it here for convenience.
C
C     Enable by: "set dft:print_mats T"
C
C     Prints the following matrices:
C
C     - overlap matrix
C     - eigenvalue list
C     - eigenvectors of the overlap matrix
C     - canonical transform matrix
C     - converged density matrix (alpha, beta for UKS)
C     - converged Fock matrix  (alpha, beta for UKS)
C     - transition dipole matrices
C
C     All quantities are in the atomic orbital basis.
C     
      subroutine dft_print_mats (g_s, svals, g_svecs,
     $     nmo, g_dens, g_focks)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "cdft.fh"
      
C     == Inputs ==
      integer g_s
      double precision svals(*)
      integer nmo(2)
      integer g_svecs(2), g_dens(2), g_focks(2)

C     == Parameters ==
      character(len=*), parameter :: fmt = "(i10, i10, 1e20.10)"

      
C     == Variables ==
      integer me
      integer ibf, jbf
      double precision val, val3(3)
      integer g_tmp_kl1, g_tmp_kl2, g_tmp_kl3

      me = ga_nodeid()

      if (.not. ga_duplicate(g_s, g_tmp_kl1, 'tmp1 matrix'))
     $     call errquit ("failed to create X", 0, GA_ERR)
      
      if (.not. ga_duplicate(g_s, g_tmp_kl2, 'tmp2 matrix'))
     $     call errquit ("failed to create X", 0, GA_ERR)
      
      if (.not. ga_duplicate(g_s, g_tmp_kl3, 'tmp3 matrix'))
     $     call errquit ("failed to create X", 0, GA_ERR)
      
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "Overlap matrix"
         write (luout, *) "--------------"
         do ibf = 1, nbf_ao
            do jbf = 1, nbf_ao
               call ga_get(g_s, ibf, ibf, jbf, jbf, val, 1)
               write (luout,fmt) ibf, jbf, val
            enddo
         enddo
         write (luout, *) ""
      endif
      
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "Eigenvalues of the overlap matrix"
         write (luout, *) "---------------------------------"
         do ibf = 1, nbf_ao
            write (luout,"(i0, 1e20.10)") ibf, svals(ibf)
         enddo
         write (luout, *) ""
      endif

      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "Eigenvectors of the overlap matrix"
         write (luout, *) "----------------------------------"
         do ibf = 1, nbf_ao
            do jbf = 1, nbf_ao
               call ga_get(g_svecs, ibf, ibf, jbf, jbf, val, 1)
               write (luout,fmt) ibf, jbf, val
            enddo
         enddo
         write (luout, *) ""
      endif

      call ga_zero (g_tmp_kl1)
      call dft_canorg (nmo(1), svals, g_svecs, g_tmp_kl1) !compute X matrix
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "Canonical transform matrix (X)"
         write (luout, *) "------------------------------"
         do ibf = 1, nbf_ao
            do jbf = 1, nbf_ao
               call ga_get(g_tmp_kl1, ibf, ibf, jbf, jbf, val, 1)
               write (luout,fmt) ibf, jbf, val
            enddo
         enddo
         write (luout, *) ""
      endif
      
      if (ipol.eq.1) then 
         if (me.eq.0) then
            write (luout, *) ""
            write (luout, *) "Ground state density matrix"
            write (luout, *) "---------------------------"
            do ibf = 1, nbf_ao
               do jbf = 1, nbf_ao
                  call ga_get(g_dens(1), ibf, ibf, jbf, jbf, val, 1)
                  write (luout,fmt) ibf, jbf, val
               enddo
            enddo
            write (luout, *) ""
         endif
         
         if (me.eq.0) then
            write (luout, *) ""
            write (luout, *) "Ground state Fock matrix"
            write (luout, *) "------------------------"
            do ibf = 1, nbf_ao
               do jbf = 1, nbf_ao
                  call ga_get(g_focks(1), ibf, ibf, jbf, jbf, val, 1)
                  write (luout,fmt) ibf, jbf, val
               enddo
            enddo
            write (luout, *) ""
         endif
         
      elseif (ipol.eq.2) then
         if (me.eq.0) then
            write (luout, *) ""
            write (luout, *) "Ground state alpha density matrix"
            write (luout, *) "---------------------------------"
            do ibf = 1, nbf_ao
               do jbf = 1, nbf_ao
                  call ga_get(g_dens(1), ibf, ibf, jbf, jbf, val, 1)
                  write (luout,fmt) ibf, jbf, val
               enddo
            enddo
            write (luout, *) ""
         endif
         
         if (me.eq.0) then
            write (luout, *) ""
            write (luout, *) "Ground state beta density matrix"
            write (luout, *) "--------------------------------"
            do ibf = 1, nbf_ao
               do jbf = 1, nbf_ao
                  call ga_get(g_dens(2), ibf, ibf, jbf, jbf, val, 1)
                  write (luout,fmt) ibf, jbf, val
               enddo
            enddo
            write (luout, *) ""
         endif
         
         if (me.eq.0) then
            write (luout, *) ""
            write (luout, *) "Ground state alpha Fock matrix"
            write (luout, *) "------------------------------"
            do ibf = 1, nbf_ao
               do jbf = 1, nbf_ao
                  call ga_get(g_focks(1), ibf, ibf, jbf, jbf, val, 1)
                  write (luout,fmt) ibf, jbf, val
               enddo
            enddo
            write (luout, *) ""
         endif
         
         if (me.eq.0) then
            write (luout, *) ""
            write (luout, *) "Ground state beta Fock matrix"
            write (luout, *) "-----------------------------"
            do ibf = 1, nbf_ao
               do jbf = 1, nbf_ao
                  call ga_get(g_focks(2), ibf, ibf, jbf, jbf, val, 1)
                  write (luout,fmt) ibf, jbf, val
               enddo
            enddo
            write (luout, *) ""
         endif
         
      else
         call errquit("Invalid ipol",0,0)
         
      endif
      
      
      call ga_zero (g_tmp_kl1)  ! dipx
      call ga_zero (g_tmp_kl2)  ! dipy
      call ga_zero (g_tmp_kl3)  ! dipz
      
      call int_dip_ga(ao_bas_han, ao_bas_han,
     $     g_tmp_kl1, g_tmp_kl2, g_tmp_kl3)
      
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "x Transition Dipole Matrix"
         write (luout, *) "--------------------------"
         do ibf = 1, nbf_ao
            do jbf = 1, nbf_ao
               call ga_get(g_tmp_kl1, ibf, ibf, jbf, jbf, val, 1)
               write (luout,fmt) ibf, jbf, val
            enddo
         enddo
         write (luout, *) ""
      endif
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "y Transition Dipole Matrix"
         write (luout, *) "--------------------------"
         do ibf = 1, nbf_ao
            do jbf = 1, nbf_ao
               call ga_get(g_tmp_kl2, ibf, ibf, jbf, jbf, val, 1)
               write (luout,fmt) ibf, jbf, val
            enddo
         enddo
         write (luout, *) ""
      endif
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "z Transition Dipole Matrix"
         write (luout, *) "--------------------------"
         do ibf = 1, nbf_ao
            do jbf = 1, nbf_ao
               call ga_get(g_tmp_kl3, ibf, ibf, jbf, jbf, val, 1)
               write (luout,fmt) ibf, jbf, val
            enddo
         enddo
         write (luout, *) ""
      endif

      if (.not. ga_destroy(g_tmp_kl1))
     $     call errquit ("destroy failed", 0, 0)
      if (.not. ga_destroy(g_tmp_kl2))
     $     call errquit ("destroy failed", 0, 0)
      if (.not. ga_destroy(g_tmp_kl3))
     $     call errquit ("destroy failed", 0, 0)
      
      end subroutine
