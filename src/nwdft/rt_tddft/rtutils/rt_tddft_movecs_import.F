C====================================================================
C
C     Determine name of initial movecs file.  If one was not specified
C     in the input deck, it falls back to SCF output.
C
      subroutine rt_tddft_movecs_fname (params, rt_movecs_fname)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params

      
C     == Outputs ==
      character(len=*), intent(out) :: rt_movecs_fname


C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_movecs_fname: "



C
C     If we found a "vectors" directive while parsing the input we use
C     that, otherwise fall back on the output vectors "movecs_in" from
C     the DFT (via common block in cdft.fh header).
C
      if (.not. rtdb_cget (params%rtdb, "rt_tddft:init_movecs",
     $     1, rt_movecs_fname)) then

         if (movecs_in .eq. "atomic")
     $        call errquit ('Cannot use atomic guess for rt_tddft'//
     $        ' initial state--first run an SCF to generate'//
     $        ' valid movecs', 0, 0)

         call rt_tddft_print_warning ("Starting movecs not specified"//
     $        "--trying SCF output: "//char(10)//"   "//
     $        trim(movecs_out))

         rt_movecs_fname = movecs_in
      endif

      call util_file_name_resolve (rt_movecs_fname, .false.)

      end subroutine



C====================================================================
C
C     Print imported movecs information to stdout.
C
      subroutine rt_tddft_movecs_print_header (params, fname,
     $     title, basis_name, scftype, nbf, mo_nsets, nmo)

      implicit none

#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      character(len=*), intent(in)  :: fname, title, basis_name, scftype
      integer, intent(in)           :: nbf, mo_nsets, nmo(mo_nsets)


C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_movecs_print_header: "


C     == Variables ==
      integer me

      
      me = ga_nodeid ()

      if (me.eq.0) then
         write (luout, "(5x,a,a)")  "File name              : ",
     $        trim (fname)
         write (luout, "(5x,a,a)")  "Job title              : ",
     $        trim (title)
         write (luout, "(5x,a,a)")  "Basis set name         : ",
     $        trim (basis_name)
         write (luout, "(5x,a,a)")  "SCF type               : ",
     $        trim (scftype)
         write (luout, "(5x,a,i0)") "Atomic orbitals        : ",
     $        nbf

         if (mo_nsets .eq. 1) then
            write (luout, "(5x,a,i0)") "Molecular orbitals     : ",
     $           nmo
         elseif (mo_nsets .eq. 2) then
            write (luout, "(5x,a,i0,a,i0)") "Molecular orbitals     : ",
     $           nmo(1), ", ", nmo(2)
         endif
            
         call util_flush (luout)
      endif
      call ga_sync ()

      end subroutine


C====================================================================
C
C     Print imported movecs eigenvalues and occupations to stdout.
C
      subroutine rt_tddft_movecs_print_evals (params, occ, evals)

      implicit none

#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: occ(*), evals(*)


C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_movecs_print_evals: "


C     == Variables ==
      integer me
      integer iv

      
      me = ga_nodeid ()

C     
C     Print eigenvalue list
C     
      if (me.eq.0) then
            
         write (luout, *) ""
            
         call util_print_centered (luout, "  Vector      "//
     $        "Occupation     Eigenvalue [au]", 0, .true.)
         
         do iv = 1, params%ns_ao  
            if (iv .le. params%ns_mo) then ! vals not removed by canorg
               write (luout, "(2x,i8,7x,1f5.2,7x,es15.8)")
     $              iv, occ(iv), evals(iv)
            else ! lindep vecs
               write (luout, "(2x,i8,4x,a)")
     $              iv, "  (linearly dependent; removed)  "
            endif
         enddo
         
         write (luout, *) ""
         call util_flush (luout)
      endif
      call ga_sync ()

      end subroutine


C====================================================================
C
C     Convert real-valued ns_ao x ns_ao into complex ns_ao x ns_ao dens
C     mat in AO basis (used for importing from the ground state code).
C
C     P_uv = \sum_i n_i C_ui C_vi^*
C
C     where n_i is the ith orbital occupation, and C is the movecs
C     coefficient matrix.
C
C     We first build mask of orbital occupations (n3 = occupation of 3rd
C     MO).
C
C     [ n1   n2  ... nN ]
C     [ n1   n2      nN ]
C     [ ...  ...     ...]
C     [ n1   n2      nN ]
C
C     We apply this mask to the movecs matrix (g_movecs_pad) then take
C     matrix multiplication with the unmasked transpose (note: full
C     matrices are multiplied, since the mask has zeroed out the
C     non-occupied orbitals anyways).
C     
      subroutine rt_tddft_movecs_zdens(params, nsets, g_zdens_ao)

      implicit none

#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: ispin  !target movecs set    : 1 for CS/SO; 1 or 2 for OS
      integer, intent(in)           :: nsets  !total num movecs sets: 1 for CS/SO; 2 for OS


C     == Outputs ==
      integer, intent(in)           :: g_zdens_ao(nsets)
C      integer, intent(in)           :: g_movecs_ao_gs(nsets)


C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_movecs_zdens: "


C     == External ==
      logical, external  :: movecs_read_header, movecs_read
      logical, external  :: dmat_from_file

      
C     == Variables ==
      character*256 rt_movecs_fname
      integer g_movecs_pad          !double, ns_ao x ns_ao w/ 0's for canorg removed lindeps
      integer dtype, dim1, dim2
      integer me
      integer g_densre_ao

C     (movecs header stuff)
      character*256 mo_title, mo_basis_name, mo_scftype
      integer mo_nbf, mo_nsets, mo_ns_mo(nsets)

      integer locc, iocc, levals, ievals  !MO occupations and eigenvalues (from movecs file)

      integer g_tmp
      integer iorb
      integer is
      double precision occ
      character*255 rt_dmat_fname

      me = ga_nodeid ()


C
C     Checks
C
C      if ((ispin .lt. 1).or.(ispin .gt. nsets))
C     $     call errquit (pname//"invalid ispin supplied", ispin, 0)

      if ((nsets .ne. 1).and.(nsets .ne. 2))
     $     call errquit (pname//"invalid nsets supplied", nsets, 0)

      if (nsets .ne. 1)
     $     call errquit (pname//"nsets > 1 needs to be implemented",0,0)

      do is = 1, nsets
         call ga_check_handle (g_zdens_ao(is),
     $        "third argument of "//pname//"not a valid GA")
         call ga_inquire (g_zdens_ao(is), dtype, dim1, dim2)
      
         if (dtype .ne. mt_dcpl) call errquit (pname//
     $        "expecting complex-valued GA as third argument", 0, 0)
         if (dim1 .ne. dim2)
     $        call errquit (pname//"dim1 must equal dim2", 0, 0)
         if (dim1 .ne. params%ns_ao)
     $        call errquit (pname//
     $        "bad size P--expecting ns_ao x ns_ao", 0, 0)
      enddo


C
C     Allocation
C
      if (.not. ga_create(mt_dbl,params%ns_ao,params%ns_ao,
     $     "movecs", 0, 0, g_movecs_pad))
     $     call errquit ("couldnt create movecs_pad", 0, GA_ERR)

C
C     Just an alias to the ouput GA; note padded with 0's if have
C     lindep)
C     
CXXX  [KAL]: ONLY FOR nsets = 1
C      g_movecs_pad = g_movecs_ao_gs(1)


      if (.not. ga_duplicate(g_movecs_pad, g_densre_ao, "real P"))
     $     call errquit ("couldnt duplicate movecs", 0, GA_ERR)

      if (.not. ga_duplicate(g_movecs_pad, g_tmp, "tmp"))
     $     call errquit ("couldnt duplicate movecs", 0, GA_ERR)

c
c     If calculation starting from density, read dmat file
c     Else process movecs file
c
      if (params%use_dmat) then
        if (.not. rtdb_cget (params%rtdb, "rt_tddft:init_dmat1",
     $       1, rt_dmat_fname))
     $    call errquit ("Failed to read dmat file name from rtdb",
     $         0,RTDB_ERR)
        call util_file_name_resolve (rt_dmat_fname, .false.)
        if (.not.dmat_from_file(g_densre_ao,rt_dmat_fname))
     $    call errquit ("Failed to read dmat file",0,RTDB_ERR)
        call convert_d2z(1d0,g_densre_ao,0d0,g_densre_ao,g_zdens_ao)
c
      else

C
C     Read in header to get file name, check, then read in movecs.  Note
C     that the g_movecs_pad is ns_ao x ns_ao, with the last few columns
C     possibly 0 (if lindep), which is how the SCF code does it, which
C     is unlike my way, which has ns_mo x ns_mo.
C
C
        call ga_zero (g_movecs_pad)

        call rt_tddft_movecs_fname (params, rt_movecs_fname)

        if (.not. movecs_read_header (rt_movecs_fname, mo_title,
     $       mo_basis_name, mo_scftype, mo_nbf, mo_nsets,
     $       mo_ns_mo, nsets))
     $       call errquit (pname//"Failed to read movecs header", 0, 0)


C
C     Check that movecs are legit.
C
        if (mo_scftype .ne. "dft")
     $       call errquit (pname//
     $       'Initial movecs should have scftype "dft"', 0, 0)

        if (mo_nbf .ne. params%ns_ao)
     $       call errquit (pname//
     $       'Initial movecs wrong size: mo_nbf /= ns_ao', mo_nbf, 0)

        do is = 1, nsets
           if (mo_ns_mo(is) .ne. params%ns_mo)
     $          call errquit (pname//
     $          'Initial movecs wrong size: mo_ns_mo /= ns_mo',
     $          mo_ns_mo(is), 0)
        enddo
        
        if (mo_nsets .ne. nsets)
     $       call errquit (pname//"Wrong number of initial movecs,",
     $       mo_nsets, 0)

        
        call rt_tddft_movecs_print_header (params, rt_movecs_fname,
     $       mo_title, mo_basis_name, mo_scftype, mo_nbf,
     $       mo_nsets, mo_ns_mo)



C
C     Allocate buffers and movecs (ns_ao x ns_ao padded with zero *not* ns_mo x ns_mo).
C
        if (.not.ma_push_get(mt_dbl, params%ns_ao, 'occ', locc, iocc))
     &       call errquit(pname//'cannot allocate occ',0, MA_ERR)

        if (.not.ma_push_get(mt_dbl, params%ns_ao, 'evals',
     $       levals, ievals))
     &       call errquit(pname//'cannot allocate evals',0, MA_ERR)


C
C     Read in movecs (note ispin).
C
        if (.not. movecs_read (rt_movecs_fname, 1, dbl_mb(iocc),  !note 1!!!
     $       dbl_mb(ievals), g_movecs_pad))
     $       call errquit (pname//"Failed to read movecs data", 0, 0)

        call rt_tddft_movecs_print_evals (params,
     $       dbl_mb(iocc), dbl_mb(ievals))


C
C     Mask the movecs and multiply into the unmasked transposed movecs
C     to make the real-valued dens mat in AO basis, then cast to
C     complex.
C
        call ga_zero (g_tmp)

        do iorb = 1, params%ns_ao
           occ = dbl_mb(iocc + iorb - 1)
           call ga_fill_patch (g_tmp, 1, params%ns_ao,
     $          iorb, iorb, occ)
        enddo

        call ga_sync () !XXX needed?

        call ga_elem_multiply (g_tmp, g_movecs_pad, g_tmp) ! g_tmp now holds masked movecs
        call ga_zero (g_densre_ao)
        call ga_dgemm ("N", "T", params%ns_ao, params%ns_ao,
     $       params%ns_ao,1d0, g_tmp, g_movecs_pad, 0d0, g_densre_ao)

        call convert_d2z(1d0, g_densre_ao, 0d0, g_densre_ao, g_zdens_ao)

C        call ga_print (g_movecs_pad)
        if (.not. ma_pop_stack (levals))
     $       call errquit (pname//"Failed to pop evals", 0, ma_err)
        
        if (.not. ma_pop_stack (locc))
     $       call errquit (pname//"Failed to pop occ", 0, ma_err)
      end if !dmat

C     Clean up
C
      if (.not. ga_destroy (g_movecs_pad))
     $     call errquit (pname//"failed to destroy movecs",0,0)

      if (.not. ga_destroy (g_densre_ao))
     $     call errquit (pname//"failed to destroy densao_tmp",0,0)

      if (.not. ga_destroy (g_tmp))
     $     call errquit ("couldnt destroy tmp", 0, GA_ERR)

      end subroutine 





CXXX  [KAL]: REMOVE ONCE UNIFIED WITH CLOSEDSHELL ROUTINE
C====================================================================
C
C     rt_tddft_os_movecs_zdens.F
C
C     Read in initial state movecs from file and convert to OPEN SHELL
C     complex dens matrix in AO basis.  Although the output dens mat is
C     complex data type, it is pure real (as the movecs are from the SCF
C     and thus pure real).
C
C
      subroutine rt_tddft_os_movecs_zdens (params, g_zdens_ao)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"
#include "matutils.fh"

      
C     == Inputs ==
      type (rt_params_t), intent(in) :: params

      
C     == Outputs ==
      integer, intent(in) :: g_zdens_ao(2)     !dble complex, ns_ao x ns_ao; alpha, beta


C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_os_movecs_zdens: "

      
C     == External ==
      logical, external  :: movecs_read_header, movecs_read
      logical, external  :: dmat_from_file

      
C     == Variables ==
      character*256 rt_movecs_fname
      integer dtype, dim1, dim2
      integer me
      integer g_movecs_pad
      integer g_densao_tmp
      logical ok

C     (movecs header stuff)
      character*256 mo_title, mo_basis_name, mo_scftype
      integer mo_nbf, mo_nsets, mo_ns_mo(2) !2 for alpha, beta spins

      integer locc, iocc, levals, ievals  !MO occupations and eigenvalues (from movecs file)
      integer is
      character*255 rt_dmat_fname

      call rt_tddft_os_confirm(params,'rt_tddft_movecs_import.F')
      
      me = ga_nodeid ()
      
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Initial state: Imported open shell MO eigenvectors",
     $        20, .true.)
      endif
         

C
C     Check GAs
C
      call ga_check_handle (g_zdens_ao(1),
     $     "third argument of "//pname//"not a valid GA")
      call ga_inquire (g_zdens_ao(1), dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as third argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1 .ne. params%ns_ao)
     $     call errquit (pname//"bad size P--expecting ns_ao x ns_ao",
     $     0, 0)

      call ga_check_handle (g_zdens_ao(2),
     $     "third argument of "//pname//"not a valid GA")
      call ga_inquire (g_zdens_ao(2), dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as third argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1 .ne. params%ns_ao)
     $     call errquit (pname//"bad size P--expecting ns_ao x ns_ao",
     $     0, 0)


C
C     Allocation
C
      if (.not. ga_create(mt_dbl,params%ns_ao,params%ns_ao,
     $     "movecs", 0, 0, g_movecs_pad))
     $     call errquit ("couldnt create movecs_pad", 0, GA_ERR)

      if (.not. ga_duplicate(g_movecs_pad, g_densao_tmp, "real P"))
     $     call errquit ("couldnt duplicate movecs", 0, GA_ERR)
c
c     If calculation starting from density, read dmat file
c     Else process movecs file
c
      if (params%use_dmat) then
        if (.not. rtdb_cget (params%rtdb, "rt_tddft:init_dmat1",
     $       1, rt_dmat_fname))
     $    call errquit ("Failed to read alpha dmat file name from rtdb",
     $         0,RTDB_ERR)
        call util_file_name_resolve (rt_dmat_fname, .false.)
        if (.not.dmat_from_file(g_densao_tmp,rt_dmat_fname))
     $    call errquit ("Failed to read dmat file",0,RTDB_ERR)
        call convert_d2z(1d0,g_densao_tmp,0d0,g_densao_tmp,
     $                   g_zdens_ao(1))
        if (.not. rtdb_cget (params%rtdb, "rt_tddft:init_dmat2",
     $       1, rt_dmat_fname))
     $    call errquit ("Failed to read beta dmat file name from rtdb",
     $         0,RTDB_ERR)
        call util_file_name_resolve (rt_dmat_fname, .false.)
        if (.not.dmat_from_file(g_densao_tmp,rt_dmat_fname))
     $    call errquit ("Failed to read dmat file",0,RTDB_ERR)
        call convert_d2z(1d0,g_densao_tmp,0d0,g_densao_tmp,
     $                   g_zdens_ao(2))
c
      else

C
C     Read in header to get file name, check, then read in movecs.  Note
C     that the g_movecs_pad_pad is ns_ao x ns_ao, with the last few columns
C     possibly 0 (if lindep), which is how the SCF code does it, which
C     is unlike my way, which has ns_mo x ns_mo.  The 2 is for open shell.
C
        call rt_tddft_movecs_fname (params, rt_movecs_fname)


        if (.not. movecs_read_header (rt_movecs_fname, mo_title,
     $       mo_basis_name, mo_scftype, mo_nbf, mo_nsets,
     $       mo_ns_mo, 2))
     $       call errquit (pname//"Failed to read movecs header", 0, 0)

C        call rt_tddft_movecs_print_header (params, rt_movecs_fname,
C     $       mo_title, mo_basis_name, mo_scftype, mo_nbf)


C
C     Check that movecs are legit.
C
        if (mo_scftype .ne. "dft")
     $       call errquit (pname//
     $       'Initial movecs should have scftype "dft"', 0, 0)

        if (mo_nbf .ne. params%ns_ao)
     $       call errquit (pname//
     $       'Initial movecs wrong size: mo_nbf /= ns_ao', 0, 0)

        if (mo_ns_mo(1) .ne. params%ns_mo)
     $       call errquit (pname//
     $       'Initial movecs wrong size: mo_ns_mo(1) /= ns_mo', 0, 0)

        if (mo_ns_mo(2) .ne. params%ns_mo)
     $       call errquit (pname//
     $       'Initial movecs wrong size: mo_ns_mo(2) /= ns_mo', 0, 0)
        
C        call rt_tddft_print_warning ("Didnt check mo_ns_mo(2) == ns_mo")

        if (mo_nsets .ne. 2)
     $       call errquit (pname//"Wrong number of initial movecs,"//
     $       " mo_nsets should be 2 for open shell", 0, 0)


C
C     Allocate tmp buffers for occupations and evals (ns_ao x ns_ao
C     padded with zero *not* ns_mo x ns_mo).
C
        if (.not.ma_push_get(mt_dbl, params%ns_ao, 'occ', locc, iocc))
     &       call errquit(pname//'cannot allocate occ 1',0, MA_ERR)

        if (.not.ma_push_get(mt_dbl, params%ns_ao, 'evals',
     $       levals, ievals))
     &       call errquit(pname//'cannot allocate evals',0, MA_ERR)


C
C     Read in movecs.  Note we loop over the two spins since closed
C     shell, and we recycle all the MA and GA space, and finally store
C     in g_zdens_ao which has two components (alpha, beta spins).
C
        call ga_zero (g_zdens_ao(1))
        call ga_zero (g_zdens_ao(2))

C       (alpha part)
        call ga_zero (g_movecs_pad)
           
        if (.not. movecs_read (rt_movecs_fname, 1, dbl_mb(iocc),
     $       dbl_mb(ievals), g_movecs_pad))
     $       call errquit (pname//"Failed to read movecs data", 0, 0)
        
        call rt_tddft_movecs_print_evals (params,
     $       dbl_mb(iocc), dbl_mb(ievals))
        
        call ga_zero (g_densao_tmp)
        call ga_dgemm ("N", "T", params%ns_ao, params%ns_ao,
     $       params%nalpha, 1d0, g_movecs_pad, g_movecs_pad,
     $       0d0, g_densao_tmp)
        
        call convert_d2z(1d0, g_densao_tmp, 0d0, g_densao_tmp,
     $       g_zdens_ao(1))


C       (beta part)
        call ga_zero (g_movecs_pad)
           
        if (.not. movecs_read (rt_movecs_fname, 2, dbl_mb(iocc),
     $       dbl_mb(ievals), g_movecs_pad))
     $       call errquit (pname//"Failed to read movecs data", 0, 0)
        
        call rt_tddft_movecs_print_evals (params,
     $       dbl_mb(iocc), dbl_mb(ievals))
        
        
        call ga_zero (g_densao_tmp)
        call ga_dgemm ("N", "T", params%ns_ao, params%ns_ao,
     $       params%nbeta, 1d0, g_movecs_pad, g_movecs_pad,
     $       0d0, g_densao_tmp)
        
        call convert_d2z(1d0, g_densao_tmp, 0d0, g_densao_tmp,
     $       g_zdens_ao(2))
        
        if (.not. ma_pop_stack (levals))
     $       call errquit (pname//"Failed to pop evals", 0, ma_err)
        
        if (.not. ma_pop_stack (locc))
     $       call errquit (pname//"Failed to pop occ", 0, ma_err)

      end if ! dmat
C
C     Clean up
C
      if (.not. ga_destroy (g_movecs_pad))
     $     call errquit (pname//"failed to destroy movecs",0,0)

      if (.not. ga_destroy (g_densao_tmp))
     $     call errquit (pname//"failed to destroy densao_tmp",0,0)


      end subroutine
c $Id$
