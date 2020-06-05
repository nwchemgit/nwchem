C
C     rt_tddft_header.F
C
C     Prints parameters to stdout.
C
      subroutine rt_tddft_header (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "matutils.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_header: "
      character(len=*), parameter :: fmtint = "(5x,a,i0)"
      character(len=*), parameter :: fmtint2vals = "(5x,a,i0,a,i0,a)"
      character(len=*), parameter :: fmtdbl = "(5x,a,es10.4)"
      character(len=*), parameter :: fmtdbl2vals =
     $     "(5x,a,es10.4,a,es10.4,a)"
      character(len=*), parameter :: fmtdbl3vals =
     $     "(5x,a,es11.4,a,es11.4,a,es11.4,a)"  !11 instead of 10 for neg sign
      character(len=*), parameter :: fmtcha = "(5x,a)"


C     == Variables ==
      integer me
      integer everyval
      character(500) popts
      logical opt_enabled
      double precision norm
      integer ifield
      integer iexcite
      type(rt_field_t) field
      type(rt_excite_t) excite
      character*32 gname, fname
      
      me = ga_nodeid ()

      if (me.eq.0) then
         
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "System parameters",20,.true.)
         write (luout, fmtint)
     $        "No. atoms              : ", params%natoms

         if (params%ipol.eq.1) then
            write (luout, fmtcha)
     $           "Wavefunction type      : restricted"
         elseif (params%ipol.eq.2) then
            write (luout, fmtcha)
     $           "Wavefunction type      : spin polarized"
         else
            call errquit (pname//"bad ipol value", 0, 0)
         endif

         write (luout, fmtint)
     $        "No. alpha electrons    : ", params%nalpha
         write (luout, fmtint)
     $        "No. beta electrons     : ", params%nbeta
         write (luout, fmtint)
     $        "Total no. electrons    : ", params%ntotel
         write (luout, fmtint)
     $        "System charge          : ", params%syscharge


C
C     2e integral info
C
      if (me.eq.0) then 
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout, "2e integrals",20,.true.)
         write (luout, fmtdbl) "tol2e for real part of Fock matrix: ",
     $        params%tol2e_re
         write (luout, fmtdbl) "tol2e for imag part of Fock matrix: ",
     $        params%tol2e_im
      endif


C     (moved to canorg)
C         write (luout, *) ""
C         write (luout, fmtint)
C     $        "No. atomic orbitals    : ", params%ns_ao
C         write (luout, fmtint)
C     $        "No. molecular orbitals : ", params%ns_mo

         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Propagation parameters",30,.true.)
C         write (luout, fmtdbl2vals) "Start time (tmin)     : ",
C     $        params%tmin, " au = ", params%tmin*au2fs, " fs"
         write (luout, fmtdbl2vals) "End time (tmax)       : ",
     $        params%tmax, " au = ", params%tmax*au2fs, " fs"
         write (luout, fmtdbl2vals) "Time step (dt)        : ",
     $        params%dt, " au = ", params%dt*au2fs, " fs"
         write (luout, fmtint)
     $        "No. time steps (nt)   : ", params%nt

         
         if (params%prop_method.eq.1) then
            write (luout, fmtcha)
     $           "Progagation method    : Euler"
         elseif (params%prop_method.eq.2) then
            write (luout, fmtcha)
     $           "Progagation method    : 4th order Runge-Kutta"
         elseif (params%prop_method.eq.3) then
            write (luout, fmtcha) "Progagation method    : "//
     $           "2nd order Magnus w/ self-consistent interpolation"
         else
            call errquit (pname//"unknown prop_method", 0, 0)
         endif

         
         if (params%exp_method.eq.1) then
            write (luout, fmtcha)
     $           "Exponentiation method : "//
     $           "Adaptive contractive power series"
         elseif (params%exp_method.eq.2) then
            write (luout, fmtcha)
     $           "Exponentiation method : "//
     $           "Diagonalization"
         elseif (params%exp_method.eq.3) then
            write (luout, fmtcha)
     $           "Exponentiation method : "//
     $           "Baker-Campbell-Hausdorff (BCH)"
            call errquit (pname//"BCH disabled",0,0)
         else
            call errquit (pname//"unknown exp_method", 0, 0)
         endif

         
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Tolerances",28,.true.)
         write(luout,fmtdbl) "Zero / comparison / check tolerance  : ",
     $        params%tol_zero
         write(luout,fmtdbl) "Series convergence tolerance         : ",
     $        params%tol_series
         write(luout,fmtdbl) "Interpolation convergence tolerance  : ",
     $        params%tol_interpol
         write(luout,fmtint) "Reqd no. zero terms in series        : ",
     $        params%terms_series
         write(luout,fmtint) "Reqd no. zero terms in interpolation : ",
     $        params%terms_interpol
         

         
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Output and checking parameters",30,.true.)


C     Echo the time-dependent quantites to be calculated and printed.
C     Note popts string has hardcoded length (shouldnt be a problem).
         popts = "Quantites to print : messages"
         opt_enabled = .false.

         
         if (params%print_dipole) then
            popts = trim(popts) // ", dipole"
            opt_enabled = .true.
         endif

         if (params%print_quad) then
            popts = trim(popts) // ", quadrupole"
            opt_enabled = .true.
         endif

         if (params%print_field) then
            popts = trim(popts) // ", field "
            opt_enabled = .true.
         endif

         if (params%print_moocc) then
            popts = trim(popts) // ", moocc"
            opt_enabled = .true.
         endif

         if (params%print_energy) then
            popts = trim(popts) // ", energy"
            opt_enabled = .true.
         endif

         if (params%print_cputime) then
            popts = trim(popts) // ", cputime"
            opt_enabled = .true.
         endif

         if (params%print_charge) then
            popts = trim(popts) // ", charge"
            opt_enabled = .true.
         endif

         if (params%print_s2) then
            popts = trim(popts) // ", <S^2>"
            opt_enabled = .true.
         endif

         if (.not. opt_enabled) then
            call errquit (pname//"no print options enabled",
     $           0          ,INPUT_ERR)
         else
            write (luout, fmtcha) trim (popts)
         endif


         if (params%nprints > 0) then
            everyval = params%nt / params%nprints
            write (luout,fmtint2vals) "No. print points   : ",
     $           params%nprints, " (once every ",
     $           everyval, " time steps)"
         else
            call errquit (pname//"must have a least 1 print point",0,0)
         endif

c$$$         if (params%nsnapshots > 0) then
c$$$            everyval = params%nt / params%nsnapshots
c$$$            write (luout,fmtint2vals) "No. dens snapshots : ",
c$$$     $           params%nsnapshots, " (once every ",
c$$$     $           everyval, " time steps)"
c$$$         else
c$$$            write (luout,fmtcha) "No. dens snapshots : None"
c$$$         endif

         if (params%nchecks > 0) then
            everyval = params%nt / params%nchecks
            write (luout, fmtint2vals) "No. checks         : " ,
     $           params%nchecks, " (once every ",
     $           everyval, " time steps)"
         else
            write (luout,fmtcha) "No. checks         : None"
         endif

         if (params%nrestarts > 0) then
            everyval = params%nt / params%nrestarts
            write (luout, fmtint2vals) "No. restart points : " ,
     $           params%nrestarts, " (once every ",
     $           everyval, " time steps)"
         else
            write (luout,fmtcha) "No. restart points : None"
         endif


C     Profiling and extra error checking.
         if (params%checklvl .eq. 1) then
            write (luout, fmtcha)
     $           "Checking level     : 1 (checkpoints only)"
         elseif (params%checklvl .eq. 2) then
            write (luout, fmtcha) "Checking level     : 2"//
     $           " (checkpoints; symmetries in each Fock build)"
         elseif (params%checklvl .eq. 3) then
            write (luout, fmtcha)
     $           "Checking level     : 3 (full--not implemented yet)"
         else
            call errquit (pname//"bad checklvl", 0, 0)
         endif

         if (params%prof) then
            write (luout, fmtcha) "Profiling          : Yes"
         else
            write (luout, fmtcha) "Profiling          : No"
         endif

         write (luout, *) ""
         write (luout, *) ""
         call util_flush (luout)


         call util_print_centered (luout,
     $        "Dipole and quadrupole matrices",15,.true.)

         write (luout, fmtdbl3vals)
     $        "Center of mass = (", params%origin(1), ", ",
     $        params%origin(2), ", ", params%origin(3), ")"
         
      endif !end proc 0


C
C     For the following, we compute norms on all procs and output on
C     proc 0.

C     (dipole mats)
      norm = mat_norm (params%g_dip%x)
      if (me.eq.0) write (luout, fmtdbl) "norm[D_x]  = ", norm

      norm = mat_norm (params%g_dip%y)
      if (me.eq.0) write (luout, fmtdbl) "norm[D_y]  = ", norm

      norm = mat_norm (params%g_dip%z)
      if (me.eq.0) write (luout, fmtdbl) "norm[D_z]  = ", norm


C     (quadrupole mats)         
      norm = mat_norm (params%g_quad%xx)
      if (me.eq.0) write (luout, fmtdbl) "norm[Q_xx] = ", norm

      norm = mat_norm (params%g_quad%xy)
      if (me.eq.0) write (luout, fmtdbl) "norm[Q_xy] = ", norm

      norm = mat_norm (params%g_quad%xz)
      if (me.eq.0) write (luout, fmtdbl) "norm[Q_xz] = ", norm

      norm = mat_norm (params%g_quad%yy)
      if (me.eq.0) write (luout, fmtdbl) "norm[Q_yy] = ", norm

      norm = mat_norm (params%g_quad%yz)
      if (me.eq.0) write (luout, fmtdbl) "norm[Q_yz] = ", norm

      norm = mat_norm (params%g_quad%zz)
      if (me.eq.0) write (luout, fmtdbl) "norm[Q_zz] = ", norm

      if (me.eq.0) then 
         write (luout, *) ""
         write (luout, *) ""
      endif


C
C     Print field information
C
      if (me.eq.0) then 
         call util_print_centered (luout, "Applied fields", 15, .true.)
      endif


      if (params%nfields .eq. 0) then
         if (me.eq.0) write (luout, fmtcha) "(none specified)"
      else
         do ifield = 1, params%nfields
            call rt_tddft_field_rtdb_get (params%rtdb, ifield, field)

C     XXX UPDATE THIS OUTPUT

            if (me.eq.0) then
               
               write (luout, *) ""
               
C               call util_print_centered (luout, trim(field%name),
C     $              15, .true.)
               write (luout, "(15x,a)") '"'//trim(field%name)//'"'
               write (luout, fmtcha)      "Type            : "//
     $              trim(field%type)
               write (luout, fmtcha)      "Polarization    : "//
     $              trim(field%polarization)
               write (luout, fmtdbl2vals) "Field maximum   : ",
     $              field%max, " au = ", field%max*au2vnm, " V/nm"
               
               if ((field%type.eq."cw").or.(field%type.eq."pulse")) then
                  write (luout, fmtdbl2vals) "Field frequency : ",
     $                 field%frequency, " au = ",
     $                 field%frequency*au2ev, " eV"
               endif

               if (field%type.eq."pulse") then
                  write (luout, fmtdbl2vals) "Pulse width     : ",
     $                 field%width, " au = ",
     $                 field%width*au2fs, " fs"
                  write (luout, fmtdbl2vals) "Pulse center    : ",
     $                 field%center, " au = ",
     $                 field%center*au2fs, " fs"
               endif
            endif               !me eq 0
         enddo
      endif

      
C
C     Print excitation rules
C
      if (me.eq.0) then 
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout, "Excitation rules",15,.true.)
      endif

      if (params%nexcites .eq. 0) then
         if (me.eq.0) write (luout, fmtcha)
     $        "(none specified--system will not be excited)"
      else
         do iexcite = 1, params%nexcites
            call rt_tddft_excite_rtdb_get (params%rtdb, iexcite, excite)

            gname = params%geom(excite%geom_indx)%name
            fname = params%field(excite%field_indx)%name

            if (me.eq.0) then
               write (luout, fmtcha) 'Excite geometry "'//trim(gname)//
     $              '" with "'//trim(fname)//'"'
            endif
         enddo
      endif
      
      if (me.eq.0) write (luout, *) ""
      if (me.eq.0) write (luout, *) ""

      end subroutine rt_tddft_header
         
      
c $Id$
