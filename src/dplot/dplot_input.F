      SubRoutine DPlot_Input(rtdb)
*
* $Id$
*
      Implicit none
#include "errquit.fh"
*
************************************************************************
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
************************************************************************
*
      Integer rtdb
      Integer ind,indUnits,indWhat,indSpin,indWhere,i,j
      Integer Lth,Lth_Tit
      Integer No_Of_Spacings(3),ValueI,iroot
      Integer nOrb,iOrb,k_Orb,l_Orb
      Double Precision Ext_From_To(3,2),ValueR
*
      Integer   Num_Dirs,Num_Units,Num_What,Num_Spin,Num_Where
      Parameter (Num_Dirs  = 18)
      Parameter (Num_Units = 3)
      Parameter (Num_What  = 2)
      Parameter (Num_Spin  = 4)
      Parameter (Num_Where = 3)
*
      Character*255 test, Title
      Character*8   Dirs(Num_Dirs),Allowed_Units(Num_Units),
     &              Allowed_Spin(Num_Spin),
     &              Allowed_Where(Num_Where)
      Character*8  Allowed_What(Num_What)
      Character*8   Units,Spin,Where
      Character*12  What
      Character*3   Charges
!      character*30 tag
      Character*256  File_Vec,File_Vec2,File_Out
      Character*256  File_Mat1,File_Mat2
      Character*256  File_civecs
      character*32  theory
      logical lgaussian,dodos,ltransden,ldiffden
      double precision tol_rho
*
      Data Dirs / 'LimitXYZ','Vectors','civecs','root',
     .            'densmat','Output','Title','Orbitals',
     &            'Spin','Where','Charges','Gaussian','print',
     A            'dos','tol_rho','transden','diffden',
     .            'End' /
c
      Data Allowed_Units / 'angstrom','au','bohr' /
      Data Allowed_What  / 'density','view' /
      Data Allowed_Spin  / 'total','spindens','alpha','beta' /
      Data Allowed_Where / 'grid','nuclei','g+n' /
*
************************************************************************
*
*---- Set up defaults
      lgaussian = .false.
      Units     = 'angstrom'
      What      = 'density'
      Spin      = 'total'
      Where     = 'grid'
      dodos     =.false.
      iroot     = 1
      ltransden = .true.
      ldiffden = .false.
      tol_rho = 1d-40
c
c     try to get a scf movecs
c
!      tag = 'scf:input vectors'
!      if (.not. rtdb_cget(rtdb, tag, 1, File_vec)) then
!         tag = 'dft:input vectors'
!         if (.not. rtdb_cget(rtdb, tag, 1, File_vec)) 
!     .        File_Vec  = 'movecs'
!      endif
!     write(0,*) ' FIle vec is ',File_Vec
c     guess dft as theory  (is this really needed???)
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     theory='dft'
      if (.not. rtdb_cget(rtdb, theory(1:inp_strlen(theory))//
     A     ':output vectors', 1, file_vec)) file_vec = ' '

      File_Vec2 = ' '
      File_civecs = ' '
      File_Mat1 = ' '
      File_Mat2 = ' '
      File_Out  = 'dplot'
      Title     = 'Unknown Title'
      Lth_Tit   = 13
      Charges   = 'no '
      nOrb=999999
      Do i = 1, 3
         No_Of_Spacings(i) = 0
         Do j = 1, 2
            Ext_From_To(i,j) = 0.0d0
         End Do
      End Do
*
*---- Read the keyword
   10 Continue
      If (.not.inp_read())
     &   Call ErrQuit('DPlot_Input: inp_read failed',0,
     &       INPUT_ERR)
      If (.not.inp_a(test))
     &   Call ErrQuit('DPlot_Input: failed to read keyword',0,
     &       INPUT_ERR)
*
      If (.not.inp_match(Num_Dirs,.false.,test,Dirs,ind))
     &   Call ErrQuit('DPlot_Input: unknown directive',0,
     &       INPUT_ERR)
      GoTo (100, 200, 225, 230, 250, 300, 400, 500, 600, 700, 800, 
     &     900, 964, 1997, 1998, 1999, 2000, 9999) ind
      Call ErrQuit('DPlot_Input: unimplemented directive',ind,
     &       INPUT_ERR)
*
*---- Read the information about the limits of the grid
  100 Continue
      If (inp_a(test)) Then
         If (.not. inp_compare(.false.,test,'units'))
     &      Call ErrQuit('DPlot_Input: LimitXYZ [units <units>]',0,
     &       INPUT_ERR)
         If (.not.inp_a(Units))
     &      Call ErrQuit('DPlot_Input: LimitXYZ [units <units>]',0,
     &       INPUT_ERR)
         If (.not. inp_match(Num_Units,.false.,Units,Allowed_Units,
     &      indUnits)) Call ErrQuit('DPlot_Input: unknown units',0,
     &       INPUT_ERR)
      End If
      Do i = 1, 3
         If (.not. inp_read())
     &      Call ErrQuit('DPlot_Input: inp_read failed',0, INPUT_ERR)
         If (.not. inp_f(ValueR))
     &      Call ErrQuit('DPlot_Input: failed to read limits',0,
     &       INPUT_ERR)
         Ext_From_To(i,1) = ValueR
         If (.not. inp_f(ValueR))
     &      Call ErrQuit('DPlot_Input: failed to read limits',0,
     &       INPUT_ERR)
         Ext_From_To(i,2) = ValueR
         If (.not. inp_i(ValueI))
     &      Call ErrQuit('DPlot_Input: failed to read spacings',0,
     &       INPUT_ERR)
         No_Of_Spacings(i) = ValueI
      End Do
      GoTo 10
*
*---- Read the information about the source of the MO vectors
  200 Continue
      If (.not.inp_a(File_Vec))
     &   Call ErrQuit('DPlot_Input: Vectors <file1> [<file2>]',0,
     &       INPUT_ERR)
      If (inp_a(test)) Then
         Lth = Inp_StrLen(test)
         File_Vec2(1:Lth) = test(1:Lth)
      End If
      GoTo 10
*
*---- Read the information about the source of the civecs
  225 Continue
      If (.not.inp_a(File_civecs))
     &   Call ErrQuit('DPlot_Input: civecs <file> ',0,
     &       INPUT_ERR)
      If (inp_a(test)) Then
       Lth = Inp_StrLen(test)
       File_civecs(1:Lth) = test(1:Lth)
      End If
      GoTo 10
*
*---- Read the information about the number of root
  230 Continue
      If (.not. inp_i(iroot))
     &  Call ErrQuit('DPlot_Input: failed to read iroot',0,
     &     INPUT_ERR)
      GoTo 10
*
*---- Read the information about the source of the one particle reduced density matrices
  250 Continue
      if (.not.inp_a(File_Mat1)) then
         call util_file_name('densmat', .false.,.false.,file_mat1)
      endif
c     1   call errquit('dplot_input: densmat <file1> [<file2>]',0,
c     1        input_err)
      if (inp_a(test)) then
         lth = inp_strlen(test)
         File_Mat2(1:lth) = test(1:lth)
      endif
      goto 10
*
*---- Read the name of the output file
  300 Continue
      If (.not.inp_a(File_Out))
     &   Call ErrQuit('DPlot_Input: Output <file>',0, INPUT_ERR)
      GoTo 10
*
*---- Read the title of the output file
  400 Continue
      test = ' '
      Title = ' '
      If (.not.Inp_a(Title))
     &   Call ErrQuit('DPlot_Input: Title <title>',0, INPUT_ERR)
      GoTo 10
*
*---- Read the orbital info
  500 Continue
      If (inp_a(test)) Then
         If (.not. inp_match(Num_What,.false.,test,Allowed_What,
     &      indWhat)) Call ErrQuit('DPlot_Input: unknown option',0,
     &       INPUT_ERR)
         Lth = Inp_StrLen(test)
         What = ' '
         What(1:Lth) = test(1:Lth)
      End If
      If (.not. inp_read())
     &   Call ErrQuit('DPlot_Input: inp_read failed',0, INPUT_ERR)
      If (.not.inp_i(nOrb))
     &   Call ErrQuit('DPlot_Input: Orbitals; <# orbs.>; <Orb1 ...>',0,
     &       INPUT_ERR)
      If (nOrb.gt.0) Then
         If (nOrb.ne.999999) Then
            If (.not. inp_read())
     &         Call ErrQuit('DPlot_Input: inp_read failed',0,
     &       INPUT_ERR)
            If (.not.MA_Push_Get(MT_Int,nOrb,'Orb',l_Orb,k_Orb))
     &         Call ErrQuit('DPlot_Input: failed to allocate Orb',0,
     &       INPUT_ERR)
            Do iOrb = 1, nOrb
               If (.not. inp_i(ValueI))
     &            Call ErrQuit('DPlot_Input: failed to read orb. no.',
     &                         iOrb, INPUT_ERR)
               Int_mb(k_Orb - 1 + iOrb) = ValueI
            End Do
         End If
      Else
         Call ErrQuit('DPlot_Input: Sorry, nOrb <= 0',0, INPUT_ERR)
      End If
      GoTo 10
*
*---- Read the spin info
  600 Continue
      If (.not.inp_a(Spin))
     &   Call ErrQuit('DPlot_Input: Spin <Spin>',0, INPUT_ERR)
      If (.not. inp_match(Num_Spin,.false.,Spin,Allowed_Spin,indSpin))
     &   Call ErrQuit('DPlot_Input: unknown spin',0, INPUT_ERR)
      GoTo 10
*
*---- Read the where info
  700 Continue
      If (.not.inp_a(Where))
     &   Call ErrQuit('DPlot_Input: Where <Where>',0, INPUT_ERR)
      If (.not. inp_match(Num_Where,.false.,Where,Allowed_Where,
     &   indWhere)) Call ErrQuit('DPlot_Input: unknown where',0,
     &       INPUT_ERR)
      GoTo 10
*
*---- Read the charges info
  800 Continue
      Charges = 'yes'
      GoTo 10
*
 900  continue
      lgaussian=.true.
      goto 10
c
 964  continue
      goto 10
      call util_print_input(rtdb, 'dft')
c
 1997 continue
      dodos=.true.
      goto 10
c
 1998 continue
      If (.not. inp_f(tol_rho))
     &  Call ErrQuit('DPlot_Input: failed to read tol_rho',0,
     &     INPUT_ERR)
      tol_rho=max(1d-99,tol_rho)
      goto 10
c
 1999 continue
      ltransden = .true.
      goto 10
c
 2000 continue
      ldiffden = .true.
      goto 10
c
 9999 Continue
*
      If (Charges.eq.'yes') What = 'charges'
*
*---- Write limits to the rtdb
      Call Inp_UCase(Units)
      if (.not. rtdb_put(rtdb, 'dplot:lgaussian', mt_log, 1,
     &   lgaussian)) call errquit('dpinput:rtdbwr failed',11,0)
      If (.not.rtdb_cput(rtdb,'dplot:Units',1,Units))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - Units',0,
     &       RTDB_ERR)
      If (.not.rtdb_put(rtdb,'dplot:Ext_From_To',mt_dbl,6,Ext_From_To))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - Ext_From_To',0,
     &       RTDB_ERR)
      If (.not.rtdb_put(rtdb,'dplot:No_Of_Spacings',mt_int,
     &   3,No_Of_Spacings))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - No_Of_Spacings',0,
     &       RTDB_ERR)
*---- Write the name of the vector file to rtdb
      If (.not.rtdb_cput(rtdb,'dplot:File_Vec',1,File_Vec))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - File_Vec',0,
     &       RTDB_ERR)
c      write(0,*) ' vec1 and 2',File_vec,File_vec2
      If (.not.rtdb_cput(rtdb,'dplot:File_Vec2',1,File_Vec2))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - File_Vec2',0,
     &       RTDB_ERR)
*---- Write the name of the civecs file to rtdb
      If (.not.rtdb_cput(rtdb,'dplot:File_civecs',1,File_civecs))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - File_civecs',0,
     &       RTDB_ERR)
*---- Write the root information into the rtdb
      If (.not.rtdb_put(rtdb,'dplot:iroot',mt_int,1,iroot))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - iroot',0,
     &       RTDB_ERR)
c---- Write the names of the density matrix files to rtdb
      if (.not.rtdb_cput(rtdb,'dplot:File_Mat1',1,File_Mat1))
     1   call errquit('Dplot_Input: rtdb_cput failed - File_Mat1',0,
     1        RTDB_ERR)
      if (.not.rtdb_cput(rtdb,'dplot:File_Mat2',1,File_Mat2))
     1   call errquit('Dplot_Input: rtdb_cput failed - File_Mat2',0,
     1        RTDB_ERR)
*---- Write the name of the output file to rtdb
      If (.not.rtdb_cput(rtdb,'dplot:File_Out',1,File_Out))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - File_Out',0,
     &       RTDB_ERR)
*---- Write the title of the output file to rtdb
      If (.not.rtdb_cput(rtdb,'dplot:Title',1,Title))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - Title',0,
     &       RTDB_ERR)
*---- Write the orbital info
      Call Inp_UCase(What)
      If (.not.rtdb_cput(rtdb,'dplot:What',1,What))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - What',0,
     &       RTDB_ERR)
      If (.not.rtdb_put(rtdb,'dplot:nOrb',mt_int,1,nOrb))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - nOrb',0,
     &       RTDB_ERR)
      If (nOrb.ne.999999) Then
         If (.not.rtdb_put(rtdb,'dplot:OrbNo',mt_int,
     &      nOrb,Int_mb(k_Orb)))
     &      Call ErrQuit('DPlot_Input: rtdb_put failed - nOrb',0,
     &       RTDB_ERR)
         If (.not.MA_Pop_Stack(l_Orb))
     &      Call ErrQuit('DPlot_Input: failed to deallocate Orb',0,
     &       MA_ERR)
      End If
*---- Write the spin info
      Call Inp_UCase(Spin)
      If (.not.rtdb_cput(rtdb,'dplot:Spin',1,Spin))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - Spin',0,
     &       RTDB_ERR)
*---- Write 'where' info
      Call Inp_UCase(Where)
      If (.not.rtdb_cput(rtdb,'dplot:Where',1,Where))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - Where',0,
     &       RTDB_ERR)
*
      If (.not.rtdb_put(rtdb,'dplot:dodos',mt_log,
     &   1,dodos))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - dodos',0,
     &       RTDB_ERR)
*
      If (.not.rtdb_put(rtdb,'dplot:tol_rho',mt_dbl,
     &   1,tol_rho))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - tol_rho',0,
     &       RTDB_ERR)
*
      If (.not.rtdb_put(rtdb,'dplot:transden',mt_log,
     &   1,ltransden))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - ltransden',0,
     &       RTDB_ERR)
*
      If (.not.rtdb_put(rtdb,'dplot:diffden',mt_log,
     &   1,ldiffden))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - ldiffden',0,
     &       RTDB_ERR)

      End
