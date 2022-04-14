      Logical Function DPlot(rtdb)
*
* $Id$
*
      Implicit None
#include "errquit.fh"
*     
************************************************************************
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
************************************************************************
*     
*---- NWChem handles
      Integer rtdb, Geom, AO_Bas_Han
*---- Basis set variables
      Integer nBF_AO
*---- MO Vectors variables
      Logical  movecs_read_header, movecs_read, Int_Normalize
      logical ao_rdm_read
      External movecs_read_header, movecs_read, Int_Normalize
      external ao_rdm_read
c
      Character*256 File_Vec,File_Vec2
      Character*256 File_Mat1,File_Mat2
      character*20 scftype, scftype2
      Character*255 Title_Vec, Title_Vec2, Basis_Name, Basis_Name2
      Integer nBF,nBF2,nSet,nSet2,iSet,nSetTri,
     &        nMO(2),nMO2(2)
      Integer g_Vec(2),g_Dns(2),g_Dns2(2)
c
      Character*256 File_civecs
      Integer g_x(2),g_y(2)
      integer iroot
c
*---- Local variables
      Integer iProc,i,j,No_Of_Spacings(3),Lth_Tit,Out_Unit,iOff
      Integer nGrid,k_Grid,l_Grid,nPoint,nAtom,iAtom,Atom_Off
      Integer k_Tags,l_Tags,k_Charge,l_Charge
      Integer k_Vls,l_Vls,g_Vls,k_xyz,l_xyz,k_AMt,l_AMt,k_BMt,l_BMt
      Integer k_Dst,l_Dst,k_rq0,l_rq0
      Integer k_Ene,l_Ene,k_Occ,l_Occ
      Integer Lth_Vec,Lth_Out
      Integer nOrb,k_Orb,l_Orb
      Double Precision Ext_From_To(3,2),Scale,Volume
      Double Precision STD,Tot_Charge,maxelt
      Character Coord(3)
      Character*8 Units,What,Spin,Where
      Character*256 File_Out
      character*255 file_out_dir
      character*30 tag
      character*16 tag16
      Character*255 Title
      Logical LResult, DeBug
      logical lgaussian,dodos
      double precision tol_rho
      logical ltransden,ldiffden
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
*
************************************************************************
*     
      iProc = ga_nodeid()
      DeBug = util_print('dplot', print_high)
      Coord(1) = 'X'
      Coord(2) = 'Y'
      Coord(3) = 'Z'
*
*---- Push context down to DPlot
      If (.not.context_push('dplot'))
     &     Call ErrQuit('DPlot: context_push failed', 0, UNKNOWN_ERR)
*
***************************************
*     Read the data from the rtdb     *
***************************************
c
c --  Read output file format
      if (.not. rtdb_get(rtdb, 'dplot:lgaussian', mt_log, 1,
     &   lgaussian)) then
        call errquit('DPlot: rtdb_get failed - lgaussian: '
     &             //'No DPLOT input block?',11, INPUT_ERR)
      endif
c
c     --  Read file names of MO vectors
      If (.not.rtdb_cget(rtdb,'dplot:File_Vec',1,File_Vec))
     &   Call ErrQuit('DPlot: rtdb_get failed - File_Vec',0,
     &       RTDB_ERR)
c
      If (.not.rtdb_cget(rtdb,'dplot:File_Vec2',1,File_Vec2))
     1   call errquit('DPlot: rtdb_get failed - File_Vec2',0,
     1        RTDB_ERR)
c
c --  Read civecs file
      If (.not.rtdb_cget(rtdb,'dplot:File_civecs',1,File_civecs))
     1   call errquit('DPlot: rtdb_get failed - File_civecs',0,
     1        RTDB_ERR)
c
c --  Read iroot
      If (.not.rtdb_get(rtdb,'dplot:iroot',mt_int,1,iroot))
     &   Call ErrQuit('DPlot: rtdb_get failed - iroot',0, RTDB_ERR)
c
c --  Read file names of one particle reduced density matrices
      If (.not.rtdb_cget(rtdb,'dplot:File_Mat1',1,File_Mat1))
     1   call errquit('DPlot: rtdb_get failed - File_Mat1',0,
     1        RTDB_ERR)
c
      If (.not.rtdb_cget(rtdb,'dplot:File_Mat2',1,File_Mat2))
     1   call errquit('DPlot: rtdb_get failed - File_Mat2',0,
     1        RTDB_ERR)
c
c --  Read dodos switch
      if (.not. rtdb_get(rtdb, 'dplot:dodos', mt_log, 1,
     &   dodos)) call errquit('dpinput:rtdb_get failed',11, RTDB_ERR)
      if(dodos) then
         call dplot_dosinit(rtdb,file_vec)
         DPlot = .TRUE.
         return
      endif
c
c --  Read tol_rho 
      if (.not. rtdb_get(rtdb, 'dplot:tol_rho', mt_dbl, 1, tol_rho)) 
     &   call errquit('dpinput:rtdb_get failed',11, RTDB_ERR)
c
c     Transition density
      If (.not.rtdb_get(rtdb,'dplot:transden',mt_log,1,ltransden))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - ltransden',0,
     &       RTDB_ERR)
c
c     Difference density
      If (.not.rtdb_get(rtdb,'dplot:diffden',mt_log,1,ldiffden))
     &   Call ErrQuit('DPlot_Input: rtdb_put failed - ldiffden',0,
     &       RTDB_ERR)
c
c --  Read units
      If (.not.rtdb_cget(rtdb,'dplot:Units',1,Units))
     &   Call ErrQuit('DPlot: rtdb_get failed - Units',0, RTDB_ERR)
      If (Debug . and . iProc.eq.0)
     &   Write(LuOut,*)' Units                : ',Units
c
c --  Read plot range
      If (.not.rtdb_get(rtdb,'dplot:Ext_From_To',mt_dbl,6,Ext_From_To))
     &   Call ErrQuit('DPlot: rtdb_get failed - Ext_From_To',0,
     &       RTDB_ERR)
c
c --  Read number of spacings for building the grids
      If (.not.rtdb_get(rtdb,'dplot:No_Of_Spacings',mt_int,
     &   3,No_Of_Spacings))
     &   Call ErrQuit('DPlot: rtdb_get failed - No_Of_Spacings',0,
     &       RTDB_ERR)
c
c --  Read output filename
      If (.not.rtdb_cget(rtdb,'dplot:File_Out',1,File_Out))
     &   Call ErrQuit('DPlot: rtdb_get failed - File_Out',0, RTDB_ERR)
c
c --  Read Title
      If (.not.rtdb_cget(rtdb,'dplot:Title',1,Title))
     &   Call ErrQuit('DPlot: rtdb_get failed - Title',0, RTDB_ERR)
c
c --  Read What
      If (.not.rtdb_cget(rtdb,'dplot:What',1,What))
     &   Call ErrQuit('DPlot: rtdb_get failed - What',0, RTDB_ERR)
c
c --  Read Spin
      If (.not.rtdb_cget(rtdb,'dplot:Spin',1,Spin))
     &   Call ErrQuit('DPlot: rtdb_get failed - Spin',0, RTDB_ERR)
c
c --  Read nOrb
      If (.not.rtdb_get(rtdb,'dplot:nOrb',mt_int,1,nOrb))
     &   Call ErrQuit('DPlot: rtdb_get failed - nOrb',0, RTDB_ERR)
c
c --  Read Where
      If (.not.rtdb_cget(rtdb,'dplot:Where',1,Where))
     &   Call ErrQuit('DPlot: rtdb_get failed - Where',0, RTDB_ERR)
c
c --  Read the total charge of the molecule
      If (.not.rtdb_get(rtdb,'charge',mt_dbl,1,Tot_Charge))
     &   Tot_Charge = 0.d0
c
c -- If there is no specified file names of MO vectors, then looking for
c    the scf or dft MO vectors names in RTDB.
      if (inp_compare(.false.,File_Vec,' ')) then
        tag = 'scf:input  vectors'
c        if (iproc.eq.0) then
c         write(luout,*) ' trying scf '
c        end if
        if (.not. rtdb_cget(rtdb, tag, 1, File_vec)) then
c        if (iproc.eq.0) then
c           write(luout,*) ' trying dft '
c        end if
           tag = 'dft:input  vectors'
           if (.not. rtdb_cget(rtdb, tag, 1, File_vec))
     .          call util_file_name('movecs',.false.,.false.,File_Vec)
           if (iproc.eq.0) write(luout,*) ' File vec is ',File_Vec
        endif
      endif
c
c ---- Reading information from RTDB is completed
c
c --- Printing out
      If (iProc.eq.0) Then
         Write(LuOut,*)
         Write(LuOut,*)' Limits (a.u.) specified for the density plot:'
         Write(LuOut,*)' ---------------------------------------------'
         Write(LuOut,*)
         Write(LuOut,*)'       From        To      # of spacings'
      End If
      Scale = 1.8897265d0
      If (Units.eq.'ANGSTROM') Call dscal(6,Scale,Ext_From_To,1)
c
c --  Calculate the total number of the specified grid points
      nGrid = 1
      Do i = 1, 3
         nGrid = nGrid*(No_Of_Spacings(i)+1)
         If (iProc.eq.0) Write(LuOut,'(a,2(2x,f10.5),5x,i5)')Coord(i),
     &      (Ext_From_To(i,j),j = 1, 2), No_Of_Spacings(i)
      End Do
      If (iProc.eq.0) Then
         Write(LuOut,*)
         Write(LuOut,*)' Total number of grid points = ',nGrid
         Write(LuOut,*)
      End If
c
c --  print out filenames of MO vectors 
      Lth_Vec = Inp_StrLen(File_Vec)
      If (iProc.eq.0 .and. File_Vec .ne. ' ') 
     &   Write(LuOut,*)' 1-st set of MOs      : ',File_Vec(1:Lth_Vec)
c
      Lth_Vec = Inp_StrLen(File_Vec2)
      If (iProc.eq.0 . and . File_Vec2.ne.' ') Then
         Write(LuOut,*)' 2-nd set of MOs      : ',File_Vec2(1:Lth_Vec)
      End If
c
c --  print out civecs filename 
      Lth_Vec = Inp_StrLen(File_civecs)
      If (iProc.eq.0 .and. File_civecs .ne. ' ') 
     &   Write(LuOut,*)' CIVECS file     : ',File_civecs(1:Lth_Vec)
c
c --  Print out the filenames of one particle reduced density matrices
      Lth_Vec = Inp_StrLen(File_Mat1)
      if (iproc.eq.0 .and. File_Mat1 .ne. ' ') Then
         Write(LuOut,*) ' 1-st One Particle Reduced Density Matrix :',
     1   File_Mat1(1:Lth_Vec)
      endif
c
      Lth_Vec = Inp_StrLen(File_Mat2)
      if (iproc.eq.0 .and. File_Mat2 .ne. ' ') Then
         write(LuOut,*) ' 2-nd One Particle Reduced Density Matrix :',
     1   File_Mat2(1:Lth_Vec)
      endif
c
c --  print out the output filename
      Lth_Out = Inp_StrLen(File_Out)
      If (iProc.eq.0) 
     &   Write(LuOut,*)' Output is written to : ',File_Out(1:Lth_Out)
c
c --  print out the Title
      Lth_Tit = Inp_StrLen(Title)
      If (Debug . and . iProc.eq.0)
     &   Write(LuOut,*)' Title                : ',Title(1:Lth_Tit)
c
c --  print out type of picture and output file format
      If (iProc.eq.0) Then
         If (What.eq.'DENSITY') Then
            Write(LuOut,*)' Type of picture      : CHARGE DENSITY'
         Else If (What.eq.'TRANSDEN') Then
            Write(LuOut,*)' Type of picture      : TRANSITION DENSITY'
         Else If (What.eq.'VIEW') Then
            Write(LuOut,*)' Type of picture      : ORBITAL VIEW'
         Else If (What.eq.'CHARGES') Then
            Write(LuOut,*)' Type of picture      : POTENTIAL'
         End If
      End If
c
      if(lgaussian) then
       if (iProc.eq.0)
     &      Write(LuOut,*)' Format used          : Gaussian9x Cube'
      else
       if (iProc.eq.0)
     &      Write(LuOut,*)' Format used          : InsightII'
      endif
c
c --  Set up the right combination of Spin and What 
      If (What.eq.'VIEW' . and . Spin.ne.'BETA')          ! ***Set Up***
     &   Spin = 'ALPHA'
      If (iProc.eq.0 . and . What.ne.'CHARGES')
     &   Write(LuOut,*)' Spin                 : ',Spin
*
c --  debuging and checking
      If (DeBug . and . iProc.eq.0)
     &   Write(LuOut,*)' nOrb                 = ', nOrb
      If (File_Vec2.ne.' ' . and . What.eq.'VIEW')         ! ***Check***
     &   Call ErrQuit('DPlot: VIEW is meaningless for 2 orbital sets',0,
     &       UNKNOWN_ERR)
      If (nOrb.ne.999999) Then
         LResult = MA_Push_Get(MT_Int,nOrb,'Orb',l_Orb,k_Orb)
         Call ErrMem(LResult,'DPlot: Failed to allocate Orb.',0)
         If (.not.rtdb_get(rtdb,'dplot:OrbNo',mt_int,
     &      nOrb,Int_mb(k_Orb)))
     &      Call ErrQuit('DPlot: rtdb_get failed - OrbNo',0, RTDB_ERR)
         If (iProc.eq.0) Then
            If (What.eq.'DENSITY') Then
               Write(LuOut,'(a,50(1x,i3))')
     &          '  The density is computed using orbitals : ',
     &          (Int_mb(k_Orb-1+j),j=1,nOrb)
            Else
               If (nOrb.gt.1 .and . What.eq.'VIEW')        ! ***Check***
     &            Call ErrQuit('DPlot: VIEW and nOrb > 1 conflict',0,
     &                 UNKNOWN_ERR)
               Write(LuOut,'(a,i3,a)')
     &         '  The orbital ',Int_mb(k_Orb),' is plotted'
            End If
         End If
      Else
         If (iProc.eq.0 . and . What.ne.'CHARGES') then
            if (File_Mat1.eq.' ') then
               Write(LuOut,*)' The density is computed using all',
     1                       ' orbitals'
            else
               Lth_Vec = Inp_StrLen(File_Mat1)
               Write(LuOut,*) ' The density is computed using ',
     1                        ' density matrices'
            endif
         endif
      End If
c      If (iProc.eq.0 . and . 
c     1   (File_Vec2.ne.' '.or.File_Mat2.ne.' ')) Then
c         Write(LuOut,*)
c     &   ' The density is computed as Rho = Rho(1) - Rho(2)'
c      End If
*
c --  Set up the grid
      If (What.eq.'VIEW') Where = 'GRID'                  ! ***Set Up***
      If (iProc.eq.0 . and . What.eq.'DENSITY') Then
         If (Where.eq.'GRID') Then
            Write(LuOut,*)' The density is computed on the specified',
     &                    ' grid'
         Else If (Where.eq.'NUCLEI') Then
            Write(LuOut,*)' The density is computed at the position',
     &                    ' of nuclei'
         Else If (Where.eq.'G+N') Then
            Write(LuOut,*)' The density is computed on the grid and',
     &                    ' at the position of nuclei'
         End If
      End If
*
      If (nOrb.ne.999999 . and . What.eq.'CHARGES')        ! ***Check***
     &   Call ErrQuit('DPlot: Charges and Orbitals conflict',0,
     &       UNKNOWN_ERR)
      If (What.eq.'CHARGES' . and . Where.ne.'GRID')       ! ***Check***
     &   Call ErrQuit('DPlot: Charges can be combined with GRID only',
     & 0, 0)
      If (What.eq.'CHARGES' . and . File_Vec2.ne.' ')      ! ***Check***
     &   Call ErrQuit('DPlot: Specify one orb. file for opt. Charges',
     & 0, 0)
*
*********************************
*     Start the calculation     *
*********************************
*
*---- Open the output file
      Out_Unit = LuOut
      If (iProc.eq.0 . and . Where.ne.'NUCLEI' . and .
     &    File_Out.ne.'OUTPUT') Then
c     prefix directory only if no fowardslash character is present
         if(file_out(1:1).eq.'/') then
            write(0,*) ' found slash'
            file_out_dir=file_out
         else
            call util_file_name_noprefix(file_out,.false.,.false.,
     F           file_out_dir)
         endif
            Out_Unit=20
            Open(unit=Out_Unit,
     &           file=file_out_dir,
     &           status='unknown',
     &           access='sequential',
     &           form='formatted',
     &           err=111)
      End If
*
*---- Write the title and the format to the output file
*
*---- Get the geometry
      If (.not.Geom_Create(Geom,'geometry'))
     &     Call ErrQuit('DPlot: Geom_Create failed',0,0)
      If (.not.Geom_RTDB_Load(rtdb,Geom,'geometry'))
     &     Call ErrQuit('DPlot: Geom_RTDB_Load failed',0,0)
*
      LResult = Geom_Ncent(Geom,nAtom)
      If (Where.eq.'GRID') Then
         nPoint   = nGrid
         Atom_Off = -1
      Else If (Where.eq.'NUCLEI') Then
         nPoint   = nAtom
         Atom_Off = 0
      Else If (Where.eq.'G+N') Then
         nPoint   = nAtom + nGrid
         Atom_Off = nGrid
      End If
*
*---- Get the basis set
      If (.not.Bas_Create(AO_Bas_Han,'ao basis'))
     &   Call ErrQuit('DPlot: Bas_Create failed',0,0)
      If (.not.Bas_RTDB_Load(rtdb,Geom,AO_Bas_Han,'ao basis'))
     &   Call ErrQuit('DPlot: Bas_RTDB_Load failed',0,0)
      If (.not.Int_Normalize(rtdb,AO_Bas_Han))
     &   Call ErrQuit('DPlot: Int_Normalize failed',0,0)
      If (.not.Bas_NumBF(AO_bas_han,nBF_AO))
     &   Call ErrQuit('DPlot: Bas_NumBF failed',0,0)
      If (Debug . and . iProc.eq.0)
     &   Write(LuOut,*)' nBF_AO     = ',nBF_AO
*
*---- Allocate memory and get tags, charges and cartesian coordiantes
      LResult = MA_Push_Get(MT_Byte,nAtom*16,'Tags',l_Tags,k_Tags)
      Call ErrMem(LResult,'C_C: Failed to allocate Tags',0)
      LResult = MA_Push_Get(MT_Dbl,nAtom,'Charge',l_Charge,k_Charge)
      Call ErrMem(LResult,'C_C: Failed to allocate Charge',0)
      LResult = MA_Push_Get(MT_Dbl,nAtom*3,'xyz',l_xyz,k_xyz)
      Call ErrMem(LResult,'C_C: Failed to allocate xyz',0)
*
      If ( .not.Geom_Cart_Get(Geom,nAtom,Byte_mb(k_Tags),
     &   Dbl_mb(k_xyz),Dbl_mb(k_Charge)) )
     &   Call ErrQuit('C_C: Geom_Cart_Get failed',0,0)
*
*---- Allocate memory for grid points and get the grid
      LResult = MA_Push_Get(MT_Dbl,3*nPoint,'Grid',l_Grid,k_Grid)
      Call ErrMem(LResult,'DPlot: Failed to allocate grid',0)
*
      If (Where.ne.'NUCLEI') Then
         Call Get_Grid(Dbl_mb(k_Grid),nGrid,
     &                 Ext_From_To(1,1),Ext_From_To(1,2),
     &                 No_Of_Spacings,Volume,Out_Unit,iProc,Where,
     &                 nAtom,Dbl_mb(k_xyz),What,
     .        lgaussian,dbl_mb(k_charge))
         If (What.eq.'CHARGES') nPoint = nGrid
         Call GA_Sync
      End If
*
      If (Where.ne.'GRID') Then
         Call dcopy(3*nAtom,Dbl_mb(k_xyz),1,Dbl_mb(k_Grid+3*Atom_Off),1)
         Call GA_Sync
      End If
*
*---- Read the info from the first MO file
      call util_file_name_resolve(File_Vec,.false.)
      If (.not.movecs_read_header(File_Vec,Title_Vec,Basis_Name, 
     &     scftype, nBF,nSet,nMO,2))
     &   Call ErrQuit('DPlot: movecs_read_header failed',1,0)
      nSetTri = nSet*(nSet+1)/2
*---- Check, if the MO Vectors correspond to the specified basis set
      If (nBF.ne.nBF_AO)
     &   Call ErrQuit('DPlot: Basis set is not consistent with MOVec',
     &   0,0)
*
      If (Debug . and . iProc.eq.0) Then
         Write(LuOut,*)
         Write(LuOut,*)' Title_Vec  = ',Title_Vec(1:60)
         Write(LuOut,*)' Basis_Name = ',Basis_Name(1:60)
         Write(LuOut,*)' nBF        = ',nBF
         Write(LuOut,*)' nSet       = ',nSet
         Write(LuOut,*)' nMO        = ',nMO(1),nMO(2)
      End If
*
*---- Create global array for MOs and the density matrix
      g_Vec(1) = ga_create_atom_blocked(Geom, AO_Bas_Han, 'Vec1')
      If (nSet.eq.2)then
         g_Vec(2) = ga_create_atom_blocked(Geom, AO_Bas_Han, 'Vec2')
      endif
      If (What.ne.'VIEW') Then
         g_Dns(1) = ga_create_atom_blocked(Geom, AO_Bas_Han, 'Dns1')
         If (nSet.eq.2)then
             g_Dns(2) = ga_create_atom_blocked(Geom, AO_Bas_Han, 'Dns2')
         endif
      End If
c
      If (What.ne.'VIEW') Then
         if (nOrb.eq.999999 .and. File_Mat1.ne.' ') then
            call ga_zero(g_dns(1))
c           Read first density matrix from file
            call util_file_name_resolve(File_Mat1, .false.)
#if 0
            if (.not.ao_rdm_read(File_Mat1,ao_bas_han,g_Dns(1)))
     1           call errquit('DPlot: ao_rdm_read failed ',0,0)
#else
            call ao_1prdm_read(nbf,g_dns(1),file_mat1)
#endif
c
c           Read second density matrix from file
            if (File_Mat2.ne.' ') then
               If (iProc.eq.0) Write(LuOut,*)' The density is',
     1                         ' computed as Rho = Rho(1) - Rho(2)'
               g_Dns2(1) = ga_create_atom_blocked(Geom, AO_Bas_Han, 
     1                         'Dn1')
               call ga_zero(g_Dns2(1))
               call util_file_name_resolve(File_Mat2, .false.)
#if 0
               if(.not.ao_rdm_read(File_Mat2,ao_bas_han,g_Dns2(1)))
     1         call errquit('DPlot: ao_rdm_read failed',1,0)
#else
            call ao_1prdm_read(nbf,g_dns(2),file_mat2)
#endif
c
c              Compute the difference between the density matrices
               call GA_dAdd(1.d0,g_Dns(1),-1.d0,g_Dns2(1),g_Dns(1))
               if (.not. ga_destroy(g_Dns2(1))) call errquit
     &            ('dplot: could not destroy g_Dns2(1)', 0,0)
            endif
c
         else if (File_civecs.ne.' ') then
c
          LResult = MA_Push_Get(MT_Dbl,nBF,'Ene',l_Ene,k_Ene)
          Call ErrMem(LResult,'DPlot: Failed to allocate orb. en.',0)
          LResult = MA_Push_Get(MT_Dbl,nBF,'Occ',l_Occ,k_Occ)
          Call ErrMem(LResult,'DPlot: Failed to allocate occ. num.',0)
          Do iSet = 1, nSet
            call util_file_name_resolve(File_Vec,.false.)
            If (.not.movecs_read(File_Vec,iSet,Dbl_mb(k_Occ),
     &       Dbl_mb(k_Ene),g_Vec(iSet)))
     &       Call ErrQuit('DPlot: movecs_read failed',0,0)
          End Do
          LResult = MA_Pop_Stack(l_Occ)
          Call ErrMem(LResult,'DPlot: Failed to deallocate occ. num.',0)
          LResult = MA_Pop_Stack(l_Ene)
          Call ErrMem(LResult,'DPlot: Failed to deallocate orb. en.',0)
c
c         Construct the transition density or difference density
          call util_file_name_resolve(File_civecs,.false.)
          call int_init(rtdb,1,AO_Bas_Han)
          if (iproc.eq.0) write(luout,*) ' Root: ', iroot
          if (ltransden) then
           call get_transden(spin,iroot,nset,File_civecs,ao_bas_han,
     &        g_Vec,g_Dns)
          else if (ldiffden) then
           call get_diffden(spin,iroot,nset,File_civecs,ao_bas_han,
     &        g_Vec,g_Dns)
          else
           call errquit('dplot: please specify transden or diffden',0,0)
          end if
          call int_terminate()
         else
*------- Get the density matrix for the first orbital file
           Call Get_Dens(Geom,AO_Bas_Han,nSet,nBF,g_Vec,g_Dns,
     &                   Spin,File_Vec,nOrb,k_Orb,iProc)
c          Calculate the difference
           if (File_Vec2.ne.' ') then
              If (iProc.eq.0) Write(LuOut,*)' The density is',
     1                   ' computed as Rho = Rho(1) - Rho(2)'
              call util_file_name_resolve(File_Vec2,.false.)
              If (.not.movecs_read_header(File_Vec2,Title_Vec2,
     1           Basis_Name2,scftype2, nBF2,nSet2,nMO2,2))
     &        Call ErrQuit('DPlot: movecs_read_header failed',2,0)
*------- Check, if the the orbital files match
              If (nBF2.ne.nBF . or . nSet2.ne.nSet)
     &        Call ErrQuit('DPlot: Orbital files do not match',0,0)
*------- Create global array for density matrix
              g_Dns2(1) = ga_create_atom_blocked(Geom, AO_Bas_Han, 
     1                    'Dn1')
              If (nSet.eq.2)then
                 g_Dns2(2) = ga_create_atom_blocked(Geom, AO_Bas_Han, 
     1                       'Dn2')
              endif
*------- Get the density matrix for the second orbital file
              Call Get_Dens(Geom,AO_Bas_Han,nSet,nBF,g_Vec,g_Dns2,
     &                      Spin,File_Vec2,nOrb,k_Orb,iProc)
*------- Form the final density
              Call GA_dAdd(1.d0,g_Dns(1),-1.d0,g_Dns2(1),g_Dns(1))
*------- Deallocate memory
              If (nSet.eq.2)then
                 if (.not. ga_destroy(g_Dns2(2))) call errquit
     &              ('dplot: could not destroy g_Dns(2)', 0,0)
              endif
              if (.not. ga_destroy(g_Dns2(1))) call errquit
     &           ('dplot: could not destroy g_Dns(1)', 0,0)
           endif
         endif
      Else
*------- Read the orbital file
         LResult = MA_Push_Get(MT_Dbl,nBF,'Ene',l_Ene,k_Ene)
         Call ErrMem(LResult,'DPlot: Failed to allocate orb. en.',0)
         LResult = MA_Push_Get(MT_Dbl,nBF,'Occ',l_Occ,k_Occ)
         Call ErrMem(LResult,'DPlot: Failed to allocate occ. num.',0)
         Do iSet = 1, nSet
            call util_file_name_resolve(File_Vec,.false.)
            If (.not.movecs_read(File_Vec,iSet,Dbl_mb(k_Occ),
     &         Dbl_mb(k_Ene),g_Vec(iSet)))
     &         Call ErrQuit('DPlot: movecs_read failed',0,0)
         End Do
         LResult = MA_Pop_Stack(l_Occ)
         Call ErrMem(LResult,'DPlot: Failed to deallocate occ. num.',0)
         LResult = MA_Pop_Stack(l_Ene)
         Call ErrMem(LResult,'DPlot: Failed to deallocate orb. en.',0)
      End If
*
*---- Get the values of the density at the grid points
      LResult = MA_Push_Get(MT_Dbl,nPoint,'Vls',l_Vls,k_Vls)
      Call ErrMem(LResult,'DPlot: Failed to allocate values',0)
      call dcopy(npoint,0d0,0,dbl_mb(k_vls),1)
*ga:1:0
      If (What.ne.'CHARGES') Then
         Call Create_Contour(rtdb,Geom,AO_Bas_Han,nBF,nPoint,
     &                       nSet,nSetTri,g_Dns(1),g_Vec(1),
     &                       Dbl_mb(k_Grid),Dbl_mb(k_Vls),
     .                       no_of_spacings, tol_rho,
     &                       What,Spin,Int_mb(k_Orb))
         maxelt=0d0
         do i=1,npoint
            maxelt=max(maxelt,abs(dbl_mb(k_vls+i-1)))
         enddo
         if(iproc.eq.0) write(luout,*) ' max element ',maxelt
      Else
         if (.not. ga_create(mt_dbl, nPoint, 1, 'Vls', 
     &                       0, 0, g_Vls)) 
     &   call errquit('dplot: error creating g_Vls',0,0)
         LResult = MA_Push_Get(MT_Dbl,(nAtom+1)**2,'AMt',l_AMt,k_AMt)
         Call ErrMem(LResult,'C_C: Failed to allocate AMt',0)
         LResult = MA_Push_Get(MT_Dbl,nAtom+1,'BMt',l_BMt,k_BMt)
         Call ErrMem(LResult,'C_C: Failed to allocate BMt',0)
         LResult = MA_Push_Get(MT_Dbl,nAtom*nPoint,'Dst',l_Dst,k_Dst)
         Call ErrMem(LResult,'C_C: Failed to allocate Dst',0)
         LResult = MA_Push_Get(MT_Dbl,nAtom,'rq0',l_rq0,k_rq0)
         Call ErrMem(LResult,'C_C: Failed to allocate rq0',0)
*
         Call qdist(Dbl_mb(k_rq0),Dbl_mb(k_Dst),
     &              Dbl_mb(k_Grid),Dbl_mb(k_xyz),
     &              nPoint,nAtom)
         LResult = MA_Pop_Stack(l_rq0)
*
         Call Get_Charges(rtdb,Geom,AO_Bas_Han,nBF,nPoint,nAtom,
     &                   g_Dns(1),Dbl_mb(k_Grid),Dbl_mb(k_Vls),g_Vls,
     &                   Dbl_mb(k_xyz),Dbl_mb(k_Charge),Byte_mb(k_Tags),
     &                   Dbl_mb(k_AMt),Dbl_mb(k_BMt),Dbl_mb(k_Dst),
     &                   STD,Tot_Charge)
*
         LResult = MA_Pop_Stack(l_Dst)
         LResult = MA_Pop_Stack(l_BMt)
         LResult = MA_Pop_Stack(l_AMt)
         Call GA_Get(g_Vls,1,nPoint,1,1,Dbl_mb(k_Vls),nPoint)
         Call GA_Sync
         if (.not. ga_destroy(g_Vls)) call errquit
     &      ('dplot: could not destroy g_Vls', 0,0)
      End If
*
*---- Print out the values of the density matrix to the output
      call dplot_dump(lgaussian,out_unit,iproc,luout,
     ,     where,what,
     ,     title(1:Lth_Tit),
     ,     natom,dbl_mb(k_xyz),dbl_mb(k_charge),volume,tol_rho,
     ,     Ext_From_To(1,1),Ext_From_To(1,2),
     ,     no_of_spacings,ngrid,dbl_mb(k_vls),Byte_mb(k_Tags))
      If (iProc.eq.0 .and. (Where.ne.'GRID'.or.What.eq.'CHARGES')) Then
         Write(LuOut,*)
         If (What.eq.'CHARGES') Then
            Write(LuOut,*)'             Tag',
     &                    '         C o o r d i n a t e s',
     &                    '            Charges'
         Else
            Write(LuOut,*)'             Tag',
     &                    '         C o o r d i n a t e s',
     &                    '            Density'
         End If
         Write(LuOut,*)'             -------------------------------',
     &                 '-------------------------'
         Do iAtom = 1, nAtom
            iOff = (iAtom - 1)*3
            If (What.ne.'CHARGES')
     &           Dbl_mb(k_Charge-1+iAtom)=Dbl_mb(k_Vls+Atom_Off-1+iAtom)
            do j=1,16
               tag16(j:j)=Byte_mb(k_Tags+ (iAtom-1)*16+j-1)
            enddo
            Write(LuOut,'(A16,3f12.6,3x,f15.10)')
     T       tag16(1:inp_strlen(tag16)),
     &           (Dbl_mb(k_xyz+iOff+j),j=0,2),
     &           Dbl_mb(k_Charge-1+iAtom)
         End Do
         Write(LuOut,*)'             -------------------------------',
     &                 '-------------------------'
         If (What.eq.'CHARGES')
     &      Write(*,'(a,a,e8.2)')'              Standard deviation',
     &                         ' for the fitted charges = ',STD
         Write(LuOut,*)
      End If
*
*---- Deallocate memory
      LResult = MA_Pop_Stack(l_Vls)
      Call ErrMem(LResult,'DPlot: Failed to deallocate values',0)
      LResult = MA_Pop_Stack(l_Grid)
      Call ErrMem(LResult,'DPlot: Failed to deallocate grid',0)
      LResult = MA_Pop_Stack(l_xyz)
      Call ErrMem(LResult,'DPlot: Failed to deallocate xyz',0)
      LResult = MA_Pop_Stack(l_Charge)
      Call ErrMem(LResult,'DPlot: Failed to deallocate charge',0)
      LResult = MA_Pop_Stack(l_Tags)
      Call ErrMem(LResult,'DPlot: Failed to deallocate tags',0)
      If (nOrb.ne.999999) Then
         LResult = MA_Pop_Stack(l_Orb)
         Call ErrMem(LResult,'DPlot: Failed to deallocate Orb.',0)
      End If
*
*---- Destroy global arrays
c      Call GADest(g_Vls)
      If (What.eq.'DENSITY') Then
c         If (nSet.eq.2) Call GADest(g_Dns(2))
         If (nSet.eq.2)then
            if (.not. ga_destroy(g_Dns(2))) call errquit
     &         ('dplot: could not destroy g_Dns(2)', 0,0)
         endif
c         Call GADest(g_Dns(1))
         if (.not. ga_destroy(g_Dns(1))) call errquit
     &      ('dplot: could not destroy g_Dns(1)', 0,0)
      End If
c      If (nSet.eq.2) Call GADest(g_Vec(2))
      If (nSet.eq.2)then
         if (.not. ga_destroy(g_Vec(2))) call errquit
     &      ('dplot: could not destroy g_Vec(2)', 0,0)
      endif
c      Call GADest(g_Vec(1))
      if (.not. ga_destroy(g_Vec(1))) call errquit
     &   ('dplot: could not destroy g_Vec(1)', 0,0)
      If (.not.Geom_destroy(Geom))
     &     Call ErrQuit('DPlot: Geom_destory failed',0,0)
*
*---- Close the output file
      If (iProc.eq.0 . and . Where.ne.'NUCLEI' . and .
     &    File_Out.ne.'OUTPUT') Close(Out_Unit)
*
*---- Pop context
      If (.not.context_pop('dplot')) 
     &     Call ErrQuit('DPlot: context_pop failed',0,0)
*
*---- Normal termination
      DPlot = .TRUE.

      If (.not.Bas_destroy(AO_Bas_Han))
     &   Call ErrQuit('DPlot: Bas_destroy failed',0,0)
      Return
*
*---- Error exit
  111 Continue
      Call ErrQuit('Create_Contour: Failed to open output file',0,0)
      End
