# Microsoft Developer Studio Project File - Name="nwchem" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=nwchem - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "nwchem.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "nwchem.mak" CFG="nwchem - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "nwchem - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "nwchem - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "nwchem - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /define:"WIN32" /define:COMPILATION_DATE='(unknown)' /define:NWCHEM_BRANCH='(unknown)' /define:COMPILATION_DIR='(unknown)' /extend_source:132 /fpp /include:"include" /include:"tools/include" /nodefine /nologo /warn:nofileopt /warn:nouncalled /fpp:"/c /m"
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 nwctask.lib ccsd.lib mcscf.lib selci.lib mp2.lib moints.lib stepper.lib driver.lib dftgrad.lib nwdft.lib gradients.lib cphf.lib esp.lib ddscf.lib guess.lib hessian.lib vib.lib util.lib rimp2.lib property.lib nwints.lib prepar.lib nwargos.lib nwmd.lib cafe.lib space.lib analyze.lib pfft.lib dplot.lib drdy.lib pario.lib global.lib ma.lib peigs.lib linalg.lib tcgmsg-mpi.lib armci.lib cvwmpi.lib wsock32.lib lapack.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib $(LINK_F90) /nologo /subsystem:console /machine:I386 /libpath:"$(NWCHEM_TOP_WIN32)/lib/win32 $(NWCHEM_TOP_WIN32)/src/tools/lib/win32" /libpath:"$(MPI_LIB)"
# Begin Custom Build
OutDir=.\Release
InputPath=.\Release\nwchem.exe
SOURCE="$(InputPath)"

"$(OutDir)/nwchem.pg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	echo local 0 > $(OutDir)/nwchem.pg

# End Custom Build

!ELSEIF  "$(CFG)" == "nwchem - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "nwchem___Win32_Debug0"
# PROP BASE Intermediate_Dir "nwchem___Win32_Debug0"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /define:"WIN32" /define:COMPILATION_DATE='(unknown)' /define:NWCHEM_BRANCH='(unknown)' /define:COMPILATION_DIR='(unknown)' /extend_source:132 /fpp /include:"include" /include:"tools/include" /nodefine /nologo /traceback /warn:argument_checking /warn:nofileopt /warn:nouncalled /fpp:"/c /m"
# SUBTRACT F90 /dbglibs
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /ML /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 nwctask.lib ccsd.lib mcscf.lib selci.lib mp2.lib moints.lib stepper.lib driver.lib dftgrad.lib nwdft.lib gradients.lib cphf.lib esp.lib ddscf.lib guess.lib hessian.lib vib.lib util.lib rimp2.lib property.lib nwints.lib prepar.lib nwargos.lib nwmd.lib cafe.lib space.lib analyze.lib pfft.lib dplot.lib drdy.lib pario.lib global.lib ma.lib peigs.lib linalg.lib tcgmsg-mpi.lib armci.lib cvwmpi.lib wsock32.lib lapack.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib $(LINK_F90) /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept /libpath:"$(NWCHEM_TOP_WIN32)/lib/win32_Debug" /libpath:"$(NWCHEM_TOP_WIN32)/src/tools/lib/win32" /libpath:"$(MPI_LIB)"
# SUBTRACT LINK32 /nodefaultlib /force
# Begin Custom Build
OutDir=.\Debug
InputPath=.\Debug\nwchem.exe
SOURCE="$(InputPath)"

"$(OutDir)/nwchem.pg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	echo local 0 > $(OutDir)/nwchem.pg

# End Custom Build

!ENDIF 

# Begin Target

# Name "nwchem - Win32 Release"
# Name "nwchem - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\nwchem.F
DEP_F90_NWCHE=\
	".\include\bgj_common.fh"\
	".\include\inp.fh"\
	".\include\printlevels.fh"\
	".\include\pstat.fh"\
	".\include\pstat_consts.fh"\
	".\include\rtdb.fh"\
	".\include\stdio.fh"\
	".\include\util.fh"\
	".\tools\include\global.fh"\
	".\tools\include\macommon.h"\
	".\tools\include\mafdecls.fh"\
	".\tools\include\tcgmsg.fh"\
	
# End Source File
# Begin Source File

SOURCE=.\stubs_win32.f
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\tools\lib\win32\armci.lib
# End Source File
# Begin Source File

SOURCE=.\tools\lib\win32\global.lib
# End Source File
# Begin Source File

SOURCE=.\tools\lib\win32\linalg.lib
# End Source File
# Begin Source File

SOURCE=.\tools\lib\win32\ma.lib
# End Source File
# Begin Source File

SOURCE=.\tools\lib\win32\pario.lib
# End Source File
# Begin Source File

SOURCE=".\tools\lib\win32\tcgmsg-mpi.lib"
# End Source File
# End Target
# End Project
