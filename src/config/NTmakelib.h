#
#	$Id$
#

LIBRARY_PATH = $(LIB_DISTRIB)\$(LIBRARY)
OBJS = $(OBJ_OPTIMIZE) $(OBJ)

!IF EXIST ($(CNFDIR)\win32\file_exists.exe)
# !!! Use extra program to see if archive exists
!IF [ $(CNFDIR)\win32\file_exists.exe $(LIBRARY_PATH) ]
LIBRARY_PATH_IF_EXISTS = $(LIBRARY_PATH)
!ELSE
LIBRARY_PATH_IF_EXISTS =
!ENDIF
!ELSE
LIBRARY_PATH_IF_EXISTS = path_error
!ENDIF

#
# Use pseudotarget for library since that is the only way to ensure
# in general that subdirs will be checked properly.
#
!IFNDEF FORCE_LIB_UPDATE
library: $(OBJDIR) $(LIB_DISTRIB) obj_opt obj_noopt $(LIBRARY_PATH)
!ELSE
library: $(OBJDIR) $(LIB_DISTRIB) obj_opt obj_noopt force_lib_update
!ENDIF
!IFDEF SUBDIRS
	@nmake -nologo foreach_subdir
!ENDIF

#
# Targets to handle separate compilation of objects specified with and
# without optimization
#
obj_opt: $(OBJ_OPTIMIZE)

obj_noopt:
!IFDEF OBJ
	@nmake -nologo NWDEBUG=1 obj_noopt_target
!ENDIF

obj_noopt_target: $(OBJ)

#
# The use of the $? macro below assumes that any object that needs to
# be placed into the library is newer than the library. This could fail
# in principle if the library file is touched (e.g. by another build
# updating the same library) after a needed object file has been built
# but before the following rule is executed. This seems very unlikely,
# but if it proves to be a problem simply replace $? with $(OBJS). The
# drawback to this is that it causes a "Replacing ..." message to be
# printed for every object file, not just the ones out of date.
#
$(LIBRARY_PATH): $(OBJS)
	$(AR) @<<
	$(ARFLAGS) $(LIBRARY_PATH_IF_EXISTS) $?
<<

#
# Force the AR command to run. This is useful for updating libraries
# such as util made from multiple dirs without having to delete and
# recompile the objects.
#
force_lib_update: $(OBJS)
	$(AR) @<<
	$(ARFLAGS) $(LIBRARY_PATH_IF_EXISTS) $(OBJS)
<<

WIN32.stamp: $(HEADERS)
!IFDEF HEADERS
	!copy $? $(INCDIR)
!ENDIF
	@if exist WIN32.stamp erase WIN32.stamp
	@echo "" > WIN32.stamp
!IFDEF SUBDIRS
	@nmake -nologo SUBDIR_TARGET=WIN32.stamp foreach_subdir
!ENDIF

"$(LIB_DISTRIB)" :
	@if not exist "$(TOPDIR)\lib\$(NULL)" mkdir "$(TOPDIR)\lib"
	@if not exist "$(LIB_DISTRIB)\$(NULL)" mkdir "$(LIB_DISTRIB)"

"$(BINDIR)" :
	@if not exist "$(TOPDIR)\bin\$(NULL)" mkdir "$(TOPDIR)\bin"
	@if not exist "$(BINDIR)\$(NULL)" mkdir "$(BINDIR)"

"$(OBJDIR)" :
	@if not exist "$(OBJDIR)\$(NULL)" mkdir "$(OBJDIR)"

clean: $(SUBDIRS)
!IFDEF SUBDIRS
	@nmake -nologo SUBDIR_TARGET=clean foreach_subdir
!ENDIF
	@if exist WIN32.stamp erase WIN32.stamp
	@if exist "*.i" erase "*.i"
#	@if exist "*.exe" erase "*.exe"
#	@if exist "*.ilk" erase "*.ilk"
#	@if exist "*.pdb" erase "*.pdb"
#	@if exist $(LIBRARY_PATH) erase $(LIBRARY_PATH)
	@if exist "$(OBJDIR)\*.obj" erase "$(OBJDIR)\*.obj"
	@if exist "$(OBJDIR)" rmdir "$(OBJDIR)"

#
# These targets allow automation of running nmake recursively in subdirs
# This is needed on Win98 since DOS syntax
#
#   cd $(TARGET_DIR) & nmake $(SUBDIR_TARGET)
#
# doesn't work because the '&' is not recognized
#
foreach_subdir: $(SUBDIRS)
	!@nmake -nologo TARGET_DIR=$** subdir_target

subdir_target:
	@echo Making $(SUBDIR_TARGET) in $(TARGET_DIR)
	@cd $(TARGET_DIR)
	@nmake -nologo $(SUBDIR_TARGET)

#
# Rules for automatic generation of MakeFile's from GNUmakefile's
# whenever possible
#

#
# The generation of MakeFile's in subdirs must be carried out in two
# steps.  We must ensure that the parent dir's MakeFile exists because
# it is used to determine its own subdirs.  Therefore the first target
# simply generates the MakeFile's for the current dir's subdirs before
# the second target steps down into the subdirs and fires up their
# Makefile's to generate MakeFile's in their subdirs recursively.
#
MakeFiles: $(SUBDIRS)
!IFDEF SUBDIRS
	!@nmake -nologo TARGET_DIR=$** MakeFile_1
	!@nmake -nologo TARGET_DIR=$** MakeFile_2
!ENDIF

MakeFile_1:
	@echo Generating MakeFile in $(TARGET_DIR)
	@cd $(TARGET_DIR)
	@nmake -nologo /F $(CNFDIR)\win32\MakeFile CNFDIR=$(CNFDIR) MakeFile

Makefile_2:
	@cd $(TARGET_DIR)
	@nmake -nologo MakeFiles

#
# To allow 'nmake MakeFile' to be invoked directly in an individual dir,
# such as after modifying the GNUmakefile by hand.
#
MakeFile: GNUmakefile
	@nmake -nologo /F $(CNFDIR)\win32\MakeFile CNFDIR=$(CNFDIR) MakeFile

#
# Clean the auto-generated MakeFile's - the dir tree must be traversed
# depth-first to do this!
#
clean_MakeFiles: $(SUBDIRS)
!IFDEF SUBDIRS
	!@nmake -nologo TARGET_DIR=$** clean_subdir_MakeFiles
!ENDIF
	@$(CNFDIR)\win32\cleanmake.exe $(LVL)

clean_subdir_MakeFiles:
	@cd $(TARGET_DIR)
	@nmake -nologo LVL="$(LVL)\$(TARGET_DIR)" clean_MakeFiles

#
# Make a Visual Studio project file
#
VS_project:
	@$(CNFDIR)\win32\convmake.exe $(TEMPLATE) $(TARGET_DIR)
