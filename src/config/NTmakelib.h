#
#	$Id: NTmakelib.h,v 1.3 2000-02-08 22:05:42 bjohnson Exp $
#

LIBRARY_PATH = $(LIB_DISTRIB)\$(LIBRARY)
OBJS = $(OBJ_OPTIMIZE) $(OBJ)

# !!! Use extra program to see if archive exists
!IF [ $(SRCDIR)\config\file_exists.exe $(LIBRARY_PATH) ]
LIBRARY_PATH_IF_EXISTS = $(LIBRARY_PATH)
!ELSE
LIBRARY_PATH_IF_EXISTS =
!ENDIF

$(LIBRARY_PATH): $(OBJDIR) $(LIB_DISTRIB) $(OBJS)
	$(AR) @<<
	$(ARFLAGS) $(LIBRARY_PATH_IF_EXISTS) $(OBJS)
<<
!IFDEF SUBDIRS
	@nmake -nologo foreach_subdir
!ENDIF

WIN32.stamp: $(HEADERS)
!IFDEF HEADERS
	!copy $** $(INCDIR)
!ENDIF
	@if exist WIN32.stamp erase WIN32.stamp
	@echo "" > WIN32.stamp
!IFDEF SUBDIRS
	@nmake -nologo SUBDIR_TARGET=WIN32.stamp foreach_subdir
!ENDIF

"$(LIB_DISTRIB)" :
	@if not exist "$(LIB_DISTRIB)/$(NULL)" mkdir "$(LIB_DISTRIB)"

"$(OBJDIR)" :
	@if not exist "$(OBJDIR)/$(NULL)" mkdir "$(OBJDIR)"

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
