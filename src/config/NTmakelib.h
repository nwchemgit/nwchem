#
#	$Id: NTmakelib.h,v 1.1 1999-11-13 03:43:56 bjohnson Exp $
#

LIBRARY_PATH = $(LIB_DISTRIB)\$(LIBRARY)
OBJS = $(OBJ_OPTIMIZE) $(OBJ)
STAMP = WIN32.stamp

# !!! Use extra program to see if archive exists
!IF [ $(SRCDIR)\config\file_exists.exe $(LIBRARY_PATH) ]
LIBRARY_PATH_IF_EXISTS = $(LIBRARY_PATH)
!ELSE
LIBRARY_PATH_IF_EXISTS =
!ENDIF

$(LIBRARY_PATH): $(STAMP) $(OBJDIR) $(LIB_DISTRIB) $(OBJS) force_lib_update
	$(AR) @<<
	$(ARFLAGS) $(LIBRARY_PATH_IF_EXISTS) $(OBJS)
<<
!IFDEF SUBDIRS
	@nmake -nologo subdirs
!ENDIF

$(STAMP): $(HEADERS)
!IFDEF HEADERS
	!copy $** $(INCDIR)
!ENDIF
	@if exist $(STAMP) erase $(STAMP)
	@echo "" > $(STAMP)

"$(LIB_DISTRIB)" :
	@if not exist "$(LIB_DISTRIB)/$(NULL)" mkdir "$(LIB_DISTRIB)"

"$(OBJDIR)" :
	@if not exist "$(OBJDIR)/$(NULL)" mkdir "$(OBJDIR)"

subdirs: $(SUBDIRS)
	!cd $** & nmake -nologo

clean: $(SUBDIRS)
!IFDEF SUBDIRS
	!cd $** & nmake -nologo clean
!ENDIF
	@if exist $(STAMP) erase $(STAMP)
	@if exist "*.i" erase "*.i"
#	@if exist "*.exe" erase "*.exe"
#	@if exist "*.ilk" erase "*.ilk"
#	@if exist "*.pdb" erase "*.pdb"
	@if exist $(LIBRARY_PATH) erase $(LIBRARY_PATH)
	@if exist "$(OBJDIR)\*" erase /q "$(OBJDIR)\*"
	@if exist "$(OBJDIR)" rmdir "$(OBJDIR)"

force_lib_update:
