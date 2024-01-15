# $Id$
# check gitlab connectivity
GITLAB_UP := $(shell ./check_gitlab.sh 2> /dev/null)
ifndef GITLAB_UP
$(info cannot download from gitlab)
$(info stubbing nwc_columbus module)
OBJ = sifs_stubs.o
LIBRARY = libnwc_columbus.a
else
$(info gitlab connectivity OK)
SUBDIRS = sifs aoints
endif

include ../config/makefile.h
include ../config/makelib.h
