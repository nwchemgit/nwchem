# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../../config/makefile.h

simint_install:
	./build_simint.sh 

LIB_TARGETS += simint_install simint-chem-simint-generator.tar.gz simint.l* *-simint-generator-???????
  MAKEFLAGS = -j 2 --no-print-directory

include ../../../config/makelib.h

clean:
	@rm -rf simint_install simint-chem-simint-generator.tar.gz simint.l* *-simint-generator-??????? dependencies
