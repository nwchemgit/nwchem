#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <ga.h>
#include <mpi.h>
#include "armci.h"


char armci_util_char_getval(char *p) { return *p; }
short int armci_util_sint_getval(short int *p) { return *p; }
/* $Id$ */
