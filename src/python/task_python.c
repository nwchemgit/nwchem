#include <Python.h>
#include <import.h>
#include <graminit.h>
#include <pythonrun.h>
#include <stdlib.h>

#include "macdecls.h"


int task_python_(int *rtdb_ptr)
{
   FILE *F;
   PyObject *phndl, *pmod, *pdict;
   char buf[20], *pbuf;
   int ret;

   
   Py_Initialize();		/* set the PYTHONPATH env   */
   initnwchem();
   PyRun_SimpleString("from nwchem import *");
				/* var to the dir where the */ 
				/* py_rtdb_wrap.so is, when */
				/* dynamicaly linking       */
				/*  */
   pbuf = buf;			/* pass the rtdb_handle to  */
   sprintf(pbuf, "pass_handle(%d)\n", *rtdb_ptr); /* the python warping mod */
   if (PyRun_SimpleString(pbuf)) {
       fprintf(stderr,"task_python: failed to pass rtdb handle\n");
       return 0;
   }

   sprintf(pbuf, "INT     = %d", MT_F_INT);       PyRun_SimpleString(pbuf);
   sprintf(pbuf, "FLOAT   = %d", MT_FLOAT);     PyRun_SimpleString(pbuf);
   sprintf(pbuf, "DBL     = %d", MT_F_DBL);       PyRun_SimpleString(pbuf);
   sprintf(pbuf, "CHAR    = %d", MT_CHAR);      PyRun_SimpleString(pbuf);
   sprintf(pbuf, "LOGICAL = %d", MT_BASE + 11); PyRun_SimpleString(pbuf);

   if (!(F = fopen("nwchem.py", "r"))) {
      printf ("\nCannot open file nwchem.py\n"); 
      exit(0);
   }
   ret = PyRun_SimpleFile(F, "nwchem.py"); 
   fclose(F);

   return !ret;
}


