#include "macdecls.h"
#include "global.h"

#include <Python.h>
#include <import.h>
#include <graminit.h>
#include <pythonrun.h>
#include <stdlib.h>


#if defined(CRAY_T3E) || defined(CRAY_T3D)
int TASK_PYTHON(int *rtdb_ptr)
#else
int task_python_(int *rtdb_ptr)
#endif
{
   FILE *F;
   PyObject *phndl, *pmod, *pdict;
   char buf[20], *pbuf;
   char filename[256];
   int ret;
   
   Py_Initialize();		/* set the PYTHONPATH env   */
   initnwchem();
   PyRun_SimpleString("from nwchem import *");
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

   if (ga_nodeid_())
       sprintf(filename,"nwchem.py-%d",ga_nodeid_());
   else
       strcpy(filename,"nwchem.py");

   util_file_parallel_copy(filename, filename);

   if (!(F = fopen(filename, "r"))) {
       fprintf(stderr,"task_python: cannot open file %s\n",filename); 
       return 0;
   }
   ret = PyRun_SimpleFile(F, filename); 
   fclose(F);
   /*   unlink(filename); */

   Py_Finalize();

   return !ret;
}


