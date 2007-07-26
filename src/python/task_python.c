/*
 $Id: task_python.c,v 1.12 2007-07-26 21:05:03 d3p852 Exp $
*/
#include "macdecls.h"
#include "global.h"

#include <Python.h>
#include <import.h>
#include <graminit.h>
#include <pythonrun.h>
#include <stdlib.h>
#include "typesf2c.h"

extern void initnwchem();
extern void util_file_parallel_copy(const char *, const char *);


#if defined(CRAY_T3E) || defined(CRAY_T3D)  || defined(WIN32)
int FATR TASK_PYTHON(Integer *rtdb_ptr)
#else
int FATR task_python_(Integer *rtdb_ptr)
#endif
{
   FILE *F;
   char buf[20], *pbuf;
   char filename[256];
   int ret;
   
   Py_SetProgramName("NWChem");
   Py_Initialize();		/* set the PYTHONPATH env   */
   initnwchem();
   if (PyRun_SimpleString("from nwchem import *")) {
     fprintf(stderr,"import of NWCHEM failed\n");
     return 0;
   }
   pbuf = buf;			/* pass the rtdb_handle to  */
   sprintf(pbuf, "pass_handle(%d)\n", *rtdb_ptr); /* the python warping mod */
   if (PyRun_SimpleString(pbuf)) {
       fprintf(stderr,"task_python: failed to pass rtdb handle\n");
       return 0;
   }

   ret = 0;
   sprintf(pbuf, "INT     = %d", MT_F_INT);      
   ret += PyRun_SimpleString(pbuf);
   sprintf(pbuf, "DBL     = %d", MT_F_DBL);       
   ret += PyRun_SimpleString(pbuf);
   sprintf(pbuf, "CHAR    = %d", MT_CHAR);      
   ret += PyRun_SimpleString(pbuf);
   sprintf(pbuf, "LOGICAL = %d", MT_BASE + 11); 
   ret += PyRun_SimpleString(pbuf);
   sprintf(pbuf, "taskid = %d", ga_nodeid_());
   ret += PyRun_SimpleString(pbuf);
   sprintf(pbuf, "np = %d", ga_nnodes_());
   ret += PyRun_SimpleString(pbuf);

   if (ret) {
     fprintf(stderr,"setting RTDB types failed\n");
     return 0;
   }

   if (ga_nodeid_())
       sprintf(filename,"nwchem.py-%d",ga_nodeid_());
   else
       strcpy(filename,"nwchem.py");

   util_file_parallel_copy(filename, filename);

   /* PyRun_SimpleFile is unreliable on windows since you're passing
      a file pointer to Python which requires that it is compiled with
      a compatible compiler ... which it most likely is not */
 
#if defined(WIN32)
   ret = PyRun_SimpleString("execfile('nwchem.py')"); 
#else
   if (!(F = fopen(filename, "r"))) {
       fprintf(stderr,"task_python: cannot open file %s\n",filename); 
       return 0;
   }
   ret = PyRun_SimpleFile(F, filename); 
#endif


   /*fclose(F);*/
   /*   unlink(filename); */

   Py_Finalize();

   return !ret;
}


