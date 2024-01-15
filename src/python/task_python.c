#include <Python.h>
#include "macdecls.h"
#include "ga.h"

#include <import.h>
#if ( PY_MAJOR_VERSION >= 3 && PY_MINOR_VERSION >= 9)
/* might require new headers and code */
#else
#include <graminit.h>
#endif
#include <pythonrun.h>
#include <stdlib.h>
#include "typesf2c.h"

#if PY_MAJOR_VERSION >= 3
extern PyMODINIT_FUNC PyInit_nwchem();
#else
extern void initnwchem();
#endif
extern void util_file_parallel_copy(const char *, const char *);


#if (defined(WIN32)) && !defined(__MINGW32__)
int FATR TASK_PYTHON(Integer *rtdb_ptr)
#else
int FATR task_python_(Integer *rtdb_ptr)
#endif
{
   FILE *F;
   char buf[20], *pbuf;
   char filename[256];
   int ret;
#if PY_VERSION_HEX >= 0x03080000
   PyStatus status;
   PyPreConfig preconfig;
#endif
   
#if PY_MAJOR_VERSION < 3
   Py_SetProgramName("NWChem");
#else
#if PY_MINOR_VERSION < 5
   wchar_t nwprogram[20];
   mbstowcs(nwprogram, "NWChem", strlen("NWChem") + 1);
#else
#if PY_VERSION_HEX < 0x03080000
   wchar_t *nwprogram = Py_DecodeLocale("NWChem", NULL);
#endif   
#endif
#if PY_VERSION_HEX < 0x03080000
   Py_SetProgramName(nwprogram);
#else
   /* Py_SetProgramName deprecated by python 3.11 */
   /* https://docs.python.org/3/c-api/init_config.html#preinitialize-python-with-pypreconfig */
   PyPreConfig_InitPythonConfig(&preconfig);
   status = Py_PreInitialize(&preconfig);
   if(PyStatus_Exception(status)) {
     fprintf(stderr,"python exception in Py_PreInitialize\n");
     return 0;
   }
#endif
#endif   
#if PY_MAJOR_VERSION >= 3
   PyImport_AppendInittab("nwchem", PyInit_nwchem);
#endif
   Py_Initialize();		/* set the PYTHONPATH env   */
#if PY_MAJOR_VERSION < 3
   initnwchem();
#endif
   if (PyRun_SimpleString("from nwchem import *")) {
     fprintf(stderr,"import of NWCHEM failed\n");
     return 0;
   }
   pbuf = buf;			/* pass the rtdb_handle to  */
   sprintf(pbuf, "pass_handle(%d)\n", (int) *rtdb_ptr); /* the python warping mod */
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
   sprintf(pbuf, "LOGICAL = %d", MT_F_LOG); 
   ret += PyRun_SimpleString(pbuf);
   sprintf(pbuf, "taskid = %d", GA_Nodeid());
   ret += PyRun_SimpleString(pbuf);
   sprintf(pbuf, "np = %d", GA_Nnodes());
   ret += PyRun_SimpleString(pbuf);

   if (ret) {
     fprintf(stderr,"setting RTDB types failed\n");
     return 0;
   }

   if (GA_Nodeid())
       sprintf(filename,"nwchem.py-%d",GA_Nodeid());
   else
       strcpy(filename,"nwchem.py");

   util_file_parallel_copy(filename, filename);

   /* PyRun_SimpleFile is unreliable on windows since you're passing
      a file pointer to Python which requires that it is compiled with
      a compatible compiler ... which it most likely is not */
 
#if defined(WIN32)
   ret = PyRun_SimpleString("exec(open('nwchem.py').read())");
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


