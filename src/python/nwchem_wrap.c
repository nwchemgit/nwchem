/*
 $Id: nwchem_wrap.c,v 1.12 2000-07-07 01:17:51 edo Exp $
*/
#if defined(DECOSF)
#include <alpha/varargs.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <Python.h>
#include <abstract.h>

#include "rtdb.h"
#include "macdecls.h"
#include "global.h"


static PyObject *NwchemError;

static int rtdb_handle;            /* handle to the rtdb */

#if defined(CRAY) || defined(CRAY_T3E)
#define task_energy_ TASK_ENERGY
#define task_gradient_ TASK_GRADIENT
#define task_optimize_ TASK_OPTIMIZE
#endif

extern int nw_inp_from_string(int, const char *);

extern int task_energy_(const int *);
extern int task_gradient_(const int *);
extern int task_optimize_(const int *);

static PyObject *
wrap_rtdb_open(PyObject *self, PyObject *args)
{
   const char *filename, *mode;

   if (PyArg_Parse(args, "(ss)", &filename, &mode)) {
       if (!rtdb_open(filename, mode, &rtdb_handle)) {
	   PyErr_SetString(NwchemError, "rtdb_open failed");
	   return NULL;
       }
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: rtdb_open(filename, mode)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}

static PyObject *wrap_rtdb_close(PyObject *self, PyObject *args)
{
   const char *mode;
   int  result;

   if (PyArg_Parse(args, "s", &mode)) {
       if (!(result = rtdb_close(rtdb_handle, mode))) {
	   PyErr_SetString(NwchemError, "rtdb_close failed");
	   return NULL;
       }
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: rtdb_close(mode)");
       return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}

static PyObject *wrap_pass_handle(PyObject *self, PyObject *args)
{
   if (!(PyArg_Parse(args, "i", &rtdb_handle))) {
      PyErr_SetString(PyExc_TypeError, "Usage: pass_handle(rtdb_handle)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}

static PyObject *wrap_rtdb_print(PyObject *self, PyObject *args)
{
   int flag;

   if (PyArg_Parse(args, "i", &flag)) {
      if (!rtdb_print(rtdb_handle, flag)) 
	   PyErr_SetString(NwchemError, "rtdb_print failed");
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: rtdb_print(flag)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}


static PyObject *wrap_rtdb_put(PyObject *self, PyObject *args)
{
    int i, list, list_len;
    int ma_type = -1;
    char *name;
    int* int_array;
    double *dbl_array;
    char *char_array;
    char cbuf[8192], *ptr;
    void *array = 0;
    PyObject *obj, *option_obj;

    if ((PyTuple_Size(args) == 2) || (PyTuple_Size(args) == 3)) {
      obj = PyTuple_GetItem(args, 0);      /* get var name */
      PyArg_Parse(obj, "s", &name);
      obj = PyTuple_GetItem(args, 1);      /* get an array or single value */
      
      if (PyList_Check(obj)) 
	list = 1; 
      else 
	list = 0;
      
      if (list) {
	list_len = PyList_Size(obj);
	if (   PyInt_Check(PyList_GetItem(obj, 0)))  ma_type = MT_F_INT;
	if ( PyFloat_Check(PyList_GetItem(obj, 0)))  ma_type = MT_F_DBL;
	if (PyString_Check(PyList_GetItem(obj, 0)))  ma_type = MT_CHAR;
      } else {
	list_len = 1;
	if (   PyInt_Check(obj))  ma_type = MT_F_INT;
	if ( PyFloat_Check(obj))  ma_type = MT_F_DBL;
	if (PyString_Check(obj))  ma_type = MT_CHAR; 
      } 
      
      if (PyTuple_Size(args) == 3) {
	option_obj = PyTuple_GetItem(args, 2);      /* get optional type */
	if (!(PyArg_Parse(option_obj, "i", &ma_type))) {
	  PyErr_SetString(PyExc_TypeError, 
			  "Usage: rtdb_put(value or values,[optional type])");
	  return NULL;
	}
      }
      
      if (ma_type != MT_CHAR) {
	if (!(array = malloc(MA_sizeof(ma_type, list_len, MT_CHAR)))) {
	  PyErr_SetString(PyExc_MemoryError,
			  "rtdb_put failed allocating work array");
	  return NULL;
	}
      }
      
      switch (ma_type) {
      case MT_INT:
      case MT_F_INT:  
      case MT_BASE + 11:	/* Logical */
	int_array = array;
	for (i = 0; i < list_len; i++) {
	  if (list) 
	    PyArg_Parse(PyList_GetItem(obj, i), "i", int_array+i);
	  else 
	    PyArg_Parse(obj, "i", int_array+i);
	}
	break;
	
      case MT_DBL:  
      case MT_F_DBL:
	dbl_array = array;
	for (i = 0; i < list_len; i++) {
	  if (list) 
	    PyArg_Parse(PyList_GetItem(obj, i), "d", dbl_array+i);
	  else 
	    PyArg_Parse(obj, "d", dbl_array+i);
	}
	break;
	
      case MT_CHAR: 
	ptr = cbuf;
	*ptr = 0;
	for (i = 0; i < list_len; i++) {
	  if (list) 
	    PyArg_Parse(PyList_GetItem(obj, i), "s", &char_array); 
	  else 
	    PyArg_Parse(obj, "s", &char_array); 
	  printf("PROCESSED '%s'\n", char_array);
	  if ((ptr+strlen(char_array)) >= (cbuf+sizeof(cbuf))) {
	     PyErr_SetString(PyExc_MemoryError,"rtdb_put too many strings");
	     return NULL;
	   }
	  strcpy(ptr,char_array);
	  ptr = ptr+strlen(char_array);
	  strcpy(ptr,"\n");
	  ptr = ptr + 1;
	}		 
	list_len = strlen(cbuf) + 1;
	array = cbuf;
	break;
	
      default:
	PyErr_SetString(NwchemError, "rtdb_put: ma_type is incorrect");
	if (array) free(array);
	return NULL;
	break;
      }                
      
      if (!(rtdb_put(rtdb_handle, name, ma_type, list_len, array))) {
	PyErr_SetString(NwchemError, "rtdb_put failed");
	if ((ma_type != MT_CHAR) && array) free(array);
	return NULL;
      }
      
    } else {
      PyErr_SetString(PyExc_TypeError, 
		      "Usage: rtdb_put(value or values,[optional type])");
      if ((ma_type != MT_CHAR) && array) free(array);
      return NULL;
    }
    Py_INCREF(Py_None);
    if ((ma_type != MT_CHAR) && array) free(array);
    return Py_None;
}

PyObject *wrap_rtdb_get(PyObject *self, PyObject *args)
{
   int i;
   int nelem, ma_type;
   char *name;
#define MAXPTRS 2048
   char *ptrs[MAXPTRS];
   char *format_str=0, format_char;
   PyObject *returnObj = 0;
   void *array=0;
   int ma_handle, ind;

   if (PyArg_Parse(args, "s", &name)) {
       if (!rtdb_ma_get(rtdb_handle, name, &ma_type, &nelem, &ma_handle)) {
	   PyErr_SetString(NwchemError, "rtdb_ma_get failed");
	   return NULL;
       }
       if (!MA_get_pointer(ma_handle, &array)) {
	   PyErr_SetString(NwchemError, "rtdb_ma_get failed");
	   return NULL;
       }
       /*printf("name=%s ma_type=%d nelem=%d ptr=%x\n",name, ma_type, 
	      nelem, array);*/

       switch (ma_type) {
       case MT_F_INT:
       case MT_INT  : 
       case MT_BASE + 11  : 
	   format_char = 'i'; break;
       case MT_F_DBL: 
       case MT_DBL  : 
	   format_char = 'd'; break;
	   break;
       case MT_CHAR : 
	   format_char = 's'; break;
       default:
	   PyErr_SetString(NwchemError, "rtdb_get: ma type incorrect");
	   (void) MA_free_heap(ma_handle);
	   return NULL;
	   break;
       }
       
       if (!(format_str = malloc(nelem+3))) {
	   PyErr_SetString(PyExc_MemoryError,
			   "rtdb_get failed allocating format string");
	   (void) MA_free_heap(ma_handle);
	   return NULL;
       }

       /* For character string need to build an array of pointers */

       if (ma_type == MT_CHAR) {
	 char *ptr, *next;
	 nelem = 0;
	 next = ptr = array;
         while (1) {		/* Replace strtok to handle consectutive separators */
	   int eos = (*ptr == 0);
	   if (*ptr == '\n') *ptr = 0;
           if (*ptr == 0) {
	     if (nelem >= MAXPTRS) {
	       PyErr_SetString(PyExc_MemoryError,"rtdb_get too many strings");
	       (void) MA_free_heap(ma_handle);
	       return NULL;
	     }
	     ptrs[nelem] = next;
	     nelem++;
	     if (!eos) next = ptr+1;
	   }
	   if (eos) break;
           ptr++;
	 }
           
	 /*
	 for (next=strtok((char *)array, "\n");
	      next;
	      next=strtok((char *) 0, "\n")) {
	   if (nelem >= MAXPTRS) {
	     PyErr_SetString(PyExc_MemoryError,"rtdb_get too many strings");
	     (void) MA_free_heap(ma_handle);
	     return NULL;
	   }
	   ptrs[nelem] = next;
	   nelem++;
	 }
	 */
       }

       ind = 0;
       if (nelem > 1) format_str[ind++] = '[';
       for (i = 0; i < nelem; i++, ind++)
	   format_str[ind] = format_char;
       if (nelem > 1) format_str[ind++] = ']';
       format_str[ind] = 0;
       
       if (ma_type == MT_CHAR)
	 returnObj = Py_VaBuildValue(format_str, ptrs);
       else
	 returnObj = Py_VaBuildValue(format_str, array);
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: value = rtdb_get(name)");
       if (format_str) free(format_str);
       return NULL;
   }
   (void) MA_free_heap(ma_handle);
   if (format_str) free(format_str);
   return returnObj;
}

PyObject *wrap_rtdb_delete(PyObject *self, PyObject *args)
{
   char *name;
   PyObject *returnObj = NULL;

   if (PyArg_Parse(args, "s", &name)) {
       if (rtdb_delete(rtdb_handle, name)) {
	 returnObj = Py_None;
         Py_INCREF(Py_None);
       }
       else {
	   PyErr_SetString(NwchemError, "rtdb_delete failed");
       }
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: value = rtdb_delete(name)");
   }
   return returnObj;
}

PyObject *wrap_rtdb_get_info(PyObject *self, PyObject *args)
{
   int nelem, ma_type;
   char *name;
   char *format_str="[iis]";
   PyObject *returnObj = 0;
   char date[26];

   if (PyArg_Parse(args, "s", &name)) {
       if (!rtdb_get_info(rtdb_handle, name, &ma_type, &nelem, date)) {
	   PyErr_SetString(NwchemError, "rtdb_get_info failed");
	   return NULL;
       }
       if (!(returnObj = PyTuple_New(3))) {
	   PyErr_SetString(NwchemError, "rtdb_get_info failed with pyobj");
	   return NULL;
       }
       PyTuple_SET_ITEM(returnObj, 0, PyInt_FromLong((long) ma_type)); 
       PyTuple_SET_ITEM(returnObj, 1, PyInt_FromLong((long) nelem)); 
       PyTuple_SET_ITEM(returnObj, 2, PyString_FromString(date)); 
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: value = rtdb_get_info(name)");
       return NULL;
   }
   return returnObj;
}


PyObject *wrap_rtdb_first(PyObject *self, PyObject *args)
{
   char name[256];
   PyObject *returnObj = NULL;

   if (rtdb_first(rtdb_handle, sizeof(name), name)) {
     returnObj = PyString_FromString(name); /*Py_BuildValue("s#", name, 1); */
   }
   else {
       PyErr_SetString(NwchemError, "rtdb_first: failed");
       return NULL;
   }
   return returnObj;
}

PyObject *wrap_rtdb_next(PyObject *self, PyObject *args)
{
   char name[256];
   PyObject *returnObj = NULL;

   if (rtdb_next(rtdb_handle, sizeof(name), name)) {
     returnObj = PyString_FromString(name); /*Py_BuildValue("s#", name, 1); */
   }
   else {
       PyErr_SetString(NwchemError, "rtdb_next: failed");
       return NULL;
   }
   return returnObj;
}

static PyObject *wrap_task_energy(PyObject *self, PyObject *args)
{
    char *theory;
    double energy;
    
    if (PyArg_Parse(args, "s", &theory)) {
	if (!rtdb_put(rtdb_handle, "task:theory", MT_CHAR, 
		      strlen(theory)+1, theory)) {
	    PyErr_SetString(NwchemError, "task_energy: putting theory failed");
	    return NULL;
	}
	if (!task_energy_(&rtdb_handle)) {
	    PyErr_SetString(NwchemError, "task_energy: failed");
	    return NULL;
	}
	if (!rtdb_get(rtdb_handle, "task:energy", MT_F_DBL, 1, &energy)) {
	    PyErr_SetString(NwchemError, "task_energy: getting energy failed");
	    return NULL;
	}
    }
    else {
	PyErr_SetString(PyExc_TypeError, "Usage: task_energy(theory)");
	return NULL;
    }
    
    return Py_BuildValue("d", energy);
}

static PyObject *wrap_task_gradient(PyObject *self, PyObject *args)
{
    char *theory;
    double energy, *gradient;
    int ma_type, nelem, ma_handle, ind, i;
    char *format_str;
    PyObject *returnObj, *eObj, *gradObj;
    
    if (PyArg_Parse(args, "s", &theory)) {
	if (!rtdb_put(rtdb_handle, "task:theory", MT_CHAR, 
		      strlen(theory)+1, theory)) {
	    PyErr_SetString(NwchemError, "task_gradient: putting theory failed");
	    return NULL;
	}
	if (!task_gradient_(&rtdb_handle)) {
	    PyErr_SetString(NwchemError, "task_gradient: failed");
	    return NULL;
	}
	if (!rtdb_get(rtdb_handle, "task:energy", MT_F_DBL, 1, &energy)) {
	    PyErr_SetString(NwchemError, "task_gradient: getting energy failed");
	    return NULL;
	}
	if (!rtdb_ma_get(rtdb_handle,"task:gradient",&ma_type,&nelem,&ma_handle)) {
	    PyErr_SetString(NwchemError, "task_gradient: getting gradient failed");
	    return NULL;
	}
	if (!MA_get_pointer(ma_handle, &gradient)) {
	    PyErr_SetString(NwchemError, "task_gradient: ma_get_ptr failed");
	    return NULL;
	}
    }
    else {
	PyErr_SetString(PyExc_TypeError, "Usage: task_gradient(theory)");
	return NULL;
    }
    
    if (!(format_str = malloc(nelem+3))) {
	PyErr_SetString(PyExc_MemoryError,
			"rtdb_get failed allocating format string");
	(void) MA_free_heap(ma_handle);
	return NULL;
    }

    ind = 0;
    format_str[ind++] = '[';
    for (i = 0; i < nelem; i++, ind++) {
	format_str[ind] = 'd';
    }
    format_str[ind++] = ']';
    format_str[ind] = 0;
    
    eObj = Py_BuildValue("d",energy);
    gradObj = Py_VaBuildValue(format_str, (void *) gradient);
    returnObj = Py_BuildValue("OO", eObj, gradObj);
    Py_DECREF(eObj);
    Py_DECREF(gradObj);
    (void) MA_free_heap(ma_handle);
    (void) free(format_str);

    return returnObj;
}

static PyObject *wrap_task_optimize(PyObject *self, PyObject *args)
{
    char *theory;
    double energy, *gradient;
    int ma_type, nelem, ma_handle, ind, i;
    char *format_str;
    PyObject *returnObj, *eObj, *gradObj;
    
    if (PyArg_Parse(args, "s", &theory)) {
	if (!rtdb_put(rtdb_handle, "task:theory", MT_CHAR, 
		      strlen(theory)+1, theory)) {
	    PyErr_SetString(NwchemError, "task_optimize: putting theory failed");
	    return NULL;
	}
	if (!task_optimize_(&rtdb_handle)) {
	    PyErr_SetString(NwchemError, "task_optimize: failed");
	    return NULL;
	}
	if (!rtdb_get(rtdb_handle, "task:energy", MT_F_DBL, 1, &energy)) {
	    PyErr_SetString(NwchemError, "task_optimize: getting energy failed");
	    return NULL;
	}
	if (!rtdb_ma_get(rtdb_handle,"task:gradient",&ma_type,&nelem,&ma_handle)) {
	    PyErr_SetString(NwchemError, "task_optimize: getting gradient failed");
	    return NULL;
	}
	if (!MA_get_pointer(ma_handle, &gradient)) {
	    PyErr_SetString(NwchemError, "task_optimize: ma_get_ptr failed");
	    return NULL;
	}
    }
    else {
	PyErr_SetString(PyExc_TypeError, "Usage: task_optimize(theory)");
	return NULL;
    }
    
    if (!(format_str = malloc(nelem+3))) {
	PyErr_SetString(PyExc_MemoryError,
			"rtdb_get failed allocating format string");
	(void) MA_free_heap(ma_handle);
	return NULL;
    }

    ind = 0;
    format_str[ind++] = '[';
    for (i = 0; i < nelem; i++, ind++) {
	format_str[ind] = 'd';
    }
    format_str[ind++] = ']';
    format_str[ind] = 0;
    
    eObj = Py_BuildValue("d",energy);
    gradObj = Py_VaBuildValue(format_str, (void *) gradient);
    returnObj = Py_BuildValue("OO", eObj, gradObj);
    Py_DECREF(eObj);
    Py_DECREF(gradObj);
    (void) MA_free_heap(ma_handle);
    (void) free(format_str);

    return returnObj;
}


static PyObject *wrap_ga_nodeid(PyObject *self, PyObject *args)
{
    int nodeid = ga_nodeid_();
    if (args) {
	PyErr_SetString(PyExc_TypeError, "Usage: nodeid()");
	return NULL;
    }

    return Py_BuildValue("i", nodeid);
}
    

static PyObject *wrap_nw_inp_from_string(PyObject *self, PyObject *args)
{
   char *pchar;

   if (PyArg_Parse(args, "s", &pchar)) {
       if (!nw_inp_from_string(rtdb_handle, pchar)) {
	   PyErr_SetString(NwchemError, "input_parse failed");
	   return NULL;
      }
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: input_parse(string)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


static struct PyMethodDef nwchem_methods[] = {
   {"rtdb_open",       wrap_rtdb_open, 0}, 
   {"rtdb_close",      wrap_rtdb_close, 0}, 
   {"pass_handle",     wrap_pass_handle, 0}, 
   {"rtdb_print",      wrap_rtdb_print, 0}, 
   {"rtdb_put",        wrap_rtdb_put, 0}, 
   {"rtdb_get",        wrap_rtdb_get, 0}, 
   {"rtdb_delete",     wrap_rtdb_delete, 0}, 
   {"rtdb_get_info",   wrap_rtdb_get_info, 0}, 
   {"rtdb_first",      wrap_rtdb_first, 0}, 
   {"rtdb_next",       wrap_rtdb_next, 0}, 
   {"task_energy",     wrap_task_energy, 0}, 
   {"task_gradient",   wrap_task_gradient, 0}, 
   {"task_optimize",   wrap_task_optimize, 0}, 
   {"input_parse",     wrap_nw_inp_from_string, 0}, 
   {"ga_nodeid",       wrap_ga_nodeid, 0}, 
   {NULL, NULL}
};

void initnwchem()
{
    PyObject *m, *d;
    m = Py_InitModule("nwchem", nwchem_methods);
    d = PyModule_GetDict(m);
    NwchemError = PyErr_NewException("nwchem.error", NULL, NULL);
    PyDict_SetItemString(d, "NWChemError", NwchemError);
}

