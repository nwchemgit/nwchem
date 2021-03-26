/* the environment variable VECLIB_MAXIMUM_THREADS is used to set the no. of threads for Apple*/
/* Accelerate framework. See man 3 Accelerate  for more details */
#include <stdio.h>
#include <stdlib.h>
#include "typesf2c.h"

void FATR veclib_set_num_threads_(Integer *n){
    char var[256];
    sprintf(var, "%ld", *n);
    setenv("VECLIB_MAXIMUM_THREADS", var, 1);
  }
Integer FATR veclib_get_num_threads_(){
  const char* value;
  int n=0;
  value = getenv("VECLIB_MAXIMUM_THREADS");
  if (NULL != value) {
    n=atol(value);
  }else{
    n=1;
  }
  return (Integer) n;
}
