/*------------------------------------------hybrid execution------------*/
/* $Id$ */
#include <assert.h>
///#define NUM_DEVICES 1
static long long device_id=-1;
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include "ga.h"
#include "typesf2c.h"

extern void FATR util_getppn_(Integer *);
//int my_smp_index();

int check_device_(long *icuda) {
	/* Check whether this process is associated with a GPU */
//	extern int armci_me, armci_master;
//	if((armci_me - armci_master)<NUM_DEVICES)
  if((my_smp_index())<*icuda) return 1;
	return 0;
}

//void device_init_(int *icuda) {
int device_init_(long *icuda,long *cuda_device_number ) {
  /* Set device_id */
  
  int dev_count_check=0;
  device_id = my_smp_index();
  cudaGetDeviceCount(&dev_count_check);
  if(dev_count_check < *icuda){
    printf("Warning: Please check whether you have %ld cuda devices per node\n",*icuda);
    fflush(stdout);
    *cuda_device_number = 30;
  }
  else {
    cudaSetDevice(device_id);
  }
  return 1;
}
int my_smp_index(){
  int ppn, my_smp_index;
  Integer* ppn_out=malloc(sizeof(Integer));
  util_getppn_(ppn_out);
  ppn= (int ) *ppn_out;
  free(ppn_out);
  return GA_Nodeid()%ppn;
}


