/*------------------------------------------hybrid execution------------*/

#include <assert.h>
///#define NUM_DEVICES 1
static long long device_id=-1;
#include <stdio.h>

int check_device_(long *icuda) {
	/* Check whether this process is associated with a GPU */
	extern int armci_me, armci_master;
//	if((armci_me - armci_master)<NUM_DEVICES)
	if((armci_me - armci_master)<*icuda)
		return 1;
	return 0;
}

//void device_init_(int *icuda) {
int device_init_(long *icuda,long *cuda_device_number ) {
	/* Set device_id */

	extern int armci_me, armci_master;
	long dev_count_check;
	dev_count_check=0;
	device_id = armci_me - armci_master;
	cudaGetDeviceCount(&dev_count_check);
	if(dev_count_check < *icuda){
	 printf("Warning: Please check whether you have %d cuda devices per node\n",*icuda);
	 *cuda_device_number = 30;
	}
	else {
	cudaSetDevice(device_id);
	}
}

/* $Id$ */
