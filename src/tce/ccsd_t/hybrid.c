/*------------------------------------------hybrid execution------------*/

#include <assert.h>
/* has to be automated 
   1. code checks the number of devices
   2. user can supply this information
*/
///#define NUM_DEVICES 1
static long long device_id=-1;
#include <stdio.h>

int check_device_(long *icuda) {
//	printf("icuda hybid %d \n",*icuda);
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
//	printf("device_id %d and dev_count_check is %d \n",device_id,dev_count_check);
//	printf("cuda_device_number 1.1 %d \n",*cuda_device_number);
	if(dev_count_check < *icuda){
	 printf("Warning: Please check whether you have %d cuda devices per node\n",*icuda);
	 *cuda_device_number = 30;
//	 printf("cuda_device_number 1.2 %d \n",*cuda_device_number);
//	 fflush(stdin);
//	 return -1;
//	 exit(0);
	}
	else {
//	printf("cuda_device_number 1.3  %d \n",*cuda_device_number);
	cudaSetDevice(device_id);
	}
}

