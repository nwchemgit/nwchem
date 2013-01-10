/*------------------------------------------hybrid execution------------*/

#include <assert.h>
/* has to be automated 
   1. code checks the number of devices
   2. user can supply this information
 */
#define NUM_DEVICES 1
static long long device_id=-1;
#include <stdio.h>

int check_device_() {
	/* Check whether this process is associated with a GPU */
	extern int armci_me, armci_master;
	if((armci_me - armci_master)<NUM_DEVICES)
		return 1;
	return 0;
}

void device_init_() {
	/* Set device_id */
	extern int armci_me, armci_master;
	device_id = armci_me - armci_master;
	cudaSetDevice(device_id);
}

