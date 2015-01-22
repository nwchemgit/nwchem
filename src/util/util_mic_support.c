/*
 $Id$
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <omp.h>
#include <offload.h>

#define DEBUG_ 1
#define DEBUG2_ 1
#define DBG_NUM_DEVS 2
#ifndef DEFAULT_OFFLOAD_THREAD_MULTIPLIER
#define DEFAULT_OFFLOAD_THREAD_MULTIPLIER 4
#endif
#ifndef RANKS_PER_DEVICE
#define RANKS_PER_DEVICE 1
#endif
#define BUFSZ 256
#define MAXGETHOSTNAME 2048
#include <mpi.h>
#include "ga.h"
#include "macdecls.h"
#include "typesf2c.h"

static short int mic_get_num_initialized=0;
static int num_mic_devs;

extern int util_cgetppn();
extern int util_my_smp_index();
extern int util_my_smp_master();

/* NWC_RANKS_PER_DEVICE values 
   positive : number of ranks/device to use
   0        : use only CPUS, do not use MIC
   negative : undefined, use predefined RANKS_PER_DEVICE */
int util_getenv_nwc_ranks_per_device_() {
	int ret = -1;
	const char *name = "NWC_RANKS_PER_DEVICE";
	char *value;
	value = getenv(name);
	if (value) {
		ret = atoi(value);
	}
	return ret;
}
int util_nwc_ranks_per_device_() {
  int ranks_per_device=util_getenv_nwc_ranks_per_device_();
  if (ranks_per_device == 0) {
    return 0;
  } else if (ranks_per_device < 0){
    ranks_per_device = RANKS_PER_DEVICE;
  }
  return ranks_per_device;
}

int offload_master_(){
  int ppn,nnn,ranks_per_device=util_nwc_ranks_per_device_();
  if (ranks_per_device == 0) return 0;
  ppn=util_cgetppn();
  nnn=0;
  if(util_my_smp_index() < ranks_per_device*util_mic_get_num_devices_()) nnn = 1;


#ifdef DEBUG2
  /* debug */
    char *myhostname = (char *) malloc (MAXGETHOSTNAME);
      gethostname(myhostname, sizeof(myhostname) );
      printf( " me = %d hostname %s ppn = %d my_smp_index = %d r_p_d %d n_d %d nnn %d\n", GA_Nodeid(), myhostname,
	      ppn, util_my_smp_index(), ranks_per_device,util_mic_get_num_devices_(), nnn);
      fflush(stdout);
      free(myhostname);

#endif
  return nnn;
}

int util_mic_get_num_devices_() {
  /* only smp master does call, then bcast */
#define SIZE_GROUP 256
  MPI_Group wgroup_handle,group_handle;
  MPI_Comm group_comm;
  int err,i,ranks[SIZE_GROUP];
  int my_smp_master=util_my_smp_master();
  int size_group=util_cgetppn();

  if(mic_get_num_initialized) {
    return num_mic_devs;
  }else{

  if(util_my_smp_index() == 0) {
   num_mic_devs=_Offload_number_of_devices();
#ifdef DEBUG
   char *myhostname = (char *) malloc (MAXGETHOSTNAME);
   if(num_mic_devs != DBG_NUM_DEVS){
   gethostname(myhostname, sizeof(myhostname) );
   printf(" me %d hostname %s num_mic_devs %d \n", GA_Nodeid(), myhostname, num_mic_devs);
   if(num_mic_devs != DBG_NUM_DEVS){
   num_mic_devs=2;
   printf(" me %d reset hostname %s set num_mic_devs %d \n", GA_Nodeid(), myhostname, num_mic_devs);
   //   GA_Error("wrong number of MIC devs", (long) num_mic_devs);
   }else{
   printf(" me %d 2nd try hostname %s correct num_mic_devs %d \n", GA_Nodeid(), myhostname, num_mic_devs);
   free(myhostname);
   }
   }
#endif
  }
  
    /*get world group handle to be used later */
    err=MPI_Comm_group(MPI_COMM_WORLD, &wgroup_handle);
    if (err != MPI_SUCCESS) {
      fprintf(stdout,"util_getppn: MPI_Comm_group failed\n");
      GA_Error("util_getppn error", 0L);
    }
    for (i=0; i< size_group; i++) ranks[i] = i + my_smp_master; 
    
    /* create new group of size size_group */
    err=MPI_Group_incl(wgroup_handle, size_group, ranks, &group_handle);
    if (err != MPI_SUCCESS) {
      fprintf(stdout,"util_micdevs: MPI_Group_incl failed\n");
      GA_Error("util_micdevs error", 0L);
      fflush(stdout);
    }
    
    /* Create new new communicator for the newly created group */
    err=MPI_Comm_create(MPI_COMM_WORLD, group_handle, &group_comm);
    if (err != MPI_SUCCESS) {
      fprintf(stdout,"util_micdevs: MPI_Comm_group failed\n");
      GA_Error("util_micdevs error", 0L);
    }
    

    
    err= MPI_Bcast(&num_mic_devs, 1, MPI_INT, 0, group_comm);
    if (err != MPI_SUCCESS) {
      fprintf(stdout,"util_mics: MPI_Bcast failed\n");
      fflush(stdout);
      GA_Error("util_mic_get_num_devices error", 0L);
    }

      /*flush group and comm*/
      err=MPI_Group_free(&group_handle);
      if (err != MPI_SUCCESS) {
	fprintf(stdout,"util_micdevs: MPI_Group_free failed\n");
	GA_Error("util_micdevs error", 0L);
      }
      
      err=MPI_Comm_free(&group_comm);
      if (err != MPI_SUCCESS) {
	fprintf(stdout,"util_micdevs: MPI_Comm_free failed\n");
	GA_Error("util_micdevs error", 0L);
      }

      mic_get_num_initialized = 1;
      return num_mic_devs;
  }
}

int util_mic_get_device_() {
        int count;
	if (!offload_master_()) {
		fprintf(stdout, "%02d: need to be offload master\n", GA_Nodeid());
		GA_Error("util_mic_get_device error", 0L);
	}
        count = util_cgetppn()/util_mic_get_num_devices_();
        if (count > 0) {
          count=util_my_smp_index()/util_nwc_ranks_per_device_();
        }
	return count;


}

void FATR util_mic_set_affinity_() {
	char affinity[BUFSZ];
	char num_threads[BUFSZ];
	int pos;

	int micdev;
	int nprocs;
	int ranks_per_dev;
	int rank_on_dev;
	int nthreads;
	int ppn;
	int ranks_per_device=util_getenv_nwc_ranks_per_device_();

	if (ranks_per_device == 0) {
	  return ;
	} else if (ranks_per_device < 0){
	  ranks_per_device = RANKS_PER_DEVICE;
	}
	
	pos=snprintf(affinity, BUFSZ, "KMP_PLACE_THREADS=");
	micdev=util_mic_get_device_();
	ppn=util_cgetppn();
#pragma offload target(mic:micdev) out(nprocs)
	{
		/* do one offload to query the coprocessor for the number of cores */
		nprocs = ((int) sysconf(_SC_NPROCESSORS_ONLN) / 4) - 1;
	}

	rank_on_dev = util_my_smp_index() % util_nwc_ranks_per_device_();
	
	nthreads = nprocs / ranks_per_device * DEFAULT_OFFLOAD_THREAD_MULTIPLIER;
	
	pos+=snprintf(affinity+pos, BUFSZ-pos, "%dc,%dt,%do",
	              nprocs / ranks_per_device, DEFAULT_OFFLOAD_THREAD_MULTIPLIER,
				  rank_on_dev * (nprocs / ranks_per_device));
	snprintf(num_threads, BUFSZ, "OMP_NUM_THREADS=%d", nthreads);
	
	printf("%02d: micdev=%d nprocs=%d rank_on_dev=%d ranks_per_device=%d affinity='%s' pos=%d\n", 
	       GA_Nodeid(), micdev, nprocs, rank_on_dev, ranks_per_device, affinity, pos);
	fflush(stdout);
#pragma offload target(mic:micdev) in(affinity) in(num_threads)
	{
		/* set the affinity masks and the number of offloaded OpenMP threads */
		kmp_set_defaults("KMP_AFFINITY=compact");
		kmp_set_defaults(affinity);
		kmp_set_defaults(num_threads);
	}
}
