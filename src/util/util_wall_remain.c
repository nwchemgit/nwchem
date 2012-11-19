/*
 $Id$
*/
#include <stdio.h>
#include <stdlib.h>
#include "typesf2c.h"

#if defined(USE_FCD)
#define util_batch_job_time_remaining_ UTIL_BATCH_JOB_TIME_REMAINING
#endif

#define NOT_AVAILABLE -1

/* util_batch_job_time_remaining returns the wall time (>=0) in seconds
   remaining for job execution, or -1 if no information is available */

#if (defined(USE_LL) && defined(JOBTIMEPATH)) || defined(LSF) || defined (PBS) || defined (SLURM)
#define DONEIT 1  

#ifdef LSF
#include <stdlib.h>
#  include <tcgmsg.h>
Integer FATR util_batch_job_time_remaining_(void)
{
  Integer wallspent, lsflimit;
  char *uval;
  uval = getenv("JOB_RLIMIT_RUN"); 
  if(uval != NULL){
    sscanf(uval,"%ld",&lsflimit);
    wallspent = (Integer) tcg_time();
#ifdef DEBUG
    (void) fprintf(stderr,"ujtr: returning time= %d -- %d \n", lsflimit,wallspent);
#endif
    /* subtract  a minute spent in initialization biz */
    return ((Integer) lsflimit - wallspent - 60);
  } else {
#ifdef DEBUG
    (void) fprintf(stderr,"ujtr: returnig time= %d  \n", -1);
#endif
    return -1;
  }
  
}
#elif defined(SLURM)
#include <slurm/slurm.h>
Integer FATR util_batch_job_time_remaining_(void)
{
  Integer wallspent=0;
  uint32_t uval;
  char *cval;
  cval = getenv("SLURM_JOBID");
  if(cval != NULL){
    sscanf(cval,"%d",&uval);
    wallspent = (Integer) slurm_get_rem_time(uval);
  }
  return ((Integer) wallspent);
}
#else
#include <unistd.h>
#include <sys/wait.h>
Integer FATR util_batch_job_time_remaining_(void)
{
  FILE *p;
  char cmd[1024];
  int t, status;

#ifdef LSFOLD
  sprintf(cmd,"%s/jobtime_lsf",JOBTIMEPATH);
#elif defined(PBS)
  sprintf(cmd,"%s/jobtime_pbs",JOBTIMEPATH);
#else
  sprintf(cmd,"%s/jobtime",JOBTIMEPATH);
#endif

  if (access(cmd,F_OK|X_OK)) {	/* If cannot access perl script */
#ifdef DEBUG
    (void) fprintf(stderr,"ujtr: cannot access %s\n",cmd);
#endif
    return NOT_AVAILABLE;
  }

  if (!(p = popen(cmd,"r"))) {
#ifdef DEBUG
    (void) fprintf(stderr,"ujtr: popen %s failed\n",cmd);
#endif
    return NOT_AVAILABLE;
  }
  
#ifdef DEBUG
    (void) fprintf(stderr,"ujtr: before read time from pipe\n");
#endif
  if (fscanf(p,"%d",&t) != 1) {
#ifdef DEBUG
    (void) fprintf(stderr,"ujtr: failed to read time from pipe\n");
    (void) fprintf(stderr,"ujtr: returnig time= %d \n", t);
#endif
    (void) fclose(p);
    (void) wait(&status);
    return NOT_AVAILABLE;
  }

  (void) fclose(p);
  (void) wait(&status);

  if (t < 0) t = 0;

#ifdef DEBUG
    (void) fprintf(stderr,"ujtr: returnig time= %d \n", t);
#endif
  return t;
}
#endif
#endif


#ifndef DONEIT
Integer FATR util_batch_job_time_remaining_(void)
{
  return NOT_AVAILABLE;
}
#endif
