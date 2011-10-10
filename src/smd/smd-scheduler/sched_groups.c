/** Client-server for dynamic process group based execution of
 * tasks. This version is group-aware and has clients and server
 * within a group. The completion of tasks is handled by an array of
 * non-blocking recv requests and MPI_Waitany.  
 *
 * @author Sriram Krishnamoorthy
 */

#include <mpi.h>
#include <ga.h>
#include <assert.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <alloca.h>
#include "armci.h"

#if __STDC_VERSION__ >= 199901L
/* inline is available*/
#else
# define inline /* none */
#endif

/**
OPTS: 
* x Increase starting size of tasks buffer (reduce reallocs)
* x Assign tasks from end of tasks buffer. Use unused space to
* initiate in-place data transfer, avoiding memory copy.
* x ma_alloc() instead of malloc()? Might be faster.
*/

#define PROF_TOT_TIME   0
#define PROF_INIT_TIME  1
#define PROF_POP_TASK   0
#define PROF_RUN_TASK   0
#define PROF_SIG_START  0
#define PROF_POLL1      0
#define PROF_POLL2      0
#define PROF_LOOP1      1
#define PROF_LOOP2      0
#define PROF_CLEANUP    1

#define PROF_CLN_TOT_TIME     0
#define PROF_CLN_LOOP         0
#define PROF_CLN_WAIT_START   0
#define PROF_CLN_PREPAR       0
#define PROF_CLN_GRP          0
#define PROF_CLN_PTASK        1
#define PROF_CLN_SYNC         0
#define PROF_CLN_COMPL        0
#define PROF_CLN_DEST         0

#define TWO_SENDS             0

#define SIG_START_PRIORITY    1

#define LEADER_BCAST 1  /*leader of group bcasts task info;svr send to leader*/
#define MULTIPLE_DISPATCH 0 /*multiple tasks can be sent per request*/
#define HIERARCHICAL_DISPATCH 1 /*Multiple tasks are dispatched through a tree*/

#if HIERARCHICAL_DISPATCH && !LEADER_BCAST
#undef LEADER_BCAST
#define LEADER_BCAST 1
#endif

void walltime(double *timeus) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  *timeus = tp.tv_sec + tp.tv_usec*1.0e-6;
}

/* #define TIMESTAMP(_cond,_val) do {if(_cond) (_val)=MPI_Wtime();}while(0) */
#define TIMESTAMP(_cond,_val) do {if(_cond) \
      walltime(&(_val));}while(0)
#define TIMEDIFF(_cond,_rval,_start,_end) do {	\
    if(_cond) (_rval)+=(_end)-(_start);}while(0)
#define TIMEPRINT(_cond,_val) do { \
    if(_cond) fprintf(stdout, "%d: prof=%s val=%lfus\n",GA_Nodeid(),#_cond,(_val)*1000000.0); } while(0)


#define util_wallsec_ MPI_Wtime
/* #define util_wallsec_() 0 */

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define TAIL_STAMP 93 /*should be <sizeof(char)*/
#define MAX_SNDRCV_MSG_SIZE 32768

#define SVR 0 /*rank of server in current GA default group*/

#define NUM_INLINE_TASKS 100 /*Determines size of buf posted in client recvs*/
#define MAX_GRPS_PER_INIT 100 /*max #rungrps in one hierrachical dispatch*/

short int armci_util_sint_getval(short int *ptr);

/*Callback & other functions. Defined her and used elsewhere*/
void FATR sched_grp_server_code_();
void FATR sched_grp_client_code_();
void FATR sched_grp_insert_task_(Integer *task_list, Integer *tskid, Integer *nprocs);
double FATR sched_grp_total_ptask_time_();
double FATR sched_grp_total_ptask_min_();
double FATR sched_grp_total_ptask_max_();
double sched_grp_init_time_();
double sched_grp_cleanup_time_();
double sched_grp_loop1_time_();

/*External functions. To be defined elsewhere*/
void smd_task_populate_(void *task_list);
void smd_process_task_(Integer *tskid, Integer *p_handle);

/*Client responds to server with this structure*/
typedef struct {
  unsigned int type   ;
  unsigned int ntasks ;
  unsigned int startid;
} task_id_t;

/*Stamp to indicate sending of message for a simplistic send-recv protocol. Keep as a multiple of sizeof(int)*/
typedef struct msg_stamp_t {
  short int len;  /*len of message -- in bytes*/
/*   short int from; /\*source of the message*\/ */
  short int stamp; /*tail stamp to indicate all bytes received*/
} msg_stamp_t;

/*Message header. Keep as multiple of sizeof(int)*/
typedef struct msg_head_t {
  short int type; /*Type of message */
  short int payload; /*payload for additional info*/
} msg_head_t;

static int sg_world_me;

/*Forward declarations*/
void server_code();
void client_code();
void signal_termination(char **client_recvbufs, msg_stamp_t *stamps);
/* void armci_msg_group_barrier(ARMCI_Group *group); */

/*-------------------Server-side data structures-------------------*/

/*Position of next task to be processsed is given by (nend-ntodo+1)*/
typedef struct task_list_t {
  int nproc; /*Group size to process each task in this list*/
  int ntasks; /*#tasks in this list*/
  int nstart, nend; /*Starting and ending positions of tasks in this list*/
  int ntodo, nrunning, ndone; /*#tasks to do, running & done in this list*/
  int buf_size; /*Size of the tasks buffer (in bytes)*/
  int *tasks; /*Array of task ids (just integers)*/
  struct task_list_t *next; /*Next task list object*/
} task_list_t; 

typedef struct proc_t {
  int id;
  struct proc_t *next; 
} proc_t;

typedef enum rgrp_state_t {INIT, SENT} rgrp_state_t;

typedef struct run_grp_t {
  task_list_t *task_list; /*Each set of run tasks can have tasks from
			    only one task list*/
  int ntasks; /*#tasks in this group*/
  int nstart; /*Starting positions of these tsks in task_list*/
  proc_t *proclist; /*List of procs assgined to this object*/
  int nproc; /*#procs in the proc list*/
#if 0
  MPI_Request req; /*Request handle to poll completion*/
  task_id_t resp; /*Buf to post Irecv from client group leader*/
#else
/*   logical first; /\*1st in list of rgrps; sendbuf alloced and nbh active*\/ */
  char *sendbuf; /*pinned buffer to send to client*/
  armci_hdl_t nbh; /*corresponding non-blocking handle*/
#endif
  struct run_grp_t *prev, *next;
  int *signalptr; /*pointer at which client signals*/
  rgrp_state_t state; /*state of this run grp*/
} run_grp_t;

typedef struct serv_state_t {
  task_list_t *tlist; /*Task lists*/
  proc_t *proc_array; /*Array of proc objects initially allocated*/
  proc_t *idle_procs; /*List of idle procs*/
  int nidle_proc;    /*size of idle_procs list*/
  run_grp_t *run_grp; /*List of running groups*/
  int ntodo, nrunning, ndone; /*Summary of info in task list objs*/

  int *client_signalbuf; /*one int per client to signal to server--in default wworld*/
  msg_stamp_t *stamps; /*pinned memory to signal clients -- in comm_world*/
  char **client_tgtbufs; /*client buffers for svr to put data to -- in comm_world*/
} serv_st_t;

#define TERM_CLIENT  -1 /*Terminate. server->client*/

void armci_send(char **client_recvbufs, msg_stamp_t *stamps,
		void *src, void *dst, int bytes,
		int proc, armci_hdl_t *nbh);

static void serv_st_init(serv_st_t *serv_st) {
  const int nproc = GA_Nnodes();
  const int me = GA_Nodeid();
  int world_nproc;
  int i, c, ga_dflt_grp;
  int **arr;
  int *vbuf;
  
  /*sanity check to ensure data structure assumptions*/
  dassert(1, sizeof(short int)==2);
  dassert(1, sizeof(int)>=4);
  dassert(1, sizeof(msg_head_t)+sizeof(msg_stamp_t)<=4*sizeof(int));

/*   if(GA_Pgroup_get_default()!= GA_Pgroup_get_world()) { */
/*     printf("Should run in world group. Using ARMCI_Malloc instead of ARMCI_Malloc_group. Also using MPI_Barrier() instead of GA_Sync().\n"); */
/*     dassert(1,GA_Pgroup_get_default()== GA_Pgroup_get_world()); */
/*   } */

  MPI_Comm_rank(MPI_COMM_WORLD, &sg_world_me);
  MPI_Comm_size(MPI_COMM_WORLD, &world_nproc);
  ga_dflt_grp = GA_Pgroup_get_default();
  serv_st->tlist = NULL;
  serv_st->run_grp = NULL;
  serv_st->ntodo = serv_st->nrunning = serv_st->ndone = 0;

  ARMCI_Group default_grp;
  ARMCI_Group_get_default(&default_grp);

  arr = (int **)malloc(nproc*sizeof(int*));
/*   ARMCI_Malloc_group((void **)arr, sizeof(int)*nproc, &default_grp); */
  ARMCI_Malloc((void **)arr, sizeof(int)*nproc);
  dassert(1, arr[me] != NULL);
  serv_st->client_signalbuf = arr[me];
  bzero(serv_st->client_signalbuf, nproc*sizeof(int));

  serv_st->client_tgtbufs = (char **)calloc(world_nproc,sizeof(char *));
  dassert(1,serv_st->client_tgtbufs); 
  for(i=0; i<nproc; i++) {
    int p = GA_Pgroup_absolute_id(ga_dflt_grp, i);
    serv_st->client_tgtbufs[p] = (char *)arr[i];
  }
  free(arr);
 
  serv_st->stamps = ARMCI_Malloc_local(world_nproc*sizeof(msg_stamp_t));
  dassert(1, serv_st->stamps != NULL);

  for(i=0; i<world_nproc; i++) {
    serv_st->stamps[i].stamp = TAIL_STAMP;
  }

  serv_st->proc_array = malloc(sizeof(proc_t)*(nproc-1));
  dassert(1,serv_st->proc_array != NULL);
  for(i=0,c=0; i<nproc; i++) {
    if(i==SVR) continue;
/*     fprintf(stderr,"%d: Adding %d to proc array\n", GA_Nodeid(), i); */
    serv_st->proc_array[c].id = GA_Pgroup_absolute_id(ga_dflt_grp, i);
    serv_st->proc_array[c].next = &serv_st->proc_array[c+1];
    c++;
  }
  serv_st->proc_array[nproc-2].next = NULL;
  serv_st->idle_procs = &serv_st->proc_array[0];
  serv_st->nidle_proc = nproc-1;  
/*   GA_Sync(); */
  
  /*warmup for microbenchmarks*/
#if 0
  for(i=0; i<world_nproc; i++) {
    serv_st->client_signalbuf[i] = 0; 
  }
  vbuf = ARMCI_Malloc_local(100*sizeof(int));
  dassert(1,vbuf);
  for(i=0; i<GA_Nnodes()-1; i++) {
    int proc = serv_st->proc_array[i].id;
    armci_send(serv_st->client_tgtbufs,serv_st->stamps,vbuf,serv_st->client_tgtbufs[proc], 10*sizeof(int),proc, NULL);
  }
  ARMCI_WaitAll();
  ARMCI_AllFence();
  for(i=0; i<100; i++) {
    vbuf[i] = 0;
  }
  ARMCI_Free_local(vbuf);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
}

static inline void recur_tlist_free(task_list_t *tlist) {
  if(tlist) {
    recur_tlist_free(tlist->next);
    free(tlist->tasks);
    free(tlist);
  }
}

static void serv_st_destroy(serv_st_t *serv_st) {
  dassert(1,serv_st);
  dassert(1,serv_st->nrunning == 0);
  dassert(1,serv_st->ntodo == 0);
  dassert(1,serv_st->run_grp == NULL);
  free(serv_st->proc_array);
  recur_tlist_free(serv_st->tlist);

  ARMCI_WaitAll();
  ARMCI_Free(serv_st->client_signalbuf);
  ARMCI_Free_local(serv_st->stamps);
  free(serv_st->client_tgtbufs);
}

static inline void serv_st_reset(serv_st_t *serv_st) {
  dassert(1,serv_st);
  recur_tlist_free(serv_st->tlist);
  serv_st->tlist = NULL;
  serv_st->run_grp = NULL;
  serv_st->ntodo = serv_st->nrunning = serv_st->ndone = 0;
  serv_st->idle_procs = &serv_st->proc_array[0];
  serv_st->nidle_proc = GA_Nnodes()-1;  
}

static inline task_list_t *tlist_alloc(int nproc) {
  task_list_t *tlist = (task_list_t*)calloc(1,sizeof(task_list_t));
  tlist->nproc = nproc;
/*   tlist->buf_size = 64*sizeof(int); */
/*   tlist->tasks = (int *)malloc(tlist->buf_size); */
  return tlist;
}

static inline task_list_t *tlist_get(serv_st_t *serv_st, int nproc) {
  task_list_t *tlist = serv_st->tlist;
  while(tlist!=NULL && tlist->nproc!=nproc) 
    tlist=tlist->next;
  return tlist;
}

static task_list_t *tlist_alloc_get(serv_st_t *serv_st, int nproc) {
  task_list_t *tl, *tlist = tlist_get(serv_st, nproc);
  if(tlist==NULL) {
    tlist = tlist_alloc(nproc);
    dassert(1,tlist != NULL);
    
    for(tl=serv_st->tlist; tl && tl->next && tl->next->nproc<nproc; 
	tl=tl->next) { /*void*/ }
    if(tl) {
      dassert(1,tl->nproc > tlist->nproc);
      tlist->next = tl->next;
      tl->next = tlist;
    }
    else {
      dassert(1,!serv_st->tlist);
      serv_st->tlist = tlist;
    }
  }
  return tlist;
}

static inline void tlist_extend_buf(task_list_t *tlist) {
  dassert(1,tlist);
  tlist->tasks = realloc(tlist->tasks, 2*(tlist->buf_size+sizeof(int)));
  dassert(1,tlist->tasks);
  tlist->buf_size = 2*(tlist->buf_size+sizeof(int));
}

static inline void task_add(serv_st_t *serv_st, int tskid, int nproc) {
  task_list_t *tlist = tlist_alloc_get(serv_st, nproc);
  dassert(1,tlist != NULL);
  if(tlist->buf_size == sizeof(int)*tlist->ntasks) {
    tlist_extend_buf(tlist);
  }
  dassert(1,sizeof(int)*tlist->ntasks < tlist->buf_size);
  dassert(1,tlist->nend <tlist->buf_size);
  tlist->tasks[tlist->nend] = tskid;
  tlist->nend ++;
  tlist->ntasks ++;
  tlist->ntodo ++;
  serv_st->ntodo++;
}

static inline void serv_st_seal(serv_st_t *serv_st) {
  /*setup state -- all tasks have been inserted*/
  /*nothing to do for now*/
}

static task_list_t *serv_st_widest_tlist(serv_st_t *serv_st) {
  task_list_t *tlist=NULL, *ptr;
  int nproc = -1;
  dassert(1,serv_st);

  for(ptr=serv_st->tlist; ptr; ptr=ptr->next) {
    if(ptr->nproc > nproc) {
      nproc = ptr->nproc;
      tlist = ptr;
    }
  }
  dassert(1,serv_st->tlist==NULL || tlist!=NULL);
  return tlist;
}

static task_list_t *task_mark_running(serv_st_t *serv_st, int ntasks, int nproc) {
  int *ptr;
  task_list_t *tlist = tlist_get(serv_st, nproc);
  dassert(1,tlist);
  dassert(1,tlist->ntodo >= ntasks);
  ptr = tlist->tasks + (tlist->nend-tlist->ntodo);
  tlist->ntodo -= ntasks;
  tlist->nrunning += ntasks;
  serv_st->ntodo -= ntasks;
  serv_st->nrunning += ntasks;
  return tlist;
}

static void task_mark_done(serv_st_t *serv_st, int ntasks, int nproc) {
  task_list_t *tlist = tlist_get(serv_st, nproc);
  dassert(1,tlist);
  dassert(1,tlist);
  dassert(1,tlist->nrunning >= ntasks);
  tlist->nrunning -= ntasks;
  tlist->ndone += ntasks;
  serv_st->nrunning -= ntasks;
  serv_st->ndone += ntasks;
}

static int is_tskid_unique(serv_st_t *serv_st, int tskid) {
  task_list_t *tlist;
  dassert(1,serv_st);
  
  for(tlist=serv_st->tlist; tlist!=NULL; tlist=tlist->next) {
    int i;
    for(i=tlist->nstart; i<tlist->nend; i++) {
      if(tlist->tasks[i] == tskid)
	return 0;
    }
  }
  return 1;
}

static inline proc_t *serv_st_assign_procs(serv_st_t *serv_st, int nproc) {
  proc_t *rval, *tmp;
  int i;
  dassert(1,serv_st);
  dassert(1,nproc>0);
  if(nproc > serv_st->nidle_proc) 
    return NULL;
  rval = serv_st->idle_procs;
  for(i=0; i<nproc-1; i++) {
/*     printf("%d: included in the group is %d\n", GA_Nodeid(),serv_st->idle_procs->id); */
    serv_st->idle_procs = serv_st->idle_procs->next;
  }
  tmp = serv_st->idle_procs->next;
  serv_st->idle_procs->next = NULL;
  serv_st->idle_procs = tmp;
  serv_st->nidle_proc -= nproc;
  return rval;
}

static void serv_st_reclaim_procs(serv_st_t *serv_st, int nproc, proc_t *plist) {
  proc_t *ptr;
  dassert(1,serv_st && (nproc>0) && plist);
  for(ptr=plist; ptr->next!=NULL; ptr=ptr->next) { /*void*/ }
  ptr->next = serv_st->idle_procs;
  serv_st->idle_procs = plist;
  serv_st->nidle_proc += nproc;
}

static run_grp_t *rgrp_create(serv_st_t *serv_st, int nstart, int ntasks, 
			      int nproc, struct proc_t *procs) {
  run_grp_t *rgrp;
  dassert(1,nproc>0);
  dassert(1,procs);
  rgrp = (run_grp_t *)malloc(sizeof(run_grp_t));
  dassert(1,rgrp);

  rgrp->task_list = task_mark_running(serv_st, ntasks, nproc);
  dassert(1,rgrp->task_list);
  rgrp->ntasks = ntasks;
  rgrp->nstart = nstart;
  rgrp->proclist = procs;
  rgrp->nproc = nproc;
#if 0
  rgrp->req = MPI_REQUEST_NULL;
#else
  ARMCI_INIT_HANDLE(&rgrp->nbh);
  rgrp->sendbuf = NULL;
#endif  
  if(serv_st->run_grp) {
    rgrp->next = serv_st->run_grp;
    rgrp->prev = serv_st->run_grp->prev;
    serv_st->run_grp->prev->next = rgrp;
    serv_st->run_grp->prev = rgrp;
  } else {
    serv_st->run_grp = rgrp;
    rgrp->next = rgrp->prev = rgrp;
  }
  rgrp->signalptr = NULL;
  rgrp->state = INIT;
  return rgrp;
}

static void rgrp_finalize(serv_st_t *serv_st, run_grp_t *rgrp) {
  dassert(1,serv_st);
  dassert(1,rgrp);

  task_mark_done(serv_st, rgrp->ntasks, rgrp->nproc);
  serv_st_reclaim_procs(serv_st, rgrp->nproc, rgrp->proclist);
  rgrp->prev->next = rgrp->next;
  rgrp->next->prev = rgrp->prev;
  if(serv_st->run_grp == rgrp) {
    if(rgrp == rgrp->next)
      serv_st->run_grp = NULL;
    else
      serv_st->run_grp = serv_st->run_grp->next;
  }  
  if(rgrp->sendbuf) {
    ARMCI_Wait(&rgrp->nbh);
    ARMCI_Free_local(rgrp->sendbuf); /*allocated when sending*/
  }
  free(rgrp);
}

/*--------------SEND/RECV type functions for PUT-----------------*/

/*always send from server*/
void armci_send(char **client_recvbufs, msg_stamp_t *stamps,
		void *src, void *dst, int bytes,
		int proc, armci_hdl_t *nbh) {
#if TWO_SENDS
  msg_stamp_t *tail = (msg_stamp_t *)(MAX_SNDRCV_MSG_SIZE+(char*)client_recvbufs[proc]);
  dassert(1, bytes<MAX_SNDRCV_MSG_SIZE);
  dassert(1, bytes>0);
  dassert(1, stamps[proc].stamp == TAIL_STAMP);
  if(nbh) ARMCI_INIT_HANDLE(nbh);
  ARMCI_NbPut(src, dst, bytes, proc, NULL);
/*   ARMCI_Put(src, dst, bytes, proc); */
  stamps[proc].len = bytes;
  stamps[proc].from = GA_Pgroup_nodeid(GA_Pgroup_get_world());
/*   printf("Sending bytes=%d to=%p @proc=%d\n", bytes,dst,proc); */
/*   fflush(stdout); */
  ARMCI_NbPut(&stamps[proc], tail, sizeof(msg_stamp_t), proc, nbh);
#else
  char *dstptr;
  msg_stamp_t *tail;
  dassert(1, bytes>0);
  dstptr = (MAX_SNDRCV_MSG_SIZE-bytes)+(char *)dst;
  tail = (msg_stamp_t*)(bytes + (char*)src);

#if 0
  tail->from = sg_world_me;
#endif
  tail->len = bytes;
  tail->stamp = TAIL_STAMP;
/*   printf("Sending bytes=%d to=%p @proc=%d remote_tail=%p\n",  */
/* 	 bytes,dstptr,proc,dstptr+bytes); */
/*   fflush(stdout); */
  if(nbh) ARMCI_INIT_HANDLE(nbh);
  ARMCI_NbPut(src, dstptr, bytes+sizeof(msg_stamp_t), proc, nbh);
#endif
}

/*always send from server*/
static inline void armci_send_stamp(char **client_recvbufs, msg_stamp_t *stamps,
				    void *src, void *dst, int bytes,
				    int proc, armci_hdl_t *nbh) {
  msg_stamp_t *tail = (msg_stamp_t *)(MAX_SNDRCV_MSG_SIZE+(char*)client_recvbufs[proc]);
  char *dstptr = (MAX_SNDRCV_MSG_SIZE-bytes)+(char *)dst;
  dassert(1, bytes<MAX_SNDRCV_MSG_SIZE);
  dassert(1, bytes>0);
  dassert(1, stamps[proc].stamp == TAIL_STAMP);
  if(nbh) ARMCI_INIT_HANDLE(nbh);
  ARMCI_NbPut(src, dstptr, bytes, proc, NULL);
/*   ARMCI_Put(src, dst, bytes, proc); */
  stamps[proc].len = bytes;
/*   printf("Sending bytes=%d to=%p @proc=%d\n", bytes,dst,proc); */
/*   fflush(stdout); */
  ARMCI_NbPut(&stamps[proc], tail, sizeof(msg_stamp_t), proc, nbh);
}



/*always recv from client*/
void *armci_recvwait(void *recvbuf, int *nbytes){
  msg_stamp_t *tail = (msg_stamp_t *)(MAX_SNDRCV_MSG_SIZE + (char*)recvbuf);
/*   printf("%d: Trying to recv at buf=%p tail=%p\n", GA_Nodeid(),recvbuf,tail); */
/*   fflush(stdout); */
  while(armci_util_sint_getval(&tail->stamp)==TAIL_STAMP+1) {}
  *nbytes = tail->len;
/*   printf("%d: recv-ed bytes=%d from=%d tail-stamp=%d\n", GA_Nodeid(),*nbytes,*from, (int)tail->stamp); */
/*   fflush(stdout); */
  tail->stamp = TAIL_STAMP+1;
#if TWO_SENDS
  return recvbuf;
#else
  return (MAX_SNDRCV_MSG_SIZE - *nbytes) + (char*)recvbuf;
#endif
}

void *armci_recvpoll(void *recvbuf, int *from, int *nbytes) {
  msg_stamp_t *tail = (msg_stamp_t *)(MAX_SNDRCV_MSG_SIZE + (char*)recvbuf);
  if(armci_util_char_getval(&tail->stamp)==TAIL_STAMP+1) {
    tail->stamp = TAIL_STAMP+1;
#if 0
    *from = tail->from;
#endif
    *nbytes = tail->len;
#if TWO_SENDS
    return recvbuf;
#else
    return (MAX_SNDRCV_MSG_SIZE - *nbytes) + (char*)recvbuf;
#endif
  }
  return NULL;
}


/*-----------------------------------------------------------------*/



/*--------------------callback function--------------------*/

/* void task_add(void *task_list, int tskid, int nproc) { */
/*   task_list_t *tlist = (task_list_t *)task_list; */

/*   task_insert(tlist, tskid, nproc); */
/* } */

void signal_rgrp_start(serv_st_t *serv_st, run_grp_t *rgrp);
void signal_rgrp_start_list(serv_st_t *serv_st, int nrungrps);
run_grp_t *rpoll_completion(serv_st_t *serv_st);
run_grp_t *rwait_completion(serv_st_t *serv_st);
void full_poll_loop(serv_st_t *serv_st);
void full_wait_loop(serv_st_t *serv_st);
void full_waitany_loop(serv_st_t *serv_st);

static double t_ptask=0.0;
static double t_ptaskmin=1.0e11, t_ptaskmax=-1;
double sched_grp_total_ptask_time_() { return t_ptask; }
double sched_grp_total_ptask_min_() { 
/*   printf("Returning min=%lf\n", t_ptaskmin); */
  return t_ptaskmin; }
double sched_grp_total_ptask_max_() { 
/*   printf("Returning max=%lf\n", t_ptaskmax); */
  return t_ptaskmax; }

static int nmultiples=0; /*#times more than one run grp can be created*/
static int nrungrps=0; /*#run groups created in all*/
static int npolls=0, nwaits=0, npoll_loops=0;

int FATR sched_grp_nmultiples_() { return nmultiples; }
int FATR sched_grp_nrungrps_() { return nrungrps; }
int FATR sched_grp_npolls_() { return npolls; }
int FATR sched_grp_nwaits_() { return nwaits;}
int FATR sched_grp_npoll_loops_() { return npoll_loops;}

static double sg_tinit=0.0, sg_tloop1=0.0, sg_tcleanup=0.0;
double sched_grp_init_time_() { return sg_tinit; }
double sched_grp_cleanup_time_() { return sg_tcleanup; }
double sched_grp_loop1_time_() {  return sg_tloop1; }

#define PROF_CNT 1000
int sg_ctr=0;
double profs[PROF_CNT];

static void prof_append(double val) { 
  if(sg_ctr<PROF_CNT) profs[sg_ctr++] = val; 
}
void prof_disp() {
  int i;
  const int default_me = GA_Nodeid();
  const int cnt = MIN(PROF_CNT,sg_ctr);
  for(i=0; i<cnt; i++) {
    printf("%3d: profs[%d]=%.2lfus\n", default_me, i, profs[i]*1.0e6);
  }
}

void FATR sched_grp_reset_times_() {
  t_ptask = 0.0;
  nmultiples=nrungrps=0;
  npolls = nwaits=npoll_loops=0;
  t_ptaskmin=1.0e11;
  t_ptaskmax = -1.0;
  sg_ctr=0;
  sg_tinit = sg_tloop1 = sg_tcleanup=0.0;
}

inline int log2base(int x) {
  int r;
  for(r=0; x>1; x>>=1, r+=1);
  if(r>0 && (x & ((1<<r)-1))) r+= 1;
  return r;
}

/** Server code. Server manages a queue of tasks to be executed in
this phase. Each task in the queue is scheduled to be executed
once. When all the tasks complete execution, the server populates the
queue again and repeats the process. The server terminates when all
populating a task queue returns an empty queue. */
void server_code() {
  serv_st_t serv_st;
  task_list_t *tlist;
  run_grp_t *rgrp;
  const char *pname = __FUNCTION__;
  const int world_me = GA_Nodeid();
  double e1, e2, e3, e4, f1, f2, f3, f4, f5, f6, f7;
  double t_total=0, t_init=0,t_pop_task=0, t_sig_start=0, t_poll=0, t_loop1=0, t_loop2=0, t_poll2=0,t_cleanup=0; 
  proc_t *procs;
  int default_grp;
  int i, nstart, ntasks;
  int ngrps_to_start; /*#groups to be started*/
/*   double util_wallsec_(); */
  t_ptask = 0.0;
  double t_loop2_2=0.0;

  default_grp = GA_Pgroup_get_default();
  fprintf(stderr, "%d: 1 %s\n", world_me, pname); 

  TIMESTAMP(PROF_TOT_TIME|PROF_INIT_TIME,e1);
   fprintf(stderr, "%d: 2 %s\n", world_me, pname); 
  serv_st_init(&serv_st);
  fprintf(stderr, "%d: 3 %s\n", world_me, pname); 
  
  TIMESTAMP(PROF_INIT_TIME|PROF_LOOP1, e2);
  TIMEDIFF(PROF_INIT_TIME, sg_tinit, e1, e2);
  int niter=0;
  while(1) {
    TIMESTAMP(PROF_POP_TASK, f1);
    fprintf(stderr, "%d: 4 %s\n", world_me, pname); 
    smd_task_populate_((Integer *)&serv_st);
    serv_st_seal(&serv_st);
/*     fprintf(stderr, "%d: 5. ntodo=%d nrunning=%d ndone=%d %s\n", */
/* 	    world_me, serv_st.ntodo, serv_st.nrunning, serv_st.ndone, pname); */

    if(serv_st.ntodo == 0) break;

    tlist = serv_st.tlist;
    int first_time=1;
    TIMESTAMP(PROF_POP_TASK|PROF_LOOP2, f2);
    TIMEDIFF(PROF_POP_TASK,t_pop_task, f1, f2);
    while(serv_st.ntodo>0) {
      f1 = MPI_Wtime();
      niter+=1;
      TIMESTAMP(PROF_SIG_START, f3);
      ngrps_to_start=0;
#if SIG_START_PRIORITY
      while
#else
	if
#endif
	  (tlist && tlist->nproc<=serv_st.nidle_proc && ngrps_to_start<MAX_GRPS_PER_INIT && serv_st.ntodo>0 && 
	    (procs=serv_st_assign_procs(&serv_st,tlist->nproc))!=NULL) {
	dassert(1,tlist);

	if(first_time) first_time=0;
	else nmultiples += 1;
	
	nstart = tlist->nend - tlist->ntodo;
#if MULTIPLE_DISPATCH
	{ 
	  int x, y, ntotal, factor;
	  ntotal = tlist->ntasks;
	  /* 	    x = (int)(log(ntotal)/log(2) + 0.5); */
	  /* 	    y = (int)(log(tlist.ntodo)/log(2) + 0.5); */
	  x = log2base(ntotal);
	  y = log2base(tlist->ntodo);
	  dassert(1,x>=y);
	  factor =  1<<(x-y+1);
	  ntasks = MAX(ntotal*tlist->nproc/((GA_Nnodes()-1)*factor),1);
	  ntasks = MIN(ntasks, GA_Nnodes()-tlist->nproc+NUM_INLINE_TASKS);
	  /* 	    ntasks = 1; */
/* 	  fprintf(stderr, "%d: ntasks=%d dispatched\n", GA_Nodeid(),ntasks); */
	}
#else
	ntasks = 1;
#endif      
/* 	  fprintf(stderr, "%d: ntasks=%d dispatched\n", GA_Nodeid(),ntasks); */
	nrungrps += 1;
	rgrp = rgrp_create(&serv_st, nstart, ntasks, tlist->nproc, procs);
/* 	printf("Created rgrp=%p nstart=%d ntasks=%d nproc=%d procs[0]=%d\n", rgrp, nstart, ntasks, tlist->nproc, procs->id); */
/* 	fflush(stdout); */
	ngrps_to_start+=1;
#if !HIERARCHICAL_DISPATCH
	signal_rgrp_start_list(&serv_st, ngrps_to_start);
	ngrps_to_start=0;
#endif
	if(tlist->ntodo == 0) 
	  tlist=tlist->next;
      }
#if HIERARCHICAL_DISPATCH
      if(ngrps_to_start>0) {
/* 	printf("Dispatching tasks to %d grps\n", ngrps_to_start); */
/* 	fflush(stdout); */
	signal_rgrp_start_list(&serv_st, ngrps_to_start);
	ngrps_to_start=0;
      }
      f2 = MPI_Wtime();
      prof_append(f2-f1);
#endif
      TIMESTAMP(PROF_SIG_START|PROF_POLL1, f4);
      TIMEDIFF(PROF_SIG_START,t_sig_start,f3,f4);
#if 1
      if(serv_st.ntodo>0) {
	if(tlist->nproc>serv_st.nidle_proc) {
/* 	  printf("Calling waitany loop ntodo=%d nproc=%d nidle_proc=%d\n", */
/* 	       serv_st.ntodo, tlist->nproc, serv_st.nidle_proc); */
/* 	  fflush(stdout); */
	  full_waitany_loop(&serv_st);
	}
	else {
/* 	  full_poll_loop(&serv_st); */
	}
      }
      f3 = MPI_Wtime();
      prof_append(f3-f2);
#endif
      TIMESTAMP(PROF_POLL1,f5);
      TIMEDIFF(PROF_POLL1,t_poll,f4,f5);
    }
    TIMESTAMP(PROF_LOOP2|PROF_POLL2,f6);
    TIMEDIFF(PROF_LOOP2,t_loop2,f2,f6);
    full_wait_loop(&serv_st);
    TIMESTAMP(PROF_POLL2,f7);
    TIMEDIFF(PROF_POLL2,t_poll2,f6,f7);
    serv_st_reset(&serv_st);
  }
  TIMESTAMP(PROF_CLEANUP|PROF_LOOP1,e3);
  TIMEDIFF(PROF_LOOP1,sg_tloop1,e2,e3);
/*   printf("%d:: Signalling termination\n", GA_Nodeid()); */
  signal_termination(serv_st.client_tgtbufs, serv_st.stamps);
/*   GA_Sync(); */
  MPI_Barrier(MPI_COMM_WORLD);
  serv_st_destroy(&serv_st);
/*   fprintf(stderr, "%d:: Done signalling termination\n", GA_Nodeid());  */
  TIMESTAMP(PROF_TOT_TIME|PROF_CLEANUP, e4);
  TIMEDIFF(PROF_CLEANUP,sg_tcleanup,e3,e4);
  TIMEDIFF(PROF_TOT_TIME,t_total,e1,e4);
/*   fprintf(stderr, "%d: SERVER nmultiples=%d nrungrps=%d\n",GA_Nodeid(),nmultiples,nrungrps); */

/*   TIMEPRINT(PROF_SIG_START, t_sig_start); */
/*   TIMEPRINT(PROF_POLL1, t_poll); */
/*   TIMEPRINT(PROF_POLL2, t_poll2); */
/*   TIMEPRINT(PROF_LOOP2, t_loop2); */
/*   TIMEPRINT(PROF_LOOP1, sg_tloop1); */
/*   TIMEPRINT(PROF_CLEANUP, t_cleanup); */
/*   TIMEPRINT(PROF_INIT_TIME, t_init); */

  prof_disp();
  
/*   for(i=0; i<ctr; i++) { */
/*     printf("poll1[%d]=%lf\n",i,profs[i]); */
/*   } */
/*   printf("niter=%d t_loop2_2=%lfus\n", niter, t_loop2_2*1000000.0); */
/*   printf("%d: npoll_loops=%d npolls=%d nwaits=%d\n", ga_nodeid_(),npoll_loops,npolls,nwaits); */
  
/*   fprintf(stderr, "%d:: SERVER TOTAL time=%lf\n", ga_nodeid_(), e4-e1); */
/*   fprintf(stderr, "%d:: SERVER INIT time=%lf\n", ga_nodeid_(), e2-e1); */
/*   fprintf(stderr, "%d:: SERVER LOOP time=%lf\n", ga_nodeid_(), e3-e2); */
/*   fprintf(stderr, "%d:: SERVER TLOOP time=%lf\n", ga_nodeid_(), t_loop); */
/*   fprintf(stderr, "%d:: SERVER TLOOP1 time=%lf\n", ga_nodeid_(), t_loop1); */
/*   fprintf(stderr, "%d:: SERVER CLEANUP time=%lf\n", ga_nodeid_(), e4-e3); */
/*   fprintf(stderr, "%d:: SERVER POLLING time=%lf\n", ga_nodeid_(), t_poll); */
/*   fprintf(stderr, "%d:: SERVER POLLING2 time=%lf\n", ga_nodeid_(), t_poll2); */
/*   fprintf(stderr, "%d:: SERVER SIG_START time=%lf\n", ga_nodeid_(), t_sig_start); */
/*   fprintf(stderr, "%d:: SERVER POP_TASK time=%lf\n", ga_nodeid_(), t_pop_task); */
/*   fflush(stdout); */
}


static int int_compare(const void *v1, const void *v2) {
  int i1 = *(int *)v1;
  int i2 = *(int *)v2;
  return i1-i2;
}


#define TASK_START   -2 /*Process a task. server->client*/
#define TASK_DONE    -3 /*Done processing task. client->server*/
#define TASKLIST_START -4 /*Process list of tasks. server->client*/
#define TASKLIST_DONE  -5 /*Done processing list of tasks. client->server*/
#define RTASK_START   -6 /*Process a group of tasks. server->client*/
#define RGRPLIST_START -7 /*Process list of run groups. server->client*/
#define SIGNALGRPS_TAG  11 /*Tag to be used for signalling*/
#define SIGNALTSKS_TAG  12 /*Tag to be used for signalling*/
#define BCAST_TAG   SIGNALTSKS_TAG /*Tag for broadcast among clients*/
#define SIGNAL_TAG  SIGNALTSKS_TAG /*Tag to be used for signalling*/

/*-------------------Broadcast implementation---------------*/

/*pid_list is supposed to be ids in the world group */
void bintree(int n, int *pid_list, int root, int *Root, int *Up, int *Left, int *Right) {
  int index, up, left, right, i;
  
  index=-1;
  for(i=0; i<n; i++) {
    if(pid_list[i]==sg_world_me) {
      index=i; break;
    }
  }
  dassert(1,index!=-1);

  *Root = *Up = *Right = *Left = -1;
  if(n<=0) return;

  up    = (index-1)/2; if( up < 0) up = -1;
  left  = 2*index + 1; if(left >= n) left = -1;
  right = 2*index + 2; if(right >= n)right = -1;

  *Root = pid_list[0];  
  if(up!=-1) *Up = pid_list[up];
  if(left!=-1) *Left = pid_list[left];
  if(right!=-1) *Right = pid_list[right];
}

/*pid_list and root are in comm world, not on any sub-group*/
void broadcast(char **client_recvbufs, msg_stamp_t *stamps, int n, int *pid_list, int root, void *buf, int bytes) {
  int Root,Up,Left,Right;
  bintree(n, pid_list, root, &Root, &Up, &Left, &Right);

  dassert(1,root==Root);
/*   dassert(1,Up != -1 || src==root); */
/*   dassert(1,Up==-1 || src==Up); */

  if (Left > -1) 
    armci_send_stamp(client_recvbufs,stamps,buf,client_recvbufs[Left],bytes,Left,NULL);
  if (Right > -1) 
    armci_send_stamp(client_recvbufs,stamps,buf,client_recvbufs[Right],bytes,Right,NULL);
}

/*----------------Client implementation---------------*/

typedef struct client_st_t {
  char *recvbuf; /*recv buffer for this client*/
  char **client_recvbufs; /*recv buffer for all procs -- in comm_world*/
  int *svr_signalptr; /*remote pointer to signal work/rgrp completion*/
  int *pinbuf; /*localbuf for client to send to server*/
  msg_stamp_t *stamps; /*pinned memory to signal other clients -- in comm_world*/
#if HIERARCHICAL_DISPATCH
  char *xbuf; /*clobber buffer used by hierarchical_dispatch*/
  armci_hdl_t nbh[2]; /*sending messages to other rgrp leaders using tehse handles*/
  int flag[2]; /*Are these handles active*/
#endif
} client_st_t;

void armci_recvreset(client_st_t *client_st) {
#if HIERARCHICAL_DISPATCH
  int i;
  dassert(1, client_st);

  msg_stamp_t *tail = (msg_stamp_t *)(MAX_SNDRCV_MSG_SIZE + (char*)client_st->recvbuf);
#if 1
  for(i=0; i<2; i++) {
    if(client_st->flag[i]) {
      client_st->flag[i] = 0;
      ARMCI_Wait(&client_st->nbh[i]);
    }
  }
#else
  if(tail->stamp != TAIL_STAMP+1)
    ARMCI_WaitAll();
#endif
  tail->stamp = TAIL_STAMP+1;
#endif
}

void client_st_init(client_st_t *client_st) {
  int **arr, world_nproc, default_nproc, default_me;
  int server, i;
  msg_stamp_t *tail;
  ARMCI_Group default_grp;
  double e1, e2, e3, e4, e5, e6, e7;

/*   e1 = MPI_Wtime(); */
  dassert(1, client_st);
  MPI_Comm_rank(MPI_COMM_WORLD, &sg_world_me);
  MPI_Comm_size(MPI_COMM_WORLD, &world_nproc);
  ARMCI_Group_get_default(&default_grp);
  default_me = GA_Nodeid();
  default_nproc = GA_Nnodes();
  server = ARMCI_Absolute_id(&default_grp, SVR);

/*   e2 = MPI_Wtime(); */
  arr = (int **)malloc(sizeof(int *)*GA_Nnodes());
  dassert(1,arr);
/*   ARMCI_Malloc_group((void **)arr, MAX_SNDRCV_MSG_SIZE+sizeof(msg_stamp_t), &default_grp); */
  ARMCI_Malloc((void **)arr, MAX_SNDRCV_MSG_SIZE+sizeof(msg_stamp_t));
  dassert(1,arr[server]);
  dassert(1, arr[default_me]);
  client_st->svr_signalptr = &arr[server][default_me];

/*   e3 = MPI_Wtime(); */
  client_st->recvbuf = (char *)arr[default_me];
  tail = (msg_stamp_t*)(MAX_SNDRCV_MSG_SIZE+(char*)client_st->recvbuf);
#if 0
  tail->from=tail->len=-1;
#endif
  tail->stamp = TAIL_STAMP+1;
/*   e4 = MPI_Wtime(); */

  client_st->client_recvbufs = (char **)calloc(world_nproc,sizeof(char*));
  dassert(1,client_st->client_recvbufs);
  for(i=0; i<default_nproc; i++) {
    int p = ARMCI_Absolute_id(&default_grp,i);
    client_st->client_recvbufs[p] = (char *)arr[i];
  }
  free(arr);

/*   e5 = MPI_Wtime(); */
  /*workaround for todo fix in ARMCI*/
/*   printf("Registering locally\n"); */
/*   armci_region_register_loc(client_st->recvbuf, */
/* 			    MAX_SNDRCV_MSG_SIZE+sizeof(msg_stamp_t)); */

/*   e6 = MPI_Wtime(); */
  client_st->stamps = ARMCI_Malloc_local(world_nproc*sizeof(msg_stamp_t));
  dassert(1, client_st->stamps != NULL);
  for(i=0; i<world_nproc; i++) {
    client_st->stamps[i].stamp = TAIL_STAMP;
  }
  client_st->pinbuf = ARMCI_Malloc_local(sizeof(int));
  client_st->pinbuf[0] = 1;
/*   e7 = MPI_Wtime(); */

/*   printf("%d: s1=%.0lf s2=%.0lf s3=%.0lf s4=%.0lf s5=%.0lf s6=%.0lf\n", */
/* 	 GA_Nodeid(), (e2-e1)*1e6,(e3-e2)*1e6,(e4-e3)*1e6,(e5-e4)*1e6, */
/* 	 (e6-e5)*1e6, (e7-e6)*1e6); */
/*   fflush(stdout); */
  /*warm up cache for microbenchmark runs*/
//  for(i=0; i<10; i++) {
//    Integer tskid =  *(int *)(i+(int*)client_st->recvbuf);
//    smd_process_task_(&tskid,NULL);
//  }
#if 0
  int nbytes;
  armci_recvwait(client_st->recvbuf,&nbytes);
#endif
#if HIERARCHICAL_DISPATCH
  client_st->xbuf = (char *)malloc(MAX_SNDRCV_MSG_SIZE);
  dassert(1,client_st->xbuf);
  client_st->flag[0] = client_st->flag[1] = 0;
/*   int x=0; */
/*   for(i=0; i<MAX_SNDRCV_MSG_SIZE/sizeof(int); i++) { */
/*     x += ((int*)client_st->xbuf)[i]; */
/*   } */
#endif
  /*     printf("%d: recvbuf=%p\n", default_me, recvbuf); */
  /*     fflush(stdout); */
  MPI_Barrier(MPI_COMM_WORLD);
}

void client_st_destroy(client_st_t *client_st) {
  dassert(1, client_st);
  ARMCI_Free_local(client_st->pinbuf);
  ARMCI_Free_local(client_st->stamps);
  ARMCI_Free(client_st->recvbuf);
  free(client_st->client_recvbufs);  
#if HIERARCHICAL_DISPATCH
  free(client_st->xbuf);
#endif
}

/*---------------Various interpretations of tasks by clients--------*/

/*xbuf -- extra clobber buffer to be used internally. Of size
  MAX_SNDRCV_MSG_SIZE*/
int hierarchical_dispatch(client_st_t *client_st,
			  char *buf, int nbytes,
			  int *ntasks, int **tsks,
			  int *nproc, int **procs,
			  char *xbuf) {
#if HIERARCHICAL_DISPATCH
  msg_head_t *hdr, *fhdr;
  int ngrps, *ibuf, i;
  int fhalf, shalf, local_bytes, fbytes, sbytes;
  int fldr,sldr, fstart, sstart, ftail;
  double e1, e2, e3, e4, f1, f2;

  dassert(1,nbytes>=sizeof(int));

  hdr = (msg_head_t *)&buf[0];

  if(hdr->type == TERM_CLIENT) {
    return TERM_CLIENT;
  }
  dassert(1,hdr->type == RGRPLIST_START);
  ngrps = hdr->payload;
  ibuf = (int *)(sizeof(msg_head_t)+(char*)buf);

#if 0
  dassert(1,ngrps>0);
  if(ngrps==1) {
    hdr->type = RTASK_START;
    hdr->payload = 0;
    return regular_dispatch(client_st,buf,nbytes,ntasks,tsks,nproc,procs);
  }
  fhalf = (ngrps-1)/2;
  shalf = ngrps-1-fhalf;

  
/*   printf("%d: ngrps=%d fhalf=%d shalf=%d\n", GA_Nodeid(),ngrps,fhalf, shalf); */

/*   e1 = MPI_Wtime(); */
  {
    msg_head_t *lhdr = (msg_head_t*)xbuf;
    local_bytes = sizeof(msg_head_t)+sizeof(int)*(2+ibuf[0]+ibuf[1]);
    memcpy(xbuf, buf, local_bytes);  
    lhdr->type = RTASK_START;
    lhdr->payload = 0;
  }
/*   e2 = MPI_Wtime(); */

  fstart = 2+ibuf[0]+ibuf[1];
  ftail = fstart;
  for(i=0; i<fhalf-1; i++) {
    ftail += 2 + ibuf[ftail] + ibuf[ftail+1];
  }
  if(fhalf>0) 
    sstart = ftail + 2+ibuf[ftail]+ibuf[ftail+1];    
  else
    sstart = ftail;
  if(fhalf>0) {
    int *fbuf = &ibuf[fstart]-(sstart-ftail);
    msg_head_t *fhdr = -1 + (msg_head_t*)fbuf;

    /*find the last rgrp in fhalf and move it*/
/*     f1 = MPI_Wtime(); */
    memcpy(fbuf, &ibuf[ftail], (sstart-ftail)*sizeof(int));
/*     f2 = MPI_Wtime(); */
/*     prof_append(f2-f1); */
    
    /*send fhalf*/
    fldr = vec_min(fbuf[1], &fbuf[2+fbuf[0]]);
    fhdr->type = RGRPLIST_START;
    fhdr->payload = fhalf;
    fbytes = sizeof(msg_head_t)+(sstart-fstart)*sizeof(int);
    armci_send(client_st->client_recvbufs,client_st->stamps,fhdr, 
	       client_st->client_recvbufs[fldr],fbytes,fldr,
	       &client_st->nbh[0]);
  }
/*   e3 = MPI_Wtime(); */
  /*handle second half*/
  if(shalf>0) {
    int sldr = vec_min(ibuf[sstart+1],&ibuf[sstart+2+ibuf[sstart]]);
    msg_head_t *shdr = -1 + (msg_head_t*)&ibuf[sstart];
    int sbytes = nbytes - sstart*sizeof(int);
    shdr->type = RGRPLIST_START;
    shdr->payload = shalf;
    armci_send(client_st->client_recvbufs,client_st->stamps,shdr, 
	       client_st->client_recvbufs[sldr],sbytes,sldr,
	       &client_st->nbh[1]);
  }
/*   e4 = MPI_Wtime(); */
/*   prof_append(e2-e1); */
/*   prof_append(e3-e2); */
/*   prof_append(e4-e3); */

  if(fhalf>0) client_st->flag[0]=1;
  if(shalf>0) client_st->flag[1] = 1;

  /*handle rgrp to be processed by this proc*/
  return regular_dispatch(client_st, xbuf, local_bytes,
			  ntasks, tsks, nproc, procs); 
#elif 1
  dassert(1,ngrps>0);
  hdr->type = RTASK_START;
  hdr->payload = 0;
  if(ngrps==1) {
/*     printf("%d: ngrps==1. moving from hierarchical to regular dispatch\n",GA_Nodeid()); */
/*     fflush(stdout); */
    return regular_dispatch(client_st,buf,nbytes,ntasks,tsks,nproc,procs);
  }
  fhalf = (ngrps-1)/2;
  shalf = ngrps-1-fhalf;

  ibuf = (int *)(sizeof(msg_head_t)+(char*)buf);
  fstart = sizeof(int)*(2+ibuf[0]+ibuf[1])+sizeof(msg_head_t)+sizeof(msg_stamp_t);
  sstart=fstart;
  for(i=0; i<fhalf; i++) {
    ibuf = (int*)(sizeof(msg_head_t)+sstart+(char*)buf);
    sstart += sizeof(int)*(2+ibuf[0]+ibuf[1])+sizeof(msg_head_t)+sizeof(msg_stamp_t);
  }
  if(fhalf>0) {
    msg_head_t *fhdr = (msg_head_t*)(fstart+(char*)buf);
    int *fbuf = (int *)(fhdr+1);
    int fldr = vec_min(fbuf[1], &fbuf[2+fbuf[0]]);
    dassert(1, fhdr->type==0);
    dassert(1, fhdr->payload==0);
    fhdr->type = RGRPLIST_START;
    fhdr->payload = fhalf;
    fbytes = sstart-fstart-sizeof(msg_stamp_t);
/*     printf("%d: hierarchical dispatch. fhalf\n", GA_Nodeid()); */
/*     fflush(stdout); */
    armci_send(client_st->client_recvbufs,client_st->stamps,fhdr, 
	       client_st->client_recvbufs[fldr],fbytes,fldr,
	       &client_st->nbh[0]);
  }
/*   e3 = MPI_Wtime(); */
  /*handle second half*/
  if(shalf>0) {
    msg_head_t *shdr = (msg_head_t*)(sstart+(char*)buf);
    int *sbuf = (int *)(shdr+1);
    int sldr = vec_min(sbuf[1],&sbuf[2+sbuf[0]]);

    dassert(1, shdr->type==0);
    dassert(1, shdr->payload==0);
    shdr->type = RGRPLIST_START;
    shdr->payload = shalf;

    int sbytes = nbytes - sstart;
/*     printf("%d: hierarchical dispatch. shalf\n", GA_Nodeid()); */
/*     fflush(stdout); */
    armci_send(client_st->client_recvbufs,client_st->stamps,shdr, 
	       client_st->client_recvbufs[sldr],sbytes,sldr,
	       &client_st->nbh[1]);
  }
/*   e4 = MPI_Wtime(); */
/*   prof_append(e2-e1); */
/*   prof_append(e3-e2); */
/*   prof_append(e4-e3); */

  if(fhalf>0) client_st->flag[0]=1;
  if(shalf>0) client_st->flag[1] = 1;

  /*handle rgrp to be processed by this proc*/
  
  return regular_dispatch(client_st, buf, fstart,
			  ntasks, tsks, nproc, procs);   
#else
  dassert(1,ngrps>0);
  hdr->type = RTASK_START;
  hdr->payload = 0;
  if(ngrps==1) {
    return regular_dispatch(client_st,buf,nbytes,ntasks,tsks,nproc,procs);
  }
  ngrps-=1;
  int nend = nbytes;
  while(ngrps>0) {
    fhalf = (ngrps)/2;
    shalf = ngrps-fhalf;

    dassert(1, shalf>0);

    ibuf = (int *)(sizeof(msg_head_t)+(char*)buf);
    fstart = sizeof(int)*(2+ibuf[0]+ibuf[1])+sizeof(msg_head_t)+sizeof(msg_stamp_t);
    sstart=fstart;
    for(i=0; i<fhalf; i++) {
      ibuf = (int*)(sizeof(msg_head_t)+sstart+(char*)buf);
      sstart += sizeof(int)*(2+ibuf[0]+ibuf[1])+sizeof(msg_head_t)+sizeof(msg_stamp_t);
    }

    msg_head_t *shdr = (msg_head_t*)(sstart+(char*)buf);
    int *sbuf = (int *)(shdr+1);
    int sldr = vec_min(sbuf[1],&sbuf[2+sbuf[0]]);
    
    dassert(1, shdr->type==0);
    dassert(1, shdr->payload==0);
    shdr->type = RGRPLIST_START;
    shdr->payload = shalf;
    
    int sbytes = nend - sstart;
    armci_send(client_st->client_recvbufs,client_st->stamps,shdr, 
	       client_st->client_recvbufs[sldr],sbytes,sldr,
	       NULL);
    
    ngrps -= shalf;
    nend = sstart - sizeof(msg_stamp_t);
  }

  /*handle rgrp to be processed by this proc*/  
  return regular_dispatch(client_st, buf, fstart,
			  ntasks, tsks, nproc, procs);   
#endif

#else
  dassert(1,0);
  return 0;
#endif
}

int regular_dispatch(client_st_t *client_st, 
		     char *buf, int nbytes,
		     int *ntasks, int **tsks,
		     int *nproc, int **procs) {
  int leader, src, i;
  int world_me, *ibuf;
  msg_head_t *hdr;
  dassert(1,nbytes>=sizeof(int));

  MPI_Comm_rank(MPI_COMM_WORLD, &world_me);
  hdr = (msg_head_t *)&buf[0];  
  if(hdr->type == TERM_CLIENT) {
    return TERM_CLIENT;
  }
  dassert1(1, hdr->type == RTASK_START,hdr->type);
  ibuf = (int *)(sizeof(msg_head_t)+(char*)buf);
  *ntasks = ibuf[0];
  *nproc = ibuf[1];
  *tsks = &ibuf[2];
  *procs = &ibuf[2 + *ntasks];

/*   printf("%d. reg dispatch ntasks=%d(%d) tsks=%p nproc=%d procs=%p\n", GA_Nodeid(), */
/* 	 *ntasks,*ntasks,*tsks,*nproc,*procs); */
/*   fflush(stdout); */

/*   printf("%d. reg dispatch nelem=%d ntasks=%d(%d) tsks=%p nproc=%d procs=%p\n", GA_Nodeid(), */
/* 	 *ntasks,nelem,buf[1],*tsks,*nproc,*procs); */
/*   fflush(stdout); */

  leader = vec_min(*nproc, *procs);
  
  if(leader == world_me) {
    qsort(*procs, *nproc, sizeof(int), int_compare);
  }
#if LEADER_BCAST
  broadcast(client_st->client_recvbufs, client_st->stamps, *nproc,*procs,**procs,buf,nbytes);
#endif  
/*   printf("%d. %s(). ntasks=%d tsks=%p nproc=%d procs=%p\n", GA_Nodeid(), */
/* 	 __FUNCTION__,*ntasks,*tsks,*nproc,*procs); */
/*   fflush(stdout); */
  return hdr->type; /*The actual task*/
}


double x = 1.0;
/** Client code. Receives signals from the server to process a task or
    terminate processing and return*/
void client_code() {
  int *buf1 = NULL, bufsize1, *buf2=NULL, bufsize2;
  int flag, i;
  Integer p_handle;
  int ntsks=0, src;
  const char *pname = "client_code";
  double e1, e2, e3, e4, e5, f1, f2, f3, f4,f5,f6,f7,f8;
  double t_prepar=0, t_wait_start=0, t_grp=0,t_sync=0,t_compl=0,t_dest=0,t_loop=0,t_total=0;
  double get_doit_time_();
  double get_esp_time_();
/*   double get_gm_crt_time_(); */
/*   double get_chrg_set_time_(); */
/*   double get_gm_push_time_(); */
  const int server = GA_Pgroup_absolute_id(ga_pgroup_get_default_(),SVR);
  const int default_grp = ga_pgroup_get_default_();; /*default GA group for this dispatcher instance*/
  const int default_me = GA_Nodeid();
  int ntasks, *tsks, nproc, *procs, task_type;
  client_st_t client_st;
  char *rbuf;
  int world_grp = GA_Pgroup_get_world();
/*   assert(sizeof(int)==sizeof(task_id_t)); */
  t_ptask = 0.0;
/*   fprintf(stderr, "%d(%d): server=%d %s\n", GA_Nodeid(), GA_Nnodes(),server,pname); */

  TIMESTAMP(PROF_CLN_TOT_TIME, e1);
/*   fprintf(stderr, "%d: 0 %s\n", GA_Nodeid(), pname); */

/*   GA_Pgroup_set_default(GA_Pgroup_get_world()); */

/*   fprintf(stderr, "%d: 1 %s\n", world_me, pname); */

  client_st_init(&client_st);
/*   fprintf(stderr, "%d: 2 %s\n", world_me, pname); */

  TIMESTAMP(PROF_CLN_LOOP, e2);
  while(1) {
    int grp_me, from, nbytes;
    Integer tskid;

/*     fprintf(stderr, "%d:: Waiting for work\n", world_me); */
    f1 = MPI_Wtime();
    TIMESTAMP(PROF_CLN_WAIT_START, f1);
    rbuf = armci_recvwait(client_st.recvbuf, &nbytes);
    TIMESTAMP(PROF_CLN_WAIT_START|PROF_CLN_PREPAR,f2);
    TIMEDIFF(PROF_CLN_WAIT_START,t_wait_start,f1,f2);
    f2 = MPI_Wtime();
    prof_append(f2-f1);
/*     printf("%d:: Client got msg\n", default_me); */
/*     fflush(stdout); */
    extern int armci_me;
#if HIERARCHICAL_DISPATCH
/*     printf("%d: calling hier dispatch\n", armci_me); */
    task_type = hierarchical_dispatch(&client_st,rbuf,nbytes,
				      &ntasks, &tsks, &nproc, &procs, 
				      client_st.xbuf); 
/*     printf("%d: returned from hier dispatch\n", armci_me); */
#else
    task_type = regular_dispatch(&client_st,rbuf, nbytes,
				 &ntasks, &tsks, &nproc, &procs);
#endif
/*     printf("%d: Done dispatch. type=%d\n",GA_Nodeid(),task_type); */
/*     fflush(stdout); */
    if(task_type == TERM_CLIENT) {
/*       printf("%d: Init-ing signal termination\n", armci_me); */
      signal_termination(client_st.client_recvbufs, client_st.stamps);
/*       printf("%d: Done signal termination\n", armci_me); */
      break;
    }
#if 0
    dassert(1, task_type == RTASK_START);
    dassert(1, ntasks>0);
    dassert(1,tsks);
    dassert(1,nproc>0);
    dassert(1,procs);
    ntsks += 1;
#endif
    /*The proc ids are in world group. So create sub-group of world group*/
/*     GA_Pgroup_set_default(GA_Pgroup_get_world()); */
    p_handle = GA_Pgroup_create(procs, nproc);
    GA_Pgroup_set_default(p_handle);
/*     GA_Pgroup_sync(p_handle); */
    TIMESTAMP(PROF_CLN_GRP|PROF_CLN_PTASK, f4);
    TIMEDIFF(PROF_CLN_GRP, t_grp, f3, f4);

    for(i=0; i<ntasks; i++) {
      Integer tskid = tsks[i];
/*       printf("%d(%d):: Invoking process task tskid=%d\n", GA_Nodeid(), default_me, tskid); */
      smd_process_task_(&tskid, &p_handle);
    }
/*     for(i=0; i<world_me*2000; i++) { */
/*       x += x*0.5; */
/*     } */
    TIMESTAMP(PROF_CLN_PTASK|PROF_CLN_SYNC, f5);
/*     printf("%d min=%lf max=%lf\n",default_me,t_ptaskmin,t_ptaskmax); */
    TIMEDIFF(PROF_CLN_PTASK, t_ptask, f4, f5);
    t_ptaskmin = MIN(t_ptaskmin, f5-f4);
    t_ptaskmax = MAX(t_ptaskmax, f5-f4);
    GA_Pgroup_sync(p_handle);
    TIMESTAMP(PROF_CLN_SYNC|PROF_CLN_COMPL, f6);
    TIMEDIFF(PROF_CLN_SYNC, t_sync, f5, f6);
/*     f6 = MPI_Wtime(); */
/*     grp_me = GA_Nodeid(); */
    f3 = MPI_Wtime();
    prof_append(f3-f2);
/*     if(grp_me == 0) { */
    if(procs[0] == default_me) {
/*       printf("%d: Signalling svr signalptr=%p\n", default_me,client_st.svr_signalptr); */
/*       fflush(stdout); */
#if HIERARCHICAL_DISPATCH
      armci_recvreset(&client_st);
#endif
/*       printf("%d: signalling svr\n", armci_me); */
      ARMCI_NbPut(client_st.pinbuf,client_st.svr_signalptr,sizeof(int),server,NULL);
/*       printf("%d: Done signalling svr\n", armci_me); */
    }
    TIMESTAMP(PROF_CLN_COMPL, f7);
    TIMEDIFF(PROF_CLN_COMPL|PROF_CLN_DEST, t_compl, f6, f7);
    /*     GA_Pgroup_sync(p_handle); */
    GA_Pgroup_set_default(world_grp);
    GA_Pgroup_destroy(p_handle);
    TIMESTAMP(PROF_CLN_DEST, f8);
    TIMEDIFF(PROF_CLN_DEST, t_dest, f7, f8);
  }
/*   GA_Sync(); */
  MPI_Barrier(MPI_COMM_WORLD);
  ARMCI_WaitAll();
  GA_Pgroup_set_default(default_grp);
  client_st_destroy(&client_st);
  TIMESTAMP(PROF_CLN_TOT_TIME, e3);
  TIMEDIFF(PROF_CLN_TOT_TIME, t_total, e1, e3);
  TIMEDIFF(PROF_CLN_LOOP, t_loop, e2, e3);

  prof_disp();
/*   fprintf(stderr, "%d:: CLIENT total time=%lf\n", ga_nodeid_(), e3-e1); */
/*   fprintf(stderr, "%d:: CLIENT ntsks=%d\n", ga_nodeid_(), ntsks); */
/*   fprintf(stderr, "%d:: CLIENT loop time=%lf\n", ga_nodeid_(), e3-e2); */
/*   fprintf(stderr, "%d:: CLIENT wait start time=%lf\n", ga_nodeid_(),t_wait_start); */
/*   fprintf(stderr, "%d:: CLIENT prepare time=%lf\n", ga_nodeid_(),t_prepar); */
/*   fprintf(stderr, "%d:: CLIENT grp crt time=%lf\n", ga_nodeid_(), t_grp); */
/*   fprintf(stderr, "%d:: CLIENT ptask time=%lf\n", ga_nodeid_(), t_ptask); */
/*   fprintf(stderr, "%d:: CLIENT sync time=%lf\n", ga_nodeid_(), t_sync); */
/*   fprintf(stderr, "%d:: CLIENT compl time=%lf\n", ga_nodeid_(), t_compl); */
/*   fprintf(stderr, "%d:: CLIENT grp dstry time=%lf\n", ga_nodeid_(), t_dest); */
/*   fflush(stdout); */
/*   fprintf(stderr, "%d:: CLIENT doit time=%lf\n",ga_nodeid_(),get_doit_time_()); */
/*   fprintf(stderr, "%d:: CLIENT esp time=%lf\n",ga_nodeid_(),get_esp_time_()); */
/*   fprintf(stderr, "%d:: CLIENT chrg_set time=%lf\n",ga_nodeid_(),get_chrg_set_time_()); */
/*   fprintf(stderr, "%d:: CLIENT gm_crt time=%lf\n",ga_nodeid_(),get_gm_crt_time_()); */
/*   fprintf(stderr, "%d:: CLIENT gm_push time=%lf\n",ga_nodeid_(),get_gm_push_time_()); */
}



/*--------------Server-Client signalling code-------------*/


int ARMCI_Absolute_id(ARMCI_Group *group,int group_rank);
/*server signals all clients to terminate*/
void signal_termination(char **client_recvbufs, msg_stamp_t *stamps) {
  int child;
  msg_head_t *vbuf;
  const int rank = GA_Nodeid();
  const int size = GA_Nnodes(); 
  ARMCI_Group grp;
  
  ARMCI_Group_get_default(&grp);
  dassert(1, SVR==0);
  vbuf = ARMCI_Malloc_local(sizeof(msg_head_t)+sizeof(msg_stamp_t));
  dassert(1,vbuf);
  vbuf->type = TERM_CLIENT;
  vbuf->payload = 0; 

  child = 2*rank+1;
  if(child<size) {
    int proc = ARMCI_Absolute_id(&grp, child);
    armci_send(client_recvbufs, stamps, vbuf, client_recvbufs[proc],sizeof(msg_head_t),proc,NULL);
  }

  child = 2*rank+2;
  if(child<size) {
    int proc = ARMCI_Absolute_id(&grp, child);
    armci_send(client_recvbufs, stamps, vbuf, client_recvbufs[child],sizeof(msg_head_t),child,NULL);
  }
  ARMCI_WaitAll();
  ARMCI_Free_local(vbuf);
}

int vec_min(int n, int *v) {
  int m, i;
  dassert(1,n>0);
  m = v[0];
  for(i=1; i<n; i++) {
    m = MIN(m, v[i]);
  }
  return m;
}

void signal_rgrp_start(serv_st_t *serv_st, run_grp_t *rgrp) {
  int *ibuf, i, bytes;
  char *buf;
  proc_t *ptr, *leader_ptr;
  int leader, ldr_in_dflt_grp;
  int proc_off;
  msg_head_t *hdr;

  dassert(1,rgrp);
  dassert(1,rgrp->ntasks+2+rgrp->nproc <= 2+NUM_INLINE_TASKS + GA_Nnodes());

  bytes= sizeof(msg_head_t)+(2+rgrp->ntasks+rgrp->nproc)*sizeof(int);
  buf=ARMCI_Malloc_local(bytes+sizeof(msg_stamp_t));
  rgrp->sendbuf = (char *)buf;
  dassert(1,buf != NULL);
  hdr = (msg_head_t*)buf;
  hdr->type = RTASK_START;
  hdr->payload=0; 
  ibuf = (int *)(sizeof(msg_head_t)+(char*)buf);
  *ibuf++ = rgrp->ntasks;
  *ibuf++ = rgrp->nproc;
  for(i=0; i<rgrp->ntasks; i++){
    *ibuf++ = rgrp->task_list->tasks[rgrp->nstart+i];
/*     fprintf(stderr, "starting tskid:%d\n", buf[2+i]); */
  }
  ptr = rgrp->proclist;
  leader = -1; 
  leader_ptr = NULL;
  for(i=0; i<rgrp->nproc; i++) {
    *ibuf++ = ptr->id;
/*     fprintf(stderr,"%d: group with proc=%d of %d procs\n",GA_Nodeid(),ptr->id,rgrp->nproc); */
    if(leader==-1 || leader>ptr->id) {
      leader = ptr->id;
      leader_ptr = ptr;
    }
    ptr = ptr->next;
  }
  dassert(1,leader_ptr);
    
/*   fprintf(stderr, "%d(s): %s(): #tsks=%d tsks=%d..%d proc[0]=%d\n",GA_Nodeid(),__FUNCTION__,rgrp->ntasks,buf[2],buf[proc_off-1],buf[proc_off]); */
/*   fprintf(stderr, "%d(s): rgrp_start. nels=%d\n",GA_Nodeid(),proc_off+rgrp->nproc); */

#if 1
  {
    ldr_in_dflt_grp = leader_ptr - serv_st->proc_array;
    if(SVR<leader) ldr_in_dflt_grp += 1;
    rgrp->signalptr = &serv_st->client_signalbuf[ldr_in_dflt_grp];
    dassert(1, *rgrp->signalptr == 0);
  }
/*   printf("Setting signal data. leader=%d rgrp=%p signalptr=%p\n", */
/* 	 leader, rgrp, rgrp->signalptr); */
/*   fflush(stdout); */
#endif

#if LEADER_BCAST
  armci_send(serv_st->client_tgtbufs,serv_st->stamps, rgrp->sendbuf, serv_st->client_tgtbufs[leader], bytes, leader, &rgrp->nbh);
#else
  for(i=0; i<rgrp->nproc; i++) {
/*     fprintf(stderr, "%d: sendimg start mesg to %d\n",GA_Nodeid(), buf[proc_off+i]); */
    int proc = buf[i+2+rgrp->ntasks];
    armci_send(serv_st->client_tgtbufs,serv_st->stamps, rgrp->sendbuf, serv_st->client_tgtbufs[proc], bytes, proc,NULL);
  }  
#endif

  rgrp->state = SENT;
}


void signal_rgrp_start_list(serv_st_t *serv_st, int nrungrps) {
  int i, j, x, bufsize, dest, off, leader, ctr;
  int *grpinfo, *tskinfo;
  char *buf;
  int *ibuf;
  run_grp_t *rgrp;

  rgrp = serv_st->run_grp->state==SENT ? 
    serv_st->run_grp->prev: serv_st->run_grp;
  dassert(1,nrungrps>0);
  dassert(1,rgrp);
  dassert(1,rgrp->state == INIT);

#if !HIERARCHICAL_DISPATCH
  dassert(1,rgrp);
  for(i=0; i<nrungrps; i++,rgrp=rgrp->prev) {
    dassert(1,rgrp);
    dassert(1,rgrp->state == INIT);
/*     printf("Starting run grp %p (serv_st.run_grp=%p)\n",rgrp,serv_st->run_grp); */
    signal_rgrp_start(serv_st, rgrp);
  }
#else
/*   if(nrungrps==1) { */
/*     signal_rgrp_start(serv_st, rgrp); */
/*     return; */
/*   } */
#if 0
  rgrp = serv_st->run_grp->state==SENT ? 
    serv_st->run_grp->prev: serv_st->run_grp;
  ctr=0;
  while(ctr<nrungrps) {
    int bytes=sizeof(msg_head_t)+sizeof(msg_stamp_t);
    int cnt=0;
    int proc_dest;
    run_grp_t *ptr = rgrp, *ptr1;
    msg_head_t *hdr;
    /*find max #rgrps that can be sent in one mesg*/
    while(ctr+cnt<nrungrps) {
      int size = sizeof(int)*(2+ptr->ntasks+ptr->nproc);
      if(bytes+size>MAX_SNDRCV_MSG_SIZE)
	break;
      bytes += size;
      cnt += 1;
      ptr = ptr->prev;
    }

/*     printf("ctr=%d cnt=%d nrungrps=%d\n", ctr, cnt, nrungrps); */

    /*allocate buf and assign to head rgrp*/
    buf = ARMCI_Malloc_local(bytes);
    dassert(1, buf);
    rgrp->sendbuf = buf;
    /*initialize buffer*/
    hdr = (msg_head_t*)buf;
    hdr->type = RGRPLIST_START;
    hdr->payload = cnt;
    ibuf = (int *)(hdr+1);
    for(i=0, ptr1=rgrp; i<cnt; i++, ptr1=ptr1->prev) {
      *ibuf++ = ptr1->ntasks;
      *ibuf++ = ptr1->nproc;
      for(j=0; j<ptr1->ntasks; j++) {
	*ibuf++ = ptr1->task_list->tasks[ptr1->nstart+j];
      }
      proc_t *plist = ptr1->proclist;
      proc_t *leader_ptr=NULL;
      for(j=0; j<ptr1->nproc; j++, plist=plist->next) {
	*ibuf++ = plist->id;
	if(leader_ptr==NULL || leader_ptr->id>plist->id) {
	  leader_ptr = plist;
	}
      }

      {
	int ldr_in_dflt_grp = leader_ptr - serv_st->proc_array;
	if(SVR<leader_ptr->id) ldr_in_dflt_grp += 1;
	ptr1->signalptr = &serv_st->client_signalbuf[ldr_in_dflt_grp];
/* 	printf("leader=%d signalptr=%p\n",leader_ptr->id,ptr1->signalptr); */
/* 	fflush(stdout); */
	dassert(1, *ptr1->signalptr == 0);
      }

      if(i==0) {
	proc_dest = leader_ptr->id;
      }
      ptr1->state=SENT;
    }

    /*send buffer*/
/*     printf("Sending hierarchical dispatch msg to %d bytes=%d\n", proc_dest,bytes); */
    armci_send(serv_st->client_tgtbufs, serv_st->stamps, rgrp->sendbuf, serv_st->client_tgtbufs[proc_dest], bytes-sizeof(msg_stamp_t), proc_dest, &rgrp->nbh);

    /*move onto process remaining rgrps*/
    rgrp = ptr;
    ctr += cnt;
  }
#else
  rgrp = serv_st->run_grp->state==SENT ? 
    serv_st->run_grp->prev: serv_st->run_grp;
  ctr=0;
  while(ctr<nrungrps) {
    int bytes=0;
    int cnt=0;
    int proc_dest;
    run_grp_t *ptr = rgrp, *ptr1;
    msg_head_t *hdr;
    /*find max #rgrps that can be sent in one mesg*/
    while(ctr+cnt<nrungrps) {
      int size = sizeof(msg_head_t)+sizeof(msg_stamp_t)+sizeof(int)*(2+ptr->ntasks+ptr->nproc);
      if(bytes+size>MAX_SNDRCV_MSG_SIZE)
	break;
      bytes += size;
      cnt += 1;
      ptr = ptr->prev;
    }

/*     printf("ctr=%d cnt=%d nrungrps=%d\n", ctr, cnt, nrungrps); */
    /*allocate buf and assign to head rgrp*/
    buf = ARMCI_Malloc_local(bytes);
    dassert(1, buf);
    rgrp->sendbuf = buf;
    /*initialize buffer*/
    for(i=0, ptr1=rgrp; i<cnt; i++, ptr1=ptr1->prev) {
      hdr = (msg_head_t*)buf;
      if(i==0) {
	hdr->type = RGRPLIST_START;
	hdr->payload = cnt;
      }
      else {
	hdr->type = hdr->payload = 0;
      }
      ibuf = (int *)(hdr+1);
      *ibuf++ = ptr1->ntasks;
      *ibuf++ = ptr1->nproc;
      for(j=0; j<ptr1->ntasks; j++) {
	*ibuf++ = ptr1->task_list->tasks[ptr1->nstart+j];
      }
      proc_t *plist = ptr1->proclist;
      proc_t *leader_ptr=NULL;
      for(j=0; j<ptr1->nproc; j++, plist=plist->next) {
	*ibuf++ = plist->id;
	if(leader_ptr==NULL || leader_ptr->id>plist->id) {
	  leader_ptr = plist;
	}
      }

      {
	int ldr_in_dflt_grp = leader_ptr - serv_st->proc_array;
	if(SVR<leader_ptr->id) ldr_in_dflt_grp += 1;
	ptr1->signalptr = &serv_st->client_signalbuf[ldr_in_dflt_grp];
/* 	printf("leader=%d signalptr=%p\n",leader_ptr->id,ptr1->signalptr); */
/* 	fflush(stdout); */
	dassert(1, *ptr1->signalptr == 0);
      }

      if(i==0) {
	proc_dest = leader_ptr->id;
      }
      ptr1->state=SENT;
      buf = sizeof(msg_head_t)+sizeof(msg_stamp_t)+
	sizeof(int)*(2+ptr1->ntasks+ptr1->nproc)+(char*)buf;
    }

    /*send buffer*/
/*     printf("Sending hierarchical dispatch msg to=%d bytes=%d\n", proc_dest,bytes); */
    armci_send(serv_st->client_tgtbufs, serv_st->stamps, rgrp->sendbuf, serv_st->client_tgtbufs[proc_dest], bytes-sizeof(msg_stamp_t), proc_dest, &rgrp->nbh);

    /*move onto process remaining rgrps*/
    rgrp = ptr;
    ctr += cnt;
  }  
#endif
#endif
}


run_grp_t *rwait_completion(serv_st_t *serv_st) {
  run_grp_t *ptr = serv_st->run_grp;
  MPI_Status status;
  int nelem, x;
  const int maxspin=100;

  dassert(1, ptr); /*need something to wait on*/
  nwaits+=1;
#if 0
  MPI_Wait(&ptr->req, &status);
#else
  x=0;
  while(!armci_util_int_getval(ptr->signalptr)) {
    if(x==maxspin) {
      x=0;
      sched_yield();
    }
  }
  *ptr->signalptr = 0;
/*   armci_util_wait_int(ptr->signalptr,1,100); */
#endif
/*   MPI_Get_elements(&status, MPI_CHAR, &nelem); */
/*   dassert(1,nelem == sizeof(task_id_t)); */
/*   dassert1(1,ptr->resp.type == TASKLIST_DONE,status.MPI_SOURCE); */
/*   dassert(1,ptr->resp.startid == ptr->task_list->tasks[ptr->nstart]); */
/*   dassert(1,ptr->resp.ntasks == ptr->ntasks); */
  return ptr;
}


void full_wait_loop(serv_st_t *serv_st) {
  while(serv_st->run_grp) {
/*     dassert(1, rwait_completion(serv_st)==serv_st->run_grp); */
    rgrp_finalize(serv_st,rwait_completion(serv_st));
  }
}

void full_poll_loop(serv_st_t *serv_st) {
  run_grp_t *ptr, *ptr1, *start;
  MPI_Status status;
  int flag=0, i;
  const int nproc = GA_Nnodes();
  
  if(!serv_st->run_grp) return;
  npoll_loops += 1;
#if 0
  start = serv_st->run_grp;
  ptr=serv_st->run_grp;
  do {
/*     fprintf(stderr, "%d(s): testing completion\n", GA_Nodeid()); */
    flag=0;
    npolls+=1;
    MPI_Test(&ptr->req, &flag, &status);
/*     fprintf(stderr, "%d(s): done testing completion\n", GA_Nodeid()); */
    ptr1 = ptr->next;
    if(flag && status.MPI_SOURCE!=MPI_ANY_SOURCE) {
/*     if(flag){ */
      int nelem;
      MPI_Get_elements(&status, MPI_CHAR, &nelem);
      dassert(1,nelem == sizeof(task_id_t));
      dassert1(1,ptr->resp.type == TASKLIST_DONE,status.MPI_SOURCE);
/*       dassert(1,ptr->v[1] == ptr->tskid); */
/*       fprintf(stderr, "%d: startid=%d ntasks=%d start_pos=%d start_task=%d\n", */
/* 	      GA_Nodeid(),ptr->resp.startid,ptr->resp.ntasks,ptr->nstart, ptr->task_list->tasks[ptr->nstart]); */
      dassert(1,ptr->resp.startid == ptr->task_list->tasks[ptr->nstart]);
      dassert(1,ptr->resp.ntasks == ptr->ntasks);
      rgrp_finalize(serv_st, ptr);
    }
    ptr = ptr1;
  } while(serv_st->run_grp && ptr!=serv_st->run_grp);
#else
#if 0
  for(i=1; i<nproc; i++) {
    if(armci_util_int_getval(&serv_st->client_signalbuf[i])) {
      dassert(1,serv_st->status[i].rgrp);
      serv_st->client_signalbuf[i]=0;
      rgrp_finalize(serv_st, serv_st->status[i].rgrp);
      serv_st->status[i].rgrp = NULL;
    }
  }
#else
  start = serv_st->run_grp;
  ptr=serv_st->run_grp;
  do {
    ptr1 = ptr->next;
    if(armci_util_int_getval(ptr->signalptr)) {
      *ptr->signalptr = 0;
      rgrp_finalize(serv_st, ptr);
      break;
    }
    ptr = ptr1;
  } while(serv_st->run_grp && ptr!=serv_st->run_grp);
#endif
#endif
}


void full_waitany_loop(serv_st_t *serv_st) {
  run_grp_t *ptr, *ptr1, *start;
  MPI_Status status;
  int flag=0, i, x;
  const int nproc = GA_Nnodes();
  int maxspin=100;
  
  if(!serv_st->run_grp) return;
#if 0
  x=0;
  while(!flag)
  for(i=1; i<nproc; i++) {
    if(armci_util_int_getval(&serv_st->client_signalbuf[i])) {
      serv_st->client_signalbuf[i]=0;
      rgrp_finalize(serv_st, serv_st->status[i].rgrp);
      serv_st->status[i].rgrp = NULL;
      flag=1;
      break;
    }
    else {
      x++;
      if(x==maxspin) {x=0; sched_yield(); }
    }
  }
#else
  start = serv_st->run_grp;
  ptr=serv_st->run_grp;
  flag=0;
  do {
    ptr1 = ptr->next;
    if(ptr->signalptr && armci_util_int_getval(ptr->signalptr)) {
      *ptr->signalptr = 0;
      rgrp_finalize(serv_st, ptr);
      flag=1;
      break;
    }
    ptr = ptr1;
  } while(!flag);
#endif
}

#if 0
task_t *wait_completion(task_list_t *tlist, proc_list_t *plist) {
  task_t *ptr;
  
  if(!tlist->running) return NULL;

#if 0
  ptr = poll_completion(tlist,plist);
  if(!ptr) {
    fprintf(stderr, "MPI_Probe tobe invoked\n");
    MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    fprintf(stderr, "Back from MPI_Probe tobe invoked\n");
    ptr = poll_completion(tlist,plist);
  }
  dassert(1,ptr != NULL);
#else
  ptr=poll_completion(tlist,plist);
  if(!ptr) {
    MPI_Status status;
    MPI_Wait(&tlist->running->req, &status);
    return tlist->running;
  }
#endif
  return ptr;
}
#endif

/*--------------------fortran wrappers---------------------*/

void FATR sched_grp_server_code_() {
  server_code();
}

void FATR sched_grp_client_code_() {
/*   fprintf(stderr, "%d: sched_grp_client_code_\n", GA_Nodeid()); */
  client_code();
}

void FATR sched_grp_insert_task_(Integer *serv_state,
				 Integer *tskid, 
				 Integer *nproc) {
  serv_st_t *serv_st = (serv_st_t *)serv_state;  
  task_add(serv_st, *tskid, *nproc);  
}


/* $Id$ */
