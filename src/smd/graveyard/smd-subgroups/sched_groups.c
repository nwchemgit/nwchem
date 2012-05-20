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
#include <alloca.h>

#define LEADER_BCAST 1

#define util_wallsec_ MPI_Wtime
/* #define util_wallsec_() 0 */

#define MIN(x,y) ((x)<(y)?(x):(y))

#define SVR 0 /*rank of server in current GA default group*/

/*Callback & other functions. Defined her and used elsewhere*/
void FATR sched_grp_server_code_();
void FATR sched_grp_client_code_();
void FATR sched_grp_insert_task_(Integer *task_list, Integer *tskid, Integer *nprocs);
double FATR sched_grp_total_ptask_time_();

/*External functions. To be defined elsewhere*/
void populate_tasks_(void *task_list);
void process_task_(Integer *tskid, Integer *p_handle);

/*Forward declarations*/
void server_code();
void client_code();
void signal_termination(int server_id, int msg_src_mpi); 
/* void armci_msg_group_barrier(ARMCI_Group *group); */

/*-------------------Server-side data structures-------------------*/


struct proc_t;

typedef struct task_t {
  struct task_t *next, *prev;   /*Next & prev tasks in the list of tasks*/
  int nproc; /*#procs in the proc group to execute this task*/
  struct proc_t *proc_list; /*List of procs to process this task*/
  int tskid;  /*interpreted by the user*/

  int pos; /*position of this task in task list. used to relate to completion request handle*/
  int v[2]; /*Buf to post Irecv from client group leader*/
} task_t;


typedef struct {
  task_t *done, *running, *todo;
  int ndone, nrunning, ntodo; 
} task_list_t;

/*-----------------Request list----------------------*/

typedef struct {
  int idle_id;
  int nrunning, nidle;
  int *next_ids, *prev_ids;
  task_t **tasks;
  MPI_Request *reqs; /*Request handles to poll completion (for all tasks)*/
  int nreqs;
} req_list_t;

#define ID_NULL (-1)

void req_list_init(req_list_t *rlist, task_list_t *tlist) {
  const int nprocs=GA_Nnodes();
  const int nreqs = MIN(nprocs, tlist->ntodo);
  int i;

  rlist->next_ids=rlist->prev_ids=NULL;
  rlist->reqs = NULL;
  rlist->tasks = NULL;
  rlist->idle_id = 0;
  rlist->nrunning = 0;
  rlist->nidle = nreqs;
  rlist->nreqs = nreqs;
  
  if(nreqs) {
    rlist->next_ids = (int *)malloc(nreqs*sizeof(int));
    rlist->prev_ids = (int *)malloc(nreqs*sizeof(int));
    rlist->reqs = (MPI_Request *)malloc(nreqs*sizeof(MPI_Request));
    rlist->tasks = (task_t **)malloc(nreqs*sizeof(task_t *));
    assert(rlist->next_ids!=NULL);
    assert(rlist->prev_ids!=NULL);
    assert(rlist->reqs!=NULL);
    assert(rlist->tasks!=NULL);
    
    rlist->prev_ids[0] = ID_NULL;
    for(i=1; i<nreqs; i++) {
      rlist->next_ids[i-1] = i;
      rlist->prev_ids[i] = i-1;
    }
    rlist->next_ids[nreqs-1] = ID_NULL;

    for(i=0; i<nreqs; i++) {
      rlist->tasks[i] = NULL;
      rlist->reqs[i] = MPI_REQUEST_NULL;
    }
  }
}

void req_list_destroy(req_list_t *rlist) {
  assert(rlist);
  assert(rlist->nreqs == rlist->nidle);
  assert(rlist->nrunning == 0);
  free(rlist->next_ids);
  free(rlist->prev_ids);
  free(rlist->reqs);
  free(rlist->tasks);
}

int req_list_get_for_running(req_list_t *rlist, task_t *task) {
  int pos;

  assert(rlist->nidle>0); /*Should always be able to satisfy this req*/
  assert(rlist->idle_id != ID_NULL);

  pos = rlist->idle_id;
  rlist->idle_id = rlist->next_ids[rlist->idle_id];
  rlist->prev_ids[rlist->idle_id] = ID_NULL;

  rlist->next_ids[pos] = ID_NULL;
  assert(rlist->prev_ids[pos] == ID_NULL);
  
  rlist->nrunning += 1;
  rlist->nidle -= 1;

  rlist->reqs[pos] = MPI_REQUEST_NULL;
  rlist->tasks[pos] = task;
  return pos;
}

void req_list_return_as_idle(req_list_t *rlist, int pos) {
  assert(rlist->nrunning > 0);
  assert(rlist->next_ids[pos] == ID_NULL);
  assert(rlist->prev_ids[pos] == ID_NULL);

  rlist->next_ids[pos] = rlist->idle_id;
  rlist->prev_ids[pos] = ID_NULL;
  rlist->prev_ids[rlist->idle_id] = pos;
  rlist->idle_id = pos;

  rlist->nidle += 1;
  rlist->nrunning -= 1;
}


/*------------------Task list--------------------*/

void task_list_init(task_list_t *tlist) {
  tlist->done = tlist->running = tlist->todo = NULL;
  tlist->ndone = tlist->nrunning = tlist->ntodo = 0;
}

static void list_destroy(task_t *tsk) {
  if(!tsk) return;
  list_destroy(tsk->next);
  free(tsk);
}

void task_list_destroy(task_list_t *tlist) {
  list_destroy(tlist->done);
/*   list_destroy(running); */
/*   list_destroy(todo); */
  assert(tlist->running == NULL);
  assert(tlist->todo == NULL);
}

void task_list_reset(task_list_t *tlist) {
  task_list_destroy(tlist);
  task_list_init(tlist);
}

int is_tskid_unique(task_list_t *tlist, int tskid) {
  task_t *ptr;
  for(ptr=tlist->done; ptr!=NULL; ptr=ptr->next) {
    if(ptr->tskid == tskid)
      return 0;
  }
  for(ptr=tlist->running; ptr!=NULL; ptr=ptr->next) {
    if(ptr->tskid == tskid)
      return 0;
  }
  for(ptr=tlist->todo; ptr!=NULL; ptr=ptr->next) {
    if(ptr->tskid == tskid)
      return 0;
  }
  return 1;
}

static void task_rmv(task_list_t *tlist, task_t *task, task_t **head) {
  if(task->prev) task->prev->next = task->next;
  if(task->next) task->next->prev = task->prev;
  if(*head == task) *head = task->next;
  task->prev = task->next = NULL;
}

static void task_ins(task_list_t *tlist, task_t *task, task_t **head) {
  task->next = *head;
  task->prev=NULL;
  if(task->next) task->next->prev = task;
  *head = task;
}

void task_insert(task_list_t *tlist, int tskid, int nproc) {
  task_t *task;
  const int maxproc = GA_Nnodes();

  assert(is_tskid_unique(tlist, tskid));
  assert(nproc > 0 && nproc < maxproc);

  task = (task_t *)malloc(sizeof(task_t));
  assert(task != NULL);
  
  task->tskid = tskid;
  task->proc_list = NULL;
  task->nproc = nproc;
  task->pos = tlist->ntodo;

  task_ins(tlist, task, &tlist->todo);
  tlist->ntodo += 1;
}

void task_mark_running(task_list_t *tlist, task_t *task) {
  assert(task->proc_list != NULL); /**It has been assigned procs*/

  assert(tlist->done != task);
  assert(tlist->running != task);
  task_rmv(tlist, task, &tlist->todo);
  task_ins(tlist, task, &tlist->running);
  tlist->nrunning += 1;
  tlist->ntodo -= 1;
}

void task_mark_done(task_list_t *tlist, task_t *task) {
  assert(task->proc_list == NULL); /**The assigned procs have been reclaimed*/

  assert(tlist->done != task);
  assert(tlist->todo != task);
  task_rmv(tlist, task, &tlist->running);
  task_ins(tlist, task, &tlist->done);

  tlist->nrunning -= 1;
  tlist->ndone += 1;
}

typedef struct proc_t {
  int procid;
  struct proc_t *next;
} proc_t;

typedef struct {
  proc_t *idle;
  int nidle;
  proc_t *buf; 
} proc_list_t;

void plist_init(proc_list_t *plist) {
  int i, ctr;
  const int me = GA_Nodeid();
  const int nproc = GA_Nnodes();
  const int default_grp = ga_pgroup_get_default_();
  const char *pname = "plist_init";

/*   fprintf(stderr, "%d:: 1 %s\n", me,pname); */

  plist->buf = (proc_t *)malloc((nproc-1)*sizeof(proc_t));
  assert(plist->buf != NULL);
  for(i=0, ctr=0; i<nproc; i++) {
    if(i != SVR) {
      plist->buf[ctr].procid = GA_Pgroup_absolute_id(default_grp,i);
      plist->buf[ctr].next = &plist->buf[ctr+1];
      ctr+=1;
    }
  }
/*   fprintf(stderr, "%d:: 2 %s\n", me,pname); */
  plist->buf[nproc-2].next = NULL;
  plist->idle = &plist->buf[0];
  plist->nidle = nproc-1;
}

void plist_destroy(proc_list_t *plist) {
  const int nproc = GA_Nnodes();
  assert(plist->nidle == nproc-1);
  free(plist->buf);
}

proc_t* plist_assign_procs(proc_list_t *plist, int nproc) {
  int i;
  proc_t *ptr, *ret;
  if(nproc > plist->nidle)
    return NULL;

  ptr=plist->idle;
  for(i=0; i<nproc-1; i++) {
    ptr=ptr->next;
  }
  ret = plist->idle;
  plist->idle = ptr->next;
  plist->nidle -= nproc;
  ptr->next = NULL;
  return ret;
}

void plist_reclaim_procs(proc_list_t *plist, int nproc, proc_t *ptr) {
#if 1
  int i;
  proc_t *p=ptr;
  for(i=0; i<nproc-1 && p->next!=NULL; i++) {
    p = p->next;
  }
  assert(i==nproc-1 && p->next==NULL);

  p->next = plist->idle;
  plist->idle = ptr;
  plist->nidle += nproc;
#else
  proc_t *p=plist->idle;
  while(p->next!=NULL) {
    p = p->next;
  }
  p->next = ptr;
  plist->nidle += nproc;
#endif
}


/*--------------------callback function--------------------*/

void task_add(void *task_list, int tskid, int nproc) {
  task_list_t *tlist = (task_list_t *)task_list;

  task_insert(tlist, tskid, nproc);
}

void signal_task_start(task_t *task, const int proc, proc_t *procs,
		       req_list_t *rlist);
task_t *poll_completion(task_list_t *tlist, proc_list_t *plist, 
			req_list_t *rlist);
task_t *wait_completion(task_list_t *tlist, proc_list_t *plist,
			req_list_t *rlist);

static double t_ptask=0.0;
double FATR sched_grp_total_ptask_time_() { return t_ptask; }


/** Server code. Server manages a queue of tasks to be executed in
this phase. Each task in the queue is scheduled to be executed
once. When all the tasks complete execution, the server populates the
queue again and repeats the process. The server terminates when all
populating a task queue returns an empty queue. */
void server_code() {
  task_list_t tlist;
  proc_list_t plist;
  req_list_t rlist;
  task_t *ptr;
  const char *pname = "server_code";
  const int world_me = GA_Nodeid();
  double e1, e2, e3, e4, f1, f2, f3, f4, f5, f6, f7;
  double t_pop_task=0, t_sig_start=0, t_poll=0, t_loop=0, t_loop1=0, t_poll2=0; 
  proc_t *procs;
  int default_grp;
/*   double util_wallsec_(); */
  t_ptask = 0.0;

  default_grp = GA_Pgroup_get_default();
/*   fprintf(stderr, "%d: 1 %s\n", world_me, pname); */

  e1 = util_wallsec_();
  plist_init(&plist);
/*   fprintf(stderr, "%d: 2 %s\n", world_me, pname); */
  task_list_init(&tlist);
/*   fprintf(stderr, "%d: 3 %s\n", world_me, pname); */
  
  e2 = util_wallsec_();
  while(1) {
    f1 = util_wallsec_();
    task_list_reset(&tlist);  
/*     fprintf(stderr, "%d: 4 %s\n", world_me, pname); */
    populate_tasks_((Integer *)&tlist);
/*     fprintf(stderr, "%d: 5. ntodo=%d nrunning=%d ndone=%d %s\n", */
/* 	    world_me, tlist.ntodo, tlist.nrunning, tlist.ndone, pname); */
    req_list_init(&rlist, &tlist);
    f2 = util_wallsec_();
    t_pop_task += (f2-f1);

    if(tlist.ntodo == 0) break;

    f2 = util_wallsec_();
    while(tlist.ntodo>0) {
#warning "Simple server implementation. To be optimized"
      
/*       int nproc = tlist.todo->nproc; */

      f3 = util_wallsec_();
      if(1 /*(tlist.nrunning==0)*/) {
	while(tlist.ntodo>0 &&
	      (procs=plist_assign_procs(&plist, tlist.todo->nproc))!=NULL) {
/* 	  	  fprintf(stderr, "%d:: Signalling start of a task\n", GA_Nodeid()); */
	  signal_task_start(tlist.todo, tlist.todo->nproc, procs, &rlist);
	  tlist.todo->proc_list = procs;
	  task_mark_running(&tlist, tlist.todo);

/* 	  	  fprintf(stderr, "%d: Started a task. ntodo=%d nrunning=%d ndone=%d %s\n", */
/* 	  		  world_me, tlist.ntodo, tlist.nrunning, tlist.ndone, pname); */
	}
      }
      f4 = util_wallsec_();
      t_sig_start += (f4-f3);
      /*       fprintf(stderr, "Polling for task completion notifies\n"); */

      if((ptr = wait_completion(&tlist,&plist,&rlist))!=NULL) {
	plist_reclaim_procs(&plist, ptr->nproc, ptr->proc_list);
	ptr->proc_list = NULL;
	task_mark_done(&tlist, ptr);
/* 	fprintf(stderr, "%d: Completed a task. ntodo=%d nrunning=%d ndone=%d %s\n", */
/* 		world_me, tlist.ntodo, tlist.nrunning, tlist.ndone, pname); */
      }
      f5 = util_wallsec_();
      t_poll += (f5-f3);
    }

    f6 = util_wallsec_();
    t_loop1 += (f6-f1);
    while(tlist.nrunning > 0) {
      if((ptr = wait_completion(&tlist, &plist,&rlist)) != NULL) {
	plist_reclaim_procs(&plist, ptr->nproc, ptr->proc_list);
	ptr->proc_list = NULL;
	task_mark_done(&tlist, ptr);
/*         fprintf(stderr, "%d: Completed a task. ntodo=%d nrunning=%d ndone=%d %s\n", */
/*                 world_me, tlist.ntodo, tlist.nrunning, tlist.ndone, pname); */
      }
    }
    req_list_destroy(&rlist);
    f7 = util_wallsec_();
    t_poll2 += (f7-f6);
    t_loop += (f7-f2);
  }
  e3 = util_wallsec_();
  plist_destroy(&plist);
  task_list_destroy(&tlist);
  /*  req_list_destroy(&rlist);*/
/*   fprintf(stderr, "%d:: Signalling termination\n", GA_Nodeid());  */
  signal_termination(SVR,GA_Pgroup_absolute_id(default_grp,SVR));
/*   fprintf(stderr, "%d:: Done signalling termination\n", GA_Nodeid());  */
  e4 = util_wallsec_();

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


#define SIGNAL_TAG  10 /*Tag to be used for signalling*/
#define TERM_CLIENT  1 /*Terminate. server->client*/
#define TASK_START   2 /*Process a task. server->client*/
#define TASK_DONE    3 /*Done processing task. client->server*/
#define BCAST_TAG   SIGNAL_TAG /*Tag for broadcast among clients*/

/*-------------------Broadcast implementation---------------*/

/*pid_list is supposed to be ids in the world group */
void bintree(int n, int *pid_list, int root, int *Root, int *Up, int *Left, int *Right) {
  int me, index, up, left, right, i;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  index=-1;
  for(i=0; i<n; i++) {
    if(pid_list[i]==me) {
      index=i; break;
    }
  }
  assert(index!=-1);

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
void broadcast(int n, int *pid_list, int root, int src, void *buf, int bytes) {
  int Root,Up,Left,Right;
  int me;
  bintree(n, pid_list, root, &Root, &Up, &Left, &Right);

  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  assert(root==Root);
  assert(Up != -1 || src==root);
  assert(Up==-1 || src==Up);

  if (Left > -1)  
    MPI_Send(buf,bytes,MPI_BYTE,Left,BCAST_TAG,MPI_COMM_WORLD);
  if (Right > -1) 
    MPI_Send(buf,bytes,MPI_BYTE,Right,BCAST_TAG,MPI_COMM_WORLD);
}

/*----------------Client implementation---------------*/

/** Client code. Receives signals from the server to process a task or
    terminate processing and return*/
void client_code() {
  int *buf = NULL, buf_size;
  int flag;
  MPI_Status status;
  Integer p_handle;
  int ntsks=0, src;
  const char *pname = "client_code";
  double e1, e2, e3, e4, e5, f1, f2, f3, f4,f5,f6,f7,f8;
  double t_prepar=0, t_wait_start=0, t_grp=0,t_sync=0,t_compl=0,t_dest=0;
/*   double get_doit_time_(); */
/*   double get_esp_time_(); */
/*   double get_gm_crt_time_(); */
/*   double get_chrg_set_time_(); */
/*   double get_gm_push_time_(); */
  const int server = GA_Pgroup_absolute_id(ga_pgroup_get_default_(),SVR);
  const int default_grp = ga_pgroup_get_default_();; /*default GA group for this dispatcher instance*/
  const int world_me = GA_Nodeid();
  const int nproc = GA_Nnodes();

  t_ptask = 0.0;
/*   fprintf(stderr, "%d: 0 server=%d %s\n", GA_Nodeid(), server,pname); */

  e1 = util_wallsec_();
/*   fprintf(stderr, "%d: 0 %s\n", GA_Nodeid(), pname); */

/*   GA_Pgroup_set_default(GA_Pgroup_get_world()); */

/*   fprintf(stderr, "%d: 1 %s\n", world_me, pname); */

  buf_size = 1+ /*action to perform*/
    1+ /*task id - if TASK_SIGNAL*/
    nproc /*process group info*/
    ;

/*   buf = (int *)malloc(buf_size*sizeof(int)); */
  buf = (int *)alloca(buf_size*sizeof(int));
  assert(buf != NULL);

/*   fprintf(stderr, "%d: 2 %s\n", world_me, pname); */

  e2 = util_wallsec_();
  while(1) {
    int nelem, grp_me;
    Integer tskid;

    f1 = util_wallsec_();
/*     fprintf(stderr, "%d:: Waiting for work\n", world_me); */
    MPI_Recv(buf, buf_size, MPI_INT, MPI_ANY_SOURCE, SIGNAL_TAG, MPI_COMM_WORLD, &status);
    f2 = util_wallsec_();
    t_wait_start += (f2-f1);
/*     fprintf(stderr, "%d:: Client got msg from %d\n", world_me, status.MPI_SOURCE); */

    MPI_Get_elements(&status, MPI_INT, &nelem);
    assert(nelem >= 1);
      
    if(buf[0] == TERM_CLIENT) {
      /*process termination and return*/
/*        fprintf(stderr, "%d:: Recv-ed term signal\n", GA_Nodeid()); */
/*       free(buf); */
/*       fprintf(stderr, "%d:: Terminating client\n", GA_Nodeid()); */
#ifdef LEADER_BCAST
      signal_termination(SVR,status.MPI_SOURCE);
#endif
      break;
    }
/*     fprintf(stderr, "%d:: got a task to process\n", world_me); */
    /*Got a task to process*/
    assert(buf[0] == TASK_START);
    ntsks += 1;

    if(status.MPI_SOURCE == server) {
      qsort(buf+2, nelem-2, sizeof(int), int_compare);
    }
    f3  = util_wallsec_();
    t_prepar += (f3-f2);

#if LEADER_BCAST
    src = (server==status.MPI_SOURCE)?buf[2]:status.MPI_SOURCE;
    broadcast(nelem-2,buf+2,buf[2],src,buf,nelem*sizeof(int));
#endif

    /*The proc ids are in world group. So create sub-group of world group*/
    GA_Pgroup_set_default(GA_Pgroup_get_world());
    p_handle = GA_Pgroup_create(&buf[2], nelem-2);
    GA_Pgroup_set_default(p_handle);
/*     GA_Pgroup_sync(p_handle); */
    f4 = MPI_Wtime();
    t_grp += (f4-f3);

    tskid = buf[1];
/*     fprintf(stderr, "%d(%d):: Invoking process task tskid=%d\n", grp_me, world_me, tskid); */
    process_task_(&tskid, &p_handle);
    f5 = MPI_Wtime();
    t_ptask += (f5-f4);
    
    GA_Pgroup_sync(p_handle);
    grp_me = GA_Nodeid();
    f6 = util_wallsec_();
    t_sync += (f6-f5);

    if(grp_me == 0) {
      int v[2] = {TASK_DONE, tskid};
/*        fprintf(stderr, "%d(%d):: Sending ack for task %d to %d\n", */
/*  	      grp_me, world_me, tskid, SERVER); */
      MPI_Send(v, 2, MPI_INT, server, SIGNAL_TAG, MPI_COMM_WORLD);
    }
    f7 = util_wallsec_();
    t_compl += (f7-f6);
/*     GA_Pgroup_sync(p_handle); */
    GA_Pgroup_destroy(p_handle);
    GA_Pgroup_set_default(default_grp);
    f8 = util_wallsec_();
    t_dest += (f8-f7);
  }
  e3 = util_wallsec_();
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


/*server signals all clients to terminate*/
void signal_termination(int server_id, int msg_src_mpi) {
  int i, v=TERM_CLIENT;
  const int rank = GA_Nodeid();
  const int size = GA_Nnodes(); 
  const int default_grp = ga_pgroup_get_default_();

/*   assert(server_id == rank);  /\*call with your id, not someone else's*\/ */
  assert(server_id == SVR); /*Only server may invoke this method*/
#ifdef LEADER_BCAST
  {
    int *pid_list = (int *)alloca(size*sizeof(int));
    int root = GA_Pgroup_absolute_id(default_grp, SVR);
/*     int src = GA_Pgroup_absolute_id(default_grp, rank); */
    for(i=0; i<size; i++) {
      pid_list[i] = GA_Pgroup_absolute_id(default_grp,i);
    }
    broadcast(size,pid_list,root,msg_src_mpi,&v,sizeof(int));
  }
#else
  for(i=0; i<size; i++) {
    if(i != SVR) {
      MPI_Send(&v, 1, MPI_INT, GA_Pgroup_absolute_id(default_grp,i), 
	       SIGNAL_TAG, MPI_COMM_WORLD);
    }
  }
#endif
}

int vec_min(int n, int *v) {
  int m, i;
  assert(n>0);
  m = v[0];
  for(i=1; i<n; i++) {
    m = MIN(m, v[i]);
  }
  return m;
}

void signal_task_start(task_t *task, const int nproc, proc_t *procs, 
		       req_list_t *rlist) {
  int *buf, i;
  proc_t *ptr;
  const int tskid = task->tskid;
  int leader, reqpos;

/*   buf = (int *)malloc((2+nproc)*sizeof(int)); */
  buf = (int *)alloca((2+nproc)*sizeof(int));
  assert(buf != NULL);
  ptr = procs;
  buf[0] = TASK_START;
  buf[1] = tskid;
  for(i=0; i<nproc; i++) {
    buf[2+i] = ptr->procid;
    ptr = ptr->next;
  }

  leader = vec_min(nproc, &buf[2]);

#if LEADER_BCAST
  MPI_Send(buf, 2+nproc, MPI_INT, leader,SIGNAL_TAG, MPI_COMM_WORLD);
#else
  for(i=0; i<nproc; i++) {
    MPI_Send(buf, 2+nproc, MPI_INT, buf[2+i],SIGNAL_TAG, MPI_COMM_WORLD);
  }
#endif
  /*   free(buf); */

  reqpos = req_list_get_for_running(rlist, task);
  MPI_Irecv(task->v, 2, MPI_INT, leader, SIGNAL_TAG, MPI_COMM_WORLD, 
	    &rlist->reqs[reqpos]);
}


task_t *poll_completion(task_list_t *tlist, proc_list_t *plist,
			req_list_t *rlist) {
  task_t *ptr;
  MPI_Status status;
  int flag, id;

#if 0
  for(id=rlist->running_id; id!=ID_NULL; id=rlist->next_ids[id]) {
    MPI_Test(&rlist->reqs[id], &flag, &status);

    if(flag) {
      int nelem;
      MPI_Get_elements(&status, MPI_INT, &nelem);
      assert(nelem == 2);
      ptr = rlist->tasks[id];
      assert(ptr->v[1] == ptr->tskid);
      return ptr;
    }    
  }
#else
  assert(0); /*No implementation as yet of polling with this data structure*/
#endif
  return NULL;
}

task_t *wait_completion(task_list_t *tlist, proc_list_t *plist,
			req_list_t *rlist) {
  task_t *ptr;
  MPI_Status status;
  int index;

  MPI_Waitany(rlist->nreqs, rlist->reqs, &index, &status);
  if(index!=MPI_UNDEFINED) {
    int nelem;
    MPI_Get_elements(&status, MPI_INT, &nelem);
    assert(nelem == 2);
    ptr = rlist->tasks[index];
    assert(ptr->v[1] == ptr->tskid);
      
    req_list_return_as_idle(rlist, index);
    return ptr;
  }
  return NULL;
}

/*--------------------fortran wrappers---------------------*/

void FATR sched_grp_server_code_() {
  server_code();
}

void FATR sched_grp_client_code_() {
/*   fprintf(stderr, "%d: sched_grp_client_code_\n", GA_Nodeid()); */
  client_code();
}

void FATR sched_grp_insert_task_(Integer *task_list, 
				 Integer *tskid, 
				 Integer *nproc) {
  task_list_t *tlist = (task_list_t *)task_list;
  
  task_insert(tlist, *tskid, *nproc);  
}


/* $Id$ */
