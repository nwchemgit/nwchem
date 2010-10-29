/* $Id$ */

/* b_ortho.c */
void b_ortho (Integer *n, DoublePrecision **colB, Integer *mapB, Integer *m, DoublePrecision **colZ, Integer *mapZ, DoublePrecision **ibuffptr, Integer *iwork, DoublePrecision *work, DoublePrecision *ort, Integer *info);

/* chol_bcst.c */
void chol_pipe_bcast (char *buf, Integer len, Integer type, Integer root, Integer k_indx, Integer *map, Integer *scratch);

/* choleski9.c */
void choleski (Integer *n, DoublePrecision **vecA, Integer *mapA, Integer *iscratch, DoublePrecision *scratch, Integer *info);
void sub_chol0 (Integer me, Integer n, DoublePrecision **col, Integer *map, Integer ncols, Integer *mycols, Integer *i_scratch, DoublePrecision *scratch, Integer *info);

/* ci_entry.c */
Integer ci_entry (Integer *me, Integer *n, Integer *i, Integer *j, Integer *map);
Integer ci_entry_ (Integer *me, Integer *n, Integer *i, Integer *j, Integer *map);
Integer ci_size_ (Integer *me, Integer *n, Integer *map);
Integer fil_mapvec_ (Integer *me, Integer *n, Integer *map, Integer *mapvec);

/* clustrf.c */
Integer clustrf_ (Integer *n, DoublePrecision *d, DoublePrecision *e, Integer *m, DoublePrecision *w, Integer *mapZ, DoublePrecision **vecZ, Integer *iblock, Integer *nsplit, Integer *isplit, DoublePrecision *ptbeval, Integer *num_clustr, Integer *clustr_info, Integer *imin, Integer *proclist, Integer *nacluster, Integer *icsplit, Integer *iscratch);

Integer clustrf4_ (Integer *n, DoublePrecision *d, DoublePrecision *e, Integer *m, DoublePrecision *w, Integer *mapZ, DoublePrecision **vecZ, Integer *iblock, Integer *nsplit, Integer *isplit, Integer *clustr_info,  Integer *nacluster, Integer *icsplit, Integer *iscratch);


/* clustrxx.c */
Integer clustrinv_(Integer *, DoublePrecision *, DoublePrecision *, DoublePrecision *, Integer *, Integer *, Integer *, Integer *, DoublePrecision **, Integer *, Integer *, Integer *, Integer *, DoublePrecision *);

/* conjug.c */
void lsl_conjugation (Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **vecB, Integer *mapB, Integer *iwork, DoublePrecision *work, DoublePrecision **buff_ptr);

/* conjug22.c */
void lsl_conjugation2 (Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **vecB, Integer *mapB, Integer *iwork, DoublePrecision *work, DoublePrecision **buff_ptr);

/* conjugation.c */
void lsl_conjugation2 (Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **vecB, Integer *mapB, Integer *iwork, DoublePrecision *work, DoublePrecision **buff_ptr);

/* de_sym.c */
void de_sym (Integer n, Integer msgtype, DoublePrecision *buffer, DoublePrecision **matrix, Integer *map, DoublePrecision **matrix2, Integer *iwork, DoublePrecision **buf_ptr);
Integer load_up (Integer nvecs, Integer *mapvec, Integer nvecs_out, Integer *map_out, DoublePrecision **matrix, DoublePrecision *buffer);
Integer un_load (Integer nvecs, Integer *mapvec, Integer nvecs_in, Integer *map_in, DoublePrecision **matrix, DoublePrecision *buffer);
Integer un_load_size (Integer nvecs, Integer *mapvec, Integer nvecs_in, Integer *map_in);

/* dsteinsch.c */
Integer dsteinsch_ (Integer *m, Integer *nclustr, Integer *max_clustr_siz, Integer *clustrinfo, Integer *schedule, Integer *mapZ, Integer *nprocsZ, Integer *pnlsiz, Integer *iscratch, Integer *freespaceZ, Integer **iptr);

/* exit.c */
void g_exit2_ (Integer *n, char *array, Integer *procmap, Integer *len, Integer *iwork);
void g_exit_ (Integer *n, char *array, Integer *procmap, Integer *len, Integer *iwork, DoublePrecision *work);
void gi_sum (Integer *buf, Integer items, Integer msgtype, Integer root, Integer snumprocs, Integer *plist, DoublePrecision *work);

/* exit2.c */
void l_exit_ (Integer *info, char *array);

/* forLL.c */
void forwardLL_ (Integer *n, Integer *mapL, Integer *mapv_L, DoublePrecision **colL, Integer *mapZ, Integer *mapvecZ, DoublePrecision **colZ, DoublePrecision *scratch, Integer *nprocs, Integer *proclist, Integer *iscratch);

/* forLU.c */
void forwardLU_ (Integer *n, Integer *mapL, Integer *mapvecL, DoublePrecision **colL, Integer *mapU, Integer *mapvecU, DoublePrecision **rowU, DoublePrecision *buffer, Integer *nprocs, Integer *proclist, Integer *ibuffer, DoublePrecision *buff, DoublePrecision **buff_ptr);

/* gmax.c */
void gmax00 (char *buf, Integer items, Integer datatype, Integer msgtype, Integer root, Integer snumprocs, Integer *plist, DoublePrecision *work);

/* inv_it2.c */
Integer inv_it (Integer *n, Integer *c1, Integer *cn, Integer *b1, Integer *bn, Integer *Zbegin, Integer *map, Integer *mapvec, DoublePrecision **vector, DoublePrecision *d, DoublePrecision *e, DoublePrecision *eval, DoublePrecision *eps, DoublePrecision *stpcrt, DoublePrecision *onenrm, Integer *iwork, DoublePrecision *work);

/* inv_it3.c */
Integer inv_it (Integer *n, Integer *c1, Integer *cn, Integer *b1, Integer *bn, Integer *Zbegin, Integer *map, Integer *mapvec, DoublePrecision **vector, DoublePrecision *d, DoublePrecision *e, DoublePrecision *eval, DoublePrecision *eps, DoublePrecision *stpcrt, DoublePrecision *onenrm, Integer *iwork, DoublePrecision *work);

/* inverse.c */
void pipe_bcst_prev_col (char *buf, Integer len, Integer type, Integer root, Integer c_indx, Integer *map, Integer *scratch);
void inverseL (Integer *msize, DoublePrecision **col, Integer *map, Integer *iwork, DoublePrecision *buffer, Integer *info);

/* lower_mxm.c */
void lu_mxm (Integer n, Integer msgtype, DoublePrecision **Lmatrix, Integer *map, DoublePrecision **Umatrix, Integer *iwork, DoublePrecision **buf_ptr, DoublePrecision *work);
Integer un_load_size1 (Integer n, Integer nvecs, Integer *mapvec, Integer nvecs_in, Integer *map_in);
Integer load_up1 (Integer n, Integer nvecs, Integer *mapvec, Integer nvecs_out, Integer *map_out, DoublePrecision **matrix, DoublePrecision **matrixU, DoublePrecision *buffer);
Integer un_load1 (Integer n, Integer nvecs, Integer *mapvec, Integer nvecs_in, Integer *map_in, DoublePrecision **matrix, DoublePrecision *buffer);

/* lu_mxm2.c */
void lu_mxm2 (Integer *n, DoublePrecision **Lmatrix, Integer *mapL, Integer *m, DoublePrecision **colU, Integer *mapU, Integer *iscratch, DoublePrecision *scratch);

/* mapdif.c */
void mapdif_ (Integer *n, Integer *mapA, Integer *mapB, Integer *iscratch, Integer *ndiff);

/* mapdif1.c */
void mapdif1_ (Integer *n, Integer *mapA, Integer *m, Integer *mapB, Integer *iscratch, Integer *ndiff);

/* mapsort.c */
Integer qsortmap (Integer *mapv, DoublePrecision **v, Integer left, Integer right);

/* mdif2b.c */
void mdif2b_ (Integer *nA, Integer *mapA, Integer *nB, Integer *mapB, Integer *nC, Integer *mapC, Integer *proclist, Integer *ndiff);

/* mdiff1.c */
void mdiff1_ (Integer *nA, Integer *mapA, Integer *nZ, Integer *mapZ, Integer *proclist, Integer *ndiff);

/* mdiff2.c */
void mdiff2_ (Integer *nA, Integer *mapA, Integer *nB, Integer *mapB, Integer *nZ, Integer *mapZ, Integer *proclist, Integer *ndiff);

/* memreq.c */
void memreq_ (Integer *type, Integer *n, Integer *mapA, Integer *mapB, Integer *mapZ, Integer *isize, Integer *rsize, Integer *ptr_size, Integer *iscratch);

/* memreq2.c */
void memreq2 (Integer *type, Integer *n, Integer *mapA, Integer *mapB, Integer *mapZ, Integer *isize, Integer *rsize, Integer *ptr_size, Integer *iscratch);

/* memreq_f.c */
void fmemreq_ (Integer *type, Integer *n, Integer *mapA, Integer *mapB, Integer *mapZ, Integer *isize, Integer *rsize, Integer *ptr_size, Integer *iscratch);

/* memreq_f2.c */
void fmemreq_ (Integer *type, Integer *n, Integer *mapA, Integer *mapB, Integer *mapZ, Integer *isize, Integer *rsize, Integer *ptr_size, Integer *iscratch);

/* mgs3.c */
void mgs3(Integer *n, DoublePrecision **colF, Integer *mapF, Integer *b1, Integer *bn, Integer *nvecsZ, Integer *first, DoublePrecision *first_buf, Integer *iscratch, DoublePrecision *scratch);

/* mgs1b.c */
void mgs_1b (Integer *n, DoublePrecision **colF, Integer *mapF, Integer *b1, Integer *bn, Integer *nvecsZ, Integer *iscratch, DoublePrecision *scratch);

/* mgs2.c */
void mgs (Integer *n, DoublePrecision **matrix, Integer *map, Integer *b1, Integer *bn, Integer *nvecsZ, Integer *ibuffer, DoublePrecision *buffer);

/* mgs22.c */
void mgs22 (Integer *n1, Integer *n2, DoublePrecision **colW, Integer *mapW, DoublePrecision **colZ, Integer *iwork, DoublePrecision *work);

/* msg22.c */
void mgs (Integer *n, DoublePrecision **matrix, Integer *map, Integer *b1, Integer *bn, Integer *nvecsZ, Integer *ibuffer, DoublePrecision *buffer);

/* mxm.c */
void mxm (Integer *n, DoublePrecision **colQ, Integer *mapQ, Integer *m, DoublePrecision **colW, Integer *mapW, Integer *iwork, DoublePrecision *work);

/* mxm2.c */
void mxm2 (Integer *n, DoublePrecision **rowQ, Integer *mapQ, Integer *m, DoublePrecision **colW, Integer *mapW, Integer *iwork, DoublePrecision *work);

/* mxm25.c */
void mxm25(Integer *n1, Integer *n2, DoublePrecision **rowQ, Integer *mapQ, Integer *m, DoublePrecision **colW, Integer *mapW, DoublePrecision **colZ, Integer *iwork, DoublePrecision *work);

/* mxm3.c */
void mxm3 (Integer *n, Integer *mapQ, DoublePrecision **colQ, Integer *m, Integer *mapW, DoublePrecision **colW, Integer *iwork, DoublePrecision *work);

/* mxm35.c */
void mxm35 (Integer *n1, Integer *n2, DoublePrecision **colQ, Integer *mapQ, Integer *m, DoublePrecision **colW, Integer *mapW, DoublePrecision **colZ, Integer *iwork, DoublePrecision *work);

/* mxm4.c */
void mxm4 (Integer *n, Integer *mapQ, DoublePrecision **colQ, Integer *m, Integer *mapW, DoublePrecision **colW, Integer *iwork, DoublePrecision *work);

/* mxm5.c */
void mxm5 (Integer *n, DoublePrecision **rowU, Integer *mapU, Integer *m, DoublePrecision **colF, Integer *mapF, Integer *iscratch, DoublePrecision *scratch);

/* mxm5x.c */
void mxm5x (Integer *n, DoublePrecision **rowU, Integer *mapU, Integer *m, DoublePrecision **colF, Integer *mapF, Integer *iscratch, DoublePrecision *scratch);

/* mxm7.c */
void mxm7 (Integer *n, Integer *mapQ, Integer *mapvecQ, DoublePrecision **colQ, Integer *m, Integer *mapW, Integer *mapvecW, DoublePrecision **colW, Integer *iwork, DoublePrecision *work);

/* mxm8.c */
void mxm8 (Integer *n, DoublePrecision **colQ, Integer *mapQ, Integer *m, DoublePrecision **colW, Integer *mapW, Integer *iwork, DoublePrecision *work);

/* mxm88.c */
void mxm88 (Integer *n, DoublePrecision **colQ, Integer *mapQ, Integer *m, DoublePrecision **colW, Integer *mapW, Integer *iwork, DoublePrecision *work, DoublePrecision **iptr);

/* mxm_ll.c */
void mxm_ll (Integer *n, DoublePrecision **colL, Integer *mapL, Integer *m, DoublePrecision **colF, Integer *mapF, Integer *iscratch, DoublePrecision *scratch);

/* mxm_ll1.c */
void mxm_llx (Integer *n, DoublePrecision **colL, Integer *mapL, Integer *m, DoublePrecision **colF, Integer *mapF, Integer *iscratch, DoublePrecision *scratch);

/* mxm_ll2.c */
void mxm_llx (Integer *n, DoublePrecision **colL, Integer *mapL, Integer *m, DoublePrecision **colF, Integer *mapF, Integer *iscratch, DoublePrecision *scratch);

/* mxm_lu.c */
void mxm_llx (Integer *n, DoublePrecision **rowL, Integer *mapL, Integer *m, DoublePrecision **colF, Integer *mapF, Integer *iscratch, DoublePrecision *scratch);

/* onenorm.c */
void one_nrm (Integer *n, Integer *m, DoublePrecision **colA, Integer *mapA, DoublePrecision *norm, Integer *iwork, DoublePrecision *work);

/* ortho.c */
void ortho (Integer *n, Integer *m, DoublePrecision **colZ, Integer *mapZ, DoublePrecision **ibuffptr, Integer *iwork, DoublePrecision *work, DoublePrecision *ort, Integer *info);

/* pdiff.c */
void pdiff( Integer *n, char *data, Integer *proclist, Integer *nprocs, Integer *iwork, char *msg1, char *msg2, Integer *info );

/* pgexit.c */
void pgexit( Integer *info, char *msg, Integer *proclist, Integer *nprocs, DoublePrecision *work );

/* pdspev_c.c */
void pdspev (Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspevx.c */
void pdspevx (Integer *ivector, Integer *irange, Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision *lb, DoublePrecision *ub, Integer *ilb, Integer *iub, DoublePrecision *abstol, Integer *meigval, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspevx2.c */
void pdspevx2 (Integer *ivector, Integer *irange, Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision *lb, DoublePrecision *ub, Integer *ilb, Integer *iub, DoublePrecision *abstol, Integer *meigval, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspevx_old.c */
void pdspevx (Integer *ivector, Integer *irange, Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision *lb, DoublePrecision *ub, Integer *ilb, Integer *iub, DoublePrecision *abstol, Integer *meigval, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspgv2_c.c */
void pdspgv2 (Integer *ifact, Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **vecB, Integer *mapB, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspgv_c.c */
void pdspgv (Integer *ifact, Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **vecB, Integer *mapB, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspgvx.c */
void pdspgvx (Integer *ifact, Integer *ivector, Integer *irange, Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **vecB, Integer *mapB, DoublePrecision *lb, DoublePrecision *ub, Integer *ilb, Integer *iub, DoublePrecision *abstol, Integer *meigval, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdsptri.c */
void pdsptri ( Integer *ivector, Integer *irange, Integer *n,  DoublePrecision *dd, DoublePrecision *ee, DoublePrecision *dplus, DoublePrecision *lplus, DoublePrecision *lb, DoublePrecision *ub, Integer *ilb, Integer *iub, DoublePrecision *abstol, Integer *meigval, DoublePrecision **vecZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision **dblptr, Integer *ibuffsize, DoublePrecision *scratch, Integer *ssize, Integer *info);


/* peigvc.c */
void peigvc_ (Integer *n, DoublePrecision *dd, DoublePrecision *ee, Integer neigval, DoublePrecision *eval, Integer *iblock, Integer *nsplit, Integer *isplit, DoublePrecision *ptbeval, Integer *one_blk, Integer *clustr_info, Integer *proclist, Integer *schedule, Integer *i_scrat, Integer *cvecs, Integer **iptr, Integer *mapZ, Integer *mapvZ, DoublePrecision **vecZ, DoublePrecision *d_scrat, Integer *iscrat, int numvec);

/* pipe_fut.c */
void pipe_bcst_fut_col (Integer n, char *buf, Integer len, Integer type, Integer root, Integer c_indx, Integer *map, Integer *scratch);

/* pmmLSL.c */
void pmmLSL (Integer *n, Integer *mapA, Integer *mapvecA, DoublePrecision **vecA, Integer *mapB, Integer *mapvecB, DoublePrecision **vecB, Integer *iscratch, DoublePrecision *dscratch, DoublePrecision **buff_ptr);
void zero_out (Integer n, DoublePrecision *array);

/* pmmLUL.c */
void pmmLUL (Integer *n, Integer *mapA, Integer *mapvecA, DoublePrecision **vecA, Integer *mapB, Integer *mapvecB, DoublePrecision **vecB, Integer *iscratch, DoublePrecision *dscratch, DoublePrecision **buff_ptr);

/* pmmLUL2.c */
void pmmLUL2 (Integer *n, Integer *mapA, Integer *mapvecA, DoublePrecision **vecA, Integer *mapB, Integer *mapvecB, DoublePrecision **vecB, Integer *iscratch, DoublePrecision *dscratch, DoublePrecision **buff_ptr);

/* pmmlsl2.c */
void pmmlsl2 (Integer *n, Integer *mapA, Integer *mapvecA, DoublePrecision **vecA, Integer *mapB, Integer *mapvecB, DoublePrecision **vecB, Integer *iscratch, DoublePrecision *dscratch, DoublePrecision **buff_ptr);

/* prev_column.c */
Integer prev_column (Integer node, Integer *map, Integer indx);

/* pstein.c */
void pstein (Integer *n, DoublePrecision *dd, DoublePrecision *ee, Integer *meigval, DoublePrecision *eval, Integer *iblock, Integer *nsplit, Integer *isplit, Integer *mapZ, DoublePrecision **vecZ, DoublePrecision *ddwork, Integer *iiwork, Integer **ppiwork, Integer *info);

/* pstein_new.c */
void pstein (Integer *n, DoublePrecision *dd, DoublePrecision *ee, Integer *meigval, DoublePrecision *eval, Integer *iblock, Integer *nsplit, Integer *isplit, Integer *mapZ, DoublePrecision **vecZ, DoublePrecision *ddwork, Integer *iiwork, Integer **ppiwork, Integer *info);

/* pstein_old.c */
void pstein (Integer *n, DoublePrecision *dd, DoublePrecision *ee, Integer *meigval, DoublePrecision *eval, Integer *iblock, Integer *nsplit, Integer *isplit, Integer *mapZ, DoublePrecision **vecZ, DoublePrecision *ddwork, Integer *iiwork, Integer **ppiwork, Integer *info);

/* pxerbla2_.c */
void pxerbla2_ (Integer *n, char *array, Integer *procmap, Integer *len, Integer *iwork, Integer *info);

/* qsort.c */
Integer iqsort (Integer **v, Integer left, Integer right);

/* qsort1.c */
void qsort1 (Integer *n, Integer *mapv, DoublePrecision **v, Integer left, Integer right, DoublePrecision *scratch);

/* randomize.c */
void i_random (Integer n, Integer *list, Integer *iscratch);

/* reducelst.c */
Integer reduce_list (Integer num_list, Integer *list, Integer *scratch);
Integer reducelst (Integer num_list, Integer *list, Integer me, Integer *scratch);

/* resid.c */
void resid (Integer *n, DoublePrecision **colA, Integer *mapA, Integer *m, DoublePrecision **colZ, Integer *mapZ, DoublePrecision *eval, DoublePrecision **ibuffptr, Integer *iwork, DoublePrecision *work, DoublePrecision *res, Integer *info);

/* residtst.c */

/* residual.c */
void residual (Integer *n, DoublePrecision **colA, Integer *mapA, DoublePrecision **colB, Integer *mapB, Integer *m, DoublePrecision **colZ, Integer *mapZ, DoublePrecision *eval, DoublePrecision **ibuffptr, Integer *iwork, DoublePrecision *work, DoublePrecision *res, Integer *info);

/* sclmatrix.c */
Integer sclmatrix (Integer *n, DoublePrecision *const, DoublePrecision **vecA, Integer *mapA);

/* shellsort.c */
void gshellsort_ (Integer *n, Integer v[]);

/* soluf.c */
void upperUF_ (Integer *n, DoublePrecision **rowU, Integer *mapU, Integer *nW, DoublePrecision **colW, Integer *mapW, Integer *iwork, DoublePrecision *work);

/* solul.c */
void solul (Integer *n, Integer *mapU, DoublePrecision **rowU, Integer *nW, Integer *mapW, DoublePrecision **colW, Integer *iwork, DoublePrecision *work);

/* sonenorm.c */
     void sonenrm (Integer *n, DoublePrecision **colA, Integer *mapA, DoublePrecision *norm, Integer *iwork, DoublePrecision *work, Integer *info);

/* soort.c */
Integer qqsort(Integer *v, Integer left, Integer right);

/* sort.c */

void sort_ (Integer *m, Integer *n, Integer *nsplit, Integer *isplit, Integer *iblock, Integer *iwork, DoublePrecision *eval, DoublePrecision *work);
     
     /* tred22.c */
     void bbcast00 (char *buf, Integer len, Integer type, Integer root, Integer snumprocs, Integer *plist);
void gsum00 (char *buf, Integer items, Integer datatype, Integer msgtype, Integer root, Integer snumprocs, Integer *plist, DoublePrecision *work);
void gsum01 (char *buf, Integer items, Integer datatype, Integer msgtype, Integer root, Integer snumprocs, Integer *plist, DoublePrecision *work);
Integer tred2 (Integer *n, DoublePrecision **vecA, Integer *mapA, DoublePrecision **Q, Integer *mapQ, DoublePrecision *diag, DoublePrecision *upperdiag, Integer *iwork, DoublePrecision *work);

/* treesort1.c */
void gtreeinsert (Integer p, Integer *list, Integer *iscrat);
Integer gtree (Integer p, Integer w, Integer *count, Integer *iscrat);
Integer reduce_list3 (Integer num_list, Integer *list, Integer *scratch);
Integer reduce_list33 (Integer num_list, Integer *list, Integer *scratch, Integer node);
Integer reduce_list2 (Integer num_list, Integer *list, Integer *scratch);
Integer reduce_list22 (Integer num_list, Integer *list, Integer *scratch, Integer node);

/* tresid.c */
void tresid (Integer *n, Integer *m, DoublePrecision *d, DoublePrecision *e, DoublePrecision **colZ, Integer *mapZ, DoublePrecision *eval, Integer *iwork, DoublePrecision *work, DoublePrecision *res, Integer *info);

/* upperxfull.c */
Integer upper_x_full (Integer *n, DoublePrecision **rowU, Integer *mapU, Integer *m, DoublePrecision **colF, Integer *mapF, Integer *iscratch, DoublePrecision *scratch);

/* util.c */
Integer count_list (Integer me, Integer *list, Integer *size);
Integer indxL (Integer k, Integer nvecs, Integer map[]);
Integer indxlf_ (Integer *k, Integer *nvecs, Integer *map);
Integer indaint (Integer k, Integer nvecs, Integer *map);
void fil_int_lst (Integer n, Integer *list, Integer item);
void fil_dbl_lst (Integer n, DoublePrecision *list, DoublePrecision item);
Integer in_list (Integer *item, Integer *list, Integer *list_len);
Integer find_proc_store (Integer indx, Integer nprocs, Integer size, Integer *sizelist);
Integer find_large_store (Integer indx, Integer nprocs, Integer *size, Integer *sizelist);
void mem_cpy (Integer *list1, Integer *list2, Integer n);

/* bortho_f.c */
void borthof_ (Integer *n, DoublePrecision *matB, Integer *mapB, Integer *m, DoublePrecision *matZ, Integer *mapZ, Integer *iwork, DoublePrecision *work, DoublePrecision *ort, Integer *info);

/* choleski_f.c */
void choleskif_ (Integer *n, DoublePrecision *colQ, Integer *mapQ, Integer *iwork, DoublePrecision *work, Integer *info);

/* inverse_f.c */
void inverself_ (Integer *msize, Integer *map, DoublePrecision *matrix, Integer *iwork, DoublePrecision *work, Integer *info);

/* mxm25_f.c */
void mxm25f_ (Integer *n1, Integer *n2, DoublePrecision *rowQ, Integer *mapQ, Integer *m, DoublePrecision *colW, Integer *mapW, DoublePrecision *colZ, Integer *iwork, DoublePrecision *work);

/* mxm2_f.c */
void mxm2_ (Integer *n, DoublePrecision *rowQ, Integer *mapQ, Integer *m, DoublePrecision *colW, Integer *mapW, Integer *iwork, DoublePrecision *work);

/* mxm5x_f.c */
void mxm5xf_ (Integer *n, DoublePrecision *matA, Integer *mapA, Integer *m, DoublePrecision *matB, Integer *mapB, Integer *iwork, DoublePrecision *work);

/* mxm88_f.c */
void mxm88f_ (Integer *n, DoublePrecision *matA, Integer *mapA, Integer *m, DoublePrecision *matB, Integer *mapB, Integer *iwork, DoublePrecision *work);

/* onenorm_f.c */
void one_nrmf_ (Integer *n, Integer *m, DoublePrecision *matA, Integer *mapA, DoublePrecision *norm, Integer *iwork, DoublePrecision *work);

/* ortho_f.c */
void orthof_ (Integer *n, Integer *m, DoublePrecision *matZ, Integer *mapZ, Integer *iwork, DoublePrecision *work, DoublePrecision *ort, Integer *info);

/* pdspev_f.c */
void pdspevf_ (Integer *n, DoublePrecision *matrixA, Integer *mapA, DoublePrecision *matZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision *dblptr, Integer *ibuffsz, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspevx_f.c */
void pdspevxf_ (Integer *ivector, Integer *irange, Integer *n, DoublePrecision *matrixA, Integer *mapA, DoublePrecision *lb, DoublePrecision *ub, Integer *ilb, Integer *iub, DoublePrecision *abstol, Integer *meigval, DoublePrecision *matZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision *dblptr, Integer *ibuffsz, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspgv_f.c */
void pdspgvf_ (Integer *ifact, Integer *n, DoublePrecision *matrixA, Integer *mapA, DoublePrecision *matrixB, Integer *mapB, DoublePrecision *matZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision *dblptr, Integer *ibuffsz, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdspgvx_f.c */
void pdspgvxf_ (Integer *ifact, Integer *ivector, Integer *irange, Integer *n, DoublePrecision *matrixA, Integer *mapA, DoublePrecision *matrixB, Integer *mapB, DoublePrecision *lb, DoublePrecision *ub, Integer *ilb, Integer *iub, DoublePrecision *abstol, Integer *meigval, DoublePrecision *matZ, Integer *mapZ, DoublePrecision *eval, Integer *iscratch, Integer *iscsize, DoublePrecision *dblptr, Integer *ibuffsz, DoublePrecision *scratch, Integer *ssize, Integer *info);

/* pdsptri_f.c */


/* resid_f.c */
void residf_ (Integer *n, DoublePrecision *matrixA, Integer *mapA, Integer *m, DoublePrecision *matrixZ, Integer *mapZ, DoublePrecision *eval, Integer *iwork, DoublePrecision *work, DoublePrecision *res, Integer *info);

/* residual_f.c */
void residualf_ (Integer *n, DoublePrecision *matrixA, Integer *mapA, DoublePrecision *matrixB, Integer *mapB, Integer *m, DoublePrecision *matrixZ, Integer *mapZ, DoublePrecision *eval, Integer *iwork, DoublePrecision *work, DoublePrecision *res, Integer *info);

void sonenrmf_ (Integer *n, DoublePrecision *matrixA, Integer *mapA, DoublePrecision *norm, Integer *iwork, DoublePrecision *work, Integer *info);

/* tresid_f.c */
void tresidf_ (Integer *n, Integer *m, DoublePrecision *d, DoublePrecision *e, DoublePrecision *matrixZ, Integer *mapZ, DoublePrecision *eval, Integer *iwork, DoublePrecision *work, DoublePrecision *res, Integer *info);


/*
   fortran stuff
*/

/*
   blas 1
*/

void dcopy_(Integer *, DoublePrecision *, Integer *, DoublePrecision *, Integer *);
void daxpy_(Integer *, DoublePrecision *, DoublePrecision *, Integer *, DoublePrecision *, Integer *);

void xerbla_ (char *ptr, Integer *m);
void mxpend_();
Integer mxmynd_();
Integer mxnprc_();

/*
void mxcombv1_(char *buf, Integer *sum(), Integer *isize, Integer *items, Integer *n, Integer *list, Integer *msgtype, char *work);
*/
