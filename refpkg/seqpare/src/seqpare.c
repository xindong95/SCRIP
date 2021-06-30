//=============================================================================
//Based on AIList (John Feng, 2018)
//by Selena Feng  12/05/2019--1/25/2020
//-----------------------------------------------------------------------------
#include "seqpare.h"
#define gdata_t_key(r) ((r).start)
#define gdata0_t_key(r) ((r).start)
KRADIX_SORT_INIT(intv, gdata_t, gdata_t_key, 4)
KRADIX_SORT_INIT(intv0, gdata0_t, gdata0_t_key, 4)
KHASH_MAP_INIT_STR(str, int32_t)
typedef khash_t(str) strhash_t;

char *parse_bed(char *s, int32_t *st_, int32_t *en_)
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atol(p);
			else if (i == 2) en = atol(p);
			++i, p = q + 1;
			if (c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 3? ctg : 0;
}

uint32_t bSearch(gdata_t* As, uint32_t idxS, uint32_t idxE, uint32_t qe)
{   //find tE: index of the first item satisfying .s<qe from right
    int tL=idxS, tR=idxE-1, tM, tE=-1;
    if(As[tR].start < qe)
        return tR;
    else if(As[tL].start >= qe)
        return -1;
    while(tL<tR-1){
        tM = (tL+tR)/2; 
        if(As[tM].start >= qe)
            tR = tM-1;
        else
            tL = tM;
    }
    if(As[tR].start < qe)
        tE = tR;
    else if(As[tL].start < qe)
        tE = tL;       
    return tE; 
}

ailist_t *ailist_init(void)
{
	ailist_t *ail = malloc(1*sizeof(ailist_t));
	ail->hc = kh_init(str);
	ail->nctg = 0;
	ail->mctg = 32;
	ail->ctg = malloc(ail->mctg*sizeof(ctg_t));
	return ail;
}

ailist0_t *ailist0_init(void)
{
	ailist0_t *ail = malloc(1*sizeof(ailist0_t));
	ail->hc = kh_init(str);
	ail->nctg = 0;
	ail->mctg = 32;
	ail->ctg = malloc(ail->mctg*sizeof(ctg0_t));
	return ail;
}

void ailist_destroy(ailist_t *ail)
{
	int32_t i;
	if (ail == 0) return;
	for (i = 0; i < ail->nctg; ++i){
		free(ail->ctg[i].name);
		free(ail->ctg[i].glist);
		free(ail->ctg[i].maxE);
	}
	free(ail->ctg);
	kh_destroy(str, (strhash_t*)ail->hc);
	free(ail);
}

void ailist0_destroy(ailist0_t *ail)
{
	int32_t i;
	if (ail == 0) return;
	for (i = 0; i < ail->nctg; ++i){
		free(ail->ctg[i].name);
		free(ail->ctg[i].glist);
		free(ail->ctg[i].maxE);
	}
	free(ail->ctg);
	kh_destroy(str, (strhash_t*)ail->hc);
	free(ail);
}

void ailist_add(ailist_t *ail, const char *chr, uint32_t s, uint32_t e, int32_t v)
{
	if(s > e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_put(str, h, chr, &absent);
	if (absent) {
		if (ail->nctg == ail->mctg)
			EXPAND(ail->ctg, ail->mctg);							
		kh_val(h, k) = ail->nctg;		
		ctg_t *p = &ail->ctg[ail->nctg++];
		p->name = strdup(chr);
		p->nr=0;	p->mr=64;
		p->glist = malloc(p->mr*sizeof(gdata_t));
		kh_key(h, k) = p->name;
	}
	int32_t kk = kh_val(h, k);
	ctg_t *q = &ail->ctg[kk];
	if (q->nr == q->mr)
		EXPAND(q->glist, q->mr);	
	gdata_t *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	p->q_idx = -1;	
	p->s_max = -1.0;	//not taken
	p->cnts = 0;
	return;
}

void ailist0_add(ailist0_t *ail, const char *chr, uint32_t s, uint32_t e, int32_t v)
{
	if(s > e)return;
	int absent;
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_put(str, h, chr, &absent);
	if (absent) {
		if (ail->nctg == ail->mctg)
			EXPAND(ail->ctg, ail->mctg);							
		kh_val(h, k) = ail->nctg;		
		ctg0_t *p = &ail->ctg[ail->nctg++];
		p->name = strdup(chr);
		p->nr=0;	p->mr=64;
		p->glist = malloc(p->mr*sizeof(gdata_t));
		kh_key(h, k) = p->name;
	}
	int32_t kk = kh_val(h, k);
	ctg0_t *q = &ail->ctg[kk];
	if (q->nr == q->mr)
		EXPAND(q->glist, q->mr);	
	gdata0_t *p = &q->glist[q->nr++];
	p->start = s;
	p->end   = e;
	return;
}

//-------------------------------------------------------------------------------
ailist_t* readBED(const char* fn)
{   //faster than strtok()
	gzFile fp;
	ailist_t *ail;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t k = 0;
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	ks = ks_init(fp);
	ail = ailist_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *ctg;
		int32_t st, en;
		ctg = parse_bed(str.s, &st, &en);
		if (ctg) ailist_add(ail, ctg, st, en, k++);
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return ail;
} 

ailist0_t* readBED0(const char* fn)
{   //faster than strtok()
	gzFile fp;
	ailist0_t *ail;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t k = 0;
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	ks = ks_init(fp);
	ail = ailist0_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *ctg;
		int32_t st, en;
		ctg = parse_bed(str.s, &st, &en);
		if (ctg) ailist0_add(ail, ctg, st, en, k++);
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return ail;
} 

void ailist_construct(ailist_t *ail, int cLen)
{   
    int cLen1=cLen/2, j1, nr, minL = MAX(64, cLen);     
    cLen += cLen1;      
    int lenT, len, iter, i, j, k, k0, t;            	
	for(i=0; i<ail->nctg; i++){
		//1. Decomposition
		ctg_t *p    = &ail->ctg[i];
		gdata_t *L1 = p->glist;							//L1: to be rebuilt
		nr 			= p->nr;
		radix_sort_intv(L1, L1+nr);                 		               
        if(nr<=minL){        
            p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;                
        }
        else{         
        	gdata_t *L0 = malloc(nr*sizeof(gdata_t)); 	//L0: serve as input list
            gdata_t *L2 = malloc(nr*sizeof(gdata_t));   //L2: extracted list 
            //----------------------------------------
        	gdata_t *D0 = malloc(nr*sizeof(gdata_t)); 	//D0:            
			int32_t *di = malloc(nr*sizeof(int32_t));	//int64_t?			  
            //----------------------------------------
            memcpy(L0, L1, nr*sizeof(gdata_t));			
            iter = 0;	k = 0;	k0 = 0;
            lenT = nr;
            while(iter<MAXC && lenT>minL){  
            	//setup di---------------------------			
		        for(j=0;j<lenT;j++){				//L0:{.start= end, .end=idx, .value=idx1}
					D0[j].start = L0[j].end;
					D0[j].end = j;
				}
				radix_sort_intv(D0, D0+lenT);
				for(j=0;j<lenT;j++){				//assign i=29 to L0[i].end=2
					t = D0[j].end;
					di[t] = j-t;					//>0 indicate containment
				}  
				//----------------------------------- 
                len = 0;
		        for(t=0;t<lenT-cLen;t++){
					if(di[t]>cLen)
				        memcpy(&L2[len++], &L0[t], sizeof(gdata_t));    			
					else
						memcpy(&L1[k++], &L0[t], sizeof(gdata_t)); 
				}             
                memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata_t));   
                k += cLen, lenT = len;                
                p->idxC[iter] = k0;
                p->lenC[iter] = k-k0;
                k0 = k, iter++;
                if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                    if(lenT>0){
                        memcpy(&L1[k], L2, lenT*sizeof(gdata_t));
                        p->idxC[iter] = k;
                        p->lenC[iter] = lenT;
                        iter++;
                        lenT = 0;						//exit!
                    }
                   	p->nc = iter;                   
                }
                else memcpy(L0, L2, lenT*sizeof(gdata_t));
            }
            free(L2),free(L0), free(D0), free(di);   
        }
        //2. Augmentation
        p->maxE = malloc(nr*sizeof(uint32_t)); 
        for(j=0; j<p->nc; j++){ 
            k0 = p->idxC[j];
            k = k0 + p->lenC[j];
            uint32_t tt = L1[k0].end;
            p->maxE[k0]=tt;
            for(t=k0+1; t<k; t++){
                if(L1[t].end > tt) tt = L1[t].end;
                p->maxE[t] = tt;  
            }             
        } 
	}	
}

void ailist0_construct(ailist0_t *ail, int cLen)
{   
    int cLen1=cLen/2, j1, nr, minL = MAX(64, cLen);     
    cLen += cLen1;      
    int lenT, len, iter, i, j, k, k0, t;            	
	for(i=0; i<ail->nctg; i++){
		//1. Decomposition
		ctg0_t *p    = &ail->ctg[i];
		gdata0_t *L1 = p->glist;							//L1: to be rebuilt
		nr 			= p->nr;
		radix_sort_intv0(L1, L1+nr);                 		               
        if(nr<=minL){        
            p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;                
        }
        else{         
        	gdata0_t *L0 = malloc(nr*sizeof(gdata0_t)); 	//L0: serve as input list
            gdata0_t *L2 = malloc(nr*sizeof(gdata0_t));   //L2: extracted list 
            //----------------------------------------
        	gdata0_t *D0 = malloc(nr*sizeof(gdata0_t)); 	//D0:            
			int32_t *di = malloc(nr*sizeof(int32_t));	//int64_t?			  
            //----------------------------------------
            memcpy(L0, L1, nr*sizeof(gdata0_t));			
            iter = 0;	k = 0;	k0 = 0;
            lenT = nr;
            while(iter<MAXC && lenT>minL){  
            	//setup di---------------------------			
		        for(j=0;j<lenT;j++){				//L0:{.start= end, .end=idx, .value=idx1}
					D0[j].start = L0[j].end;
					D0[j].end = j;
				}
				radix_sort_intv0(D0, D0+lenT);
				for(j=0;j<lenT;j++){				//assign i=29 to L0[i].end=2
					t = D0[j].end;
					di[t] = j-t;					//>0 indicate containment
				}  
				//----------------------------------- 
                len = 0;
		        for(t=0;t<lenT-cLen;t++){
					if(di[t]>cLen)
				        memcpy(&L2[len++], &L0[t], sizeof(gdata0_t));    			
					else
						memcpy(&L1[k++], &L0[t], sizeof(gdata0_t)); 
				}             
                memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata0_t));   
                k += cLen, lenT = len;                
                p->idxC[iter] = k0;
                p->lenC[iter] = k-k0;
                k0 = k, iter++;
                if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                    if(lenT>0){
                        memcpy(&L1[k], L2, lenT*sizeof(gdata0_t));
                        p->idxC[iter] = k;
                        p->lenC[iter] = lenT;
                        iter++;
                        lenT = 0;						//exit!
                    }
                   	p->nc = iter;                   
                }
                else memcpy(L0, L2, lenT*sizeof(gdata0_t));
            }
            free(L2),free(L0), free(D0), free(di);   
        }
        //2. Augmentation
        p->maxE = malloc(nr*sizeof(uint32_t)); 
        for(j=0; j<p->nc; j++){ 
            k0 = p->idxC[j];
            k = k0 + p->lenC[j];
            uint32_t tt = L1[k0].end;
            p->maxE[k0]=tt;
            for(t=k0+1; t<k; t++){
                if(L1[t].end > tt) tt = L1[t].end;
                p->maxE[t] = tt;  
            }             
        } 
	}	
}

void ailist_construct0(ailist_t *ail, int cLen)
{   //New continueous memory?   
    int cLen1=cLen/2, j1, nr, minL = MAX(64, cLen);     
    cLen += cLen1;      
    int lenT, len, iter, i, j, k, k0, t;            	
	for(i=0; i<ail->nctg; i++){
		//1. Decomposition
		ctg_t *p    = &ail->ctg[i];
		gdata_t *L1 = p->glist;							//L1: to be rebuilt
		nr 			= p->nr;
		radix_sort_intv(L1, L1+nr);                 		               
        if(nr<=minL){        
            p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;                
        }
        else{         
        	gdata_t *L0 = malloc(nr*sizeof(gdata_t)); 	//L0: serve as input list
            gdata_t *L2 = malloc(nr*sizeof(gdata_t));   //L2: extracted list 
            memcpy(L0, L1, nr*sizeof(gdata_t));			
            iter = 0;	k = 0;	k0 = 0;
            lenT = nr;
            while(iter<MAXC && lenT>minL){   
                len = 0;            
                for(t=0; t<lenT-cLen; t++){
                    uint32_t tt = L0[t].end;
                    j=1;    j1=1;
                    while(j<cLen && j1<cLen1){
                        if(L0[j+t].end>=tt) j1++;
                        j++;
                    }
                    if(j1<cLen1) memcpy(&L2[len++], &L0[t], sizeof(gdata_t));
                    else memcpy(&L1[k++], &L0[t], sizeof(gdata_t));                 
                } 
                memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata_t));   
                k += cLen, lenT = len;                
                p->idxC[iter] = k0;
                p->lenC[iter] = k-k0;
                k0 = k, iter++;
                if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                    if(lenT>0){
                        memcpy(&L1[k], L2, lenT*sizeof(gdata_t));
                        p->idxC[iter] = k;
                        p->lenC[iter] = lenT;
                        iter++;
                        lenT = 0;	//exit!
                    }
                   	p->nc = iter;                   
                }
                else memcpy(L0, L2, lenT*sizeof(gdata_t));
            }
            free(L2),free(L0);     
        }
        //2. Augmentation
        p->maxE = malloc(nr*sizeof(uint32_t)); 
        for(j=0; j<p->nc; j++){ 
            k0 = p->idxC[j];
            k = k0 + p->lenC[j];
            uint32_t tt = L1[k0].end;
            p->maxE[k0]=tt;
            for(t=k0+1; t<k; t++){
                if(L1[t].end > tt) tt = L1[t].end;
                p->maxE[t] = tt;  
            }             
        } 
	}
}

int32_t get_ctg(const ailist_t *ail, const char *chr)
{
	khint_t k;
	strhash_t *h = (strhash_t*)ail->hc;
	k = kh_get(str, h, chr);
	return k == kh_end(h)? -1 : kh_val(h, k);
}


uint32_t ailist_query(ailist_t *ail, char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir)
{   
    uint32_t nr = 0, m = *mr, *r = *ir;
    int32_t gid = get_ctg(ail, chr);
    if(gid>=ail->nctg || gid<0)return 0;
    ctg_t *p = &ail->ctg[gid];	
    for(int k=0; k<p->nc; k++){					//search each component
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];			
        int32_t t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe); 	//rs<qe: inline not better 
            if(t>=cs){
		        if(nr+t-cs>=m){
		        	m = nr+t-cs + 1024;
		        	r = realloc(r, m*sizeof(uint32_t));
		        }
		        while(t>=cs && p->maxE[t]>qs){
		            if(p->glist[t].end>qs)               	
		                r[nr++] = t;
		            t--;
		        }
            }
        }
        else{
        	if(nr+ce-cs>=m){
        		m = nr+ce-cs + 1024;
        		r = realloc(r, m*sizeof(uint32_t));
        	}
            for(t=cs; t<ce; t++)
                if(p->glist[t].start<qe && p->glist[t].end>qs)
                    r[nr++] = t;                           
        }
    }    
    *ir = r, *mr = m;                  
    return nr;                              
}

void seq_compare_nn(char *iPath, int32_t *nf, seqpare_t ***sps)
{   
    //1. Get the files  
    glob_t gResult;
    int rtn = glob(iPath, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath);
        return;
    }
    char** file_ids = gResult.gl_pathv;
    *nf = gResult.gl_pathc; 
    int i, j, n = *nf;
    if(n<1)   
        printf("Too few files (add to path /*): %i\n", n);  
    seqpare_t **sm = malloc(n*sizeof(seqpare_t *));
    for(i=0;i<n; i++)
    	sm[i] = malloc(n*sizeof(seqpare_t));
    	
	//2. Calculate the s index
	for(j=0;j<n;j++){
	    printf("%i\t %s\n", j, file_ids[j]);
		ailist0_t *ail2 = readBED0(file_ids[j]);	//as query
		ailist0_construct(ail2, 20);
		//-------------------------
		for(i=j;i<n;i++){
			ailist_t *ail1 = readBED(file_ids[i]);	//database: refresh
			ailist_construct(ail1, 20);
			seq_compare(ail1, ail2, &sm[j][i]);
			ailist_destroy(ail1);
			if(i!=j)sm[i][j] = sm[j][i];
		}
		//-------------------------
		ailist0_destroy(ail2); 
	}    
	//printf("nn: %i\t %12.8f\n", n, sm[0][1]);  
	globfree(&gResult);
	*sps = sm;
}

void seq_compare_mn(char *iPath1, char *iPath2, int32_t *nf, int32_t *mf, seqpare_t ***sps)
{   
   //1. Get the files  
    glob_t gResult1;
    int rtn = glob(iPath1, 0, NULL, &gResult1);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath1);
        return;
    }
    char** file_ids1 = gResult1.gl_pathv;
    *nf = gResult1.gl_pathc;//as query
    //-------------------------
    glob_t gResult2;
    rtn = glob(iPath2, 0, NULL, &gResult2);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath2);
        return;
    }
    char** file_ids2 = gResult2.gl_pathv;
    *mf = gResult2.gl_pathc;//as database    
    //-------------------------
    int i, j, n = *nf, m = *mf;
    if(n<1 || m<1){   
        printf("Too few files (add to path /*): %i\t %i\n", n, m);  
        return;
    }
    
    seqpare_t **sm = malloc(n*sizeof(seqpare_t *));
    for(i=0;i<n; i++)
    	sm[i] = malloc(m*sizeof(seqpare_t));
    
	//2. Calculate the s index
	for(j=0;j<n;j++){
		ailist0_t *ail2 = readBED0(file_ids1[j]);	//as query
		ailist0_construct(ail2, 20);
		printf("%i\t %s\n", j, file_ids1[j]);
		//-------------------------
		for(i=0;i<m;i++){
			ailist_t *ail1 = readBED(file_ids2[i]);	//as database: fresh
			ailist_construct(ail1, 20);
			seq_compare(ail1, ail2, &sm[j][i]);
			ailist_destroy(ail1);
		}
		//-------------------------
		ailist0_destroy(ail2); 
	}  
	globfree(&gResult1);
	globfree(&gResult2);	
	*sps = sm;
}

void seq_compare_1n(char *iPath, char *qfile, int32_t *nf, seqpare_t **sps)
{   
    //1. Get the files  
    glob_t gResult;
    int rtn = glob(iPath, 0, NULL, &gResult);     
    if(rtn!=0){
        printf("wrong dir path: %s", iPath);
        return;
    }
    char** file_ids = gResult.gl_pathv;
    *nf = gResult.gl_pathc; 
    if(*nf<1)   
        printf("Too few files (add to path /*): %i\n", *nf);  
    seqpare_t *sm = malloc(*nf*sizeof(seqpare_t));
    
	//2. Calculate the s index
	ailist0_t *ail2 = readBED0(qfile);			//as query
	ailist0_construct(ail2, 20);
	//-------------------------
	for(int i=0;i<*nf;i++){
		ailist_t *ail1 = readBED(file_ids[i]);	//as database: fresh q_idx, s_max
		ailist_construct(ail1, 20);
		seq_compare(ail1, ail2, &sm[i]);
		ailist_destroy(ail1);
	}
	//-------------------------
	ailist0_destroy(ail2);   
	globfree(&gResult);
	*sps = sm;
}

void seq_compare_11(char *file1, char *file2, seqpare_t *sp)
{ 
	ailist_t *ail1 = readBED(file1);	//as Database
	ailist_construct(ail1, 20);
	//---------------------------------
	ailist0_t *ail2 = readBED0(file2);	//file2 as query
	ailist0_construct(ail2, 20);
	//---------------------------------
	seq_compare(ail1, ail2, sp);
	//printf("s index = %12.8f\n", *sm);
	ailist_destroy(ail1);
	ailist0_destroy(ail2);
}

void seq_compare(ailist_t *ail1, ailist0_t *ail2, seqpare_t *sp)
{   //need to refresh ail1: database--q_idx, s_max
	float ERR = 0.0000001;
	int32_t i, j, k;
	int32_t N1=0, N2=0;
	for(i=0;i<ail1->nctg;i++)
		N1 += (&ail1->ctg[i])->nr;
	for(i=0;i<ail2->nctg;i++)
		N2 += (&ail2->ctg[i])->nr;	
	sp->N1 = N1;
	sp->N2 = N2;		
	//-------------------------------------------------------------------------
	int32_t qs, qe, rs, re, gid, cs, ce, t, iq, tmax;
	float qlen, s, rlen, smax, st;
	ctg_t *p1;
	ctg0_t *p2;
	for(j=0;j<ail2->nctg;j++){						//ail2: query	
		p2 = &ail2->ctg[j];
		gid = get_ctg(ail1, p2->name);
		if(gid>=ail1->nctg || gid<0)continue;			//get gid in ail1: db					
		p1 = &ail1->ctg[gid];
		for(i=0;i<p2->nr;i++){	
			qs = p2->glist[i].start;
			qe = p2->glist[i].end;				    	
			qlen = qe-qs; 
			tmax=-1;
			smax=0.0;		
			//----------------------------------------
			for(k=0; k<p1->nc; k++){						//search each component
				cs = p1->idxC[k];
				ce = cs + p1->lenC[k];			
			    t = bSearch(p1->glist, cs, ce, qe); 		//rs<qe: inline not better 
				if(t>=cs){
					while(t>=cs && p1->maxE[t]>qs){
					    if((re=p1->glist[t].end)>qs){ 
			    			p1->glist[t].cnts++;					     
					    	rs = p1->glist[t].start;             	
					        s = MIN(qe, re)-MAX(qs, rs);
					        rlen = re-rs;					//this is necessary
					        s = s/(qlen+rlen-s);
					        //------------------------------skip: matched but smaller
					        if(s>smax){						//not matched or matched & larger
					        	st = p1->glist[t].s_max;
					        	if(st<ERR || (st>ERR && s>st)){
							    	smax = s;
							    	tmax = t;
					        	}
					        }
						}
					    t--;
					}
			    }
			}  
			//if(smax<ERR)continue;                 		 
			if(smax>ERR){
				st = p1->glist[tmax].s_max;
				if(st<ERR){								//glist[tmax] not taken, mark it
					p1->glist[tmax].s_max = smax;
					p1->glist[tmax].q_idx = i;
				}
				else if(st<smax){						//mutual max[i, tmax]/rerun .q_idx			
					iq = p1->glist[tmax].q_idx;
					p1->glist[tmax].s_max = smax;
					p1->glist[tmax].q_idx = i;				
					//----rerun on iq to get max other than [tmax]: will be skipped					
					seqpare0(ail1, ail2, gid, j, iq);	
				}
			}
			//else{									//search i again without [tmax]
			//}
		}// glist[i]
	}// ctg[j]	   
	//-------------------------------------------------------------------------
	//calculate si: similarity index
	float sm = 0.0;
	int64_t cnt = 0;
	for(j=0;j<ail1->nctg;j++){
		p1 = &ail1->ctg[j];
		for(i=0;i<p1->nr;i++){
			cnt += p1->glist[i].cnts;
			if((st=p1->glist[i].s_max)>ERR)
				sm += st;
		}
	}	
	//printf("%i\t %i\t %lld\t %12.8f\n", N1, N2, (long long)cnt, sm);
	sp->teo = sm;    
	sp->tc = cnt;
	return;                   
}

void seqpare0(ailist_t *ail1, ailist0_t *ail2, int gid1, int gid2, int iq)
{	//search iq for max other than ir: may recursive
	float ERR = 0.0000001;
	ctg_t *p1 = &ail1->ctg[gid1];
	ctg0_t *p2 = &ail2->ctg[gid2];			
	//------------------------------------
	int32_t qs = p2->glist[iq].start;
	int32_t qe = p2->glist[iq].end;				    	
	int32_t t, cs, ce, rs, re, iq0=iq, tmax=-1; 
	float smax=0.0, s, st, qlen = qe-qs, rlen;
	for(int k=0; k<p1->nc; k++){						//search each component
		cs = p1->idxC[k];
		ce = cs + p1->lenC[k];			
	    t = bSearch(p1->glist, cs, ce, qe); 			//rs<qe: inline not better 
	    if(t>=cs){
			while(t>=cs && p1->maxE[t]>qs){
			    if((re=p1->glist[t].end)>qs){  
			    	rs = p1->glist[t].start;             	
			        s = MIN(qe, re)-MAX(qs, rs);
			        rlen = re-rs;
			        s = s/(qlen+rlen-s);
		        	st = p1->glist[t].s_max;
		        	if(s>smax && (st<ERR || (st>ERR && s>st))){
						smax = s;
						tmax = t;
		        	}
				}
			    t--;
			}
	    }
	} 
	if(smax<ERR)return;                  		 
	//------------------------------------
	st = p1->glist[tmax].s_max;
	if(st<ERR){								//glist[tmax] not taken, mark it
		p1->glist[tmax].s_max = smax;
		p1->glist[tmax].q_idx = iq0;
	}
	else if(st<smax){						//mutual max[i, tmax]/rerun .q_idx			
		iq = p1->glist[tmax].q_idx;
		p1->glist[tmax].s_max = smax;
		p1->glist[tmax].q_idx = iq0;				
		//----rerun on iq to get max other than [tmax]	
		//printf("recursive: %i\t %i\t %i\n", gid1, gid2, iq);				
		seqpare0(ail1, ail2, gid1, gid2, iq);	
	} 
} 

/*
float seq_compare_11(char *file1, char *file2)
{   
	float sm = 0.0;					//similarity index
	int32_t i, j, k, N1=0, N2=0;
	ailist_t *ail1 = readBED(file1);	//as Database, file2 as query
	ailist_construct(ail1, 20);
	for(i=0;i<ail1->nctg;i++)
		N1 += (ail1->ctg[i])->nr;
	//------------------------------	
	ailist_t *ail2 = readBED(file2);	//as Database, file2 as query
	ailist_construct(ail2, 20);
	for(i=0;i<ail2->nctg;i++)
		N2 += (ail2->ctg[i])->nr;			
	//-------------------------------------------------------------------------
	kstream_t *ks;
	kstring_t str = {0,0,0};
	gzFile fp = gzopen(file2, "r");
	assert(fp);
	ks = ks_init(fp);
	int32_t st0=0, en0=0, qs, qe, re, gid, cs, ce, t;
	float qlen, s, smax;	
	int32_t kmax, tmax;	
	char *ctg;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		N2++;		
		ctg = parse_bed(str.s, &qs, &qe);
		if (ctg == 0) continue;	
		//---------------------
		gid = get_ctg(ail, ctg);
		if(gid>=ail->nctg || gid<0)continue;
		ctg_t *p = &ail->ctg[gid];			
		kmax=-1, tmax=-1;
		smax=0.0;
		for(k=0; k<p->nc; k++){						//search each component
		    cs = p->idxC[k];
		    ce = cs + p->lenC[k];			
	        t = bSearch(p->glist, cs, ce, qe); 		//rs<qe: inline not better 
	        if(t>=cs){
	        	qlen = qe-qs;              			//length of q
			    while(t>=cs && p->maxE[t]>qs){
			        if((re=p->glist[t].end)>qs){  
			        	rs = p->glist[t].start;             	
			            s = MIN(qe, re)-MAX(qs, rs);
			            s = s/(qlen+re-rs-s);
			            if(s>smax){
			            	smax = s;
			            	kmax = k;
			            	tmax = t;
			            }
					}
			        t--;
			    }
	        }
		}                    		 
		//---------------------
		if(glist[tmax].s_max<0.0){					//glist[tmax] not taken
			glist[tmax].s_max = smax;
			st+=smax;
		}
		else if(glist[tmax].s_max>smax){			//replacement: to get the mutual max
		
		}
		else{										//search again without [tmax]
		
		}
	}	  
  	free(str.s);
	gzclose(fp); 
	ks_destroy(ks);   
	ailist_destroy(ail);  	
	//-------------------------------------------------------------------------
    return sm/(N1+N2-sm);                        
}
*/
