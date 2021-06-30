//=============================================================================
//Based on AIList (John Feng, 2018)
//by Selena Feng  12/05/2019--1/25/2020
//-----------------------------------------------------------------------------
#include "seqpare.h"
#define PROGRAM_NAME  "seqpare"
#define MAJOR_VERSION "0"
#define MINOR_VERSION "1"
#define REVISION_VERSION "1"
#define BUILD_VERSION "0"
#define MODIFY_VERSION "XinD"
#define VERSION MAJOR_VERSION "." MINOR_VERSION "." REVISION_VERSION "." MODIFY_VERSION

int seqpare_help(int exit_code);

int main(int argc, char **argv)
{
	if (argc < 3) return seqpare_help(0);
	int i, j, m=0;   
    char out[1024] = "out";
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-m")==0 && i+1<argc){
			m = atoi(argv[i+1]);
		}
		if(strcmp(argv[i], "-o")==0 && i+1<argc){
			strcpy(out,argv[i+1]);
		}		
	}
	
	clock_t start, end1;  
    start = clock();
	
	if(m==0 && argc==3){
		//compare two interval sets
		seqpare_t *sm = malloc(1*sizeof(seqpare_t));
		seq_compare_11(argv[1], argv[2], sm);	   
		printf("%i\t %i\t %12.8f\t %lld\n", sm->N1, sm->N2, sm->teo, (long long)sm->tc);
		free(sm); 
	}
	else if(m==1 && argc>=5){
		//searching 
		//1. Get the files
		int nf = 0, n1, n2;  
		glob_t gResult;
		int rtn = glob(argv[1], 0, NULL, &gResult);     
		if(rtn!=0){
		    printf("wrong dir path: %s", argv[1]);
		    return 0;
		}
		char** file_ids = gResult.gl_pathv;
		nf = gResult.gl_pathc; 
		if(nf<1){   
		    printf("No files (add to path /*): %i\n", nf);
		    return 0;
		}  
		seqpare_t *sm = malloc(nf*sizeof(seqpare_t));
		
		//2. Calculate the s index
		ailist0_t *ail2 = readBED0(argv[2]);			//as query
		ailist0_construct(ail2, 20);
		for(int i=0; i<nf; i++){
			ailist_t *ail1 = readBED(file_ids[i]);		//as database: fresh q_idx, s_max
			ailist_construct(ail1, 20);
			seq_compare(ail1, ail2, &sm[i]);
			ailist_destroy(ail1);
		}
		ailist0_destroy(ail2);   		
		//-----------------------------------------------------------
		if(strlen(out)>1){
		    FILE *fp = fopen(out, "w");
		    if(fp==NULL)
		        printf("Can't create file %s\n", out);
		    else{
		        fprintf(fp, "%s\t%s\t%s\t%s\t%s\n", "N1", "N2", "teo", "tc", "sm");
		        for(i=0;i<nf;i++){
		        	n1 = sm[i].N1;
		        	n2 = sm[i].N2;
		        	double smi = sm[i].teo/((double)(n1+n2)-sm[i].teo);
		        	fprintf(fp, "%i\t%i\t%12.8f\t%lld\t%12.8f\t%s\n", n1, n2, sm[i].teo, (long long)sm[i].tc, smi, file_ids[i]);
				} 
		        fclose(fp);
		    }     
		}		
		globfree(&gResult);		
		free(sm);	
	}
	else if(m==2 && argc>=5){
		//mapping of two collections
		seqpare_t **sm = NULL;
		int32_t nf, mf;
		seq_compare_mn(argv[1], argv[2], &nf, &mf, &sm);
		
		if(strlen(out)>1){
		    FILE *fp = fopen(out, "w");
		    if(fp==NULL)
		        printf("Can't create file %s\n", out);
		    else{
		        fprintf(fp, "%i\t %i\n", nf, mf);
		        for(i=0;i<nf;i++){
		            for(j=0;j<mf;j++)
		                fprintf(fp, "%12.8f\t", sm[i][j].teo); 
		            fprintf(fp, "\n");
		        } 
		        fclose(fp);
		        //------------------------------------------
		        fp = fopen(strcat(out, ".sm"), "w");
		        fprintf(fp, "%i\t %i\n", nf, mf);
		        for(i=0;i<nf;i++){
		            for(j=0;j<mf;j++)
		                fprintf(fp, "%12.8f\t", sm[i][j].teo/((float)(sm[i][j].N1+sm[i][j].N2)-sm[i][j].teo)); 
		            fprintf(fp, "\n");
		        } 
		        fclose(fp);
		    }     
		}		
		
		if(sm!=NULL){
			free(sm);
		}	
	}
	else if(m==3 && argc>=4){
		//self mapping
		seqpare_t **sm = NULL;
		int32_t nf, mf;
		seq_compare_nn(argv[1], &nf, &sm);
		
		if(strlen(out)>1){
		    FILE *fp = fopen(out, "w");
		    if(fp==NULL)
		        printf("Can't create file %s\n", out);
		    else{
		        fprintf(fp, "%i\t %i\n", nf, nf);
		        for(j=0;j<nf;j++){	        	
		            for(i=0;i<nf;i++)
		                fprintf(fp, "%12.8f\t", sm[i][j].teo);	                
		            fprintf(fp, "\n");
		        } 
		        fclose(fp);
		        //----------------------------------------------
		        fp = fopen(strcat(out, ".sm"), "w");
		        fprintf(fp, "%i\t %i\n", nf, nf);
		        for(j=0;j<nf;j++){
		        	for(i=0;i<nf;i++)
		            	fprintf(fp, "%12.8f\t", sm[j][i].teo/((float)(sm[j][i].N1+sm[j][i].N2)-sm[j][i].teo));
		            fprintf(fp, "\n");
		        } 
		        fclose(fp);		        
		    }     
		}
		
		if(sm!=NULL){
			free(sm);
		}		
	}
	else 
		return seqpare_help(0);

    end1 = clock();    
    printf("execution time: %f\n", ((double)(end1-start))/CLOCKS_PER_SEC);         
  
    return 0;
}

int seqpare_help(int exit_code)
{
    fprintf(stderr,
"%s, v%s\n"
"usage: %s <seq-file 1 (.bed or .bed.gz)> <seq-file 2> \n"
"            <path to seqs folder> <query seq file> -m 1 \n"
"            <path to seqs folder 1> <path to seqs folder 2> -m 2 \n"
"            <path to seqs folder> -m 3 \n"
"         options:\n"
"             -o <output file>\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

/*
		//searching 
		seqpare_t *sm = NULL;
		int32_t nf, n1, n2;
		seq_compare_1n(argv[1], argv[2], &nf, &sm);
		if(strlen(out)>1){
		    FILE *fp = fopen(out, "w");
		    if(fp==NULL)
		        printf("Can't create file %s\n", out);
		    else{
		        fprintf(fp, "%s\t%s\t%s\t%s\t%s\n", "N1", "N2", "teo", "tc", "sm");
		        for(i=0;i<nf;i++){
		        	n1 = sm[i].N1;
		        	n2 = sm[i].N2;
		        	float smi = sm[i].teo/((float)(n1+n2)-sm[i].teo);
		        	fprintf(fp, "%i\t%i\t%12.8f\t%lld\t%12.8f\n", n1, n2, sm[i].teo, (long long)sm[i].tc, smi);
				} 
		        fclose(fp);
		    }     
		}		
		if(sm!=NULL)
			free(sm);	
*/
