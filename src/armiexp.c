/*==============================================================================
armiexp.c : analytical residual Mutual Information for expressions
(C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <gsl/gsl_errno.h>

#include "armiexp.h"

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

/*____________________________________________________________________________*/
/* random number for target integer range [a, b] */
/*____________________________________________________________________________*/
__inline__ static int get_rand_int(unsigned int min, unsigned int max)
{
	/* in [0, RAND_MAX] */
    int base_random = rand();
    if (base_random == RAND_MAX)
        return get_rand_int(min, max);

    /* now guaranteed to be in [0, RAND_MAX) */
    int range = (max + 1) - min;
    int remainder = RAND_MAX % range;
    int bin = RAND_MAX / range;

    if (base_random < (RAND_MAX - remainder)) {
        return (min + (base_random / bin));
    } else {
        return get_rand_int(min, max);
    }   
}

/*___________________________________________________________________________*/
/* compute factorial either using GSL or Stirling's */
/* log(n!) ≈ (n + 1/2)log(n) − n + 1/2log(2π) */
__inline__ static long double factorial(unsigned int n)
{
	long double nd = (long double)n;
	if (n <= GSL_SF_FACT_NMAX) {
		return (long double)gsl_sf_fact(n);
	} else {
		return expl((nd + 0.5) * logl(nd) - nd + 0.5 * logl(2.0 * M_PI));
	}
} 

/*____________________________________________________________________________*/
/* rMIxy */
/*____________________________________________________________________________*/
__inline__ static void armi(gsl_matrix *MI_mat_, unsigned int N,
						unsigned int L,	gsl_matrix *p_E_, unsigned int E, FILE *outFile_)
{
	unsigned int e, ex, ey, n, x, y, z, zz;

	/* formula terms */
	long double rMI = 0.; /* residual column mutual information */
	long double rMI_1 = 0.; /* 1. log(N) */
	long double rMI_2 = 0.; /* 2. p_x p_y */
	long double rMI_3 = 0.; /* 3. compute sum over N-1 */
	long double rMI_23 = 0.; /* rMI_2 * rMI_3 */
	long double rMI_3_1 = 0.; /* 3.1 binomial term */
	long double rMI_3_2 = 0.; /* 3.2 log(1 + n) */
	long double rMI_3_3 = 0.; /* 3.3 polynomial terms */
	long double rMI_3_3_1 = 0.; /* 3.3.1 (p_x p_y)^n */
	long double rMI_3_3_2 = 0.; /* 3.3.2 (1 - p_x p_y)^N-1-n */
	long double rMI_3_3_3 = 0.; /* 3.3.3 p_x^n (1 - p_x)^N-1-n */
	long double rMI_3_3_4 = 0.; /* 3.3.4 p_y^n (1 - p_y)^N-1-n */

	/* base set probabilities */
	long double p_x = 0.; /* probability of alphabet symbol x */
	long double p_y = 0.; /* probability of alphabet symbol y */
	long double p_z = 0.; /* to precompute probability */

	/* binomial coefficient */
	long double binomial_v[N - 1]; /* vector of binomial coefficients */

	/* pre-computable terms */
	long double *log1n = 0;
	long double ***rMI_3_3_34_mat = 0;

	/* loop completion */
	long double completion = 0.;
	unsigned int completion_i = 0;

	/*____________________________________________________________________________*/
	/* precompute some terms to speed up calculation */
	/* precompute binomial coefficients */
	for (n = (N - 1); n > 0; -- n) {
		binomial_v[n] = factorial(N - 1) / (factorial(n) * factorial(N - 1 - n)); 
		assert(binomial_v[n] <= LDBL_MAX && "value overflow: reduce sample size to max. 1755");
	}

	/* precompute log(1 + n) over n */
	log1n = safe_malloc(N * sizeof(long double));
	log1n[0] = 0.;
	for (n = 1; n <= (N - 1); ++ n) {
			log1n[n] = log(1 + n);
	}

	/* precompute exponentials of p(x) and p(y), expressed here as p(z) */
	rMI_3_3_34_mat = alloc_mat3D_longdouble(rMI_3_3_34_mat, (int)L, (int)E, (int)N);
	init_mat3D_longdouble(rMI_3_3_34_mat, (int)L, (int)E, (int)N, 0.);

	for (z = 0; z < L; ++ z) {
		for (e = 1; e < E; ++ e) {
			for (n = 1; n < N; ++ n) {
				p_z = (long double)gsl_matrix_get(p_E_, e, z);
				if (p_z > 0) {
					rMI_3_3_34_mat[z][e][n] = powl(p_z, n) * powl((1 - p_z), (N - 1 - n));
				}
			}
		}
	}

	/*____________________________________________________________________________*/
	/* 1. log(N) */
	rMI_1 = logl(N);

	/* from here onwards, column pairs determine p_x and p_y */
	/* for all first columns x */
	fprintf(stdout, "\tcompletion\tcol1\tcol2\n");
	for (x = 0, zz = 0; x < L - 1; ++ x) {
		/* ... and all second columns y */
		for (y = x + 1; y < L; ++ y) {
			/* print progress */
			++ zz;
			completion = (long double)zz / ((long double)(L*(L-1)) / 2) * 100;
			if ((int)completion > completion_i) {
				completion_i = (int)completion;
				fprintf(stdout, "\t%3d%%\t\t%d/%d\t%d/%d\r", completion_i, x, L-2, y, L-1);
				fflush(stdout);
			}

			/*____________________________________________________________________________*/
			/* 2. double sum over xeS yeS */
			/* for all elements x: ex */
			rMI_23 = 0.;
			for (ex = 1; ex < E; ++ ex) {
				/* for all elements y: ey */
				for (ey = 1; ey < E; ++ ey) {
					/* get probabilities p_x,p_y of all elements in columns x,y */
					p_x = (long double)gsl_matrix_get(p_E_, ex, x);
					p_y = (long double)gsl_matrix_get(p_E_, ey, y);

					if ((p_x > 0) && (p_y > 0)) {
						/* p_x * p_y */
						rMI_2 = p_x * p_y;

						/*____________________________________________________________________________*/
						/* 3. sum over N-1 */
						rMI_3 = 0.;
						for (n = (N - 1); n > 0; -- n) {
							/* 3.1 binomial term */
							rMI_3_1 = binomial_v[n];

							/* 3.2 log(1 + n) */
							rMI_3_2 = log1n[n];

							/* 3.3.1 (p_x p_y)^n */
							rMI_3_3_1 = powl((p_x * p_y), n);
							/* 3.3.2 (1 - p_x p_y)^N-1-n */ 
							rMI_3_3_2 = powl((1 - (p_x * p_y)), (N - 1 - n));
							/* 3.3.3 p_x^n (1 - p_x)^N-1-n */
							rMI_3_3_3 = rMI_3_3_34_mat[x][ex][n];
							/* 3.3.4 p_y^n (1 - p_y)^N-1-n */
							rMI_3_3_4 = rMI_3_3_34_mat[y][ey][n];
							/* 3.3 compute polynomial terms */
							rMI_3_3 = (rMI_3_3_1 * rMI_3_3_2) - rMI_3_3_3 - rMI_3_3_4;	

							/* 3. compute sum over N-1 */
							rMI_3 += (rMI_3_1 * rMI_3_2 * rMI_3_3);
						}
						/* 2*3. compute sum over all element pairs */
						rMI_23 += (rMI_2 * rMI_3);
					}
				}
			}

			rMI = rMI_1 + rMI_23;
			gsl_matrix_set(MI_mat_, x, y, (double)rMI);
			gsl_matrix_set(MI_mat_, y, x, (double)rMI);

			fprintf (outFile_, "%d\t%d\t%Lf\n", x+1, y+1, rMI);
		}
	}
	free(log1n);
	free_mat3D_longdouble(rMI_3_3_34_mat, L, E);
}

/*____________________________________________________________________________*/
/* MIxy */
/*____________________________________________________________________________*/
__inline__ static void columnpair_mutual_information(gsl_matrix *level_, unsigned int N,
						unsigned int L,	gsl_matrix *p_E_, unsigned int E,
						FILE *outFile_, gsl_vector *MIvec_)
{
	unsigned int e_x, e_y, n, x, y;
	double MIxy = 0.;
	unsigned int nMI = 0;
	double p_e_x, p_e_y, p_ee_xy;
	gsl_matrix *p_EE = gsl_matrix_calloc(E, E);
	double w_ee;
	double dN = (double)N;

	/* for all column pairs */
	/* for first columns x */
	for (x = 0; x < L - 1; ++ x) {
		/* for second columns y */
		for (y = x + 1; y < L; ++ y) {
			/* reset p_EE for each column pair */
			gsl_matrix_set_zero(p_EE);
			/* compute pair weight */
			w_ee = 1 / dN;
			/* compute p_EE for column pair x,y */
			/* for n sequences */
			for (n = 0; n < N; ++ n) {
				/* get probability of alignment elements (characters) at positions [n,x] and [n,y] */
				e_x = (int)gsl_matrix_get(level_, n, x);
				e_y = (int)gsl_matrix_get(level_, n, y);

				assert((e_x < E) && (e_y < E));

				/* compute element pair probabilities */
				gsl_matrix_set(p_EE, e_x, e_y, gsl_matrix_get(p_EE, e_x, e_y) + w_ee);
			}
			/* reset for each column pair */
			MIxy = 0.;
			/* compute MIxy for column pair x,y */
			for (e_x = 1; e_x < E; ++ e_x) {
				for (e_y = 1; e_y < E; ++ e_y) {
					p_ee_xy = gsl_matrix_get(p_EE, e_x, e_y);
					p_e_x = gsl_matrix_get(p_E_, e_x, x);
					p_e_y = gsl_matrix_get(p_E_, e_y, y);
					if (p_ee_xy > 0 && p_e_x > 0 && p_e_y > 0) {
						MIxy += p_ee_xy * log(p_ee_xy / (p_e_x * p_e_y));
					}
				}
			}
			fprintf (outFile_, "%d\t%d\t%lf\n", x+1, y+1, MIxy);
			gsl_vector_set(MIvec_, nMI, MIxy);
			++ nMI;
		}
	}
	gsl_matrix_free(p_EE);
}

/*____________________________________________________________________________*/
/* column shuffle by random permutation */
/*____________________________________________________________________________*/
__inline__ static double randomise_matrix(gsl_matrix *level_rand, gsl_matrix *level, unsigned int N, unsigned int L)
{
	unsigned int l, n;

	const gsl_rng_type *T = 0; /* random number generator type */
	gsl_rng *r = 0; /* random number generator */
	int c[N]; /* column array holding sequence letters */

	/* initialise random number generator */
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	/* shuffle column values */
	for (l = 0; l < L; ++ l) {
		for (n = 0; n < N; ++ n) {
			c[n] = (int)gsl_matrix_get(level, n, l);
		}

		gsl_ran_shuffle(r, c, N, sizeof(int));
		for (n = 0; n < N; ++ n) {
			gsl_matrix_set(level_rand, n, l, (double)c[n]);
		}
	}

	/* free memory */
	gsl_rng_free(r);

	return 0;
}

/*____________________________________________________________________________*/
/* arMI main routine */
/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	int e;
	unsigned int i, l, n, x, y, zz;
	char outFileName_arMI[128];
	char outFileName_MI[128];
	char outFileName_nrMI[128];
	char outFileName_rmsd[128];
	char outFileName_time[128];
	char outDirName[128];
	FILE *outFile_arMI = 0;
	FILE *outFile_MI = 0;
	FILE *outFile_nrMIit = 0;
	FILE *outFile_nrMI = 0;
	FILE *outFile_rmsd = 0;
	FILE *outFile_time = 0;

	Arg arg; /* command line arguments */
	Expr expr; /* expression values */

	unsigned int N = 0; /* number of genes (rows) */
	unsigned int L = 0; /* length of sample dimension (columns) */
	unsigned int E = 0; /* number of expression levels */

	char nrMIit[64];
	const unsigned int nIter = 100;
	unsigned int nMI = 0;
	double rmsd = 0.;

	clock_t begin, end;
	double time_spent_arMI;
	double time_spent_MI;
	double time_spent_nrMI;

	/* loop completion */
	long double completion = 0.;
	unsigned int completion_i = 0;

	/* scale factor to map expression values to 20 bins */
	float scafa = 0.;
	double l2min = 0.; /* log2 of readmin */

	/*____________________________________________________________________________*/
	/* random seed */
	srand(time(NULL));

    /*____________________________________________________________________________*/
    /* input: parse command line arguments */
	parse_args(argc, &(argv[0]), &arg);

	/*____________________________________________________________________________*/
	/* input row names */
	if (arg.rowfilename) {
		fprintf(stdout, "\nReading row names\n");

		arg.rowfile = safe_open(arg.rowfilename, "r");
		N = get_rownames(arg.rowfile, &expr);
		fclose(arg.rowfile);

		fprintf(stdout, "\tnumber of rows: %d\n", expr.nrow);
	}

	/*____________________________________________________________________________*/
	/* input column names */
	if (arg.colfilename) {
		fprintf(stdout, "\nReading column names\n");

		arg.colfile = safe_open(arg.colfilename, "r");
		L = get_colnames(arg.colfile, &expr);
		fclose(arg.colfile);

		fprintf(stdout, "\tnumber of columns: %d\n", expr.ncol);
	}

	/*____________________________________________________________________________*/
	/* input expression values */
	if (arg.exprfilename) {
		fprintf(stdout, "\nReading expression values\n");

		/* read specified expression values */
		arg.exprfile = safe_open(arg.exprfilename, "r");
		read_expression(arg.exprfile, &expr);
		fclose(arg.exprfile);

		fprintf(stdout, "\tnumber of expression values: %d\n", expr.ndat);
	}

    /*____________________________________________________________________________*/
	/* determine min/max expression values and number of expression levels */
	/* create expression level matrix */
	/* we shift the levels by +1 to keep the routines identical to those for
		sequences, where '0' indicates a gap (which is ignored) */
	gsl_matrix *level = gsl_matrix_calloc(N, L); /* expression level matrix */

	expr.maxLevel = 1;
	scafa = 19. / (log2(expr.readmax) - log2(expr.readmin));
	l2min = log2(expr.readmin);
#ifdef DEBUG
	fprintf(stderr, "readmin %f,  readmax %f, scafa %f\n",
			expr.readmin, expr.readmax, scafa);
#endif
	for (n = 0; n < expr.nrow; ++ n) {
		for (l = 0; l < expr.ncol; ++ l) {
			assert((expr.read[n][l] >= 0.) && "expression values should be positive");
			if (expr.read[n][l] > 0.) {
				gsl_matrix_set(level, n, l, (double)(roundf((log2(expr.read[n][l]) - l2min) * scafa) + 1.));
				expr.maxLevel = max((int)gsl_matrix_get(level, n, l), expr.maxLevel);
			} else {
				gsl_matrix_set(level, n, l, (double)1.);
			}

		}
	}

	E = expr.maxLevel + 1; /* number of expression levels (in bits) */
	fprintf(stdout, "\tmax: %d\n\texpression levels: %d bit\n",
			expr.maxLevel, E);

	/*____________________________________________________________________________*/
	/* compute column-wise element probability matrix p_E */
	gsl_matrix *p_E = gsl_matrix_calloc(E + 1, L); /* probabilities of elements */

	/* for all columns (samples) */
	for (l = 0; l < L; ++ l) {
		/* for all in column 'l':
			probablities of expression levels in column 'l' */
		for (n = 0; n < N; ++ n) {
			/* set current element (0 - E) */
			e = (int)gsl_matrix_get(level, n, l);
#ifdef DEBUG
			fprintf(stderr, "row %d, column %d, e %d\n", n, l, e);
#endif
			assert(e <= E);
			if ((e < 0)  || (e > E)) {
				printf("Symbol %c (ASCII %d) out of range [0-E]\n", e, e);
				printf("expression level: row %d, column %d\n", n, l);
				exit(1);
			}
			if (e > 0) {
				gsl_matrix_set(p_E, e, l, (gsl_matrix_get(p_E, e, l) + (1 / (double)N)));
			}
		}
	}

	/*____________________________________________________________________________*/
	/* analytical residual MI > 'arMI.dat' */
	fprintf(stdout, "\nComputing arMI\n");

	gsl_matrix *MI_mat = gsl_matrix_calloc(L, L); /* MI matrix */
	const unsigned int nPairs = L * (L - 1) / 2;
	gsl_vector *arMIvec = gsl_vector_calloc(nPairs); /* to compute rmsd below */

	sprintf(outFileName_arMI, "%s%s", arg.prefix, "arMI.dat");
	outFile_arMI = fopen(outFileName_arMI, "w");
	fprintf(outFile_arMI, "col1\tcol2\tarMI\n");
	begin = clock();
	armi(MI_mat, N, L, p_E, E, outFile_arMI);
	end = clock();
	time_spent_arMI = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\n\tt(arMI): %lf s\n", time_spent_arMI);
	fclose(outFile_arMI);

	/*____________________________________________________________________________*/
	/* MI > 'MI.dat': */
	fprintf(stdout, "\nComputing MI\n");

	gsl_vector *cMIvec = gsl_vector_calloc(nPairs);

	sprintf(outFileName_MI, "%s%s", arg.prefix, "MI.dat");
	outFile_MI = fopen(outFileName_MI, "w");
	fprintf(outFile_MI, "col1\tcol2\tMI\n");
	begin = clock();
	columnpair_mutual_information(level, N, L, p_E, E, outFile_MI, cMIvec);
	end = clock();
	time_spent_MI = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\tt(MI): %lf s\n", time_spent_MI);
	fclose(outFile_MI);

	/*____________________________________________________________________________*/
	/* numerical residual MI > 'nrMI.dat': */
	/* create shuffled 100-fold randomised alignment matrix */
	/*   and compute randomised column-wise element probability matrix p_E_rand */
	fprintf(stdout, "\nComputing nrMI\n");
	fprintf(stdout, "\trandomising input alignment %d times\n", nIter);

	/* randomised expression level matrix */
	gsl_matrix *level_rand = gsl_matrix_calloc(N, L);
	/* probabilities of elements */
	gsl_matrix *p_E_rand = gsl_matrix_calloc(E, L);
	/* nrMI values of repeated randomisation */
	gsl_vector *nrMIvec = gsl_vector_calloc(nPairs);
	/* nrMI mean values of repeated randomisation */
	gsl_vector *nrMImeanVec = gsl_vector_calloc(nPairs);
	/* copy of the previous */
	gsl_vector *nrMImeanVecCp = gsl_vector_calloc(nPairs);
	/* nrMI var values of repeated randomisation */
	gsl_vector *nrMIvarVec = gsl_vector_calloc(nPairs);

	sprintf(outDirName, "%s%s", arg.prefix, "nrMI_data");
	mkdir(outDirName, 0777);

	begin = clock();

	fprintf(stdout, "\tcompletion\titer\n");
	for (i = 0, zz = 0; i < nIter; ++ i) {
		/* print progress */
		++ zz;
		completion = (long double)zz / ((long double)nIter) * 100;
		if ((int)completion > completion_i) {
			completion_i = (int)completion;
			fprintf(stdout, "\t%3d%%\t\t%d/%d\r", completion_i, zz, nIter);
			fflush(stdout);
		}

		sprintf(nrMIit, "%s%s%d%s", arg.prefix, "nrMI_data/nrMI.", i, ".dat");
		outFile_nrMIit = fopen(nrMIit, "w");
		fprintf(outFile_nrMIit, "col1\tcol2\tnrMI\n");

		if (i == 0) {
			randomise_matrix(level_rand, level, N, L);
		} else {
			randomise_matrix(level_rand, level_rand, N, L);
		}
		gsl_matrix_set_zero(p_E_rand);
		gsl_vector_set_zero(nrMIvec);

		/* for the alignment length */
		for (l = 0; l < L; ++ l) {
			/* for all sequences in column 'l':
				probablities of alphabet elements in column 'l' */
			for (n = 0; n < N; ++ n) {
				/* set current element (0 - 25) */
				e = (int)gsl_matrix_get(level_rand, n, l);
				if (e > 0) {
					gsl_matrix_set(p_E_rand, e, l, (gsl_matrix_get(p_E_rand, e, l)
												+ (1 / (double)(N))));
				}
			}
		}

		columnpair_mutual_information(level_rand, N, L, p_E_rand, E, outFile_nrMIit, nrMIvec);
		/* sum(MI) of each column pair over iterations */
		gsl_vector_add(nrMImeanVec, nrMIvec);
		/* sum(MI^2) of each column pair over iterations */
		gsl_vector_mul(nrMIvec, nrMIvec);
		gsl_vector_add(nrMIvarVec, nrMIvec);

		fclose(outFile_nrMIit);
	}

	/* normalise to mean values */
	gsl_blas_dscal((1 / (double)nIter), nrMImeanVec); /* E/N */
	gsl_blas_dscal((1 / (double)nIter), nrMIvarVec); /* E^2/N */

	/* compute E^2/N - (E/N)^2 */
	gsl_vector_memcpy(nrMImeanVecCp, nrMImeanVec);
	gsl_vector_mul(nrMImeanVecCp, nrMImeanVecCp);
	gsl_vector_sub(nrMIvarVec, nrMImeanVecCp);

	/* write vectors */
	sprintf(outFileName_nrMI, "%s%s", arg.prefix, "nrMI.dat");
	outFile_nrMI = fopen(outFileName_nrMI, "w");
	fprintf(outFile_nrMI, "col1\tcol2\tmean(nrMI)\tsd(nrMI)\n");
	for (x = 0, nMI = 0; x < L - 1; ++ x) {
		for (y = x + 1; y < L; ++ y) {
			fprintf(outFile_nrMI, "%d\t%d\t%lf\t%lg\n", x+1, y+1,
						gsl_vector_get(nrMImeanVec, nMI),
						sqrt(gsl_vector_get(nrMIvarVec, nMI)));	
			++ nMI;
		}
	}
	end = clock();
	time_spent_nrMI = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\n\tt(nrMI): %lf s\n", time_spent_nrMI);
	fclose(outFile_nrMI);

	/*____________________________________________________________________________*/
	/* write times spent */
	sprintf(outFileName_time , "%s%s", arg.prefix, "time.dat");
	outFile_time = fopen(outFileName_time, "w");
	fprintf(outFile_time, "%d\t%lf\t%lf\n", arg.nsubset, time_spent_arMI, time_spent_nrMI);
	fclose(outFile_time);

	/*____________________________________________________________________________*/
	/* compute rmsd between all pairs of arMI and nrMI */
	gsl_vector *rmsdvec = gsl_vector_calloc(nPairs);
	gsl_vector_memcpy(rmsdvec, arMIvec);
	/* difference */
	gsl_vector_sub(rmsdvec,nrMIvec);
	/* Euclidean norm */
	rmsd = gsl_blas_dnrm2(rmsdvec) / nPairs;
	
	sprintf(outFileName_rmsd, "%s%s", arg.prefix, "rmsd.dat");
	outFile_rmsd = fopen(outFileName_rmsd, "w");
	fprintf(outFile_rmsd, "%d\t%lf\n", arg.nsubset, rmsd);
	fclose(outFile_rmsd);

	/*____________________________________________________________________________*/
	/* free memory */
	gsl_matrix_free(MI_mat);
	gsl_matrix_free(p_E);
	gsl_matrix_free(p_E_rand);
	gsl_matrix_free(level);
	gsl_matrix_free(level_rand);
	gsl_vector_free(arMIvec);
	gsl_vector_free(cMIvec);
	gsl_vector_free(nrMIvec);
	gsl_vector_free(nrMImeanVec);
	gsl_vector_free(nrMImeanVecCp);
	gsl_vector_free(nrMIvarVec);
	gsl_vector_free(rmsdvec);

	free_mat2D_float(expr.read, expr.nrow);
	free(expr.colname);
	free(expr.rowname);

	/*____________________________________________________________________________*/
	/* termination*/
	fprintf(stdout, "\nSuccessful program termination\n\n");
	return 0;
}

