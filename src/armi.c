/*==============================================================================
armi.c : analytical residual Mutual Information
(C) 2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <gsl/gsl_errno.h>

#include "armi.h"

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
		dNgapx = (double)Ngap[x];
		/* for second columns y */
		for (y = x + 1; y < L; ++ y) {
			dNgapy = (double)Ngap[y];
			/* reset p_EE for each column pair */
			gsl_matrix_set_zero(p_EE);
			/* compute pair weight */
			w_ee = 1 / (((dN - dNgapx) / dN) * (dN - dNgapy) / dN * dN);
			/* compute p_EE for column pair x,y */
			/* for n rows */
			for (n = 0; n < N; ++ n) {
				/* get probability of alignment elements (characters) at positions [n,x] and [n,y] */
				e_x = (int)gsl_matrix_get(level_, n, x);
				e_y = (int)gsl_matrix_get(level_, n, y);
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
__inline__ static double randomise_matrix(gsl_matrix *mali_rand, gsl_matrix *mali, unsigned int N, unsigned int L)
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
			c[n] = (int)gsl_matrix_get(mali, n, l);
		}

		gsl_ran_shuffle(r, c, N, sizeof(int));
		for (n = 0; n < N; ++ n) {
			gsl_matrix_set(mali_rand, n, l, (double)c[n]);
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
	char outFileName_seq[128];
	char outDirName[128];
	FILE *outFile_arMI = 0;
	FILE *outFile_MI = 0;
	FILE *outFile_nrMIit = 0;
	FILE *outFile_nrMI = 0;
	FILE *outFile_rmsd = 0;
	FILE *outFile_time = 0;
	FILE *outFile_seq = 0;

	Arg arg; /* command line arguments */
	Seq *sequence = 0; /* sequence alignment */

	unsigned int N = 0; /* number of sequences */
	unsigned int L = 0; /* length of sequences */
	unsigned int E = 0; /* number of alphabet elements */

	unsigned int allocated = 64;

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

	/*____________________________________________________________________________*/
	/* random seed */
	srand(time(NULL));

    /*____________________________________________________________________________*/
    /* input: parse command line arguments */
	parse_args(argc, &(argv[0]), &arg);

    /*____________________________________________________________________________*/
	E = arg.nele; /* number of alphabet elements */

	/*____________________________________________________________________________*/
	/* input alignment: random or from file */
	/* if no input alignment, create random alignment */
	if (arg.random) {
		fprintf(stdout, "\nCreating random sequence alignment\n");

		N = arg.nseq; /* number of sequences */
		L = arg.lseq; /* length of sequences */
		/* create random alignment*/
		sequence = (Seq *)safe_malloc(N * sizeof(Seq));
		sprintf(outFileName_seq, "in.fasta");
		outFile_seq = safe_open(outFileName_seq, "w");
		for (n = 0; n < N; ++ n) {
			/* sequence name */
			sequence[n].name = (char *)safe_malloc(8 * sizeof(char));
			sprintf(sequence[n].name, ">%d", n);
			/* sequence residues */
			sequence[n].residue = (char *)safe_malloc((L + 1) * sizeof(char));
			for (l = 0; l < L; ++ l) {
				sequence[n].residue[l] = get_rand_int(65, 90);
			}
			sequence[n].residue[l] = '\0';

			/* print sequences */
			fprintf(outFile_seq, "%s\n%s\n", sequence[n].name, sequence[n].residue);
		}
		fclose(outFile_seq);
	/* ... or read alignment from file */
	} else if (arg.infilename) {
		fprintf(stdout, "\nReading alignment\n");

		/* read specified alignment*/
		sequence = (Seq *)safe_malloc(allocated * sizeof(Seq));
		arg.infile = safe_open(arg.infilename, "r");
		n = 0;
		while (read_sequence(arg.infile, &(sequence[n]))) {
			++ n;
			if (n == arg.nsubset) {
				break;
			}
			if (n == allocated) {
				allocated += 64;
				sequence = safe_realloc(sequence, allocated * sizeof(Seq));
			}
		}
		fclose(arg.infile);
		N = n;
		L = strlen(sequence[0].residue);
	}

	fprintf(stdout, "\tnumber of sequences: %d\n"
					"\tsequence length: %d\n", N, L);

	/*____________________________________________________________________________*/
	/* create alignment matrix */
	/*   and count gap frequencies */
	unsigned int Ngap[L];
	for (l = 0; l < L; ++ l) {
		Ngap[l] = 0;
	}
	
	gsl_matrix *mali = gsl_matrix_calloc(N, L); /* mali matrix */

	for (n = 0; n < N; ++ n) {
		for (l = 0; l < L; ++ l) {
			gsl_matrix_set(mali, n, l, (double)(sequence[n].residue[l] - 64));
			if ((int)gsl_matrix_get(mali, n, l) == 0) {
				++ Ngap[l];
			}
		}
	}

	/*____________________________________________________________________________*/
	/* compute column-wise element probability matrix p_E */
	gsl_matrix *p_E = gsl_matrix_calloc(E, L); /* probabilities of elements */

	/* for the alignment length */
	for (l = 0; l < L; ++ l) {
		/* for all sequences in column 'l':
			probablities of alphabet elements in column 'l' */
		for (n = 0; n < N; ++ n) {
			/* set current element (0 - 26) */
			e = (int)gsl_matrix_get(mali, n, l);
			//assert(e <= E);
			if ((e < 0)  || (e > E)) {
				printf("Symbol %c (ASCII %d) out of range [@-Z]\n", e, e);
				printf("ASCII %d '%c'\n", sequence[n].residue[l], sequence[n].residue[l]);
				exit(1);
			}
			if (e > 0) {
				gsl_matrix_set(p_E, e, l, (gsl_matrix_get(p_E, e, l) + (1 / (double)(N - Ngap[l]))));
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
	columnpair_mutual_information(mali, N, L, p_E, E, Ngap, outFile_MI, cMIvec);
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

	/* randomised mali matrix */
	gsl_matrix *mali_rand = gsl_matrix_calloc(N, L);
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
			randomise_matrix(mali_rand, mali, N, L);
		} else {
			randomise_matrix(mali_rand, mali_rand, N, L);
		}
		gsl_matrix_set_zero(p_E_rand);
		gsl_vector_set_zero(nrMIvec);

		/* for the alignment length */
		for (l = 0; l < L; ++ l) {
			/* for all sequences in column 'l':
				probablities of alphabet elements in column 'l' */
			for (n = 0; n < N; ++ n) {
				/* set current element (0 - 25) */
				e = (int)gsl_matrix_get(mali_rand, n, l);
				if (e > 0) {
					gsl_matrix_set(p_E_rand, e, l, (gsl_matrix_get(p_E_rand, e, l)
												+ (1 / (double)(N - Ngap[l]))));
				}
			}
		}

		columnpair_mutual_information(mali_rand, N, L, p_E_rand, E, Ngap, outFile_nrMIit, nrMIvec);
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
	gsl_matrix_free(MI_mat);
	gsl_matrix_free(p_E);
	gsl_matrix_free(p_E_rand);
	gsl_matrix_free(mali);
	gsl_matrix_free(mali_rand);
	gsl_vector_free(arMIvec);
	gsl_vector_free(cMIvec);
	gsl_vector_free(nrMIvec);
	gsl_vector_free(nrMImeanVec);
	gsl_vector_free(nrMImeanVecCp);
	gsl_vector_free(nrMIvarVec);
	gsl_vector_free(rmsdvec);

	for (n = 0; n < N; ++ n) {
		free(sequence[n].name);
		free(sequence[n].residue);
	}

	free(sequence);

	/*____________________________________________________________________________*/
	/* termination*/
	fprintf(stdout, "\nSuccessful program termination\n\n");
	return 0;
}

