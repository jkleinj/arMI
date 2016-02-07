/*==============================================================================
armi_mpi.c : analytical residual Mutual Information
(C) 2014-2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "armi.h"

/*____________________________________________________________________________*/
/* global variables */
#ifdef MPI
#include <mpi.h>
    int my_rank; /* rank of 'this' node */
    int nodes; /* number of nodes */
#else
    int my_rank = 0;  
    int nodes = 1;  
#endif

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

/*____________________________________________________________________________*/
/* arMIxy */
/*____________________________________________________________________________*/
__inline__ static void residual_mutual_information(gsl_matrix *MI_mat_, unsigned int N,
						unsigned int L,	gsl_matrix *p_E_, unsigned int E, FILE *outFile_)
{
	unsigned int e, ex, ey, n, x, y, z;

	/* formulat terms */
	long double arMI = 0.; /* residual column mutual information */
	long double arMI_1 = 0.; /* 1. log(N) */
	long double arMI_2 = 0.; /* 2. p_x p_y */
	long double arMI_3 = 0.; /* 3. compute sum over N-1 */
	long double arMI_23 = 0.; /* arMI_2 * arMI_3 */
	long double arMI_3_1 = 0.; /* 3.1 binomial term */
	long double arMI_3_2 = 0.; /* 3.2 log(1 + n) */
	long double arMI_3_3 = 0.; /* 3.3 polynomial terms */
	long double arMI_3_3_1 = 0.; /* 3.3.1 (p_x p_y)^n */
	long double arMI_3_3_2 = 0.; /* 3.3.2 (1 - p_x p_y)^N-1-n */
	long double arMI_3_3_3 = 0.; /* 3.3.3 p_x^n (1 - p_x)^N-1-n */
	long double arMI_3_3_4 = 0.; /* 3.3.4 p_y^n (1 - p_y)^N-1-n */

	/* amino acid (code element) probabilities */
	long double p_x = 0.; /* probability of code element x */
	long double p_y = 0.; /* probability of altcode element y */
	long double p_z = 0.; /* to precompute probability */

	/* binomial coefficient */
	unsigned int binomial_gsl; /* choice of binomial coefficient function */
	long double binomial_v[N - 1]; /* vector of binomial coefficients */

	/* progress monitor */
	long double completion = 0.;
	unsigned int completion_i = 0;

	/* pre-computable terms */
	long double *log1n = 0;
	long double ***arMI_3_3_34_mat = 0;

#ifdef MPI
	unsigned int i, xmpi, ympi;
	int dimz = (int)ceil((float)L / (float)nodes); /* MPI dimension */
    long double arMI_local[dimz]; /* results of this node */
	int x_local[dimz]; /* index array for mapping */
	int y_local[dimz]; /* index array for mapping */
    long double arMI_global[dimz * nodes]; /* results of all nodes */
	int nz = 0; /* index of dimz jobs */

	for (i = 0; i < dimz; ++ i) {
		arMI_local[i] = 0.;
	}
#endif

	/*____________________________________________________________________________*/
	/* overflow avoidance */
	if ((N > 0) && (N <= (int)GSL_SF_DOUBLEFACT_NMAX)) {
		if (my_rank == 0) fprintf(stdout, "\tbinomial coefficient from GSL function 'gsl_sf_choose'\n\n");
		binomial_gsl = 1;
	} else if ((N > (int)GSL_SF_DOUBLEFACT_NMAX) && (N <= 10000)) {
		if (my_rank == 0) fprintf(stdout, "\tbinomial coefficient from Stirling approximation\n\n");
		binomial_gsl = 0;
	} else {
		if (my_rank == 0) fprintf(stderr, "\tExiting: number of sequences %d out of valid range [1,10000]\n\n", N);
		exit(1);
	}

	/*____________________________________________________________________________*/
	/* precompute some terms to speed up calculation */
	/* precompute binomial coefficients */
	for (n = (N - 1); n > 0; -- n) {
		/* GSL */
		if (binomial_gsl) {
			binomial_v[n] = (long double)gsl_sf_choose((N - 1), n);
		/* Stirling */
		} else {
			if (n == (N - 1)) {
				binomial_v[n] = 0;
			} else {
				binomial_v[n] = (long double)((N - 1 + .5) * logl(N - 1)
						- (n + .5) * logl(n)
						- (N - 1 - n + .5) * logl(N - 1 - n)
						- .5 * logl(2 * M_PI));
			}
		}
	}

	/* precompute log(n - 1) */
	log1n = safe_malloc(N * sizeof(long double));
	log1n[0] = 0.;
	for (n = 1; n <= (N - 1); ++ n) {
			log1n[n] = log(1 + n);
	}

	/* precompute exponentials of p(x) and p(y) */
	arMI_3_3_34_mat = alloc_mat3D_longdouble(arMI_3_3_34_mat, (int)L, (int)E, (int)N);
	init_mat3D_longdouble(arMI_3_3_34_mat, (int)L, (int)E, (int)N, 0.);

	for (z = 0; z < L; ++ z) {
		for (e = 1; e < E; ++ e) {
			for (n = 1; n < N; ++ n) {
				p_z = (long double)gsl_matrix_get(p_E_, e, z);
				if (p_z > 0) {
					arMI_3_3_34_mat[z][e][n] = powl(p_z, n) * powl((1 - p_z), (N - 1 - n));
				}
			}
		}
	}

	/*____________________________________________________________________________*/
	/* 1. log(N) */
	arMI_1 = logl(N);

	/* from here on, column pairs determine p_x and p_y */
	if (my_rank == 0) fprintf(stdout, "\tcompletion\tcol1\tcol2\n");
	unsigned int zz = 0;

	/* for all first columns x */
	for (x = 0, zz = 0; x < L - 1; ++ x) {
		/* ... and all second columns y */
		for (y = x + 1; y < L; ++ y) {
			/* print progress */
			++ zz;
			completion = (long double)zz / ((long double)(L * (L-1)) / 2) * 100;
			if ((int)completion  > completion_i) {
				completion_i = (int)completion;
				if (my_rank == 0) fprintf(stdout, "\t%3d%%\t\t%d/%d\t%d/%d\n", completion_i, x, L-2, y, L-1);
			}

#ifdef MPI
			if (zz % nodes == my_rank) { 
#endif
			/*____________________________________________________________________________*/
			/* 2. double sum over xeS yeS */
			/* for all elements x: ex */
			arMI_23 = 0.;
			for (ex = 1; ex < E; ++ ex) {
				/* for all elements y: ey */
				for (ey = 1; ey < E; ++ ey) {
					/* get probabilities p_x,p_y of all elements in columns x,y */
					p_x = (long double)gsl_matrix_get(p_E_, ex, x);
					p_y = (long double)gsl_matrix_get(p_E_, ey, y);

					if ((p_x > 0) && (p_y > 0)) {
						/* p_x * p_y */
						arMI_2 = p_x * p_y;

						/*____________________________________________________________________________*/
						/* 3. sum over N-1 */
						arMI_3 = 0.;
						for (n = (N - 1); n > 0; -- n) {
							/* 3.1 binomial term */
							arMI_3_1 = binomial_v[n];

							/* 3.2 log(1 + n) */
							arMI_3_2 = log1n[n];

							/* 3.3.1 (p_x p_y)^n */
							arMI_3_3_1 = powl((p_x * p_y), n);
							/* 3.3.2 (1 - p_x p_y)^N-1-n */ 
							arMI_3_3_2 = powl((1 - (p_x * p_y)), (N - 1 - n));
							/* 3.3.3 p_x^n (1 - p_x)^N-1-n */
							arMI_3_3_3 = arMI_3_3_34_mat[x][ex][n];
							/* 3.3.4 p_y^n (1 - p_y)^N-1-n */
							arMI_3_3_4 = arMI_3_3_34_mat[y][ey][n];
							/* 3.3 compute polynomial terms */
							arMI_3_3 = (arMI_3_3_1 * arMI_3_3_2) - arMI_3_3_3 - arMI_3_3_4;	

							/* 3. compute sum over N-1 */
							if (binomial_gsl) {
								arMI_3 += arMI_3_1 * arMI_3_2 * arMI_3_3;
							} else {
								arMI_3 += expl(arMI_3_1) * arMI_3_2 * arMI_3_3;
							}
						}
						/* 2*3. compute sum over all element pairs */
						arMI_23 += arMI_2 * arMI_3;
					}
				}
			}

#ifdef MPI
			x_local[nz] = x;
			y_local[nz] = y;
			arMI_local[nz] = arMI_1 + arMI_23;
			++ nz;
			}
#else
			arMI = arMI_1 + arMI_23;
			gsl_matrix_set(MI_mat_, x, y, (double)arMI);
			gsl_matrix_set(MI_mat_, y, x, (double)arMI);
#endif
		}
	}

#ifdef MPI
    /* communicate data between all nodes */
    MPI_Allgather(arMI_local, dimx, MPI_LONG_DOUBLE, arMI_global, dimx, MPI_LONG_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(x_local, dimz, MPI_INT, x_global, dimz, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(y_local, dimz, MPI_INT, y_global, dimz, MPI_INT, MPI_COMM_WORLD);

	/* copy global data to target matrix */
	for(xmpi = 0; xmpi < L; ++ xmpi) {
		for (ympi = xmpi + 1; ympi < L; ++ ympi) {
			gsl_matrix_set(MI_mat_, xmpi, ympi, (double)arMI_global[xmpi * L + ympi]);
			gsl_matrix_set(MI_mat_, ympi, xmpi, (double)arMI_global[xmpi * L + ympi]);
		}
	}
#endif

	for (x = 0; x < L - 1; ++ x) {
		for (y = x + 1; y < L; ++ y) {
			if (my_rank == 0) fprintf(outFile_, "%d\t%d\t%Lf\n", x+1, y+1, arMI);
		}
	}

	free(log1n);
	free_mat3D_longdouble(arMI_3_3_34_mat, L, E);
}

/*____________________________________________________________________________*/
/* cMIxy */
/*____________________________________________________________________________*/
__inline__ static void columnpair_mutual_information(gsl_matrix *mali_, unsigned int N,
						unsigned int L,	gsl_matrix *p_E_, unsigned int E, unsigned int *Ngap, FILE *outFile_)
{
	/*unsigned int e;
	double p_sum;*/
	unsigned int e_x, e_y, n, x, y;
	double MIxy = 0.;
	double p_e_x, p_e_y, p_ee_xy;
	gsl_matrix *p_EE = gsl_matrix_calloc(E, E);
	double w_ee;
	double dN = (double)N;
	double dNgapx, dNgapy;

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
			/* for n sequences */
			for (n = 0; n < N; ++ n) {
				/* get probability of alignment elements (characters) at positions [n,x] and [n,y] */
				e_x = (int)gsl_matrix_get(mali_, n, x);
				e_y = (int)gsl_matrix_get(mali_, n, y);
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
			if (my_rank == 0) fprintf(outFile_, "%d\t%d\t%lf\n", x+1, y+1, MIxy);
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
	unsigned int l, n;
	FILE *outFile = 0;
	FILE *seqOutFile = 0;

	Arg arg; /* command line arguments */
	Seq *sequence = 0; /* sequence alignment */

	unsigned int N = 0; /* number of sequences */
	unsigned int L = 0; /* length of sequences */
	unsigned int E = 0; /* number of alphabet elements */

	unsigned int allocated = 64;

    /*________________________________________________________________________*/
    /* MPI */
#ifdef MPI
    /* initialize MPI routines */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

	/*____________________________________________________________________________*/
	/* random seed */
	srand(time(NULL));

    /*____________________________________________________________________________*/
    /* input: parse command line arguments */
	parse_args(argc, &(argv[0]), &arg);

    /*____________________________________________________________________________*/
	E = arg.nEle; /* number of alphabet elements */

	/*____________________________________________________________________________*/
	/* input alignment: random or from file */
	/* if no input alignment, create random alignment */
	if (arg.random) {
		if (my_rank == 0) fprintf(stdout, "Creating random sequence alignment\n\n");

		N = arg.nSeq; /* number of sequences */
		L = arg.lSeq; /* length of sequences */
		/* create random alignment*/
		sequence = (Seq *)safe_malloc(N * sizeof(Seq));
		seqOutFile = safe_open("randseq.fasta", "w");
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

			if (my_rank == 0) fprintf(stdout, "\t%s\n\t%s\n", sequence[n].name, sequence[n].residue);
			if (my_rank == 0) fprintf(seqOutFile, "%s\n%s\n", sequence[n].name, sequence[n].residue);
		}
		fclose(seqOutFile);
	/* ... or read alignment from file */
	} else if (arg.inFileName) {
		if (my_rank == 0) fprintf(stdout, "Reading alignment\n");

		/* read specified alignment*/
		sequence = (Seq *)safe_malloc(allocated * sizeof(Seq));
		arg.inFile = safe_open(arg.inFileName, "r");
		n = 0;
		while (read_sequence(arg.inFile, &(sequence[n]))) {
			++ n;
			if (n == allocated) {
				allocated += 64;
				sequence = safe_realloc(sequence, allocated * sizeof(Seq));
			}
		}
		fclose(arg.inFile);
		N = n;
		L = strlen(sequence[0].residue);
	}

	if (my_rank == 0) fprintf(stdout, "\n\tnumber of sequences: %d\n"
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
				if (my_rank == 0) printf("Symbol %c (ASCII %d) out of range [@-Z]\n", e, e);
				if (my_rank == 0) printf("ASCII %d '%c'\n", sequence[n].residue[l], sequence[n].residue[l]);
				exit(1);
			}
			if (e > 0) {
				gsl_matrix_set(p_E, e, l, (gsl_matrix_get(p_E, e, l) + (1 / (double)(N - Ngap[l]))));
			}
		}
	}

	/*____________________________________________________________________________*/
	/* create randomised alignment matrix */
	if (my_rank == 0) fprintf(stdout, "\nRandomising input alignment\n");

	gsl_matrix *mali_rand = gsl_matrix_calloc(N, L); /* randomised mali matrix */
	randomise_matrix(mali_rand, mali, N, L);

	/*____________________________________________________________________________*/
	/* compute randomised column-wise element probability matrix p_E_rand */
	gsl_matrix *p_E_rand = gsl_matrix_calloc(E, L); /* probabilities of elements */

	/* for the alignment length */
	for (l = 0; l < L; ++ l) {
		/* for all sequences in column 'l':
			probablities of alphabet elements in column 'l' */
		for (n = 0; n < N; ++ n) {
			/* set current element (0 - 25) */
			e = (int)gsl_matrix_get(mali_rand, n, l);
			if (e > 0) {
				gsl_matrix_set(p_E_rand, e, l, (gsl_matrix_get(p_E_rand, e, l) + (1 / (double)(N - Ngap[l]))));
			}
		}
	}

	/*____________________________________________________________________________*/
	/* analytical formula for residual MI, results go to 'arMI.dat': */
	/*   compute the analytical residual column MI (arMI) */
	/*   ordinal numbers in comments refer to terms of equation (11) 
         in accompanying theory */
	if (my_rank == 0) fprintf(stdout, "\nComputing arMI\n");

	gsl_matrix *MI_mat = gsl_matrix_calloc(L, L); /* MI matrix */

	outFile = fopen("arMI.dat", "w");
	residual_mutual_information(MI_mat, N, L, p_E, E, outFile);
	fclose(outFile);

	/*____________________________________________________________________________*/
	/* classical MI formula, results go to 'cMI.dat': */
	/*   compute column MI and randomised residual column MI */

	/* column-pair MI */
	if (my_rank == 0) fprintf(stdout, "\nComputing cMI\n");

	outFile = fopen("cMI.dat", "w");
	columnpair_mutual_information(mali, N, L, p_E, E, Ngap, outFile);
	fclose(outFile);

	/*____________________________________________________________________________*/
	/* randomised residual column-pair MI, results go to 'srMI.dat': */
	if (my_rank == 0) fprintf(stdout, "\nComputing srMI\n");

	outFile = fopen("srMI.dat", "w");
	columnpair_mutual_information(mali_rand, N, L, p_E_rand, E, Ngap, outFile);
	fclose(outFile);

	/*____________________________________________________________________________*/
	gsl_matrix_free(MI_mat);
	gsl_matrix_free(p_E);
	gsl_matrix_free(p_E_rand);
	gsl_matrix_free(mali);
	gsl_matrix_free(mali_rand);

	for (n = 0; n < N; ++ n) {
		free(sequence[n].name);
		free(sequence[n].residue);
	}

	free(sequence);

    /*________________________________________________________________________*/
    /* MPI */
#ifdef MPI
    /* stop MPI processes */
    MPI_Finalize();
#endif

	/*____________________________________________________________________________*/
	/* termination*/
	if (my_rank == 0) fprintf(stdout, "\nSuccessful program termination\n\n");
	return 0;
}
