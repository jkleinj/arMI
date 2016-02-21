/*==============================================================================
arg.c : parse command line arguments
Copyright (C) 2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "arg.h"

/*____________________________________________________________________________*/
/** print version */
static void print_version(Arg *arg)
{
	if (! arg->silent) fprintf(stdout, "\nversion %s\n", VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header(Arg *arg)
{
    if (! arg->silent) fprintf(stdout, "\narMI: analytical residual MI between pairs of alignment columns\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license(Arg *arg)
{
    if (! arg->silent) fprintf(stdout, "\nCopyright (C) 2015 Jens Kleinjung and Ton Coolen\n"
			"This program is free software and comes with ABSOLUTELY NO WARRANTY.\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the COPYING file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation(Arg *arg)
{
	if (! arg->silent) fprintf(stdout, "\nunpublished\n\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg)
{
    arg->infilename = "";
	arg->nseq = 8; /* number of sequences */
	arg->lseq = 16; /* sequence length */
	arg->nele = 27; /* number of alphabet elements */
	arg->random = 0; /* create random alignment */
	arg->nsubset = 0; /* number of sequences in subset */
	arg->prefix = ""; /* output prefix */
	arg->silent = 0; /* stdout */
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg)
{
	assert(arg->nseq > 1 && "Need at least 2 sequences");
	assert(arg->lseq > 0 && "Sequence length must be >0");
	assert(arg->nele > 0 && "Number of alphabet symbols must be >0");
	assert(((strlen(arg->infilename) > 0) || (arg->random == 1)) && "Need an alignment");
	assert(arg->nsubset >= 0 && "Number of sequences in subset must be >0");
}

/*____________________________________________________________________________*/
/** print command line arguments */
static void print_args(Arg *arg)
{
    time_t now;
    time(&now);

    if (! arg->silent) fprintf(stdout, "\ndate: %s", ctime(&now));
    if (! arg->silent) fprintf(stdout, \
                    "input file: %s\n"
                    "random: %d\n"
					"number of sequences in subset: %d\n"
                    "prefix: %s\n"
                    "number of alphabet elements: %d\n"
                    "GSL factorial maximum: %d\n"
                    "long double maximum: %Le\n",
        arg->infilename, arg->random,
		arg->nsubset, arg->prefix, arg->nele,
		GSL_SF_FACT_NMAX, LDBL_MAX);
	if (! arg->silent) if (arg->random) {
		fprintf(stdout, "number of random sequences: %d\n"
						"length of random sequences: %d\n"
						"input alignment stored in: %sin.fasta\n",
		arg->nseq, arg->lseq, arg->prefix);
	}
    if (! arg->silent) fflush(stdout);
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg)
{
	int c;
	const char usage[] = "\narmi [--mali ... || --random ] [OPTIONS ...]\n\
	INPUT SEQUENCE\n\
	   --mali <input alignment>\t(mode: optional, type: char  , default: void)\n\
	   --nsubset <number of seqs.>\t(mode: optional, type: int   , default: 0)\n\
	RANDOM SEQUENCE\n\
	   --random \t\t\t(mode: optional, type: no_arg, default: off)\n\
	   --nseq <number of seqs.>\t(mode: optional, type: int   , default: 8)\n\
	   --lseq <length of seqs.>\t(mode: optional, type: int   , default: 16)\n\
	OUTPUT\n\
	   --prefix <output prefix>\t(mode: optional, type: char  , default: void)\n\
		-- silent\n\
	HELP\n\
	   --cite\t\t\t(mode: optional, type: no_arg, default: off)\n\
	   --version\t\t\t(mode: optional, type: no_arg, default: off)\n\
	   --help\n";

    set_defaults(arg);

    if (argc < 2) {
		print_header(arg);
        fprintf(stderr, "%s", usage);
		print_license(arg);
        exit(1);
    }

    /** long option definition */
    static struct option long_options[] = {
        {"mali", required_argument, 0, 1},
        {"random", no_argument, 0, 2},
        {"nseq", required_argument, 0, 3},
        {"lseq", required_argument, 0, 4},
        {"nsubset", required_argument, 0, 5},
        {"prefix", required_argument, 0, 6},
        {"nele", required_argument, 0, 7},
        {"silent", no_argument, 0, 19},
        {"cite", no_argument, 0, 20},
        {"version", no_argument, 0, 21},
        {"help", no_argument, 0, 22},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3:4:5:6:7:19 20 21 22", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->infilename = optarg;
                break;
            case 2:
                arg->random = 1;
                break;
            case 3:
                arg->nseq = atoi(optarg);
				break;
            case 4:
                arg->lseq = atoi(optarg);
                break;
            case 5:
                arg->nsubset = atoi(optarg);
				break;
            case 6:
                arg->prefix = optarg;
                break;
            case 7:
                arg->nele = atoi(optarg);
                break;
            case 19:
                arg->silent = 1;
                break;
            case 20:
                print_citation(arg);
                exit(0);
            case 21:
				print_version(arg);
				print_license(arg);
                exit(0);
            case 22:
                fprintf(stderr, "%s", usage);
				print_license(arg);
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license(arg);
                exit(1);
        }
    }

	check_input(arg);
    print_header(arg);
    print_args(arg);

    return 0;
}

