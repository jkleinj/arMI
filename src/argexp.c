/*==============================================================================
argexp.c : parse command line arguments
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "argexp.h"

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
	fprintf(stdout, "\nversion %s\n", VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header()
{
    fprintf(stdout, "\narmiexp: analytical residual MI between pairs of gene expression samples\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stdout, "\nCopyright (C) 2015-2016 Jens Kleinjung and Ton Coolen\n"
			"This program is free software and comes with ABSOLUTELY NO WARRANTY.\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the COPYING file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
	fprintf(stdout, "\nunpublished\n\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg)
{
    arg->exprfilename = "";
    arg->rowfilename = "";
    arg->colfilename = "";
	arg->nsubset = 0; /* number of genes in subset */
	arg->prefix = ""; /* output prefix */
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg)
{
	assert((strlen(arg->exprfilename) > 0) && "Need an expression file");
	assert((strlen(arg->rowfilename) > 0) && "Need a row name file");
	assert((strlen(arg->colfilename) > 0) && "Need a column name file");
	assert(arg->nsubset >= 0 && "Number of genes in subset must be >0");
}

/*____________________________________________________________________________*/
/** print command line arguments */
static void print_args(Arg *arg)
{
    time_t now;
    time(&now);

    fprintf(stdout, "\ndate: %s", ctime(&now));
    fprintf(stdout, \
                    "expression file: %s\n"
                    "row name file: %s\n"
                    "column name file: %s\n"
                    "nsubset: %d\n"
                    "prefix: %s\n"
                    "GSL factorial maximum: %d\n"
                    "long double maximum: %Le\n",
        arg->exprfilename, arg->rowfilename, arg->colfilename,
		arg->nsubset, arg->prefix,
		GSL_SF_FACT_NMAX, LDBL_MAX);
    fflush(stdout);
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg)
{
	int c;
	const char usage[] = "\narmi [--expr ... --row ... --col ... ] [OPTIONS ...]\n\
	INPUT SEQUENCE\n\
	   --expr <expression values>\t(mode: mandatory, type: char  , default: void)\n\
	   --row <column names>\t(mode: mandatory, type: char  , default: void)\n\
	   --col <row names>\t(mode: mandatory, type: char  , default: void)\n\
	   --nsubset <number of genes>\t(mode: optional, type: int   , default: 0)\n\
	OUTPUT\n\
	   --prefix <output prefix>\t(mode: optional, type: char  , default: void)\n\
	HELP\n\
	   --cite\t\t\t(mode: optional, type: no_arg, default: off)\n\
	   --version\t\t\t(mode: optional, type: no_arg, default: off)\n\
	   --help\n";

    set_defaults(arg);

    if (argc < 4) {
		print_header();
        fprintf(stderr, "%s", usage);
		print_license();
        exit(1);
    }

    /** long option definition */
    static struct option long_options[] = {
        {"expr", required_argument, 0, 1},
        {"row", required_argument, 0, 2},
        {"col", required_argument, 0, 3},
        {"nsubset", required_argument, 0, 5},
        {"prefix", required_argument, 0, 6},
        {"cite", no_argument, 0, 20},
        {"version", no_argument, 0, 21},
        {"help", no_argument, 0, 22},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3:5:6:20 21 22", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->exprfilename = optarg;
                break;
            case 2:
                arg->rowfilename = optarg;
                break;
            case 3:
                arg->colfilename = optarg;
				break;
            case 5:
                arg->nsubset = atoi(optarg);
				break;
            case 6:
                arg->prefix = optarg;
                break;
            case 20:
                print_citation();
                exit(0);
            case 21:
				print_version();
				print_license();
                exit(0);
            case 22:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(1);
        }
    }

	check_input(arg);
    print_header();
    print_args(arg);

    return 0;
}

