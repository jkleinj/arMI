/*===============================================================================
getseqs.h : Read structural sequence from FASTA file format
(C) 2004 John Romein
(C) 2007 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
================================================================================*/

#ifndef GETSEQS_H
#define GETSEQS_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "safe.h"
#include "seq.h"

/*____________________________________________________________________________*/
/* prototypes */
int read_sequence(FILE *file, Seq *sequence);

#endif
