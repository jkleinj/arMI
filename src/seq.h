/*===============================================================================
seq.h : sequence data structure
Copyright (C) 2008 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
================================================================================*/

#ifndef SEQ_H
#define SEQ_H

/*___________________________________________________________________________*/
/** data structures */

/* sequence */
typedef struct  
{
    char *name; /* sequence name */
    char *residue; /* array of residues = sequence */
    int length; /* length of sequence */
} Seq;

#endif

