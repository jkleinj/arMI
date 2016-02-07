/*==============================================================================
getseqs.c : read structural sequence from FASTA file format
(C) 2004 John Romein
(C) 2007 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#include "getseqs.h"

/*____________________________________________________________________________*/
static int read_sequence_name(FILE *file, Seq *sequence)
{
    /* Reads a line starting with '>' and followed by the sequence name */
    /* Leading white space is skipped.  Reading proceeds until the end  */
    /* of the line is reached.  The name is read into sequence->name.   */

    int ch, length = 0, allocated = 64;

    while ((ch = getc(file)) != EOF && isspace(ch))
	;

    if (ch != '>')
		return 0;

    sequence->name = safe_malloc(allocated);

    do {
		sequence ->name[length ++] = ch;

		if (length == allocated)
			sequence->name = safe_realloc(sequence->name, allocated += 64);
    } while ((ch = getc(file)) != EOF && ch != '\n' && isprint(ch));

    sequence->name[length] = '\0';

	return 1;
}

/*____________________________________________________________________________*/
static void read_sequence_residues(FILE *file, Seq *sequence)
{
    /* Reads the residues in a sequence, up to (but not including) the
       next sequence header (starting with '>'), or up to end of file.
       Residues may span multiple lines.  White space and gaps are skipped.
       Nonalpha characters are rejected, resulting in an error message.
       Alpha characters are NOT converted to upper case.  The string is read
       into sequence->residues and zero-terminated; the length is stored
       into sequence->length. */

    int ch, length = 0, allocated = 64;

    sequence->residue = safe_malloc(allocated);

    while ((ch = getc(file)) != EOF && ch != '>')
		if (isalpha(ch) || ch == '@' || ch == '.') {
			/* use '@' as gap character internally */
			if (ch == '.') {
				ch = '@';
			}
			sequence->residue[length ++] = ch;

			if (length == allocated) {
				sequence->residue = safe_realloc(sequence->residue, allocated += 64);
			}
		} else if (!isspace(ch)) {
			fprintf(stderr, "Exiting: illegal character '%c' in protein sequence\n", ch);
			exit(1);
		}

    if (ch == '>')
		ungetc('>', file);

    if (length == 0) {
		printf("zero-sized sequence\n");
		exit(1);
    }

    sequence->residue[length] = '\0';
    sequence->length = length;
}

/*____________________________________________________________________________*/
int read_sequence(FILE *file, Seq *sequence)
{
    if (read_sequence_name(file, sequence)) {
		read_sequence_residues(file, sequence);
		return 1;
    }
	else
		return 0;
}

