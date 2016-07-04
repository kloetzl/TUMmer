/**
 * @file
 * @brief This file contains the definitions for various io methods.
 */
#define _GNU_SOURCE
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <pfasta.h>
#include <compat-string.h>

#include "global.h"
#include "io.h"

/**
 * @brief Joins all sequences from a file into a single long sequence.
 *
 * Apart from reading all sequences from a file, this function also
 * merges them into one long sequence.
 *
 * "I didn't learn joined up handwriting for nothing, you know."
 * ~ Gilderoy Lockhart
 *
 * @param file_name - The name of the file to be used for reading. The name is
 *  also used to infer the sequence name.
 * @param dsa - (output parameter) An array that holds found sequences.
 */
void read_fasta_join(const char *file_name, dsa_t *dsa) {
	if (!file_name || !dsa) return;

	dsa_t single;
	dsa_init(&single);
	read_fasta(file_name, &single);

	if (dsa_size(&single) == 0) {
		return;
	}

	seq_t joined = dsa_join(&single);

	/* In join mode we try to be clever about the sequence name. Given the file
	 * path we extract just the file name. ie. path/file.ext -> file
	 * This obviously fails on Windows.
	 */

	const char *left = strrchr(file_name, '/'); // find the last path separator
	left = (left == NULL) ? file_name : left + 1;
	// left is the position one of to the right of the path separator

	const char *dot = strchrnul(left, '.'); // find the extension

	// copy only the file name, not its path or extension
	joined.name = strndup(left, dot - left);
	CHECK_MALLOC(joined.name);

	dsa_push(dsa, joined);
	dsa_free(&single);
}

/**
 * @brief This function reads sequences from a file.
 * @param file_name - The file to read.
 * @param dsa - (output parameter) An array that holds found sequences.
 */
void read_fasta(const char *file_name, dsa_t *dsa) {
	if (!file_name || !dsa) return;

	int file_descriptor =
		strcmp(file_name, "-") ? open(file_name, O_RDONLY) : STDIN_FILENO;

	if (file_descriptor < 0) {
		warn("%s", file_name);
		return;
	}

	int l;
	int check;

	seq_t top = {};
	pfasta_file pf;

	if ((l = pfasta_parse(&pf, file_descriptor)) != 0) {
		warnx("%s: %s", file_name, pfasta_strerror(&pf));
		goto fail;
	}

	pfasta_seq ps;
	while ((l = pfasta_read(&pf, &ps)) == 0) {
		check = seq_init(&top, ps.seq, ps.name);

		// skip broken sequences
		if (check != 0) continue;

		dsa_push(dsa, top);
		pfasta_seq_free(&ps);
	}

	if (l < 0) {
		warnx("%s: %s", file_name, pfasta_strerror(&pf));
		pfasta_seq_free(&ps);
	}

fail:
	pfasta_free(&pf);
	close(file_descriptor);
}
