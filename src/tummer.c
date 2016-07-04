/**
 * @file
 *
 * This is the main file. It contains functions to parse the commandline
 * arguments, read files etc.
 *
 * @brief The main file
 * @author Fabian Klötzl
 *
 * @section License
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "global.h"
#include "process.h"
#include "io.h"
#include "sequence.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Global variables */
int FLAGS = F_FORWARD;
int THREADS = 1;
double RANDOM_ANCHOR_PROP = 0.05;
int MIN_LENGTH = 0;

void usage(void);
void version(void);

/**
 * @brief The main function.
 *
 * The main function reads and parses the commandline arguments. Depending on
 * the set flags it reads the input files and forwards the contained sequences
 * to processing.
 */
int main(int argc, char *argv[]) {
	int c;
	int version_flag = 0;

	struct option long_options[] = {
		{"version", no_argument, &version_flag, 1},
		{"help", no_argument, NULL, 'h'},
		{"verbose", no_argument, NULL, 'v'},
		{"join", no_argument, NULL, 'j'},
		{"min-length", required_argument, NULL, 'l'},
		// {"threads", required_argument, NULL, 't'},
		{0, 0, 0, 0}};

#ifdef _OPENMP
	// Use all available processors by default.
	THREADS = omp_get_num_procs();
#endif

	// parse arguments
	while (1) {

		int option_index = 0;

		c = getopt_long(argc, argv, "bhjrvp:l:", long_options, &option_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: break;
			case 'b': FLAGS |= F_FORWARD | F_REVCOMP; break;
			case 'h': usage(); break;
			case 'j': FLAGS |= F_JOIN; break;
			case 'r':
				FLAGS &= ~F_FORWARD;
				FLAGS |= F_REVCOMP;
				break;
			case 'v':
				FLAGS |= FLAGS & F_VERBOSE ? F_EXTRA_VERBOSE : F_VERBOSE;
				break;
			case 'p': {
				errno = 0;
				char *end;
				double prop = strtod(optarg, &end);

				if (errno || end == optarg || *end != '\0') {
					warnx(
						"Expected a floating point number for -p argument, but "
						"'%s' was given. Skipping argument.",
						optarg);
					break;
				}

				if (prop < 0.0 || prop > 1.0) {
					warnx("A probability should be a value between 0 and 1; "
						  "Ignoring -p %f argument.",
						  prop);
					break;
				}

				RANDOM_ANCHOR_PROP = prop;
				break;
			}
			case 'l': {
				errno = 0;
				char *end;
				long unsigned int length = strtoul(optarg, &end, 10);

				if (errno || end == optarg || *end != '\0') {
					warnx("Expected a number for -l argument, but '%s' was "
						  "given. Ignoring -l argument.",
						  optarg);
					break;
				}

				MIN_LENGTH = length;
				break;
			}
			case '?': /* intentional fall-through */
			default: usage(); break;
		}
	}

	if (version_flag) {
		version();
	}

	argc -= optind;
	argv += optind;

	// at least one file name must be given
	if (FLAGS & F_JOIN && argc == 0) {
		errx(1, "In join mode at least one filename needs to be supplied.");
	}

	dsa_t dsa;
	if (dsa_init(&dsa)) {
		errx(errno, "Out of memory.");
	}

	const char *file_name;

	// parse all files
	int minfiles = FLAGS & F_JOIN ? 2 : 1;
	for (;; minfiles--) {
		if (!*argv) {
			if (minfiles <= 0) break;

			// if no files are supplied, read from stdin
			file_name = "-";
		} else {
			file_name = *argv++;
		}

		if (FLAGS & F_JOIN) {
			read_fasta_join(file_name, &dsa);
		} else {
			read_fasta(file_name, &dsa);
		}
	}

	size_t n = dsa_size(&dsa);

	if (n < 2) {
		errx(1,
			 "I am truly sorry, but with less than two sequences (%zu given) "
			 "there is nothing to compare.",
			 n);
	}

	// Warn about non ACGT residues.
	if (FLAGS & F_NON_ACGT) {
		warnx("The input sequences contained characters other than acgtACGT. "
			  "These were mapped to N to ensure correct results.");
	}

	// validate sequence correctness
	const seq_t *seq = dsa_data(&dsa);
	for (size_t i = 0; i < n; ++i, ++seq) {

		// The length limit should only apply to the reference
		const size_t LENGTH_LIMIT = (INT_MAX - 1) / 2;
		if (seq->len > LENGTH_LIMIT) {
			errx(1, "The sequence %s is too long. The technical limit is %zu.",
				 seq->name, LENGTH_LIMIT);
		}

		if (seq->len == 0) {
			errx(1, "The sequence %s is empty.", seq->name);
		}
	}

	if (FLAGS & F_VERBOSE) {
		fprintf(stderr, "Comparing %zu sequences\n", n);
	}

	run(dsa_data(&dsa), n);

	dsa_free(&dsa);
	return 0;
}

/**
 * Prints the usage to stdout and then exits successfully.
 */
void usage(void) {
	const char str[] = {
		"Usage: tummer [-bjvr] [-p FLOAT] [-l INT] FILES...\n"
		"\tFILES... can be any sequence of FASTA files. If no files are "
		"supplied, stdin is used instead. The first provided sequence is used "
		"as the reference.\n"
		"Options:\n"
		"  -b                Compute forward and revere complement matches; "
		"default: forward only\n"
		"  -j, --join        Treat all sequences from one file as a single "
		"genome\n"
		"  -l, --min-length <INT>  Minimum length of a MUM; uses p-value by "
		"default\n"
		"  -p <FLOAT>        Significance of a MUM; default: 0.05\n"
		"  -r                Compute only reverse complement matches; default: "
		"forward only\n"
		"  -v, --verbose     Prints additional information\n"
		"  -h, --help        Display this help and exit\n"
		"      --version     Output version information\n"};

	printf("%s", str);
	exit(EXIT_SUCCESS);
}

/**
 * This function just prints the version string and then aborts
 * the program. It conforms to the [GNU Coding
 * Standard](http://www.gnu.org/prep/standards/html_node/_002d_002dversion.html#g_t_002d_002dversion).
 */
void version(void) {
	const char str[] = {
		"tummer " VERSION "\n"
		"Copyright (C) 2016 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
