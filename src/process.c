

/**
 * @file
 * @brief This file contains various distance methods.
 */
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "esa.h"
#include "global.h"
#include "io.h"
#include "process.h"
#include "sequence.h"

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double shuprop(size_t x, double g, size_t l);

/**
 * @brief Calculates the minimum anchor length.
 *
 * Given some parameters calculate the minimum length for anchors according
 * to the distribution from Haubold et al. (2009).
 *
 * @param p - The probability with which an anchor is allowed to be random.
 * @param g - The the relative amount of GC in the subject.
 * @param l - The length of the subject.
 * @returns The minimum length of an anchor.
 */
size_t minAnchorLength(double p, double g, size_t l) {
	size_t x = 1;

	double prop = 0.0;
	while (prop < 1 - p) {
		prop = shuprop(x, g / 2, l);
		x++;
	}

	return x;
}

/**
 * @brief Calculates the binomial coefficient of n and k.
 *
 * We used to use gsl_sf_lnchoose(xx,kk) for this functionality.
 * After all, why implement something that has already been done?
 * Well, the reason is simplicity: GSL is used for only this one
 * function and the input (n<=20) is not even considered big.
 * Hence its much easier to have our own implementation and ditch
 * the GSL dependency even if that means our code is a tiny bit
 * less optimized and slower.
 *
 * @param n - The n part of the binomial coefficient.
 * @param k - analog.
 * @returns (n choose k)
 */
size_t binomial_coefficient(size_t n, size_t k) {
	if (n <= 0 || k > n) {
		return 0;
	}

	if (k == 0 || k == n) {
		return 1;
	}

	if (k > n - k) {
		k = n - k;
	}

	size_t res = 1;

	for (size_t i = 1; i <= k; i++) {
		res *= n - k + i;
		res /= i;
	}

	return res;
}

/**
 * @brief Given `x` this function calculates the probability of a shustring
 * with a length less than `x`.
 *
 * Let X be the longest shortest unique substring (shustring) at any position.
 * Then this function computes P{X <= x} with respect to the given parameter
 * set. See Haubold et al. (2009).
 *
 * @param x - The maximum length of a shustring.
 * @param g - The the half of the relative amount of GC in the DNA.
 * @param l - The length of the subject.
 * @returns The probability of a certain shustring length.
 */
double shuprop(size_t x, double p, size_t l) {
	double xx = (double)x;
	double ll = (double)l;
	size_t k;

	double s = 0.0;

	for (k = 0; k <= x; k++) {
		double kk = (double)k;
		double t = pow(p, kk) * pow(0.5 - p, xx - kk);

		s += pow(2, xx) * (t * pow(1 - t, ll)) *
			 (double)binomial_coefficient(x, k);
		if (s >= 1.0) {
			s = 1.0;
			break;
		}
	}

	return s;
}

/**
 * @param C - The enhanced suffix array of the subject.
 * @param query - The actual query string.
 * @param query_length - The length of the query string. Needed for speed
 * reasons.
 * @param gc - The gc-content of the subject.
 */
void dist_anchor(const esa_s *C, const char *query, size_t query_length,
				 double gc) {
	lcp_inter_t inter;

	size_t last_pos_Q = 0;
	size_t last_pos_S = 0;
	size_t last_length = 0;
	// This variable indicates that the last anchor was the right anchor of a
	// pair.
	size_t last_was_right_anchor = 0;

	size_t this_pos_Q = 0;
	size_t this_pos_S;
	size_t this_length;

	size_t num_right_anchors = 0;

	size_t threshold;
	if (MIN_LENGTH != 0) {
		threshold = MIN_LENGTH;
	} else {
		threshold = minAnchorLength(RANDOM_ANCHOR_PROP, gc, C->len);
	}

	// Iterate over the complete query.
	while (this_pos_Q < query_length) {
		inter =
			get_match_cached(C, query + this_pos_Q, query_length - this_pos_Q);

		this_length = inter.l <= 0 ? 0 : inter.l;

		this_pos_S = C->SA[inter.i];
		while (this_pos_Q > 0 &&
			   query[this_pos_Q - 1] == C->S[this_pos_S - 1]) {
			this_pos_S--;
			this_pos_Q--;
			this_length++;
		}

		if (inter.i == inter.j && this_length >= threshold) {
			printf("%8zu  %8zu  %8zu\n", this_pos_S + 1, this_pos_Q + 1,
				   this_length);
		}

		// Advance
		this_pos_Q += this_length + 1;
	}

	// Very special case: The sequences are identical
	if (last_length >= query_length) {
		printf("%zu %zu %zu\n", last_pos_S, (size_t)0, query_length);
	}
}

/**
 * @param sequences - An array of pointers to the sequences.
 * @param n - The number of sequences.
 */
void run(seq_t *sequences, size_t n) {
	seq_t *subject = &sequences[0];
	esa_s E;

	if (seq_subject_init(subject) || esa_init(&E, subject)) {
		errx(1, "Failed to create index for %s.", subject->name);
	}

	size_t i = 0;
	// now compare every other sequence to the subject
	for (size_t j = 0; j < n; j++) {
		if (j == i) {
			continue;
		}

		// TODO: Provide a nicer progress indicator.
		if (FLAGS & F_EXTRA_VERBOSE) {
#pragma omp critical
			{ fprintf(stderr, "comparing %zu and %zu\n", i, j); }
		}

		size_t ql = sequences[j].len;

		if (FLAGS & F_FORWARD) {
			printf("> %s\n", sequences[j].name);
			dist_anchor(&E, sequences[j].S, ql, subject->gc);
		}

		if (FLAGS & F_REVCOMP) {
			char *R = revcomp(sequences[j].S, ql);

			printf("> %s Reverse\n", sequences[j].name);
			dist_anchor(&E, R, ql, subject->gc);
			free(R);
		}
	}

	esa_free(&E);
	seq_subject_free(subject);
}
