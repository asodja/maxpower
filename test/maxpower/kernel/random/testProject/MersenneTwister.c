#include "random_mt.h"
#include <limits.h>

extern max_file_t* MersenneTwister_init();

#define N 1000000

int main () {
	unsigned int* mt = malloc(N * sizeof(unsigned int));
	if (mt == NULL) {
		fprintf(stderr, "Was not able to allocate enough memory");
		exit(-1);
	}

	max_file_t    *maxfile = MersenneTwister_init();
	max_actions_t *actions = max_actions_init( maxfile, NULL );
	max_engine_t  *engine  = max_load( maxfile, "*" );

	max_set_ticks   ( actions, "MersenneTwisterKernel", N );
	max_queue_output( actions, "mt", mt, N * sizeof(float) );
	random_mt_init(maxfile, actions, "MersenneTwisterKernel", "mt", time(NULL));
	max_run ( engine, actions );

	max_actions_free( actions );
	max_unload( engine );

	// Pearson's chi squared test
	// Critical value for 99%
	const int n = 100;
	const double critical_value = 118.498;
	int* bucketsmt = calloc(n, sizeof(int));
	if (bucketsmt == NULL) {
		fprintf(stderr, "Was not able to allocate enough memory");
		exit(-1);
	}

	for (int i = 0; i < N; i++) {
		bucketsmt[(int)(((double)mt[i] / ((long)INT_MAX * 2)) * n)]++;
	}

	const int E = N / n;

	double summt = 0;
	for (int i = 0; i < 100; i++) {
		summt += (double) ((bucketsmt[i] - E) * (bucketsmt[i] - E)) / (double) E;
	}

	if (summt <= critical_value) {
		printf("Test Passed!\n");
	} else {
		printf("Not uniformly distributed!\n");
	}

	free(bucketsmt);
	free(mt);

	return 0;
}
