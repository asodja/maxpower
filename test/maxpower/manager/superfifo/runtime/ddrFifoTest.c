/*
 * tcpFifoTest.c
 *
 *  Created on: 16 May 2013
 *      Author: itay
 */


#include <MaxSLiCInterface.h>

#define _CONCAT(x, y) x ## y
#define CONCAT(x, y) _CONCAT(x, y)
#define INIT_NAME CONCAT(DESIGN_NAME, _init)

#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>

typedef struct {
	uint8_t data[384];
} single_burst_t;


#define NUM_BURSTS_HW  (1024 * 1024 * 16)
#define NUM_BURSTS_SIM (1024 * 16)

int main(int argc, char *argv[])
{
	(void)argc; (void)argv;

	max_file_t *maxfile = INIT_NAME();
	if(!maxfile) {
		fprintf(stderr, "Failed to init MAX file\n");
		exit(-1);
	}

	max_config_set_bool(MAX_CONFIG_PRINTF_TO_STDOUT, true);

	const int isSimulation = max_get_constant_uint64t(maxfile, "IS_SIMULATION") == 1ul;
	const size_t numBursts = isSimulation ? NUM_BURSTS_SIM : NUM_BURSTS_HW;

	single_burst_t *inputData = malloc(sizeof(single_burst_t) * numBursts);
	single_burst_t *outputData = calloc(numBursts, sizeof(single_burst_t));

	printf("Building data...\n"); fflush(stdout);
	for (size_t burst=0; burst < numBursts; burst++) {
		uint64_t *d = (uint64_t *)inputData[burst].data;
		size_t quadsPerBurst = sizeof(single_burst_t) / sizeof(uint64_t);

		for (size_t q = 0; q < quadsPerBurst; q++) {
			d[q] = burst * quadsPerBurst + q;
		}
	}

	printf("Loading engine...\n"); fflush(stdout);
	max_engine_t *engine = max_load(maxfile, "*");
	if(!engine) {
		fprintf(stderr, "Failed to open Max device\n");
		exit(-1);
	}

	printf("Running test...\n"); fflush(stdout);
	max_actions_t *action = max_actions_init(maxfile, NULL);
	if (!action) {
		fprintf(stderr, "Failed to create action set\n");
		exit(-1);
	}
	max_queue_input(action, "write_fifo", inputData, sizeof(single_burst_t) * numBursts);
	max_queue_output(action, "read_fifo", outputData, sizeof(single_burst_t) * numBursts);

	max_run(engine, action);

	max_actions_free(action);
	max_unload(engine);
	max_file_free(maxfile);

	printf("Comparing...\n"); fflush(stdout);
	uint8_t fail = 0;
	for (size_t burst=0; burst < numBursts; burst++) {
		uint64_t *input = (uint64_t *)inputData[burst].data;
		uint64_t *output = (uint64_t *)outputData[burst].data;
		size_t quadsPerBurst = sizeof(single_burst_t) / sizeof(uint64_t);

		for (size_t q = 0; q < quadsPerBurst; q++) {
			if (input[q] != output[q]) {
				fail = 1;
				printf("[Burst: %zd, Quad: %zd] Mismatch: input 0x%lx, output 0x%lx\n", burst, q, input[q], output[q]);
			}
		}
	}
	free(inputData);
	free(outputData);

	printf("%s\n", fail ? "FAILED!" : "Success");
	return fail;
}
