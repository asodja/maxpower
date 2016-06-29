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

#define PACKED __attribute__((packed))

typedef struct PACKED configWord_s {
	uint64_t  base;
	uint64_t  wordCount;
} configWord_t;

typedef struct PACKED {
	uint64_t  data[4];
} single_entry_t;

#define PAGE_ALIGNMENT 4096
#define READ_LEN_HW (1024 * 1024)
#define READ_LEN_SIM 1024
#define NUM_READS 128

uint8_t compare(single_entry_t* outputData, uint64_t configBase, size_t len)
{
	uint8_t fail = 0;
	for (size_t entryIx=0; entryIx < len; entryIx++) {
		uint64_t *output = (uint64_t *)outputData[entryIx].data;
		size_t quadsPerEntry = sizeof(single_entry_t) / sizeof(uint64_t);

		uint64_t expected = (configBase + entryIx);
		if (expected != output[0]) {
			fail = 1;
			fprintf(stderr, "[Entry: %zd, Quad: %zd] Mismatch: input 0x%lx, output 0x%lx\n", entryIx, 0L, expected, output[0]);
		}

		for (size_t q = 1; !fail && q < quadsPerEntry; q++) {
			if (0 != output[q]) {
				fail = 1;
				fprintf(stderr, "[Entry: %zd, Quad: %zd] Mismatch: input 0x%lx, output 0x%lx\n", entryIx, q, 0L, output[q]);
			}
		}
	}

	return fail;
}

int main(int argc, char *argv[])
{
	(void) argc;
	(void) argv;

	// init
	max_file_t *maxfile = INIT_NAME();
	if(!maxfile) {
		fprintf(stderr, "Failed to init MAX file\n");
		return -1;
	}

	max_config_set_bool(MAX_CONFIG_PRINTF_TO_STDOUT, true);

	const int isSimulation = max_get_constant_uint64t(maxfile, "IS_SIMULATION") == 1ul;
	const size_t readLen = isSimulation ? READ_LEN_SIM : READ_LEN_HW;

	const char *device_name = "*";
	printf("Opening device: %s\n", device_name);

	max_engine_t *engine = max_load(maxfile, device_name);
	if(!engine) {
		fprintf(stderr, "Failed to open Max device\n");
		exit(-1);
	}

	uint64_t configBase = 0;

	max_actions_t *action = max_actions_init(maxfile, NULL);
	configWord_t configWord;
	configWord.wordCount = NUM_READS * readLen;
	configWord.base = configBase;
	max_queue_input(action, "configWord", &configWord, sizeof(configWord));
	max_run(engine, action);
	max_actions_free(action);
	action = NULL;

	// prompt user
	sleep(1);
	printf("Press ENTER to begin test...\n");
	getchar();

	// test
	printf("Streaming 'read_fifo'...\n");

	single_entry_t *outputData = calloc(readLen, sizeof(single_entry_t));
	assert(outputData != NULL);

	action = max_actions_init(maxfile, NULL);
	max_queue_output(action, "read_fifo", outputData, sizeof(single_entry_t) * readLen);
	max_disable_reset(action);
	max_disable_validation(action); //FIXME Work around for validation issue when reset disabled - ticket #7528

	uint8_t fail = 0;
	printf("Checkout output ");
	for(size_t i = 0; i < NUM_READS; i++ ){
		printf("."); fflush(stdout);
		max_run(engine, action);

		size_t start = i * readLen;
		size_t base = configBase + start;
		fail |= compare(outputData, base, readLen);
	}
	printf("\n");

	max_actions_free(action);
	action = NULL;

	free(outputData);
	outputData = NULL;

	max_unload(engine);
	engine = NULL;

	max_file_free(maxfile);
	maxfile = NULL;

	// report result
	printf("%s\n", fail ? "FAILED!" : "Success");

	return fail;
}
