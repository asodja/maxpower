#include <MaxSLiCInterface.h>
#include <maxhash.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>

#define _CONCAT(x, y) x ## y
#define CONCAT(x, y) _CONCAT(x, y)
#define INIT_NAME CONCAT(DESIGN_NAME, _init)

int main(int argc, char *argv[]) {
	max_file_t *maxfile = INIT_NAME();
	max_engine_t * engine = max_load(maxfile, "*");

	maxhash_engine_state_t *engine_state = malloc(sizeof(maxhash_engine_state_t));
	engine_state->maxfile = maxfile;
	engine_state->engine = engine;
	max_config_set_bool(MAX_CONFIG_PRINTF_TO_STDOUT, true);

	maxhash_table_t *table;

	if (maxhash_hw_table_init(&table, "maxHashTestKernel", "maxHashTestTable", engine_state) != MAXHASH_ERR_OK) {
		fprintf(stderr, "Failed to initialize hash-table\n");
	}

	unsigned int num_values = 256;
	uint32_t keys[num_values];
	uint32_t values[num_values], results[num_values];

	for(unsigned int n = 0; n < num_values; n++) {
		keys[n] = (uint32_t) rand();
		values[n] = (uint32_t) rand();
		if(maxhash_put(table, &keys[n], sizeof(keys[n]), &values[n], sizeof(values[n])) != MAXHASH_ERR_OK) {
			fprintf(stderr, "Failed to put entry %u (key %x, value %x)\n", n, keys[n], values[n]);
			return 1;
		}
	}

	if(maxhash_commit(table) != MAXHASH_ERR_OK) {
		fprintf(stderr, "Failed to commit table\n");
		return 1;
	}

	sleep(2);

	max_actions_t *action = max_actions_init(maxfile, NULL);
	max_disable_validation(action);
	max_set_ticks(action, "maxHashTestKernel", num_values);
	max_queue_input(action, "keyFromCPU", keys, sizeof(keys[0])*num_values);
	max_queue_output(action, "valueToCPU", results, sizeof(results[0])*num_values);
	max_run(engine, action);
	max_actions_free(action);

	sleep(2);

	max_unload(engine);
	max_file_free(maxfile);
	maxhash_free(table);

	unsigned int errors_found = 0;
	for(unsigned int n = 0; n < num_values; n++) {
		//printf("%3u %8x %8x %8x\n", n, keys[n], values[n], results[n]);
		if(values[n] != results[n]) {
			errors_found++;
			fprintf(stderr, "Error at entry %u with key %x: Returned value %x does not match expected value %x\n",
					n, keys[n], results[n], values[n]);
		}
	}
	if(errors_found > 0) {
		printf("%u errors found (in %u keys)\n", errors_found, num_values);
	} else {
		printf("Test successful\n");
	}

	return errors_found != 0;
}
