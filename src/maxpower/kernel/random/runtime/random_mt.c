#include <stdio.h>
#include "random_mt.h"

typedef struct {
	char name[512];
	double* data;
	size_t size;
} mapped_rom_contents;

static void* malloc_wrapper(size_t size) {
	void * pointer = malloc(size);
	if (pointer == NULL) {
		fprintf(stderr, "Was not able to allocate enough memory");
		exit(-1);
	}
	return pointer;
}

static long* get_seeds(long seed, int mt_size) {
	long* mt = malloc_wrapper(mt_size * sizeof(long));
	mt[0] = abs(seed);
	for (int i = 1; i < mt_size; i++) {
		mt[i] = abs(((1812433253 * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i) & 0xFFFFFFFFL));
	}

	return mt;
}

static double* get_mapped_rom_content(long seed, int num, int mt_size, int mt_seed_offset) {
	long* seeds = get_seeds(seed + num * mt_seed_offset, mt_size);
	double* seeds_D = malloc_wrapper(mt_size * sizeof(double));
	for (int i = 0; i < mt_size; i++)
		seeds_D[i] = seeds[i];

	free(seeds);
	return seeds_D;
}

/**
 * Creates all seeds for the mapped ROMs
 * @param maxfile Pointer to the maxfile
 * @param rom_data The addressed point is used to store the pointer to the Array with all data about the different ROMs (has to be freed by caller).
 * @param kernel_name Name of the kernel containing the random number generator
 * @param name Name of the random number generator. If gaussianRng gets used add "_mt" to the name
 * @param seed Seed for the random number generator
 *
 * @return Number of generated ROMs (size of the array)
 */
static int get_mapped_roms_mersenne_twister(max_file_t* maxfile, mapped_rom_contents** rom_data, const char* kernel_name, const char* name, const long seed) {
	int mt_size = max_get_constant_uint64t(maxfile, "MersenneTwister_n");
	int mt_seed_offset = max_get_constant_uint64t(maxfile, "MersenneTwister_seedOffset");

	int string_length = snprintf(NULL, 0, "%s%s", kernel_name, name);
	char* full_name = malloc_wrapper((string_length + 1) * sizeof(char));
	snprintf(full_name, string_length + 1, "%s%s", kernel_name, name);

	int count = max_get_constant_uint64t(maxfile, full_name);

	*rom_data = (mapped_rom_contents*) malloc_wrapper(count * sizeof(mapped_rom_contents));

	for (int i = 0; i < count; i++) {
		snprintf((*rom_data)[i].name, sizeof((*rom_data)[i].name), "%s_seeds%d", full_name, i);
		(*rom_data)[i].data = get_mapped_rom_content(seed, i, mt_size, mt_seed_offset);
		(*rom_data)[i].size = mt_size;
	}

	free(full_name);

	return count;
}

void random_mt_init(max_file_t* maxfile, max_actions_t* actions, const char* kernel_name, const char* name, const long seed) {
	mapped_rom_contents* rom_data = NULL;
	int count = get_mapped_roms_mersenne_twister(maxfile, &rom_data, kernel_name, name, seed);

	for (int i = 0; i < count; i++) {
		max_set_mem_range_double(actions, kernel_name, rom_data[i].name, 0, rom_data[i].size, rom_data[i].data);
		free(rom_data[i].data);
	}

	free(rom_data);
}
