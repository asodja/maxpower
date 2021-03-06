/*
 * lmem_cpu_access.c
 *
 *  Created on: 17 Apr 2015
 *      Author: itay
 */

#include <MaxSLiCInterface.h>
#include "lmem_cpu_access.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>

#include "../../runtime/lmem.h"


#define MAX_NAME_LENGTH 64
#define MAX_STREAM_NAME_LENGTH (MAX_NAME_LENGTH + 64)
#define PAGE_SIZE 4096
#define MAX_NUM_SLOTS 512
#define MAX_BURSTS_PER_CMD 127

#define LMEM_CMD_STREAM_NAME "CpuAccessLMemCommands"

#define LMEM_WRITE_CPU_STREAM_NAME "CpuToLMem"
#define LMEM_READ_CPU_STREAM_NAME "LMemToCpu"

#define LMEM_WRITE_STREAM_NAME "LMemCpuWrite"
#define LMEM_READ_STREAM_NAME "LMemCpuRead"

#define TIMEOUT_SECONDS 5

struct lmem_cpu_access_s {
	max_file_t *maxfile;
	max_engine_t *engine;
	size_t cmd_buffer_size;
	void *cmd_buffer;

	uint32_t burst_size_bytes;

	max_llstream_t *cmd_stream;
	uint16_t to_lmem_stream_id;
	uint16_t from_lmem_stream_id;
	char *write_cpu_stream_name;
	char *read_cpu_stream_name;
};

typedef struct ATTRIB_PACKED {
	lmem_cmd_t cmd1;
	lmem_cmd_t cmd2;
} mem_cmd_stream_slot_t;


size_t lmem_get_burst_size_bytes(lmem_cpu_access_t *handle)
{
	return handle->burst_size_bytes;
}

lmem_cpu_access_t *lmem_init_cpu_access_named(max_file_t *maxfile, max_engine_t *engine,
		const char *cmd_stream_name, const char *write_cpu_stream_name, const char *read_cpu_stream_name,
		const char *write_stream_name, const char *read_stream_name, const char *interface_name)
{
	assert(maxfile != NULL);
	assert(engine != NULL);
	assert(cmd_stream_name != NULL);
	// Must enable at least one of reads or writes
	assert(write_cpu_stream_name != NULL || read_cpu_stream_name != NULL);
	// The two streams must be either both enabled or both disabled
	assert((write_cpu_stream_name != NULL) == (write_stream_name != NULL));
	assert((read_cpu_stream_name != NULL) == (read_stream_name != NULL));

	lmem_cpu_access_t *handle = malloc(sizeof(lmem_cpu_access_t));
	if (handle == NULL) {
		printf("%s: Failed to allocate memory for handle.\n", __func__);
		abort();
	}

	handle->engine = engine;
	handle->maxfile = maxfile;
	handle->burst_size_bytes = max_get_burst_size(maxfile, interface_name);

	handle->cmd_buffer_size = MAX_NUM_SLOTS * sizeof(mem_cmd_stream_slot_t);
	posix_memalign(&handle->cmd_buffer, PAGE_SIZE, handle->cmd_buffer_size);
	handle->cmd_stream = max_llstream_setup(handle->engine, cmd_stream_name,
			MAX_NUM_SLOTS, sizeof(mem_cmd_stream_slot_t), handle->cmd_buffer);

	if(write_stream_name) {
		handle->to_lmem_stream_id = 1 << max_memctl_get_id_within_group(maxfile, interface_name, write_stream_name);
		handle->write_cpu_stream_name = strdup(write_cpu_stream_name);
	} else {
		handle->write_cpu_stream_name = NULL;
	}
	if(read_stream_name) {
		handle->from_lmem_stream_id = 1 << max_memctl_get_id_within_group(maxfile, interface_name, read_stream_name);
		handle->read_cpu_stream_name = strdup(read_cpu_stream_name);
	} else {
		handle->read_cpu_stream_name = NULL;
	}

	return handle;
}

lmem_cpu_access_t *lmem_init_cpu_access(max_file_t *maxfile, max_engine_t *engine)
{
	return lmem_init_cpu_access_named(maxfile, engine, LMEM_CMD_STREAM_NAME, LMEM_WRITE_CPU_STREAM_NAME,
			LMEM_READ_CPU_STREAM_NAME, LMEM_WRITE_STREAM_NAME, LMEM_READ_STREAM_NAME, NULL);
}

static void *acquire_memory_command_slot(lmem_cpu_access_t *handle, int timeout_seconds)
{
	struct timeval timeout = { .tv_sec = timeout_seconds, .tv_usec = 0 };
	struct timeval time_start, time_now, time_delta;

	void *cmd_slot;

	gettimeofday(&time_start, NULL);
	while (max_llstream_write_acquire(handle->cmd_stream, 1, &cmd_slot) != 1) {
		gettimeofday(&time_now, NULL);
		timersub(&time_now, &time_start, &time_delta);

		if (timercmp(&time_delta, &timeout, >=)) {
			printf("%s: Timed-out while waiting for a memory command stream\n", __func__);
			abort();
		} else {
			usleep(1);
		}
	}

	return cmd_slot;
}

static void commit_memory_command_slot(lmem_cpu_access_t *handle)
{
	max_llstream_write(handle->cmd_stream, 1);
}

#define MIN(x, y) ((x) < (y) ? (x) : (y))

static void send_mem_commands(
		lmem_cpu_access_t *handle,
		uint16_t stream_id,
		uint32_t base_address_bursts,
		size_t data_size_bursts)
{
	size_t address   = base_address_bursts;
	size_t remaining = data_size_bursts;

	while (remaining > 0) {
		mem_cmd_stream_slot_t *cmdSlot = acquire_memory_command_slot(handle, TIMEOUT_SECONDS);

		size_t now = MIN(remaining, MAX_BURSTS_PER_CMD);
		// First command in the slot
		cmdSlot->cmd1 = lmem_cmd_data(stream_id, address, now, false, false);

		address   += now;
		remaining -= now;

		now = MIN(remaining, MAX_BURSTS_PER_CMD);
		// Second command in the slot
		cmdSlot->cmd2 = (remaining == 0) ? lmem_cmd_padding() : lmem_cmd_data(stream_id, address, now, false, false);

		address   += now;
		remaining -= now;

		commit_memory_command_slot(handle);
	}
}

void lmem_write(lmem_cpu_access_t *handle, uint32_t address_bursts, const void *data, size_t data_size_bursts)
{
	if(handle->write_cpu_stream_name == NULL) {
		printf("%s: Writing disabled on this handle.\n", __func__);
		abort();
	}
	max_actions_t *actions = max_actions_init(handle->maxfile, NULL);
	max_disable_validation(actions);
	max_disable_reset(actions);
	max_queue_input(actions, handle->write_cpu_stream_name, data, data_size_bursts * handle->burst_size_bytes);
	max_run_t *runContext = max_run_nonblock(handle->engine, actions);

	send_mem_commands(handle, handle->to_lmem_stream_id, address_bursts, data_size_bursts);

	max_wait(runContext);
	max_actions_free(actions);
}

void lmem_read(lmem_cpu_access_t *handle, uint32_t address_bursts, void *data, size_t data_size_bursts)
{
	if(handle->read_cpu_stream_name == NULL) {
		printf("%s: Reading disabled on this handle.\n", __func__);
		abort();
	}
	max_actions_t *actions = max_actions_init(handle->maxfile, NULL);
	max_disable_validation(actions);
	max_disable_reset(actions);
	max_queue_output(actions, handle->read_cpu_stream_name, data, data_size_bursts * handle->burst_size_bytes);
	max_run_t *runContext = max_run_nonblock(handle->engine, actions);

	send_mem_commands(handle, handle->from_lmem_stream_id, address_bursts, data_size_bursts);

	max_wait(runContext);
	max_actions_free(actions);
}
